/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include <cmath>
#include "stdlib.h"
#include "bond_gran.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "fix_property_atom.h"
#include "error.h"
#include "update.h"
#include "vector_liggghts.h"

using namespace LAMMPS_NS;

/*NP
large TODO list for granular bonds:  (could be a diploma thesis?)

+ need a better dissipative formulation than the hardcoded
  'dissipate' value which produces plastic deformation
  need some vel-dependant damping
+ need to carefully debug and validate this bond style
  valiation against fix rigid
+ check whether run this bond style w/ or w/o gran pair style active,
  (neigh_modify command)
+ need to store bond radii per particle, not per type
+ parallel debugging and testing not yet done
+ need evtally implemetation
*/

/* Matt Schramm edits for bond dampening --> MS */
/* Yu Guo edits for bond dampening --> YG */

enum{
     BREAKSTYLE_SIMPLE,
     BREAKSTYLE_STRESS,
     BREAKSTYLE_STRESS_TEMP
    };
enum {
  DAMPSTYLE_NONE = 0,
  DAMPSTYLE_YU_GUO,
  DAMPSTYLE_NON_LINEAR,
  DAMPSTYLE_NON_LINEAR_REDUCED,
  DAMPSTYLE_NON_LINEAR_REDUCED_SIMPLIFIED,
  DAMPSTYLE_NON_LINEAR_REL
};
/* ---------------------------------------------------------------------- */

BondGran::BondGran(LAMMPS *lmp) : Bond(lmp)
{
    // we need 13 history values - the 6 for the forces, 6 for torques from the last time-step and 1 for initial bond length
    n_granhistory(13);
    /*	NP
       number of entries in bondhistlist. bondhistlist[number of bond][number of value (from 0 to number given here)]
       so with this number you can modify how many pieces of information you savae with every bond
       following dependencies and processes for saving,copying,growing the bondhistlist:
    */
     
    /* NP
       gibt groesse der gespeicherten werte  pro bond wieder 
       neighbor.cpp:       memory->create(bondhistlist,maxbond,atom->n_bondhist,"neigh:bondhistlist");
       neigh_bond.cpp:     memory->grow(bondhistlist,maxbond,atom->n_bondhist,"neighbor:bondhistlist");
       bond.cpp: void Bond::n_granhistory(int nhist) {ngranhistory = nhist;     atom->n_bondhist = ngranhistory; if(){FLERR}}
       atom_vec_bond_gran.cpp:  memory->grow(atom->bond_hist,nmax,atom->bond_per_atom,atom->n_bondhist,"atom:bond_hist");
     */
    if(!atom->style_match("bond/gran"))
      error->all(FLERR,"A granular bond style can only be used together with atom style bond/gran");
    if(comm->me == 0)
        error->warning(FLERR,"Bond granular: This is a beta version - be careful!");
    fix_Temp = NULL;
}

/* ---------------------------------------------------------------------- */

BondGran::~BondGran()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(ro);
    memory->destroy(ri);
    memory->destroy(Sn);
    memory->destroy(St);
    memory->destroy(damp);
    memory->destroy(beta0);
    memory->destroy(beta1);
    memory->destroy(r_break);
    memory->destroy(sigma_break);
    memory->destroy(tau_break);
    memory->destroy(T_break);
  }
}

/* ---------------------------------------------------------------------- */

void  BondGran::init_style()
{
    if(breakmode == BREAKSTYLE_STRESS_TEMP)
       fix_Temp = static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp","property/atom","scalar",1,0,"bond gran"));
}

/* ---------------------------------------------------------------------- */
#define SIGNUM_DOUBLE(x) (((x) > 0.0) ? 1.0 : (((x) < 0.0) ? -1.0 : 0.0))
void BondGran::compute(int eflag, int vflag)
{

  double rsq,r,rinv,rsqinv;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3,vtr1,vtr2,vtr3,tor1,tor2,tor3; //,vrel
  double wnnr,wn1,wn2,wn3,wt1,wt2,wt3;
  // double fs1,fs2,fs3;

  int i1,i2,n,type;
  double delx,dely,delz; // ,ebond;
  // double dnforce[3];
  double dtforce[3];
  double dntorque[3],dttorque[3];
  double force_damp_n[3] = {0.0};
  double force_damp_t[3] = {0.0};
  double torque_damp_n[3] = {0.0};
  double torque_damp_t[3] = {0.0};
  double rot;
  double A,J;
  
  double Ip, Me, I, Js; // MS
  double *density = atom->density; //MS
  double Kn,Kt,K_ben,K_tor; //MS
  double sndt, stdt, K_tor_dt, K_ben_dt; //MS
  double rin,rout,m1,m2; //MS
  double J1, J2, bondLength; //YG
  double ft_bond_total[3], fn_bond[3], vel_temp[3]; //YG
  double vel_norm, f_norm; //YG
  double nforce_mag, tforce_mag, ntorque_mag, ttorque_mag;
  bool nstress, tstress, toohot;

  // ebond = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double *radius = atom->radius;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;

  double **omega = atom->omega;
  double **torque = atom->torque;

  // int *tag = atom->tag; // tag of atom is their ID number
  int **bondlist = neighbor->bondlist;
  double **bondhistlist = neighbor->bondhistlist;

  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  double dt = update->dt;
  double cutoff=neighbor->skin;

  if(breakmode == BREAKSTYLE_STRESS_TEMP)
  {
    if(!fix_Temp) error->all(FLERR, "Internal error in BondGran");
    Temp = fix_Temp->vector_atom;
  }
  
  for (n = 0; n < nbondlist; n++) { // Loop through Bond list
    //1st check if bond is broken,
    if(bondlist[n][3])
    {
      //printf("bond %d allready broken\n",n);
      continue;
    }

    i1 = bondlist[n][0]; // Sphere 1
    i2 = bondlist[n][1]; // Sphere 2

    //2nd check if bond overlap the box-borders
    if (x[i1][0]<(domain->boxlo[0]+cutoff)) {
      bondlist[n][3]=1;
      continue;
    } else if (x[i1][0]>(domain->boxhi[0]-cutoff)) {
      bondlist[n][3]=1;
      continue;
    } else if (x[i1][1]<(domain->boxlo[1]+cutoff)) {
      bondlist[n][3]=1;
      continue;
    } else if (x[i1][1]>(domain->boxhi[1]-cutoff)) {
      bondlist[n][3]=1;
      continue;
    } else if (x[i1][2]<(domain->boxlo[2]+cutoff)) {
      bondlist[n][3]=1;
      continue;
    } else if (x[i1][2]>(domain->boxhi[2]-cutoff)) {
      bondlist[n][3]=1;
      continue;
    }

    if (x[i2][0]<(domain->boxlo[0]+cutoff)) {
      bondlist[n][3]=1;
      continue;
    } else if (x[i2][0]>(domain->boxhi[0]-cutoff)) {
      bondlist[n][3]=1;
      continue;
    } else if (x[i2][1]<(domain->boxlo[1]+cutoff)) {
      bondlist[n][3]=1;
      continue;
    } else if (x[i2][1]>(domain->boxhi[1]-cutoff)) {
      bondlist[n][3]=1;
      continue;
    } else if (x[i2][2]<(domain->boxlo[2]+cutoff)) {
      bondlist[n][3]=1;
      continue;
    } else if (x[i2][2]>(domain->boxhi[2]-cutoff)) {
      bondlist[n][3]=1;
      continue;
    }

    type = bondlist[n][2]; // Get current bond type properties      

    rin = ri[type]*fmin(radius[i1],radius[i2]); 
    rout= ro[type]*fmin(radius[i1],radius[i2]);

    A = M_PI * (rout*rout - rin*rin); // Bond Area
    J = 0.5*M_PI*(rout*rout*rout*rout - rin*rin*rin*rin);

    m1 = 4.1887902047863909846168578443*density[i1]*radius[i1]*radius[i1]*radius[i1];
    m2 = 4.1887902047863909846168578443*density[i2]*radius[i2]*radius[i2]*radius[i2];
    Me = m1*m2/(m1+m2);

    Ip = 0.5*M_PI*(rout*rout*rout*rout - rin*rin*rin*rin); // MS
    I  = 0.5*Ip;

    J1 = 0.4 * m1 * radius[i1]*radius[i1];
    J2 = 0.4 * m2 * radius[i2]*radius[i2];
    Js = J1*J2/(J1+J2);

    delx = x[i1][0] - x[i2][0]; // x-directional seperation
    dely = x[i1][1] - x[i2][1]; // y-directional seperation
    delz = x[i1][2] - x[i2][2]; // z-directional seperation 
    domain->minimum_image(delx,dely,delz); // Not 100% sure what this is... 

    rsq = delx*delx + dely*dely + delz*delz; 
    rsqinv = 1./rsq;
    r = sqrt(rsq);
    rinv = 1./r;

    // Check if bond just formed and set eq distance
    if (bondhistlist[n][12] == 0.0) {
      bondhistlist[n][12] = r;
#     ifdef LIGGGHTS_DEBUG
        fprintf(screen, "INFO: Setting bond length between %i and %i at %g\n", i1, i2, bondhistlist[n][12]);
#     endif
    }

    // Set Bond Length
    bondLength = fabs(bondhistlist[n][12]);
    double bondLengthInv = 1.0/bondLength;
    if (bondLength < 1.0e-10) {
      fprintf(screen,"BondLength = %g\n",bondLength);
      error->all(FLERR,"bondlength too small\n");
    }

    // Set Stiffness Values
    Kn = Sn[type]*A*bondLengthInv;
    Kt = St[type]*A*bondLengthInv;
    K_tor = St[type]*Ip*bondLengthInv;
    K_ben = Sn[type]*I*bondLengthInv;

    // relative translational velocity
    vr1 = v[i1][0] - v[i2][0]; // relative velocity between sphere1 and sphere2 in the x-direction
    vr2 = v[i1][1] - v[i2][1]; // relative velocity between sphere1 and sphere2 in the y-direction
    vr3 = v[i1][2] - v[i2][2]; // relative velocity between sphere1 and sphere2 in the z-direction

    // normal component
    vnnr = vr1*delx + vr2*dely + vr3*delz;
    vnnr *= rsqinv;
    vn1 = delx*vnnr;
    vn2 = dely*vnnr;
    vn3 = delz*vnnr;

    // tangential component at the center of the sphere
    vt1 = vr1 - vn1;
    vt2 = vr2 - vn2;
    vt3 = vr3 - vn3;

    // relative rotational velocity for shear
    wr1 = (radius[i1]*omega[i1][0] + radius[i2]*omega[i2][0]) * rinv;
    wr2 = (radius[i1]*omega[i1][1] + radius[i2]*omega[i2][1]) * rinv;
    wr3 = (radius[i1]*omega[i1][2] + radius[i2]*omega[i2][2]) * rinv;

    // relative velocities for shear at contact
    vtr1 = vt1 - (delz*wr2-dely*wr3);
    vtr2 = vt2 - (delx*wr3-delz*wr1);
    vtr3 = vt3 - (dely*wr1-delx*wr2);

    // relative angular velocity
    wr1 = omega[i1][0] - omega[i2][0];
    wr2 = omega[i1][1] - omega[i2][1];
    wr3 = omega[i1][2] - omega[i2][2];

    // normal component
    wnnr = wr1*delx + wr2*dely + wr3*delz;
    wnnr *= rsqinv;
    wn1 = delx*wnnr;
    wn2 = dely*wnnr;
    wn3 = delz*wnnr;

    // tangential component
    wt1 = wr1 - wn1;
    wt2 = wr2 - wn2;
    wt3 = wr3 - wn3;

    // calc change in bond normal forces (Can be directly calculated does not need history!!!)
    sndt = Kn * (r-bondLength)*rinv;
    fn_bond[0] = - sndt*delx;
    fn_bond[1] = - sndt*dely;
    fn_bond[2] = - sndt*delz; 

    // calc change in shear forces
    stdt = Kt*dt; 
    dtforce[0] = -vtr1*stdt;
    dtforce[1] = -vtr2*stdt;
    dtforce[2] = -vtr3*stdt;

    // calc change in normal torque
    K_tor_dt = K_tor*dt;
    dntorque[0] = -wn1*K_tor_dt;
    dntorque[1] = -wn2*K_tor_dt;
    dntorque[2] = -wn3*K_tor_dt;

    // calc change in tang torque
    K_ben_dt = K_ben*dt; // K_ben will become an input parameter
    dttorque[0] = -wt1*K_ben_dt;
    dttorque[1] = -wt2*K_ben_dt;
    dttorque[2] = -wt3*K_ben_dt;

    // rotate forces from previous time step

    //rotate tangential force
    rot = bondhistlist[n][3]*delx + bondhistlist[n][4]*dely + bondhistlist[n][5]*delz;
    rot *= rsqinv;
    vel_temp[0] = bondhistlist[n][3] - rot*delx;
    vel_temp[1] = bondhistlist[n][4] - rot*dely;
    vel_temp[2] = bondhistlist[n][5] - rot*delz;
    vel_norm = sqrt (vel_temp[0]*vel_temp[0]+vel_temp[1]*vel_temp[1]+vel_temp[2]*vel_temp[2]);
    f_norm = bondhistlist[n][3]*bondhistlist[n][3] + bondhistlist[n][4]*bondhistlist[n][4] + bondhistlist[n][5]*bondhistlist[n][5];
    if (vel_norm == 0)
      f_norm =0;
    else 
      f_norm = sqrt(f_norm) /vel_norm;

    bondhistlist[n][3] = f_norm*vel_temp[0];
    bondhistlist[n][4] = f_norm*vel_temp[1];
    bondhistlist[n][5] = f_norm*vel_temp[2];

    //rotate normal torque
    rot = bondhistlist[n][6]*delx + bondhistlist[n][7]*dely + bondhistlist[n][8]*delz;
    rot *= rsqinv;
    vel_temp[0] = rot*delx;
    vel_temp[1] = rot*dely;
    vel_temp[2] = rot*delz;
    vel_norm = sqrt (vel_temp[0]*vel_temp[0]+vel_temp[1]*vel_temp[1]+vel_temp[2]*vel_temp[2]);
    f_norm = bondhistlist[n][6]*bondhistlist[n][6] + bondhistlist[n][7]*bondhistlist[n][7] + bondhistlist[n][8]*bondhistlist[n][8];
    if (vel_norm == 0)
      f_norm =0;
    else
      f_norm = sqrt (f_norm) /vel_norm;

    bondhistlist[n][6] = f_norm*vel_temp[0];
    bondhistlist[n][7] = f_norm*vel_temp[1];
    bondhistlist[n][8] = f_norm*vel_temp[2];

    //rotate tangential torque
    rot = bondhistlist[n][9]*delx + bondhistlist[n][10]*dely + bondhistlist[n][11]*delz;
    rot *= rsqinv;
    vel_temp[0] = bondhistlist[n][9] - rot*delx;
    vel_temp[1] = bondhistlist[n][10] - rot*dely;
    vel_temp[2] = bondhistlist[n][11] - rot*delz;
    vel_norm = sqrt (vel_temp[0]*vel_temp[0]+vel_temp[1]*vel_temp[1]+vel_temp[2]*vel_temp[2]);
    f_norm = bondhistlist[n][9]*bondhistlist[n][9] + bondhistlist[n][10]*bondhistlist[n][10] + bondhistlist[n][11]*bondhistlist[n][11];
    if (vel_norm == 0)
      f_norm =0;
    else
      f_norm = sqrt (f_norm) /vel_norm;

    bondhistlist[n][ 9] = f_norm*vel_temp[0];
    bondhistlist[n][10] = f_norm*vel_temp[1];
    bondhistlist[n][11] = f_norm*vel_temp[2];

    //increment normal and tangential force and torque 
    bondhistlist[n][ 0] = fn_bond[0];
    bondhistlist[n][ 1] = fn_bond[1];
    bondhistlist[n][ 2] = fn_bond[2];
    bondhistlist[n][ 3] = bondhistlist[n][ 3] +  dtforce[0] ;
    bondhistlist[n][ 4] = bondhistlist[n][ 4] +  dtforce[1] ;
    bondhistlist[n][ 5] = bondhistlist[n][ 5] +  dtforce[2] ;
    bondhistlist[n][ 6] = bondhistlist[n][ 6] + dntorque[0] ;
    bondhistlist[n][ 7] = bondhistlist[n][ 7] + dntorque[1] ;
    bondhistlist[n][ 8] = bondhistlist[n][ 8] + dntorque[2] ;
    bondhistlist[n][ 9] = bondhistlist[n][ 9] + dttorque[0] ;
    bondhistlist[n][10] = bondhistlist[n][10] + dttorque[1] ;
    bondhistlist[n][11] = bondhistlist[n][11] + dttorque[2] ;

    double tot_force_x = bondhistlist[n][ 0] + bondhistlist[n][ 3];
    double tot_force_y = bondhistlist[n][ 1] + bondhistlist[n][ 4];
    double tot_force_z = bondhistlist[n][ 2] + bondhistlist[n][ 5];
    double tot_torqu_x = bondhistlist[n][ 6] + bondhistlist[n][ 9];
    double tot_torqu_y = bondhistlist[n][ 7] + bondhistlist[n][10];
    double tot_torqu_z = bondhistlist[n][ 8] + bondhistlist[n][11];

//---torque due to tangential bond force 
    ft_bond_total[0] = bondhistlist[n][3];
    ft_bond_total[1] = bondhistlist[n][4];
    ft_bond_total[2] = bondhistlist[n][5];

    // Apply Damping
    if (dampmode == DAMPSTYLE_NONE) {

        force_damp_n[0] = force_damp_n[1] = force_damp_n[2] = 0.0;
        force_damp_t[0] = force_damp_t[1] = force_damp_t[2] = 0.0;
        torque_damp_n[0] = torque_damp_n[1] = torque_damp_n[2] = 0.0;
        torque_damp_t[0] = torque_damp_t[1] = torque_damp_t[2] = 0.0;

    } else if (dampmode == DAMPSTYLE_YU_GUO) { // of the form fd = 2*b*sqrt(k*m)*v

        double damp2 = 2.0 * damp[type];
        // normal force dampening
        double d_fn_sqrt_2_Me_Sn = -damp2 * sqrt(Me*Kn);
        force_damp_n[0] = d_fn_sqrt_2_Me_Sn*vn1;
        force_damp_n[1] = d_fn_sqrt_2_Me_Sn*vn2;
        force_damp_n[2] = d_fn_sqrt_2_Me_Sn*vn3;

        // tangential force dampening
        double d_ft_sqrt_2_Me_St = -damp2 * sqrt(Me*Kt);
        force_damp_t[0] = d_ft_sqrt_2_Me_St*vtr1;
        force_damp_t[1] = d_ft_sqrt_2_Me_St*vtr2;
        force_damp_t[2] = d_ft_sqrt_2_Me_St*vtr3;

        // normal moment dampening
        double d_mn_sqrt_2_Js_Ktor = -damp2 * sqrt(Js*K_tor);
        torque_damp_n[0] = d_mn_sqrt_2_Js_Ktor*wn1;
        torque_damp_n[1] = d_mn_sqrt_2_Js_Ktor*wn2;
        torque_damp_n[2] = d_mn_sqrt_2_Js_Ktor*wn3;

        // tangential moment dampening
        double d_mt_sqrt_2_Js_Kben = -damp2 * sqrt(Js*K_ben);
        torque_damp_t[0] = d_mt_sqrt_2_Js_Kben*wt1;
        torque_damp_t[1] = d_mt_sqrt_2_Js_Kben*wt2;
        torque_damp_t[2] = d_mt_sqrt_2_Js_Kben*wt3;

    } else if (dampmode == DAMPSTYLE_NON_LINEAR) {
      // Get energy in bond F(gamma) = c*(2.0*E/k)^gamma
      double I_n = 0.5*Me*(rout*rout+rin*rin);
      double I_t = Me*(3.0*(rout*rout+rin*rin)+bondLength*bondLength)/12.0;
      double Kn_inv = 1.0/Kn;
      double Kt_inv = 1.0/Kt;
      double K_tor_inv = 1.0/K_tor;
      double K_ben_inv = 1.0/K_ben;

      double beta0dt = beta0[type]*dt;
      double beta2dtInv = damp[type]/dt;

      // Get 2E/k for each dirrection
      // E = 0.5*Kx^2 + 0.5*mv^2
      // E = 0.5(K^2x^2/K + mv^2)
      // E = 0.5(bondhist^2/K + mv^2)
      // 2E/k == (bondhist^2/k)/k == bondhist^2/K^2
      double E_force_norm_x = Kn_inv*(Kn_inv*bondhistlist[n][ 0]*bondhistlist[n][ 0] + Me*vn1*vn1);
      double E_force_norm_y = Kn_inv*(Kn_inv*bondhistlist[n][ 1]*bondhistlist[n][ 1] + Me*vn2*vn2);
      double E_force_norm_z = Kn_inv*(Kn_inv*bondhistlist[n][ 2]*bondhistlist[n][ 2] + Me*vn3*vn3);

      double E_force_tang_x = Kt_inv*(Kt_inv*bondhistlist[n][ 3]*bondhistlist[n][ 3] + Me*vtr1*vtr1);
      double E_force_tang_y = Kt_inv*(Kt_inv*bondhistlist[n][ 4]*bondhistlist[n][ 4] + Me*vtr2*vtr2);
      double E_force_tang_z = Kt_inv*(Kt_inv*bondhistlist[n][ 5]*bondhistlist[n][ 5] + Me*vtr3*vtr3);

      double E_torqu_norm_x = K_tor_inv*(K_tor_inv*bondhistlist[n][ 6]*bondhistlist[n][ 6] + I_n*wn1*wn1);
      double E_torqu_norm_y = K_tor_inv*(K_tor_inv*bondhistlist[n][ 7]*bondhistlist[n][ 7] + I_n*wn2*wn2);
      double E_torqu_norm_z = K_tor_inv*(K_tor_inv*bondhistlist[n][ 8]*bondhistlist[n][ 8] + I_n*wn3*wn3);

      double E_torqu_tang_x = K_ben_inv*(K_ben_inv*bondhistlist[n][ 9]*bondhistlist[n][ 9] + I_t*wt1*wt1);
      double E_torqu_tang_y = K_ben_inv*(K_ben_inv*bondhistlist[n][10]*bondhistlist[n][10] + I_t*wt2*wt2);
      double E_torqu_tang_z = K_ben_inv*(K_ben_inv*bondhistlist[n][11]*bondhistlist[n][11] + I_t*wt3*wt3);

      force_damp_n[0] = -Kn*(beta0dt + beta1[type]*Me*sqrt(E_force_norm_x) + beta2dtInv*Me*E_force_norm_x)/Me;
      force_damp_n[1] = -Kn*(beta0dt + beta1[type]*Me*sqrt(E_force_norm_y) + beta2dtInv*Me*E_force_norm_y)/Me;
      force_damp_n[2] = -Kn*(beta0dt + beta1[type]*Me*sqrt(E_force_norm_z) + beta2dtInv*Me*E_force_norm_z)/Me;

      force_damp_t[0] = -Kt*(beta0dt + beta1[type]*Me*sqrt(E_force_tang_x) + beta2dtInv*Me*E_force_tang_x)/Me;
      force_damp_t[1] = -Kt*(beta0dt + beta1[type]*Me*sqrt(E_force_tang_y) + beta2dtInv*Me*E_force_tang_y)/Me;
      force_damp_t[2] = -Kt*(beta0dt + beta1[type]*Me*sqrt(E_force_tang_z) + beta2dtInv*Me*E_force_tang_z)/Me;

      torque_damp_n[0] = -K_tor*(beta0dt + beta1[type]*Js*sqrt(E_torqu_norm_x) + beta2dtInv*Js*E_torqu_norm_x)/Js;
      torque_damp_n[1] = -K_tor*(beta0dt + beta1[type]*Js*sqrt(E_torqu_norm_y) + beta2dtInv*Js*E_torqu_norm_y)/Js;
      torque_damp_n[2] = -K_tor*(beta0dt + beta1[type]*Js*sqrt(E_torqu_norm_z) + beta2dtInv*Js*E_torqu_norm_z)/Js;

      torque_damp_t[0] = -K_ben*(beta0dt + beta1[type]*Js*sqrt(E_torqu_tang_x) + beta2dtInv*Js*E_torqu_tang_x)/Js;
      torque_damp_t[1] = -K_ben*(beta0dt + beta1[type]*Js*sqrt(E_torqu_tang_y) + beta2dtInv*Js*E_torqu_tang_y)/Js;
      torque_damp_t[2] = -K_ben*(beta0dt + beta1[type]*Js*sqrt(E_torqu_tang_z) + beta2dtInv*Js*E_torqu_tang_z)/Js;

    } else if (dampmode == DAMPSTYLE_NON_LINEAR_REDUCED) {

      // Get energy in bond F(gamma) = c*(2.0*E/k)^gamma
      double I_n = 0.5*Me*(rout*rout+rin*rin);
      double I_t = Me*(3.0*(rout*rout+rin*rin)+bondLength*bondLength)/12.0;
      double Kn_inv = 1.0/Kn;
      double Kt_inv = 1.0/Kt;
      double K_tor_inv = 1.0/K_tor;
      double K_ben_inv = 1.0/K_ben;

      // Get 2E/k for each dirrection
      // E = 0.5*Kx^2 + 0.5*mv^2
      // E = 0.5(K^2x^2/K + mv^2)
      // E = 0.5(bondhist^2/K + mv^2)
      // 2E/k == (bondhist^2/k)/k == bondhist^2/K^2
      double E_force_norm_x = Kn_inv*(Kn_inv*bondhistlist[n][ 0]*bondhistlist[n][ 0] + Me*vn1*vn1);
      double E_force_norm_y = Kn_inv*(Kn_inv*bondhistlist[n][ 1]*bondhistlist[n][ 1] + Me*vn2*vn2);
      double E_force_norm_z = Kn_inv*(Kn_inv*bondhistlist[n][ 2]*bondhistlist[n][ 2] + Me*vn3*vn3);

      double E_force_tang_x = Kt_inv*(Kt_inv*bondhistlist[n][ 3]*bondhistlist[n][ 3] + Me*vtr1*vtr1);
      double E_force_tang_y = Kt_inv*(Kt_inv*bondhistlist[n][ 4]*bondhistlist[n][ 4] + Me*vtr2*vtr2);
      double E_force_tang_z = Kt_inv*(Kt_inv*bondhistlist[n][ 5]*bondhistlist[n][ 5] + Me*vtr3*vtr3);

      double E_torqu_norm_x = K_tor_inv*(K_tor_inv*bondhistlist[n][ 6]*bondhistlist[n][ 6] + I_n*wn1*wn1);
      double E_torqu_norm_y = K_tor_inv*(K_tor_inv*bondhistlist[n][ 7]*bondhistlist[n][ 7] + I_n*wn2*wn2);
      double E_torqu_norm_z = K_tor_inv*(K_tor_inv*bondhistlist[n][ 8]*bondhistlist[n][ 8] + I_n*wn3*wn3);

      double E_torqu_tang_x = K_ben_inv*(K_ben_inv*bondhistlist[n][ 9]*bondhistlist[n][ 9] + I_t*wt1*wt1);
      double E_torqu_tang_y = K_ben_inv*(K_ben_inv*bondhistlist[n][10]*bondhistlist[n][10] + I_t*wt2*wt2);
      double E_torqu_tang_z = K_ben_inv*(K_ben_inv*bondhistlist[n][11]*bondhistlist[n][11] + I_t*wt3*wt3);

      force_damp_n[0] = -Kn*damp[type]*sqrt(E_force_norm_x);
      force_damp_n[1] = -Kn*damp[type]*sqrt(E_force_norm_y);
      force_damp_n[2] = -Kn*damp[type]*sqrt(E_force_norm_z);

      force_damp_t[0] = -Kt*damp[type]*sqrt(E_force_tang_x);
      force_damp_t[1] = -Kt*damp[type]*sqrt(E_force_tang_y);
      force_damp_t[2] = -Kt*damp[type]*sqrt(E_force_tang_z);

      torque_damp_n[0] = -K_tor*damp[type]*sqrt(E_torqu_norm_x);
      torque_damp_n[1] = -K_tor*damp[type]*sqrt(E_torqu_norm_y);
      torque_damp_n[2] = -K_tor*damp[type]*sqrt(E_torqu_norm_z);

      torque_damp_t[0] = -K_ben*damp[type]*sqrt(E_torqu_tang_x);
      torque_damp_t[1] = -K_ben*damp[type]*sqrt(E_torqu_tang_y);
      torque_damp_t[2] = -K_ben*damp[type]*sqrt(E_torqu_tang_z);

    } else if (dampmode == DAMPSTYLE_NON_LINEAR_REDUCED_SIMPLIFIED) {

      // Get energy in bond F(gamma) = c*(2.0*E/k)^gamma
      // Get 2E/k for each dirrection
      // E = 0.5*Kx^2 + 0.5*mv^2
      // E = 0.5(K^2x^2/K + mv^2)
      // E = 0.5(bondhist^2/K + mv^2)
      // 2E/k == (bondhist^2/k)/k == bondhist^2/K^2
      double E_force_norm_x = fabs(bondhistlist[n][ 0]);
      double E_force_norm_y = fabs(bondhistlist[n][ 1]);
      double E_force_norm_z = fabs(bondhistlist[n][ 2]);

      double E_force_tang_x = fabs(bondhistlist[n][ 3]);
      double E_force_tang_y = fabs(bondhistlist[n][ 4]);
      double E_force_tang_z = fabs(bondhistlist[n][ 5]);

      double E_torqu_norm_x = fabs(bondhistlist[n][ 6]);
      double E_torqu_norm_y = fabs(bondhistlist[n][ 7]);
      double E_torqu_norm_z = fabs(bondhistlist[n][ 8]);

      double E_torqu_tang_x = fabs(bondhistlist[n][ 9]);
      double E_torqu_tang_y = fabs(bondhistlist[n][10]);
      double E_torqu_tang_z = fabs(bondhistlist[n][11]);

      force_damp_n[0] = -damp[type]*E_force_norm_x;
      force_damp_n[1] = -damp[type]*E_force_norm_y;
      force_damp_n[2] = -damp[type]*E_force_norm_z;

      force_damp_t[0] = -damp[type]*E_force_tang_x;
      force_damp_t[1] = -damp[type]*E_force_tang_y;
      force_damp_t[2] = -damp[type]*E_force_tang_z;

      torque_damp_n[0] = -damp[type]*E_torqu_norm_x;
      torque_damp_n[1] = -damp[type]*E_torqu_norm_y;
      torque_damp_n[2] = -damp[type]*E_torqu_norm_z;

      torque_damp_t[0] = -damp[type]*E_torqu_tang_x;
      torque_damp_t[1] = -damp[type]*E_torqu_tang_y;
      torque_damp_t[2] = -damp[type]*E_torqu_tang_z;

    } else if (dampmode == DAMPSTYLE_NON_LINEAR_REL) {
      // Get energy in bond F(gamma) = c*(2.0*E/k)^gamma
      double I_n = 0.5*Me*(rout*rout+rin*rin);
      double I_t = Me*(3.0*(rout*rout+rin*rin)+bondLength*bondLength)/12.0;
      double Kn_inv = 1.0/Kn;
      double Kt_inv = 1.0/Kt;
      double K_tor_inv = 1.0/K_tor;
      double K_ben_inv = 1.0/K_ben;
      double Me_inv = 1.0/Me;
      double Js_inv = 1.0/Js;

      double beta0dt = dt*beta0[type];
      double beta2dtInv = damp[type]/dt;

      // Get 2E/k for each dirrection
      // E = 0.5*Kx^2 + 0.5*mv^2
      // E = 0.5(K^2x^2/K + mv^2)
      // E = 0.5(bondhist^2/K + mv^2)
      // 2E/k == (bondhist^2/k + mv^2)/k

      // friction values do not change... pre compute them.
      double fdn = -Kn*beta0dt*Me_inv;
      double fdt = -Kt*beta0dt*Me_inv;
      double mdn = -K_tor*beta0dt*Js_inv;
      double mdt = -K_ben*beta0dt*Js_inv;

      // Force Damping in x direction ---------------------------------------------------------------------------------
      // Internal Friction
      force_damp_n[0]  = fdn;
      force_damp_t[0]  = fdt;

      // Fluidic Damping
      force_damp_n[0] += -Kn*beta1[type]*sqrt(Kn_inv*(bondhistlist[n][0]*bondhistlist[n][0]*Kn_inv + Me*vn1*vn1));
      force_damp_t[0] += -Kt*beta1[type]*sqrt(Kt_inv*(bondhistlist[n][3]*bondhistlist[n][3]*Kt_inv + Me*vtr1*vtr1));

      // Quadratic Damping
      force_damp_n[0] += -beta2dtInv*(bondhistlist[n][0]*bondhistlist[n][0]*Kn_inv + Me*vn1*vn1);
      force_damp_t[0] += -beta2dtInv*(bondhistlist[n][3]*bondhistlist[n][3]*Kt_inv + Me*vtr1*vtr1);
      
      // Apply damping direction
      force_damp_n[0] *= SIGNUM_DOUBLE(vr1);
      force_damp_t[0] *= SIGNUM_DOUBLE(vr1);

      // Force Damping in y direction ---------------------------------------------------------------------------------
      // Internal Friction
      force_damp_n[1]  = fdn;
      force_damp_t[1]  = fdt;

      // Fluidic Damping
      force_damp_n[1] += -Kn*beta1[type]*sqrt(Kn_inv*(bondhistlist[n][1]*bondhistlist[n][1]*Kn_inv + Me*vn2*vn2));
      force_damp_t[1] += -Kt*beta1[type]*sqrt(Kt_inv*(bondhistlist[n][4]*bondhistlist[n][4]*Kt_inv + Me*vtr2*vtr2));
      
      // Quadratic Damping
      force_damp_n[1] += -beta2dtInv*(bondhistlist[n][1]*bondhistlist[n][1]*Kn_inv + Me*vn2*vn2);
      force_damp_t[1] += -beta2dtInv*(bondhistlist[n][4]*bondhistlist[n][4]*Kt_inv + Me*vtr2*vtr2);

      // Apply damping direction
      force_damp_n[1] *= SIGNUM_DOUBLE(vr2);
      force_damp_t[1] *= SIGNUM_DOUBLE(vr2);

      // Force Damping in z direction ---------------------------------------------------------------------------------
      // Internal Friction
      force_damp_n[2]  = fdn;
      force_damp_t[2]  = fdt;
      
      // Fluidic Damping
      force_damp_n[2] += -Kn*beta1[type]*sqrt(Kn_inv*(bondhistlist[n][2]*bondhistlist[n][2]*Kn_inv + Me*vn3*vn3));
      force_damp_t[2] += -Kt*beta1[type]*sqrt(Kt_inv*(bondhistlist[n][5]*bondhistlist[n][5]*Kt_inv + Me*vtr3*vtr3));
      
      // Quadratic Damping
      force_damp_n[2] += -beta2dtInv*(bondhistlist[n][2]*bondhistlist[n][2]*Kn_inv + Me*vn3*vn3);
      force_damp_t[2] += -beta2dtInv*(bondhistlist[n][5]*bondhistlist[n][5]*Kt_inv + Me*vtr3*vtr3);

      // Apply damping direction
      force_damp_n[2] *= SIGNUM_DOUBLE(vr3);
      force_damp_t[2] *= SIGNUM_DOUBLE(vr3);

      // Torque Damping in x direction --------------------------------------------------------------------------------
      // Internal Friction
      torque_damp_n[0]  = mdn;
      torque_damp_t[0]  = mdt;
      
      // Fluidic Damping
      torque_damp_n[0] += -K_tor*beta1[type]*sqrt(K_tor_inv*(bondhistlist[n][6]*bondhistlist[n][6]*K_tor_inv + I_n*wn1*wn1));
      torque_damp_t[0] += -K_ben*beta1[type]*sqrt(K_ben_inv*(bondhistlist[n][9]*bondhistlist[n][9]*K_ben_inv + I_t*wt1*wt1));
      
      // Quadratic Damping
      torque_damp_n[0] += -beta2dtInv*(bondhistlist[n][6]*bondhistlist[n][6]*K_tor_inv + I_n*wn1*wn1);
      torque_damp_t[0] += -beta2dtInv*(bondhistlist[n][9]*bondhistlist[n][9]*K_ben_inv + I_t*wt1*wt1);

      // Apply damping direction
      torque_damp_n[0] *= SIGNUM_DOUBLE(wr1);
      torque_damp_t[0] *= SIGNUM_DOUBLE(wr1);

      // Torque Damping in y direction --------------------------------------------------------------------------------
      // Internal Friction
      torque_damp_n[1]  = mdn;
      torque_damp_t[1]  = mdt;

      // Fluidic Damping
      torque_damp_n[1] += -K_tor*beta1[type]*sqrt(K_tor_inv*(bondhistlist[n][7]*bondhistlist[n][7]*K_tor_inv + I_n*wn2*wn2));
      torque_damp_t[1] += -K_ben*beta1[type]*sqrt(K_ben_inv*(bondhistlist[n][10]*bondhistlist[n][10]*K_ben_inv + I_t*wt2*wt2));
      
      // Quadratic Damping
      torque_damp_n[1] += -beta2dtInv*(bondhistlist[n][7]*bondhistlist[n][7]*K_tor_inv + I_n*wn2*wn2);
      torque_damp_t[1] += -beta2dtInv*(bondhistlist[n][10]*bondhistlist[n][10]*K_ben_inv + I_t*wt2*wt2);

      // Apply damping direction
      torque_damp_n[1] *= SIGNUM_DOUBLE(wr2);
      torque_damp_t[1] *= SIGNUM_DOUBLE(wr2);

      // Torque Damping in z direction --------------------------------------------------------------------------------
      // Internal Friction
      torque_damp_n[2]  = mdn;
      torque_damp_t[2]  = mdt;
      
      // Fluidic Damping
      torque_damp_n[2] += -K_tor*beta1[type]*sqrt(K_tor_inv*(bondhistlist[n][8]*bondhistlist[n][8]*K_tor_inv + I_n*wn3*wn3));
      torque_damp_t[2] += -K_ben*beta1[type]*sqrt(K_ben_inv*(bondhistlist[n][11]*bondhistlist[n][11]*K_ben_inv + I_t*wt3*wt3));
      
      // Quadratic Damping
      torque_damp_n[2] += -beta2dtInv*(bondhistlist[n][8]*bondhistlist[n][8]*K_tor_inv + I_n*wn3*wn3);
      torque_damp_t[2] += -beta2dtInv*(bondhistlist[n][11]*bondhistlist[n][11]*K_ben_inv + I_t*wt3*wt3);

      // Apply damping direction
      torque_damp_n[2] *= SIGNUM_DOUBLE(wr3);
      torque_damp_t[2] *= SIGNUM_DOUBLE(wr3);
      
    } else {
        // You shouldn't have gotten here...
        error->all(FLERR,"Damp style does not exist\n");

    }

    double tot_force_damp_x = force_damp_n[0]+force_damp_t[0];
    double tot_force_damp_y = force_damp_n[1]+force_damp_t[1];
    double tot_force_damp_z = force_damp_n[2]+force_damp_t[2];

    double tot_torqu_damp_x = torque_damp_n[0]+torque_damp_t[0];
    double tot_torqu_damp_y = torque_damp_n[1]+torque_damp_t[1];
    double tot_torqu_damp_z = torque_damp_n[2]+torque_damp_t[2];

    tor1 = - rinv * (dely*ft_bond_total[2] - delz*ft_bond_total[1]);
    tor2 = - rinv * (delz*ft_bond_total[0] - delx*ft_bond_total[2]);
    tor3 = - rinv * (delx*ft_bond_total[1] - dely*ft_bond_total[0]);


    //flag breaking of bond if criterion met
    if(breakmode == BREAKSTYLE_SIMPLE)
    {
        if(r > 2. * r_break[type])
        {
            //NP fprintf(screen,"r %f, 2. * r_break[type] %f \n",r,2. * r_break[type]);
            bondlist[n][3] = 1;
            //NP error->all(FLERR,"broken");
            error->all(FLERR,"broken");
        }
    }
    else //NP stress or stress_temp
    {
        nforce_mag  = sqrt(fn_bond[0]*fn_bond[0] + fn_bond[1]*fn_bond[1] + fn_bond[2]*fn_bond[2]);
        tforce_mag  = sqrt(bondhistlist[n][3]*bondhistlist[n][3] + bondhistlist[n][ 4]*bondhistlist[n][ 4] + bondhistlist[n][ 5]*bondhistlist[n][ 5]);
        ntorque_mag = sqrt(bondhistlist[n][6]*bondhistlist[n][6] + bondhistlist[n][ 7]*bondhistlist[n][ 7] + bondhistlist[n][ 8]*bondhistlist[n][ 8]);
        ttorque_mag = sqrt(bondhistlist[n][9]*bondhistlist[n][9] + bondhistlist[n][10]*bondhistlist[n][10] + bondhistlist[n][11]*bondhistlist[n][11]);

        nstress = sigma_break[type] < (nforce_mag/A + 2.*ttorque_mag/J*(rout-rin));
        tstress = tau_break[type]   < (tforce_mag/A +    ntorque_mag/J*(rout-rin));
        toohot = false;

        if(breakmode == BREAKSTYLE_STRESS_TEMP)
        {
            toohot = 0.5 * (Temp[i1] + Temp[i2]) > T_break[type];
            fprintf(screen,"Temp[i1] %f Temp[i2] %f, T_break[type] %f\n",Temp[i1],Temp[i2],T_break[type]);
        }

        if(nstress || tstress || toohot)
        {
            bondlist[n][3] = 1; // set back to 1...
            fprintf(screen,"broken bond %d at step %ld\n",n,update->ntimestep);
            if(toohot)fprintf(screen, "   it was too hot\n");
            if(nstress)fprintf(screen,"   it was nstress\n");
            if(tstress)fprintf(screen,"   it was tstress\n");
            
            fprintf(screen,"   sigma_break == %e\n      mag_force == %e\n",sigma_break[type],(nforce_mag/A + 2.*ttorque_mag/J*(rout-rin)));
            fprintf(screen,"     tau_break == %e\n      mag_force == %e\n",tau_break[type]  ,(tforce_mag/A +    ntorque_mag/J*(rout-rin)));
        }
    }

    // apply force to each of 2 atoms
    if (newton_bond || i1 < nlocal) {
      if (!isSymmetricUpdate) {
        f[i1][0] += (tot_force_x) + (tot_force_damp_x*SIGNUM_DOUBLE(v[i1][0]));
        f[i1][1] += (tot_force_y) + (tot_force_damp_y*SIGNUM_DOUBLE(v[i1][1]));
        f[i1][2] += (tot_force_z) + (tot_force_damp_z*SIGNUM_DOUBLE(v[i1][2]));

        torque[i1][0] += radius[i1]*tor1 + (tot_torqu_x)+(tot_torqu_damp_x*SIGNUM_DOUBLE(omega[i1][0]));
        torque[i1][1] += radius[i1]*tor2 + (tot_torqu_y)+(tot_torqu_damp_y*SIGNUM_DOUBLE(omega[i1][1]));
        torque[i1][2] += radius[i1]*tor3 + (tot_torqu_z)+(tot_torqu_damp_z*SIGNUM_DOUBLE(omega[i1][2]));
      } else {
        f[i1][0] += (tot_force_x) + (tot_force_damp_x);
        f[i1][1] += (tot_force_y) + (tot_force_damp_y);
        f[i1][2] += (tot_force_z) + (tot_force_damp_z);

        torque[i1][0] += radius[i1]*tor1 + (tot_torqu_x)+(tot_torqu_damp_x);
        torque[i1][1] += radius[i1]*tor2 + (tot_torqu_y)+(tot_torqu_damp_y);
        torque[i1][2] += radius[i1]*tor3 + (tot_torqu_z)+(tot_torqu_damp_z);
      }
    }

    if (newton_bond || i2 < nlocal) {
      if (!isSymmetricUpdate) {
        f[i2][0] += -(tot_force_x) + (tot_force_damp_x*SIGNUM_DOUBLE(v[i2][0]));
        f[i2][1] += -(tot_force_y) + (tot_force_damp_y*SIGNUM_DOUBLE(v[i2][1]));
        f[i2][2] += -(tot_force_z) + (tot_force_damp_z*SIGNUM_DOUBLE(v[i2][2]));

        torque[i2][0] += radius[i2]*tor1 - (tot_torqu_x) + (tot_torqu_damp_x*SIGNUM_DOUBLE(omega[i2][0]));
        torque[i2][1] += radius[i2]*tor2 - (tot_torqu_y) + (tot_torqu_damp_y*SIGNUM_DOUBLE(omega[i2][1]));
        torque[i2][2] += radius[i2]*tor3 - (tot_torqu_z) + (tot_torqu_damp_z*SIGNUM_DOUBLE(omega[i2][2]));
      } else {
        f[i2][0] -= (tot_force_x) + (tot_force_damp_x);
        f[i2][1] -= (tot_force_y) + (tot_force_damp_y);
        f[i2][2] -= (tot_force_z) + (tot_force_damp_z);

        torque[i2][0] += radius[i2]*tor1 - (tot_torqu_x) - (tot_torqu_damp_x);
        torque[i2][1] += radius[i2]*tor2 - (tot_torqu_y) - (tot_torqu_damp_y);
        torque[i2][2] += radius[i2]*tor3 - (tot_torqu_z) - (tot_torqu_damp_z);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void BondGran::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;
  
  // Create bond property variables
  memory->create(ro,n+1,"bond:ro");
  memory->create(ri,n+1,"bond:ri");
  memory->create(Sn,n+1,"bond:Sn");
  memory->create(St,n+1,"bond:St");
  memory->create(damp,n+1,"bond:damp");
  memory->create(beta0,n+1,"bond:beta0");
  memory->create(beta1,n+1,"bond:beta1");
  
  // Create bond break variables
  memory->create(r_break,n+1,"bond:r_break");
  memory->create(sigma_break,n+1,"bond:sigma_break");
  memory->create(tau_break,n+1,"bond:tau_break");
  memory->create(T_break,n+1,"bond:T_break");
  
  memory->create(setflag,(n+1),"bond:setflag");
  for (int i = 1; i <= n; i++)
    setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondGran::coeff(int narg, char **arg)
{
  if(narg < 8) error->all(FLERR,"Incorrect args for bond coefficients (ro, ri, sn, st, damp_style, damp_perams, break_style, break_perams)"); // Matt Schramm

  int arg_id = 1;
  int extra_args = 0;
  double ro_one = force->numeric(FLERR,arg[arg_id++]);
  double ri_one = force->numeric(FLERR,arg[arg_id++]);
  double Sn_one = force->numeric(FLERR,arg[arg_id++]);
  double St_one = force->numeric(FLERR,arg[arg_id++]);
  double damp_one = 0.0;
  double beta0_one = 0.0;
  double beta1_one = 0.0;

  isSymmetricUpdate = true;

  if (ro_one <= ri_one)
    error->all(FLERR,"ro must be greater than ri");

  // Get damping style options and set damping coefficients.
  if(force->numeric(FLERR,arg[arg_id]) == 0.0) {
    dampmode = DAMPSTYLE_NONE;
    damp_one = 0.0;
    ++arg_id;
  } else if (force->numeric(FLERR,arg[arg_id]) == 1.0) {
    dampmode = DAMPSTYLE_YU_GUO;
    damp_one = force->numeric(FLERR,arg[++arg_id]);
  } else if (force->numeric(FLERR, arg[arg_id]) == 2.0) { 
    dampmode = DAMPSTYLE_NON_LINEAR;
    isSymmetricUpdate = false;
    beta0_one = force->numeric(FLERR, arg[++arg_id]);
    beta1_one = force->numeric(FLERR, arg[++arg_id]);
    damp_one = force->numeric(FLERR, arg[++arg_id]);
    extra_args = 2;
  } else if (force->numeric(FLERR, arg[arg_id]) == 3.0) { 
    dampmode = DAMPSTYLE_NON_LINEAR_REDUCED;
    isSymmetricUpdate = false;
    damp_one = force->numeric(FLERR, arg[++arg_id]);
  } else if (force->numeric(FLERR, arg[arg_id]) == 4.0) {
    dampmode = DAMPSTYLE_NON_LINEAR_REDUCED_SIMPLIFIED;
    isSymmetricUpdate = false;
    damp_one = force->numeric(FLERR, arg[++arg_id]);
  } else if (force->numeric(FLERR, arg[arg_id]) == 5.0) {
    dampmode = DAMPSTYLE_NON_LINEAR_REL;
    beta0_one = force->numeric(FLERR, arg[++arg_id]);
    beta1_one = force->numeric(FLERR, arg[++arg_id]);
    damp_one = force->numeric(FLERR, arg[++arg_id]);
    extra_args = 2;
  } else {
    error->all(FLERR,"Damping style does not exist");
  }

  // Get break style options
  if(force->numeric(FLERR,arg[++arg_id]) == 0.0) { // arg_id == 7 for extra == 0
    breakmode = BREAKSTYLE_SIMPLE;
    if (narg != 9+extra_args)
      error->all(FLERR,"Incorrect args for bond coefficients");
  } else if(force->numeric(FLERR,arg[arg_id]) == 1.0) {
    breakmode = BREAKSTYLE_STRESS;
    if (narg != 10+extra_args)
      error->all(FLERR,"Incorrect args for bond coefficients");
  } else if(force->numeric(FLERR,arg[arg_id]) == 2.0) {
    breakmode = BREAKSTYLE_STRESS_TEMP;
    if (narg != 11+extra_args)
      error->all(FLERR,"Incorrect args for bond coefficients");
  } 
  else
    error->all(FLERR,"Incorrect args for bond coefficients");

  if (!allocated)
    allocate();

  double r_break_one = 0.0;
  double sigma_break_one = 0.0;
  double tau_break_one = 0.0;
  double T_break_one = 0.0;

  if(breakmode == BREAKSTYLE_SIMPLE)
    r_break_one = force->numeric(FLERR,arg[++arg_id]);
  else {
    sigma_break_one = force->numeric(FLERR,arg[++arg_id]);
    tau_break_one = force->numeric(FLERR,arg[++arg_id]);

    if(breakmode == BREAKSTYLE_STRESS_TEMP)
      T_break_one = force->numeric(FLERR,arg[++arg_id]);
  }

  if(comm->me == 0) {
    fprintf(screen,"\n--- Bond Parameters Being Set ---\n");
    fprintf(screen,"   ro == %g\n",ro_one);
    fprintf(screen,"   ri == %g\n",ri_one);
    fprintf(screen,"   Sn == %g\n",Sn_one);
    fprintf(screen,"   St == %g\n",St_one);
    fprintf(screen,"   damp_type == %i\n",dampmode);
    if (extra_args > 0) {
      fprintf(screen,"   beta0 = %g\n   beta1 = %g\n   beta2 = %g\n",beta0_one, beta1_one, damp_one);
    } else {
      fprintf(screen,"   damp_val == %g\n",damp_one);
    }
    fprintf(screen,"   breakmode == %i\n",breakmode);
    if(breakmode == BREAKSTYLE_SIMPLE)
      fprintf(screen,"   r_break == %g\n\n",r_break_one);
    else {
      fprintf(screen,"   sigma_break == %g\n",sigma_break_one);
      fprintf(screen,"   tau_break == %g\n",tau_break_one);
      if(breakmode == BREAKSTYLE_STRESS_TEMP)
        fprintf(screen,"   T_break == %g\n",T_break_one);
    }
    fprintf(screen,"--- End Bond Parameters ---\n\n");
  }

  int ilo, ihi;
  force->bounds(arg[0],atom->nbondtypes,ilo,ihi);
  int count = 0;
  for (int i = ilo; i <= ihi; ++i) {
    ro[i] = ro_one;
    ri[i] = ri_one;
    Sn[i] = Sn_one;
    St[i] = St_one;
    damp[i] = damp_one;
    beta0[i] = beta0_one;
    beta1[i] = beta1_one;

    if(breakmode == BREAKSTYLE_SIMPLE)
      r_break[i] = r_break_one;
    else {
      sigma_break[i] = sigma_break_one;
      tau_break[i] = tau_break_one;
      if(breakmode == BREAKSTYLE_STRESS_TEMP)
        T_break[i] = T_break_one;
    }
    setflag[i] = 1;
    ++count;
  }

  if (count == 0)
    error->all(FLERR,"Incorrect args for bond coefficients - or the bonds are not initialized in create_atoms");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondGran::equilibrium_distance(int i)
{
  //NP ATTENTION: this is _not_ correct - and has implications on fix shake, pair_lj_cut_coul_long and pppm
  //NP it is not possible to define a general equilibrium distance for this bond model
  //NP as rotational degree of freedom is present
  return 0.;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondGran::write_restart(FILE *fp)
{
  fwrite(&ro[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&ri[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&Sn[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&St[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&damp[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&beta0[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&beta1[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondGran::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&ro[1],sizeof(double),atom->nbondtypes,fp);
    fread(&ri[1],sizeof(double),atom->nbondtypes,fp);
    fread(&Sn[1],sizeof(double),atom->nbondtypes,fp);
    fread(&St[1],sizeof(double),atom->nbondtypes,fp);
    fread(&damp[1],sizeof(double),atom->nbondtypes,fp); //MS
    fread(&beta0[1],sizeof(double),atom->nbondtypes,fp); //MS
    fread(&beta1[1],sizeof(double),atom->nbondtypes,fp); //MS
  }
  MPI_Bcast(&ro[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ri[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&Sn[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&St[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&damp[1],atom->nbondtypes,MPI_DOUBLE,0,world); //MS
  MPI_Bcast(&beta0[1],atom->nbondtypes,MPI_DOUBLE,0,world); //MS
  MPI_Bcast(&beta1[1],atom->nbondtypes,MPI_DOUBLE,0,world); //MS


  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ---------------------------------------------------------------------- */

double BondGran::single(int type, double rsq, int i, int j,
                          double &fforce)
{
  error->all(FLERR,"Bond granular does not support this feature");
  /*double r = sqrt(rsq);
  double dr = r - r0[type];
  double rk = k[type] * dr;
  return rk*dr;*/
  return 0.;
}

// Use method given by Guo et al. (2013) to determine stable time step for a bonded particle
double BondGran::getMinDt()
{
  int nbondlist = neighbor->nbondlist;
  int **bondlist = neighbor->bondlist;

  double minDtSq = 1.0;
  double *radius = atom->radius;
  double *density = atom->density;
  double **bondhistlist = neighbor->bondhistlist;

  for (int k = 0; k < nbondlist; k++) {
    if (bondlist[k][3]) continue;
    int i1 = bondlist[k][0];
    int i2 = bondlist[k][1];
    int type = bondlist[k][2];
    double rin = ri[type]*fmin(radius[i1],radius[i2]); 
    double rout= ro[type]*fmin(radius[i1],radius[i2]);
    double A = M_PI * (rout*rout - rin*rin);
    double K = Sn[type]*A/fabs(bondhistlist[k][12]);
    double m1 = 4.1887902047863909846168578443*density[i1]*radius[i1]*radius[i1]*radius[i1];
    double m2 = 4.1887902047863909846168578443*density[i2]*radius[i2]*radius[i2]*radius[i2];
    double Me = m1*m2/(m1+m2);
    double curDtSq = Me/K;
    if (curDtSq < minDtSq) minDtSq = curDtSq;
  }
  return sqrt(minDtSq);
}
