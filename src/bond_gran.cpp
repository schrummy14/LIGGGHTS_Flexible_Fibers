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
#include "global_properties.h"

using namespace LAMMPS_NS;

/*NP
large TODO list for granular bonds:  (could be a diploma thesis?)

+ need a better dissipative formulation than the hardcoded
  'dissipate' value which produces plastic deformation
  need some vel-dependant damping
MS This has been done with a velocity based version and a force based version

+ need to carefully debug and validate this bond style
  valiation against fix rigid
MS Validation has been done against the cantilever beam and beam theory.

+ check whether run this bond style w/ or w/o gran pair style active,
  (neigh_modify command)
Needs grain to handle torque

+ need to store bond radii per particle, not per type
+ parallel debugging and testing not yet done
+ need evtally implemetation
*/

/* Matt Schramm edits for bond dampening --> MS */
/* Yu Guo edits for bond dampening --> YG */
/* D. Kramolis edits for domain end detection and other bug fixes --> KRAMOLIS */

enum
{
    BREAKSTYLE_SIMPLE = 0,
    BREAKSTYLE_STRESS,
    BREAKSTYLE_STRESS_TEMP,
    BREAKSTYLE_SOFT_STRESS,
    BREAKSTYLE_SOFT_CONTACT_STRESS,
    BREAKSTYLE_SOFT_CONTACT_STRESS_HERTZ
};
enum
{
    DAMPSTYLE_NONE = 0,
    DAMPSTYLE_LINEAR,
    DAMPSTYLE_NON_LINEAR
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
    if (!atom->style_match("bond/gran"))
        error->all(FLERR, "A granular bond style can only be used together with atom style bond/gran");
    if (comm->me == 0)
        error->warning(FLERR, "Bond granular: This is a beta version - be careful!");
    fix_Temp = NULL;
}

/* ---------------------------------------------------------------------- */

BondGran::~BondGran()
{
    if (allocated)
    {
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

void BondGran::init_style()
{
    if (breakmode == BREAKSTYLE_STRESS_TEMP)
        fix_Temp = static_cast<FixPropertyAtom *>(modify->find_fix_property("Temp", "property/atom", "scalar", 1, 0, "bond gran"));
}

/* ---------------------------------------------------------------------- */
#define SIGNUM_DOUBLE(x) (((x) > 0.0) ? 1.0 : (((x) < 0.0) ? -1.0 : 0.0))
void BondGran::compute(int eflag, int vflag)
{
    if (eflag || vflag)
        ev_setup(eflag, vflag);
    else
        evflag = 0;

    double *const radius = atom->radius;
    double *const density = atom->density;

    double *const *const x = atom->x;
    double *const *const v = atom->v;
    double *const *const omega = atom->omega;

    double **f = atom->f;
    double **torque = atom->torque;

    int **bondlist = neighbor->bondlist;
    double **bondhistlist = neighbor->bondhistlist;

    const int nbondlist = neighbor->nbondlist;
    const int nlocal = atom->nlocal;
    const int newton_bond = force->newton_bond;
    const double dt = update->dt;

    if (breakmode == BREAKSTYLE_STRESS_TEMP)
    {
        if (!fix_Temp)
            error->all(FLERR, "Internal error in BondGran");
        Temp = fix_Temp->vector_atom;
    }

    double totBondEnergy = 0.0;
    for (int n = 0; n < nbondlist; n++)
    { // Loop through Bond list
        //1st check if bond is broken,
        if (bondlist[n][3] == 1)
        {
            // printf("bond %d allready broken\n",n);
            continue;
        }

        const int i1 = bondlist[n][0]; // Sphere 1
        const int i2 = bondlist[n][1]; // Sphere 2

        // KRAMOLIS do not use neighbor cutoff - it is usually set as 4 * atom radius or more
        // (bugfix - correct sign - bonds were breaking inside domain far away from wall, added info about broken bonds)
        // should be rather bond_skin or bond_length or radius, enable periodic detection
        // 2nd check if bond overlap the box-borders
        double maxoverlap = 0.75 * radius[i1]; // max overlap is diameter - but position is in the middle of atom
        if ((x[i1][0] < (domain->boxlo[0] + maxoverlap)) && (domain->xperiodic == 0))
        {
            bondlist[n][3] = 1;
#           ifdef LIGGGHTS_BOND_DEBUG
                fprintf(screen, "broken bond %d at step %ld\n", n, update->ntimestep);
                fprintf(screen, "   bond overlaped domain border at xmin\n");
#           endif
            continue;
        }
        else if ((x[i1][0] > (domain->boxhi[0] - maxoverlap)) && (domain->xperiodic == 0))
        {
            bondlist[n][3] = 1;
#           ifdef LIGGGHTS_BOND_DEBUG
                fprintf(screen, "broken bond %d at step %ld\n", n, update->ntimestep);
                fprintf(screen, "   bond overlaped domain border at xmax\n");
#           endif
            continue;
        }
        else if ((x[i1][1] < (domain->boxlo[1] + maxoverlap)) && (domain->yperiodic == 0))
        {
            bondlist[n][3] = 1;
#           ifdef LIGGGHTS_BOND_DEBUG
                fprintf(screen, "broken bond %d at step %ld\n", n, update->ntimestep);
                fprintf(screen, "   bond overlaped domain border at ymin\n");
#           endif
            continue;
        }
        else if ((x[i1][1] > (domain->boxhi[1] - maxoverlap)) && (domain->yperiodic == 0))
        {
            bondlist[n][3] = 1;
#           ifdef LIGGGHTS_BOND_DEBUG
                fprintf(screen, "broken bond %d at step %ld\n", n, update->ntimestep);
                fprintf(screen, "   bond overlaped domain border at ymax\n");
#           endif
            continue;
        }
        else if ((x[i1][2] < (domain->boxlo[2] + maxoverlap)) && (domain->zperiodic == 0))
        {
            bondlist[n][3] = 1;
#           ifdef LIGGGHTS_BOND_DEBUG
                fprintf(screen, "broken bond %d at step %ld\n", n, update->ntimestep);
                fprintf(screen, "   bond overlaped domain border at zmin\n");
#           endif
            continue;
        }
        else if ((x[i1][2] > (domain->boxhi[2] - maxoverlap)) && (domain->zperiodic == 0))
        {
            bondlist[n][3] = 1;
#           ifdef LIGGGHTS_BOND_DEBUG
                fprintf(screen, "broken bond %d at step %ld\n", n, update->ntimestep);
                fprintf(screen, "   bond overlaped domain border at zmax\n");
#           endif
            continue;
        }
        maxoverlap = 0.75 * radius[i2]; // max overlap is diameter - but position is in the middle of atom
        if ((x[i2][0] < (domain->boxlo[0] + maxoverlap)) && (domain->xperiodic == 0))
        {
            bondlist[n][3] = 1;
#           ifdef LIGGGHTS_BOND_DEBUG
                fprintf(screen, "broken bond %d at step %ld\n", n, update->ntimestep);
                fprintf(screen, "   bond overlaped domain border at xmin\n");
#           endif
            continue;
        }
        else if ((x[i2][0] > (domain->boxhi[0] - maxoverlap)) && (domain->xperiodic == 0))
        {
            bondlist[n][3] = 1;
#           ifdef LIGGGHTS_BOND_DEBUG
                fprintf(screen, "broken bond %d at step %ld\n", n, update->ntimestep);
                fprintf(screen, "   bond overlaped domain border at xmax\n");
#           endif
            continue;
        }
        else if ((x[i2][1] < (domain->boxlo[1] + maxoverlap)) && (domain->yperiodic == 0))
        {
            bondlist[n][3] = 1;
#           ifdef LIGGGHTS_BOND_DEBUG
                fprintf(screen, "broken bond %d at step %ld\n", n, update->ntimestep);
                fprintf(screen, "   bond overlaped domain border at ymin\n");
#           endif
            continue;
        }
        else if ((x[i2][1] > (domain->boxhi[1] - maxoverlap)) && (domain->yperiodic == 0))
        {
            bondlist[n][3] = 1;
#           ifdef LIGGGHTS_BOND_DEBUG
                fprintf(screen, "broken bond %d at step %ld\n", n, update->ntimestep);
                fprintf(screen, "   bond overlaped domain border at ymax\n");
#           endif
            continue;
        }
        else if ((x[i2][2] < (domain->boxlo[2] + maxoverlap)) && (domain->zperiodic == 0))
        {
            bondlist[n][3] = 1;
#           ifdef LIGGGHTS_BOND_DEBUG
                fprintf(screen, "broken bond %d at step %ld\n", n, update->ntimestep);
                fprintf(screen, "   bond overlaped domain border at zmin\n");
#           endif
            continue;
        }
        else if ((x[i2][2] > (domain->boxhi[2] - maxoverlap)) && (domain->zperiodic == 0))
        {
            bondlist[n][3] = 1;
#           ifdef LIGGGHTS_BOND_DEBUG
                fprintf(screen, "broken bond %d at step %ld\n", n, update->ntimestep);
                fprintf(screen, "   bond overlaped domain border at zmax\n");
#           endif
            continue;
        }

        const int type = abs(bondlist[n][2]); // Get current bond type properties. These could be negative so take abs()

        const double rin = ri[type] * fmin(radius[i1], radius[i2]);
        const double rout = ro[type] * fmin(radius[i1], radius[i2]);

        const double A = M_PI * (rout * rout - rin * rin); // Bond Area
        const double J = 0.5 * M_PI * (rout * rout * rout * rout - rin * rin * rin * rin);

        const double m1 = 4.1887902047863909846168578443 * density[i1] * radius[i1] * radius[i1] * radius[i1];
        const double m2 = 4.1887902047863909846168578443 * density[i2] * radius[i2] * radius[i2] * radius[i2];
        const double Me = m1 * m2 / (m1 + m2);

        const double Ip = 0.5 * M_PI * (rout * rout * rout * rout - rin * rin * rin * rin); // MS
        const double I = 0.5 * Ip;

        const double J1 = 0.4 * m1 * radius[i1] * radius[i1];
        const double J2 = 0.4 * m2 * radius[i2] * radius[i2];
        const double Js = J1 * J2 / (J1 + J2);

        double delx = x[i1][0] - x[i2][0];       // x-directional seperation
        double dely = x[i1][1] - x[i2][1];       // y-directional seperation
        double delz = x[i1][2] - x[i2][2];       // z-directional seperation
        domain->minimum_image(delx, dely, delz); // Not 100% sure what this is... Allows for periodic simulations

        const double rsq = delx * delx + dely * dely + delz * delz;
        const double rsqinv = 1. / rsq;
        const double r = sqrt(rsq);
        const double rinv = 1. / r;

        // relative translational velocity
        const double vr1 = v[i1][0] - v[i2][0]; // relative velocity between sphere1 and sphere2 in the x-direction
        const double vr2 = v[i1][1] - v[i2][1]; // relative velocity between sphere1 and sphere2 in the y-direction
        const double vr3 = v[i1][2] - v[i2][2]; // relative velocity between sphere1 and sphere2 in the z-direction

        // normal component
        double vnnr = vr1 * delx + vr2 * dely + vr3 * delz;
        vnnr *= rsqinv;
        const double vn1 = delx * vnnr;
        const double vn2 = dely * vnnr;
        const double vn3 = delz * vnnr;

        // tangential component at the center of the sphere
        const double vt1 = vr1 - vn1;
        const double vt2 = vr2 - vn2;
        const double vt3 = vr3 - vn3;

        // Check if bond just formed and set eq distance
        if (bondhistlist[n][12] == 0.0)
        {
            // const double tmp1 = Sn[type] * A /r;
            // const double tmp2 = St[type] * A /r;
            // const double tmp3 = St[type] * Ip /r;
            // const double tmp4 = Sn[type] * I /r;
            // if (n==0) fprintf(screen, "Kfn = %f, Kft = %f, Kmn = %f, Kmt = %f\n", tmp1, tmp2, tmp3, tmp4);
            if (r < 1.0e-10)
            {
                fprintf(screen, "BondLength = %g\n", r);
                error->all(FLERR, "bondlength too small\n");
            }
            bondhistlist[n][12] = r;
#           ifdef LIGGGHTS_BOND_DEBUG
                fprintf(screen, "INFO: Setting bond length between %i and %i at %g\n", atom->tag[i1], atom->tag[i2], bondhistlist[n][12]);
#           endif
        }

        // Set Bond Length
        const double bondLength = fabs(bondhistlist[n][12]);
        const double bondLengthInv = 1.0 / bondLength; // Try r1r2/(r1+r2)  
        // const double bondLengthInv = (radius[i1]+radius[i2])/(radius[i1]*radius[i2]);

        // Check if the bond IS broken but the atoms need to seperate before contact physics start to occur
        if (breakmode >= BREAKSTYLE_SOFT_STRESS && bondlist[n][2] < 1)
        { // Bondlist[n][3] does not work as durring re-neighboring, the value is set to 0
            // Here we could do a more complicated method to "push" the particles away from each other...
            const double sep = r - (radius[i1] + radius[i2]);
            if (sep > 0)
            {
#               ifdef LIGGGHTS_BOND_DEBUG
                    fprintf(screen, "Spheres %i and %i have seperated enough, set bond to broken and apply contact forces.\n", atom->tag[i1], atom->tag[i2]);
#               endif
                bondlist[n][3] = 1;
                if (i1 < nlocal)
                {
                    for (int k1 = 0; k1 < atom->num_bond[i1]; k1++)
                    {
                        int j2 = atom->map(atom->bond_atom[i1][k1]); //mapped index of bond-partner
                        if (i2 == j2)
                        {
                            atom->bond_type[i1][k1] = 0;
                            break;
                        }
                    }
                }
                if (i2 < nlocal)
                {
                    for (int k1 = 0; k1 < atom->num_bond[i2]; k1++)
                    {
                        int j1 = atom->map(atom->bond_atom[i2][k1]); //mapped index of bond-partner
                        if (i1 == j1)
                        {
                            atom->bond_type[i2][k1] = 0;
                            break;
                        }
                    }
                }
            }
            else
            {
                double fn;
                switch (breakmode)
                {
                case BREAKSTYLE_SOFT_STRESS:
                {
                    fn = -0.001 * Sn[type] * A * bondLengthInv * sep * rinv;
                    break;
                }
                case BREAKSTYLE_SOFT_CONTACT_STRESS:
                {
                    fn = MAX(-Sn[type] * A * bondLengthInv * (r - bondLength) * rinv, 0.0);
                    break;
                }
                case BREAKSTYLE_SOFT_CONTACT_STRESS_HERTZ:
                {
                    const double poi = 0.25;
                    const double one_minus_p2_inv = 1. / (1.0 - poi * poi);
                    const double deltan = r - bondLength;
                    const double reff = 0.5 * rout;
                    const double Yeff = 0.5 * Sn[type] * one_minus_p2_inv; // Need Poisson's ratio here also... set to 0.25...
                    const double sqrtval = sqrt(-reff * deltan);
                    const double kn = 4. / 3. * Yeff * sqrtval;
                    fn = MAX(-kn * deltan * rinv, 0.0);
                    break;
                }
                default:
                {
                    fprintf(screen, "I should not be here\n");
                    fn = 0.0;
                }
                }

                if (newton_bond || i1 < nlocal)
                {
                    f[i1][0] += fn * delx;
                    f[i1][1] += fn * dely;
                    f[i1][2] += fn * delz;
                }
                if (newton_bond || i2 < nlocal)
                {
                    f[i2][0] -= fn * delx;
                    f[i2][1] -= fn * dely;
                    f[i2][2] -= fn * delz;
                }
            }
            continue;
        }

        // Set Stiffness Values
        // This is a new value to try and fix the stiffness value when using a smooth mega-particle
        const double newVal1 = 2.0*radius[i1]*radius[i2]/(radius[i1]+radius[i2]);
        const double newVal2 = bondLength/(radius[i1]+radius[i2]);
        const double newVal = newVal1*newVal2*newVal2*newVal2;
        const double newValInv = 1.0/newVal;
        // const double newVal = newVal1*pow(newVal2,3); // a = (r1+r2)/Lb) ~~ 3  // a = gamma*F(r1,r2) ~~ 3 //  /*a(bondLength, radius[i1], radius[i2])*/


        // const double newVal = radius[i1]*radius[i2]*bondLength/((radius[i1]+radius[i2])*(radius[i1]+radius[i2]));
        const double Kn = Sn[type] * A * newValInv; // * bondLengthInv;
        const double Kt = St[type] * A * newValInv; // * bondLengthInv;
        const double K_tor = St[type] * Ip * newValInv; // * bondLengthInv;
        const double K_ben = Sn[type] * I * newValInv; // * bondLengthInv;

        // relative rotational velocity for shear
        double wr1 = (radius[i1] * omega[i1][0] + radius[i2] * omega[i2][0]) * rinv;
        double wr2 = (radius[i1] * omega[i1][1] + radius[i2] * omega[i2][1]) * rinv;
        double wr3 = (radius[i1] * omega[i1][2] + radius[i2] * omega[i2][2]) * rinv;

        // relative velocities for shear at contact
        const double vtr1 = vt1 - (delz * wr2 - dely * wr3);
        const double vtr2 = vt2 - (delx * wr3 - delz * wr1);
        const double vtr3 = vt3 - (dely * wr1 - delx * wr2);

        // relative angular velocity
        wr1 = omega[i1][0] - omega[i2][0];
        wr2 = omega[i1][1] - omega[i2][1];
        wr3 = omega[i1][2] - omega[i2][2];

        // normal component
        double wnnr = wr1 * delx + wr2 * dely + wr3 * delz;
        wnnr *= rsqinv;
        const double wn1 = delx * wnnr;
        const double wn2 = dely * wnnr;
        const double wn3 = delz * wnnr;

        // tangential component
        const double wt1 = wr1 - wn1;
        const double wt2 = wr2 - wn2;
        const double wt3 = wr3 - wn3;

        // calc change in bond normal forces (Can be directly calculated does not need history!!!)
        const double sndt = Kn * (r - bondLength) * rinv;
        const double fnx = -sndt * delx;
        const double fny = -sndt * dely;
        const double fnz = -sndt * delz;

        // calc change in shear forces
        const double stdt = Kt * dt;
        const double dftx = -vtr1 * stdt;
        const double dfty = -vtr2 * stdt;
        const double dftz = -vtr3 * stdt;

        // calc change in normal torque
        const double K_tor_dt = K_tor * dt;
        const double dtnx = -wn1 * K_tor_dt;
        const double dtny = -wn2 * K_tor_dt;
        const double dtnz = -wn3 * K_tor_dt;

        // calc change in tang torque
        const double K_ben_dt = K_ben * dt; // K_ben will become an input parameter
        const double dttx = -wt1 * K_ben_dt;
        const double dtty = -wt2 * K_ben_dt;
        const double dttz = -wt3 * K_ben_dt;

        // rotate forces from previous time step

        //rotate tangential force
        double rot = bondhistlist[n][3] * delx + bondhistlist[n][4] * dely + bondhistlist[n][5] * delz;
        rot *= rsqinv;
        double vel_temp_x = bondhistlist[n][3] - rot * delx;
        double vel_temp_y = bondhistlist[n][4] - rot * dely;
        double vel_temp_z = bondhistlist[n][5] - rot * delz;
        double vel_norm = sqrt(vel_temp_x * vel_temp_x + vel_temp_y * vel_temp_y + vel_temp_z * vel_temp_z);
        double f_norm = bondhistlist[n][3] * bondhistlist[n][3] + bondhistlist[n][4] * bondhistlist[n][4] + bondhistlist[n][5] * bondhistlist[n][5];
        if (vel_norm == 0)
            f_norm = 0;
        else
            f_norm = sqrt(f_norm) / vel_norm;

        bondhistlist[n][3] = f_norm * vel_temp_x;
        bondhistlist[n][4] = f_norm * vel_temp_y;
        bondhistlist[n][5] = f_norm * vel_temp_z;

        //rotate normal torque
        rot = bondhistlist[n][6] * delx + bondhistlist[n][7] * dely + bondhistlist[n][8] * delz;
        rot *= rsqinv;
        vel_temp_x = rot * delx;
        vel_temp_y = rot * dely;
        vel_temp_z = rot * delz;
        vel_norm = sqrt(vel_temp_x * vel_temp_x + vel_temp_y * vel_temp_y + vel_temp_z * vel_temp_z);
        f_norm = bondhistlist[n][6] * bondhistlist[n][6] + bondhistlist[n][7] * bondhistlist[n][7] + bondhistlist[n][8] * bondhistlist[n][8];
        if (vel_norm == 0)
            f_norm = 0;
        else
            f_norm = sqrt(f_norm) / vel_norm;

        bondhistlist[n][6] = f_norm * vel_temp_x;
        bondhistlist[n][7] = f_norm * vel_temp_y;
        bondhistlist[n][8] = f_norm * vel_temp_z;

        // rotate tangential torque
        rot = bondhistlist[n][9] * delx + bondhistlist[n][10] * dely + bondhistlist[n][11] * delz;
        rot *= rsqinv;
        vel_temp_x = bondhistlist[n][9] - rot * delx;
        vel_temp_y = bondhistlist[n][10] - rot * dely;
        vel_temp_z = bondhistlist[n][11] - rot * delz;
        vel_norm = sqrt(vel_temp_x * vel_temp_x + vel_temp_y * vel_temp_y + vel_temp_z * vel_temp_z);
        f_norm = bondhistlist[n][9] * bondhistlist[n][9] + bondhistlist[n][10] * bondhistlist[n][10] + bondhistlist[n][11] * bondhistlist[n][11];
        if (vel_norm == 0)
            f_norm = 0;
        else
            f_norm = sqrt(f_norm) / vel_norm;

        bondhistlist[n][9] = f_norm * vel_temp_x;
        bondhistlist[n][10] = f_norm * vel_temp_y;
        bondhistlist[n][11] = f_norm * vel_temp_z;

        // increment normal and tangential force and torque
        bondhistlist[n][0] = fnx;
        bondhistlist[n][1] = fny;
        bondhistlist[n][2] = fnz;
        bondhistlist[n][3] += dftx;
        bondhistlist[n][4] += dfty;
        bondhistlist[n][5] += dftz;
        bondhistlist[n][6] += dtnx;
        bondhistlist[n][7] += dtny;
        bondhistlist[n][8] += dtnz;
        bondhistlist[n][9] += dttx;
        bondhistlist[n][10] += dtty;
        bondhistlist[n][11] += dttz;

        const double tot_force_x = bondhistlist[n][0] + bondhistlist[n][3];
        const double tot_force_y = bondhistlist[n][1] + bondhistlist[n][4];
        const double tot_force_z = bondhistlist[n][2] + bondhistlist[n][5];
        const double tot_torqu_x = bondhistlist[n][6] + bondhistlist[n][9];
        const double tot_torqu_y = bondhistlist[n][7] + bondhistlist[n][10];
        const double tot_torqu_z = bondhistlist[n][8] + bondhistlist[n][11];

        // torque due to tangential bond force
        const double ft_bond_total_x = bondhistlist[n][3];
        const double ft_bond_total_y = bondhistlist[n][4];
        const double ft_bond_total_z = bondhistlist[n][5];

        // Get Total Bond Energy
        const double Kn_inv = 1.0 / Kn;
        const double Kt_inv = 1.0 / Kt;
        const double K_tor_inv = 1.0 / K_tor;
        const double K_ben_inv = 1.0 / K_ben;

        const double Efnx = 0.5 * (bondhistlist[n][0] * bondhistlist[n][0] * Kn_inv);
        const double Efny = 0.5 * (bondhistlist[n][1] * bondhistlist[n][1] * Kn_inv);
        const double Efnz = 0.5 * (bondhistlist[n][2] * bondhistlist[n][2] * Kn_inv);

        const double Eftx = 0.5 * (bondhistlist[n][3] * bondhistlist[n][3] * Kt_inv);
        const double Efty = 0.5 * (bondhistlist[n][4] * bondhistlist[n][4] * Kt_inv);
        const double Eftz = 0.5 * (bondhistlist[n][5] * bondhistlist[n][5] * Kt_inv);

        const double Emnx = 0.5 * (bondhistlist[n][6] * bondhistlist[n][6] * K_tor_inv);
        const double Emny = 0.5 * (bondhistlist[n][7] * bondhistlist[n][7] * K_tor_inv);
        const double Emnz = 0.5 * (bondhistlist[n][8] * bondhistlist[n][8] * K_tor_inv);

        const double Emtx = 0.5 * (bondhistlist[n][9] * bondhistlist[n][9] * K_ben_inv);
        const double Emty = 0.5 * (bondhistlist[n][10] * bondhistlist[n][10] * K_ben_inv);
        const double Emtz = 0.5 * (bondhistlist[n][11] * bondhistlist[n][11] * K_ben_inv);

        totBondEnergy += Efnx + Efny + Efnz + Eftx + Efty + Eftz + Emnx + Emny + Emnz + Emtx + Emty + Emtz;

        // Apply Damping
        double force_damp_n[3] = {0.0};
        double force_damp_t[3] = {0.0};
        double torque_damp_n[3] = {0.0};
        double torque_damp_t[3] = {0.0};
        if (dampmode == DAMPSTYLE_NONE)
        {
            // Don't need to do anything. Arrays were already set to zero
            // force_damp_n[0] = force_damp_n[1] = force_damp_n[2] = 0.0;
            // force_damp_t[0] = force_damp_t[1] = force_damp_t[2] = 0.0;
            // torque_damp_n[0] = torque_damp_n[1] = torque_damp_n[2] = 0.0;
            // torque_damp_t[0] = torque_damp_t[1] = torque_damp_t[2] = 0.0;
        }
        else if (dampmode == DAMPSTYLE_LINEAR)
        { // of the form fd = 2*b*sqrt(k*m)*v

            const double damp2 = 2.0 * damp[type];
            // normal force dampening
            const double d_fn_sqrt_2_Me_Sn = -damp2 * sqrt(Me * Kn);
            force_damp_n[0] = d_fn_sqrt_2_Me_Sn * vn1;
            force_damp_n[1] = d_fn_sqrt_2_Me_Sn * vn2;
            force_damp_n[2] = d_fn_sqrt_2_Me_Sn * vn3;

            // tangential force dampening
            const double d_ft_sqrt_2_Me_St = -damp2 * sqrt(Me * Kt);
            force_damp_t[0] = d_ft_sqrt_2_Me_St * vtr1;
            force_damp_t[1] = d_ft_sqrt_2_Me_St * vtr2;
            force_damp_t[2] = d_ft_sqrt_2_Me_St * vtr3;

            // normal moment dampening
            const double d_mn_sqrt_2_Js_Ktor = -damp2 * sqrt(Js * K_tor);
            torque_damp_n[0] = d_mn_sqrt_2_Js_Ktor * wn1;
            torque_damp_n[1] = d_mn_sqrt_2_Js_Ktor * wn2;
            torque_damp_n[2] = d_mn_sqrt_2_Js_Ktor * wn3;

            // tangential moment dampening
            const double d_mt_sqrt_2_Js_Kben = -damp2 * sqrt(Js * K_ben);
            torque_damp_t[0] = d_mt_sqrt_2_Js_Kben * wt1;
            torque_damp_t[1] = d_mt_sqrt_2_Js_Kben * wt2;
            torque_damp_t[2] = d_mt_sqrt_2_Js_Kben * wt3;
        }
        else if (dampmode == DAMPSTYLE_NON_LINEAR)
        {
            // Get energy in bond F(gamma) = c*(2.0*E/k)^gamma
            // double Kn_inv = 1.0/Kn;
            // double Kt_inv = 1.0/Kt;
            // double K_tor_inv = 1.0/K_tor;
            // double K_ben_inv = 1.0/K_ben;

            // Energies were calculated earlier...
            // 0.5*k*x^2 -> 0.5*k*k*x^2/k -> 0.5*f^2/k
            // Factor out the 0.5 as the energy is multiplied by 2

            force_damp_n[0] = -Kn * (beta0[type] + beta1[type] * sqrt(2.0 * Efnx * Kn_inv) + 2.0 * damp[type] * Efnx * Kn_inv) * SIGNUM_DOUBLE(vr1);
            force_damp_n[1] = -Kn * (beta0[type] + beta1[type] * sqrt(2.0 * Efny * Kn_inv) + 2.0 * damp[type] * Efny * Kn_inv) * SIGNUM_DOUBLE(vr2);
            force_damp_n[2] = -Kn * (beta0[type] + beta1[type] * sqrt(2.0 * Efnz * Kn_inv) + 2.0 * damp[type] * Efnz * Kn_inv) * SIGNUM_DOUBLE(vr3);

            force_damp_t[0] = -Kt * (beta0[type] + beta1[type] * sqrt(2.0 * Eftx * Kt_inv) + 2.0 * damp[type] * Eftx * Kt_inv) * SIGNUM_DOUBLE(vr1);
            force_damp_t[1] = -Kt * (beta0[type] + beta1[type] * sqrt(2.0 * Efty * Kt_inv) + 2.0 * damp[type] * Efty * Kt_inv) * SIGNUM_DOUBLE(vr2);
            force_damp_t[2] = -Kt * (beta0[type] + beta1[type] * sqrt(2.0 * Eftz * Kt_inv) + 2.0 * damp[type] * Eftz * Kt_inv) * SIGNUM_DOUBLE(vr3);

            torque_damp_n[0] = -K_tor * (beta0[type] + beta1[type] * sqrt(2.0 * Emnx * K_tor_inv) + 2.0 * damp[type] * Emnx * K_tor_inv) * SIGNUM_DOUBLE(wr1);
            torque_damp_n[1] = -K_tor * (beta0[type] + beta1[type] * sqrt(2.0 * Emny * K_tor_inv) + 2.0 * damp[type] * Emny * K_tor_inv) * SIGNUM_DOUBLE(wr2);
            torque_damp_n[2] = -K_tor * (beta0[type] + beta1[type] * sqrt(2.0 * Emnz * K_tor_inv) + 2.0 * damp[type] * Emnz * K_tor_inv) * SIGNUM_DOUBLE(wr3);

            torque_damp_t[0] = -K_ben * (beta0[type] + beta1[type] * sqrt(2.0 * Emtx * K_ben_inv) + 2.0 * damp[type] * Emtx * K_ben_inv) * SIGNUM_DOUBLE(wr1);
            torque_damp_t[1] = -K_ben * (beta0[type] + beta1[type] * sqrt(2.0 * Emty * K_ben_inv) + 2.0 * damp[type] * Emty * K_ben_inv) * SIGNUM_DOUBLE(wr2);
            torque_damp_t[2] = -K_ben * (beta0[type] + beta1[type] * sqrt(2.0 * Emtz * K_ben_inv) + 2.0 * damp[type] * Emtz * K_ben_inv) * SIGNUM_DOUBLE(wr3);
        }
        else
        {
            // You shouldn't have gotten here...
            error->all(FLERR, "Damp style does not exist\n");
        }

        const double tot_force_damp_x = force_damp_n[0] + force_damp_t[0];
        const double tot_force_damp_y = force_damp_n[1] + force_damp_t[1];
        const double tot_force_damp_z = force_damp_n[2] + force_damp_t[2];

        const double tot_torqu_damp_x = torque_damp_n[0] + torque_damp_t[0];
        const double tot_torqu_damp_y = torque_damp_n[1] + torque_damp_t[1];
        const double tot_torqu_damp_z = torque_damp_n[2] + torque_damp_t[2];

        const double tor1 = -rinv * (dely * ft_bond_total_z - delz * ft_bond_total_y);
        const double tor2 = -rinv * (delz * ft_bond_total_x - delx * ft_bond_total_z);
        const double tor3 = -rinv * (delx * ft_bond_total_y - dely * ft_bond_total_x);

        //flag breaking of bond if criterion met
        if (breakmode == BREAKSTYLE_SIMPLE)
        {
            if (r > 2. * r_break[type])
            {
                //NP fprintf(screen,"r %f, 2. * r_break[type] %f \n",r,2. * r_break[type]);
                bondlist[n][3] = 1;
                //NP error->all(FLERR,"broken");
                error->all(FLERR, "broken");
            }
        }
        else
        { //NP stress or stress_temp
            const double nforce_mag = sqrt(bondhistlist[n][0] * bondhistlist[n][0] + bondhistlist[n][1] * bondhistlist[n][1] + bondhistlist[n][2] * bondhistlist[n][2]);
            const double tforce_mag = sqrt(bondhistlist[n][3] * bondhistlist[n][3] + bondhistlist[n][4] * bondhistlist[n][4] + bondhistlist[n][5] * bondhistlist[n][5]);
            const double ntorque_mag = sqrt(bondhistlist[n][6] * bondhistlist[n][6] + bondhistlist[n][7] * bondhistlist[n][7] + bondhistlist[n][8] * bondhistlist[n][8]);
            const double ttorque_mag = sqrt(bondhistlist[n][9] * bondhistlist[n][9] + bondhistlist[n][10] * bondhistlist[n][10] + bondhistlist[n][11] * bondhistlist[n][11]);

            const bool nstress = sigma_break[type] < (nforce_mag / A + 2. * ttorque_mag / J * (rout - rin));
            const bool tstress = tau_break[type] < (tforce_mag / A + ntorque_mag / J * (rout - rin));
            const bool toohot = breakmode == BREAKSTYLE_STRESS_TEMP ? 0.5 * (Temp[i1] + Temp[i2]) > T_break[type] : false;

            if (nstress || tstress || toohot)
            {
                bondlist[n][3] = 1; // set back to 1...
#               ifdef LIGGGHTS_BOND_DEBUG
                    fprintf(screen, "\nBroken bond between spheres %i and %i at step %ld\n", atom->tag[i1], atom->tag[i2], update->ntimestep);
                    if (toohot)
                        fprintf(screen, "   it was too hot\n");
                    if (nstress)
                        fprintf(screen, "   it was nstress\n");
                    if (tstress)
                        fprintf(screen, "   it was tstress\n");

                    fprintf(screen, "   sigma_break == %e\n      mag_force == %e\n", sigma_break[type], (nforce_mag / A + 2. * ttorque_mag / J * (rout - rin)));
                    fprintf(screen, "     tau_break == %e\n      mag_force == %e\n", tau_break[type], (tforce_mag / A + ntorque_mag / J * (rout - rin)));
                    if (breakmode == BREAKSTYLE_STRESS_TEMP)
                    {
                        fprintf(screen, "Temp[i1] %f Temp[i2] %f, T_break[type] %f\n", Temp[i1], Temp[i2], T_break[type]);
                    }
#               endif
            }
        }

        // When a bond breaks, check that the two spheres are not overlapping each other, if they are
        // we will need to let the two spheres seperate some more before contact equations are used.
        int newBondType = 0;
        if (breakmode >= BREAKSTYLE_SOFT_STRESS && bondlist[n][3] == 1)
        {
            if (r < radius[i1] + radius[i2])
            {
#               ifdef LIGGGHTS_BOND_DEBUG
                    fprintf(screen, "The spheres %i and %i are overlapping, set type to negative.\n", atom->tag[i1], atom->tag[i2]);
#               endif
                bondlist[n][3] = 0;
                bondlist[n][2] *= -1;
                newBondType = -1;

                bondhistlist[n][0] = 0.0;
                bondhistlist[n][1] = 0.0;
                bondhistlist[n][2] = 0.0;
                bondhistlist[n][3] = 0.0;
                bondhistlist[n][4] = 0.0;
                bondhistlist[n][5] = 0.0;
                bondhistlist[n][6] = 0.0;
                bondhistlist[n][7] = 0.0;
                bondhistlist[n][8] = 0.0;
                bondhistlist[n][9] = 0.0;
                bondhistlist[n][10] = 0.0;
                bondhistlist[n][11] = 0.0;
                if (breakmode == BREAKSTYLE_SOFT_CONTACT_STRESS || breakmode == BREAKSTYLE_SOFT_CONTACT_STRESS_HERTZ)
                {
#ifdef LIGGGHTS_BOND_DEBUG
                    fprintf(screen, "Spheres %i and %i bond length is now %e\n", atom->tag[i1], atom->tag[i2], r);
#endif
                    bondhistlist[n][12] = r;
                }
            }
            if (i1 < nlocal)
            {
                for (int k1 = 0; k1 < atom->num_bond[i1]; k1++)
                {
                    int j2 = atom->map(atom->bond_atom[i1][k1]); //mapped index of bond-partner
                    if (i2 == j2)
                    {
                        atom->bond_type[i1][k1] *= newBondType;
                        break;
                    }
                }
            }
            if (i2 < nlocal)
            {
                for (int k1 = 0; k1 < atom->num_bond[i2]; k1++)
                {
                    int j1 = atom->map(atom->bond_atom[i2][k1]); //mapped index of bond-partner
                    if (i1 == j1)
                    {
                        atom->bond_type[i2][k1] *= newBondType;
                        break;
                    }
                }
            }
        }

        // apply force to each of 2 atoms
        if (newton_bond || i1 < nlocal)
        {
            f[i1][0] += (tot_force_x) + (tot_force_damp_x);
            f[i1][1] += (tot_force_y) + (tot_force_damp_y);
            f[i1][2] += (tot_force_z) + (tot_force_damp_z);

            torque[i1][0] += radius[i1] * tor1 + (tot_torqu_x) + (tot_torqu_damp_x);
            torque[i1][1] += radius[i1] * tor2 + (tot_torqu_y) + (tot_torqu_damp_y);
            torque[i1][2] += radius[i1] * tor3 + (tot_torqu_z) + (tot_torqu_damp_z);
        }

        if (newton_bond || i2 < nlocal)
        {
            f[i2][0] -= (tot_force_x) + (tot_force_damp_x);
            f[i2][1] -= (tot_force_y) + (tot_force_damp_y);
            f[i2][2] -= (tot_force_z) + (tot_force_damp_z);

            torque[i2][0] += radius[i2] * tor1 - (tot_torqu_x) - (tot_torqu_damp_x);
            torque[i2][1] += radius[i2] * tor2 - (tot_torqu_y) - (tot_torqu_damp_y);
            torque[i2][2] += radius[i2] * tor3 - (tot_torqu_z) - (tot_torqu_damp_z);
        }
    }
    // neigh2atom(0.5*totBondEnergy);
    atom->setBondEnergy(0.5 * totBondEnergy);
}

/* ---------------------------------------------------------------------- */

void BondGran::allocate()
{
    allocated = 1;
    const int n = atom->nbondtypes;
    const int np1 = n + 1;

    // Create bond property variables
    memory->create(ro, np1, "bond:ro");
    memory->create(ri, np1, "bond:ri");
    memory->create(Sn, np1, "bond:Sn");
    memory->create(St, np1, "bond:St");
    memory->create(damp, np1, "bond:damp");
    memory->create(beta0, np1, "bond:beta0");
    memory->create(beta1, np1, "bond:beta1");

    // Create bond break variables
    memory->create(r_break, np1, "bond:r_break");
    memory->create(sigma_break, np1, "bond:sigma_break");
    memory->create(tau_break, np1, "bond:tau_break");
    memory->create(T_break, np1, "bond:T_break");

    memory->create(setflag, np1, "bond:setflag");
    for (int i = 1; i <= n; i++)
        setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondGran::coeff(int narg, char **arg)
{
    if (narg < 8)
        error->all(FLERR, "Incorrect args for bond coefficients (ro, ri, sn, st, damp_style, damp_perams, break_style, break_perams)"); // Matt Schramm

    int arg_id = 1;
    int extra_args = 0;
    const double ro_one = force->numeric(FLERR, arg[arg_id++]);
    const double ri_one = force->numeric(FLERR, arg[arg_id++]);
    const double Sn_one = force->numeric(FLERR, arg[arg_id++]);
    const double St_one = force->numeric(FLERR, arg[arg_id++]);
    double damp_one = 0.0;
    double beta0_one = 0.0;
    double beta1_one = 0.0;

    if (ro_one <= ri_one)
        error->all(FLERR, "ro must be greater than ri");

    // Get damping style options and set damping coefficients.
    if (force->numeric(FLERR, arg[arg_id]) == 0.0)
    {
        dampmode = DAMPSTYLE_NONE;
        damp_one = 0.0;
        ++arg_id;
    }
    else if (force->numeric(FLERR, arg[arg_id]) == 1.0)
    {
        dampmode = DAMPSTYLE_LINEAR;
        damp_one = force->numeric(FLERR, arg[++arg_id]);
    }
    else if (force->numeric(FLERR, arg[arg_id]) == 2.0)
    {
        dampmode = DAMPSTYLE_NON_LINEAR;
        beta0_one = force->numeric(FLERR, arg[++arg_id]);
        beta1_one = force->numeric(FLERR, arg[++arg_id]);
        damp_one = force->numeric(FLERR, arg[++arg_id]);
        extra_args = 2;
    }
    else
    {
        error->all(FLERR, "Damping style does not exist");
    }

    // Get break style options
    if (force->numeric(FLERR, arg[++arg_id]) == 0.0)
    { // arg_id == 7 for extra == 0
        breakmode = BREAKSTYLE_SIMPLE;
        if (narg != 9 + extra_args)
            error->all(FLERR, "Incorrect args for bond coefficients");
    }
    else if (force->numeric(FLERR, arg[arg_id]) == 1.0)
    {
        breakmode = BREAKSTYLE_STRESS;
        if (narg != 10 + extra_args)
            error->all(FLERR, "Incorrect args for bond coefficients");
    }
    else if (force->numeric(FLERR, arg[arg_id]) == 2.0)
    {
        breakmode = BREAKSTYLE_STRESS_TEMP;
        if (narg != 11 + extra_args)
            error->all(FLERR, "Incorrect args for bond coefficients");
    }
    else if (force->numeric(FLERR, arg[arg_id]) == 3.0)
    {
        breakmode = BREAKSTYLE_SOFT_STRESS;
        if (narg != 10 + extra_args)
            error->all(FLERR, "Incorrect args for bond coefficients");
    }
    else if (force->numeric(FLERR, arg[arg_id]) == 4.0)
    {
        breakmode = BREAKSTYLE_SOFT_CONTACT_STRESS;
        if (narg != 10 + extra_args)
            error->all(FLERR, "Incorrect args for bond coefficients");
    }
    else if (force->numeric(FLERR, arg[arg_id]) == 5.0)
    {
        breakmode = BREAKSTYLE_SOFT_CONTACT_STRESS_HERTZ;
        if (narg != 10 + extra_args)
            error->all(FLERR, "Incorrect args for bond coefficients");
    }
    else
        error->all(FLERR, "Incorrect args for bond coefficients");

    if (!allocated)
        allocate();

    double r_break_one = 0.0;
    double sigma_break_one = 0.0;
    double tau_break_one = 0.0;
    double T_break_one = 0.0;

    if (breakmode == BREAKSTYLE_SIMPLE)
        r_break_one = force->numeric(FLERR, arg[++arg_id]);
    else
    {
        sigma_break_one = force->numeric(FLERR, arg[++arg_id]);
        tau_break_one = force->numeric(FLERR, arg[++arg_id]);

        if (breakmode == BREAKSTYLE_STRESS_TEMP)
            T_break_one = force->numeric(FLERR, arg[++arg_id]);
    }

    if (comm->me == 0 && screen)
    {
        fprintf(screen, "\n--- Bond Parameters Being Set ---\n");
        fprintf(screen, "   ro == %g\n", ro_one);
        fprintf(screen, "   ri == %g\n", ri_one);
        fprintf(screen, "   Sn == %g\n", Sn_one);
        fprintf(screen, "   St == %g\n", St_one);
        fprintf(screen, "   damp_type == %i\n", dampmode);
        if (extra_args > 0)
        {
            fprintf(screen, "   beta0 = %g\n   beta1 = %g\n   beta2 = %g\n", beta0_one, beta1_one, damp_one);
        }
        else
        {
            fprintf(screen, "   damp_val == %g\n", damp_one);
        }
        fprintf(screen, "   breakmode == %i\n", breakmode);
        if (breakmode == BREAKSTYLE_SIMPLE)
            fprintf(screen, "   r_break == %g\n\n", r_break_one);
        else
        {
            fprintf(screen, "   sigma_break == %g\n", sigma_break_one);
            fprintf(screen, "   tau_break == %g\n", tau_break_one);
            if (breakmode == BREAKSTYLE_STRESS_TEMP)
                fprintf(screen, "   T_break == %g\n", T_break_one);
        }
        fprintf(screen, "--- End Bond Parameters ---\n\n");
    }

    int ilo, ihi;
    force->bounds(arg[0], atom->nbondtypes, ilo, ihi);
    int count = 0;
    for (int i = ilo; i <= ihi; ++i)
    {
        ro[i] = ro_one;
        ri[i] = ri_one;
        Sn[i] = Sn_one;
        St[i] = St_one;
        damp[i] = damp_one;
        beta0[i] = beta0_one;
        beta1[i] = beta1_one;

        if (breakmode == BREAKSTYLE_SIMPLE)
            r_break[i] = r_break_one;
        else
        {
            sigma_break[i] = sigma_break_one;
            tau_break[i] = tau_break_one;
            if (breakmode == BREAKSTYLE_STRESS_TEMP)
                T_break[i] = T_break_one;
        }
        setflag[i] = 1;
        ++count;
    }

    if (count == 0)
        error->all(FLERR, "Incorrect args for bond coefficients - or the bonds are not initialized in create_atoms");
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
    fwrite(&ro[1], sizeof(double), atom->nbondtypes, fp);
    fwrite(&ri[1], sizeof(double), atom->nbondtypes, fp);
    fwrite(&Sn[1], sizeof(double), atom->nbondtypes, fp);
    fwrite(&St[1], sizeof(double), atom->nbondtypes, fp);
    fwrite(&damp[1], sizeof(double), atom->nbondtypes, fp);
    fwrite(&beta0[1], sizeof(double), atom->nbondtypes, fp);
    fwrite(&beta1[1], sizeof(double), atom->nbondtypes, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondGran::read_restart(FILE *fp)
{
    allocate();

    if (comm->me == 0)
    {
        fread(&ro[1], sizeof(double), atom->nbondtypes, fp);
        fread(&ri[1], sizeof(double), atom->nbondtypes, fp);
        fread(&Sn[1], sizeof(double), atom->nbondtypes, fp);
        fread(&St[1], sizeof(double), atom->nbondtypes, fp);
        fread(&damp[1], sizeof(double), atom->nbondtypes, fp);  //MS
        fread(&beta0[1], sizeof(double), atom->nbondtypes, fp); //MS
        fread(&beta1[1], sizeof(double), atom->nbondtypes, fp); //MS
    }
    MPI_Bcast(&ro[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
    MPI_Bcast(&ri[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
    MPI_Bcast(&Sn[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
    MPI_Bcast(&St[1], atom->nbondtypes, MPI_DOUBLE, 0, world);
    MPI_Bcast(&damp[1], atom->nbondtypes, MPI_DOUBLE, 0, world);  //MS
    MPI_Bcast(&beta0[1], atom->nbondtypes, MPI_DOUBLE, 0, world); //MS
    MPI_Bcast(&beta1[1], atom->nbondtypes, MPI_DOUBLE, 0, world); //MS

    for (int i = 1; i <= atom->nbondtypes; i++)
        setflag[i] = 1;
}

/* ---------------------------------------------------------------------- */

double BondGran::single(int type, double rsq, int i, int j, double &fforce)
{
    error->all(FLERR, "Bond granular does not support this feature");
    return 0.;
}

// Use method given by Guo et al. (2013) to determine stable time step for a bonded linearly damped particle
double BondGran::getMinDt()
{
    const int nbondlist = neighbor->nbondlist;
    int *const *const bondlist = neighbor->bondlist;

    double curDt = 0.0;
    double minDt = 1.0;
    double *const radius = atom->radius;
    double *const density = atom->density;
    double *const *const bondhistlist = neighbor->bondhistlist;

    for (int k = 0; k < nbondlist; k++)
    {
        if (bondlist[k][3] > 0)
            continue;
        const int i1 = bondlist[k][0];
        const int i2 = bondlist[k][1];
        const int type = bondlist[k][2];

        const double rin = ri[type] * fmin(radius[i1], radius[i2]);
        const double rout = ro[type] * fmin(radius[i1], radius[i2]);
        const double A = M_PI * (rout * rout - rin * rin);

        // Set Stiffness Values
        // This is a new value to try and fix the stiffness value when using a smooth mega-particle
        const double newVal1 = 2.0*radius[i1]*radius[i2]/(radius[i1]+radius[i2]);
        const double newVal2 = fabs(bondhistlist[k][12])/(radius[i1]+radius[i2]);
        const double newVal = newVal1*pow(newVal2,3.0);
        const double K = Sn[type] * A / newVal;

        // const double K = Sn[type] * A / fabs(bondhistlist[k][12]);
        const double m1 = 4.1887902047863909846168578443 * density[i1] * radius[i1] * radius[i1] * radius[i1];
        const double m2 = 4.1887902047863909846168578443 * density[i2] * radius[i2] * radius[i2] * radius[i2];
        const double Me = m1 * m2 / (m1 + m2);

        // curDt = sqrt(Me/K);
        if (dampmode == DAMPSTYLE_LINEAR)
        {
            curDt = sqrt(Me / K) / (1.0 + 2.93 * damp[type]);
        }
        else
        {
            // A fiber bond should contain two atoms with the same poisson's ratio, take the average for other systems
            // Assume a value of 0.25 for now...
            const double poi = 0.25;
            if (damp[type] == 0.0 && beta0[type] == 0.0)
            {
                curDt = sqrt(Me / K) * sqrt(2.0) * sqrt(1.0 + poi) / (1.6 + 0.76 * beta1[type]);
            }
            else
            {
                curDt = sqrt(Me / K) * sqrt(2.0) * sqrt(1.0 + poi) / (2.29 - 5.49 * beta1[type] + 0.0023 * damp[type]);
            }
        }
        if (curDt < minDt)
            minDt = curDt;
    }
    return minDt;
}
