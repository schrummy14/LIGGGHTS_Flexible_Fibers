/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fix_nve_sphere_limit.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "modify.h"
#include "comm.h"
#include "error.h"
#include "domain.h" 

using namespace LAMMPS_NS;
using namespace FixConst;

#define INERTIA 0.4          // moment of inertia prefactor for sphere
enum{NONE,DIPOLE};

/* ---------------------------------------------------------------------- */

FixNVESphereLimit::FixNVESphereLimit(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal fix nve/limit/sphere command"); 

  time_integrate = 1;
  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;

  if(strcmp(arg[3],"radius_ratio") == 0)
      relflag = 1;
  else if(strcmp(arg[3],"absolute") == 0)
      relflag = 0;
  else
      error->fix_error(FLERR,this,"expecting keyword 'absolute' or 'radius_ratio'"); 

  xlimit = atof(arg[4]);

  ncount = 0;
}

/* ---------------------------------------------------------------------- */

int FixNVESphereLimit::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVESphereLimit::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  vlimitsq = (xlimit/dtv) * (xlimit/dtv);
  ncount = 0;

  if(relflag == 1 && (!atom->radius_flag || !atom->rmass_flag))
        error->fix_error(FLERR,this,"using 'radius_ratio' needs per-atom radius and mass");

  if (strstr(update->integrate_style,"respa"))
    step_respa = ((Respa *) update->integrate)->step;
  // warn if using fix shake, which will lead to invalid constraint forces

  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"shake") == 0) {
      if (comm->me == 0)
        error->warning(FLERR,"Should not use fix nve/limit with fix shake");
    }
    
  // check that all particles are finite-size spheres
  // no point particles allowed

  double *radius = atom->radius;
  int *mask_ = atom->mask;
  int nlocal_ = atom->nlocal;

  for (int i = 0; i < nlocal_; i++)
    if (mask_[i] & groupbit)
      if (radius[i] == 0.0)
        error->one(FLERR,"Fix nve/sphere requires extended particles");
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixNVESphereLimit::initial_integrate(int vflag)
{
  double dtfm,vsq,scale, rsq; 

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  double *radius = atom->radius; 
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double dtfrotate, dtirotate; 
  if (domain->dimension == 2) dtfrotate = dtf / 0.5; // for discs the formula is I=0.5*Mass*Radius^2
  else dtfrotate  = dtf / INERTIA;
  if (relflag == 1) 
  {
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];

          vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
          rsq = radius[i]*radius[i];
          if (vsq > vlimitsq*rsq) {
            ncount++;
            scale = sqrt(vlimitsq*rsq/vsq);
            v[i][0] *= scale;
            v[i][1] *= scale;
            v[i][2] *= scale;
          }
          
          x[i][0] += dtv * v[i][0];
          x[i][1] += dtv * v[i][1];
          x[i][2] += dtv * v[i][2];
          
          
          // rotation update
          dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
          omega[i][0] += dtirotate * torque[i][0];
          omega[i][1] += dtirotate * torque[i][1];
          omega[i][2] += dtirotate * torque[i][2];
          
          vsq = omega[i][0]*omega[i][0];
          vsq+= omega[i][1]*omega[i][1];
          vsq+= omega[i][2]*omega[i][2];
          vsq*= radius[i]*radius[i];
          
          if(vsq>vlimitsq*rsq){
            ncount++;
            scale = sqrt(vlimitsq*rsq/vsq);
            omega[i][0] *= scale/radius[i];
            omega[i][1] *= scale/radius[i];
            omega[i][2] *= scale/radius[i];
          }
        }
      }
  }
  else
  {
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];

          vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
          if (vsq > vlimitsq) {
            ncount++;
            scale = sqrt(vlimitsq/vsq);
            v[i][0] *= scale;
            v[i][1] *= scale;
            v[i][2] *= scale;
          }

          x[i][0] += dtv * v[i][0];
          x[i][1] += dtv * v[i][1];
          x[i][2] += dtv * v[i][2];
          
          // rotation update
          dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
          omega[i][0] += dtirotate * torque[i][0];
          omega[i][1] += dtirotate * torque[i][1];
          omega[i][2] += dtirotate * torque[i][2];
          
          vsq = omega[i][0]*omega[i][0];
          vsq+= omega[i][1]*omega[i][1];
          vsq+= omega[i][2]*omega[i][2];
          vsq*= radius[i]*radius[i];
          
          if(vsq>vlimitsq*rsq){
            ncount++;
            scale = sqrt(vlimitsq/vsq);
            omega[i][0] *= scale/radius[i];
            omega[i][1] *= scale/radius[i];
            omega[i][2] *= scale/radius[i];
          }
        }
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVESphereLimit::final_integrate()
{
  double dtfm,vsq,scale,rsq; 

  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    if (relflag == 1) 
    {
        for (int i = 0; i < nlocal; i++) {
          if (mask[i] & groupbit) {
            dtfm = dtf / rmass[i];
            v[i][0] += dtfm * f[i][0];
            v[i][1] += dtfm * f[i][1];
            v[i][2] += dtfm * f[i][2];

            vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
            rsq = radius[i]*radius[i];
            if (vsq > vlimitsq*rsq) {
              ncount++;
              scale = sqrt(vlimitsq*rsq/vsq);
              v[i][0] *= scale;
              v[i][1] *= scale;
              v[i][2] *= scale;
            }
          }
        }
    }
    else
    {
        for (int i = 0; i < nlocal; i++) {
          if (mask[i] & groupbit) {
            dtfm = dtf / rmass[i];
            v[i][0] += dtfm * f[i][0];
            v[i][1] += dtfm * f[i][1];
            v[i][2] += dtfm * f[i][2];

            vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
            if (vsq > vlimitsq) {
              ncount++;
              scale = sqrt(vlimitsq/vsq);
              v[i][0] *= scale;
              v[i][1] *= scale;
              v[i][2] *= scale;
            }
          }
        }
    }

  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];

        vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
        if (vsq > vlimitsq) {
          ncount++;
          scale = sqrt(vlimitsq/vsq);
          v[i][0] *= scale;
          v[i][1] *= scale;
          v[i][2] *= scale;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVESphereLimit::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVESphereLimit::final_integrate_respa(int ilevel, int iloop)
{
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVESphereLimit::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  vlimitsq = (xlimit/dtv) * (xlimit/dtv);
}

/* ----------------------------------------------------------------------
   energy of indenter interaction
------------------------------------------------------------------------- */

double FixNVESphereLimit::compute_scalar()
{
  double one = ncount;
  double all;
  MPI_Allreduce(&one,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}
