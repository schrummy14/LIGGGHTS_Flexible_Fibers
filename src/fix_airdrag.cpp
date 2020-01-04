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
    This file is from LAMMPS
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_airdrag.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixAirDrag::FixAirDrag(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 8) error->all(FLERR,"Illegal fix air drag command\nair_viscosity, air_density, wx, wy, wz");

  air_viscosity = force->numeric(FLERR,arg[3]);
  air_density   = force->numeric(FLERR,arg[4]);
  wx =  force->numeric(FLERR,arg[5]);
  wy =  force->numeric(FLERR,arg[6]);
  wz =  force->numeric(FLERR,arg[7]);
  
}

/* ---------------------------------------------------------------------- */

FixAirDrag::~FixAirDrag()
{
  
}

/* ---------------------------------------------------------------------- */

int FixAirDrag::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAirDrag::init()
{
  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixAirDrag::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }

  // error checks on coarsegraining
  if(force->cg_active())
    error->cg(FLERR,this->style);
}

/* ---------------------------------------------------------------------- */

void FixAirDrag::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAirDrag::post_force(int vflag)
{
  // apply drag force to atoms in group
  // direction is opposed to velocity vector
  // magnitude depends on atom type

  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **T = atom->torque;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  const double drag_PI = 3.14159265358979323846264338328;

  double b,c,d; // Coefficients for Drag
  double temp;
  
  double vrx, vry, vrz;
  double vel_norm; // norm of the velocity

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
    
      // FORCES ----------------------------------------------------------
      // Calculate b and c coefficients for forces
      // Equations from "Classical Mechanics" by John R. Taylor
      // Linear term b = 3*pi*nu*D
      b = 6.0 * drag_PI * air_viscosity * radius[i];
      // Quadratic term c = k*p*A
      c = 0.25 * air_density * drag_PI*radius[i]*radius[i];

      vrx = wx - v[i][0];
      vry = wy - v[i][1];
      vrz = wz - v[i][2];
      vel_norm = sqrt(vrx*vrx + vry*vry + vrz*vrz);
      
      temp = b + c*vel_norm;
      f[i][0] += vrx*temp;
      f[i][1] += vry*temp;
      f[i][2] += vrz*temp;
      
      // TORQUES ---------------------------------------------------------
      // Equations from "Viscous torque on a sphere under arbitrary rotation" by U. Lei et all
      d = 8.0*drag_PI*air_viscosity*radius[i]*radius[i]*radius[i];
      T[i][0] -= d*omega[i][0];
      T[i][1] -= d*omega[i][1];
      T[i][2] -= d*omega[i][2];

    }
}

/* ---------------------------------------------------------------------- */

void FixAirDrag::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAirDrag::min_post_force(int vflag)
{
  post_force(vflag);
}
