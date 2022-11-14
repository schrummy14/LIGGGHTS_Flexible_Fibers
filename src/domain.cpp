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
    This file is from LAMMPS, but has been modified. Copyright for
    modification:

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz

    Copyright of original file:
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author (triclinic) : Pieter in 't Veld (SNL)
------------------------------------------------------------------------- */

#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <cmath>
#include "domain.h"
#include "style_region.h"
#include "atom.h"
#include "force.h"
#include "kspace.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "fix_deform.h"
#include "region.h"
#include "lattice.h"
#include "comm.h"
#include "universe.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "neighbor.h"

// include last to ensure correct macros
#include "domain_definitions.h"

#include "lbalance.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 1.0e-4
#define DELTAREGION 4
#define BONDSTRETCH 1.1

enum{NO_REMAP,X_REMAP,V_REMAP};                   // same as fix_deform.cpp

/* ----------------------------------------------------------------------
   default is periodic
------------------------------------------------------------------------- */

Domain::Domain(LAMMPS *lmp) :
    Pointers(lmp),
    box_exist(0),
    dimension(3),
    nonperiodic(0),
    xperiodic(1),
    yperiodic(1),
    zperiodic(1),
    triclinic(0),
    tiltsmall(1),
    xprd(0.0),
    yprd(0.0),
    zprd(0.0),
    xprd_half(0.0),
    yprd_half(0.0),
    zprd_half(0.0),
    minxlo(0.0),
    minxhi(0.0),
    minylo(0.0),
    minyhi(0.0),
    minzlo(0.0),
    minzhi(0.0),
    xy(0.0),
    xz(0.0),
    yz(0.0),
    box_change(0),
    box_change_size(0),
    box_change_shape(0),
    box_change_domain(0),
    deform_flag(0),
    deform_vremap(0),
    deform_groupbit(0),
    lattice(NULL),
    nregion(0),
    maxregion(0),
    regions(NULL),

    is_wedge(false)
{
    periodicity[0] = xperiodic;
    periodicity[1] = yperiodic;
    periodicity[2] = zperiodic;

    boundary[0][0] = boundary[0][1] = 0;
    boundary[1][0] = boundary[1][1] = 0;
    boundary[2][0] = boundary[2][1] = 0;

    prd[0] = prd[1] = prd[2] = 1.0;
    prd_half[0] = prd_half[1] = prd_half[2] = 0.5;
    prd_lamda[0] = prd_lamda[1] = prd_lamda[2] = 1.0;
    prd_half_lamda[0] = prd_half_lamda[1] = prd_half_lamda[2] = 0.5;

    boxlo[0] = boxlo[1] = boxlo[2] = -0.5;
    boxhi[0] = boxhi[1] = boxhi[2] = 0.5;
    boxlo_lamda[0] = boxlo_lamda[1] = boxlo_lamda[2] = 0.0;
    boxhi_lamda[0] = boxhi_lamda[1] = boxhi_lamda[2] = 1.0;
    boxlo_bound[0] = boxlo_bound[1] = boxlo_bound[2] = -0.5;
    boxhi_bound[0] = boxhi_bound[1] = boxhi_bound[2] = 0.5;

    corners[0][0] = corners[0][1] = corners[0][2] = 0.0;
    corners[1][0] = corners[1][1] = corners[1][2] = 0.0;
    corners[2][0] = corners[2][1] = corners[2][2] = 0.0;
    corners[3][0] = corners[3][1] = corners[3][2] = 0.0;
    corners[4][0] = corners[4][1] = corners[4][2] = 0.0;
    corners[5][0] = corners[5][1] = corners[5][2] = 0.0;
    corners[6][0] = corners[6][1] = corners[6][2] = 0.0;
    corners[7][0] = corners[7][1] = corners[7][2] = 0.0;

    sublo[0] = sublo[1] = sublo[2] = -0.5;
    subhi[0] = subhi[1] = subhi[2] = 0.5;
    sublo_lamda[0] = sublo_lamda[1] = sublo_lamda[2] = 0.0;
    subhi_lamda[0] = subhi_lamda[1] = subhi_lamda[2] = 1.0;

    h[3] = h[4] = h[5] = 0.0;
    h_inv[3] = h_inv[4] = h_inv[5] = 0.0;
    h_rate[0] = h_rate[1] = h_rate[2] =
    h_rate[3] = h_rate[4] = h_rate[5] = 0.0;
    h_ratelo[0] = h_ratelo[1] = h_ratelo[2] = 0.0;

    char **args = new char*[2];
    args[0] = (char *) "none";
    args[1] = (char *) "1.0";
    set_lattice(2,args);
    delete [] args;
}

/* ---------------------------------------------------------------------- */

Domain::~Domain()
{
  delete lattice;
  for (int i = 0; i < nregion; i++) delete regions[i];
  memory->sfree(regions);
}

/* ---------------------------------------------------------------------- */

void Domain::init()
{
  // set box_change flags if box size/shape/sub-domains ever change
  // due to shrink-wrapping or fixes that change box size/shape/sub-domains

  box_change_size = box_change_shape = box_change_domain = 0;

  if (nonperiodic == 2) box_change_size = 1;
  for (int i = 0; i < modify->nfix; i++) {
    if (modify->fix[i]->box_change_size) box_change_size = 1;
    if (modify->fix[i]->box_change_shape) box_change_shape = 1;
    if (modify->fix[i]->box_change_domain) box_change_domain = 1;
  }

  box_change = 0;
  if (box_change_size || box_change_shape || box_change_domain) box_change = 1;

  // check for fix deform

  deform_flag = deform_vremap = deform_groupbit = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"deform") == 0) {
      deform_flag = 1;
      if (((FixDeform *) modify->fix[i])->remapflag == V_REMAP) {
        deform_vremap = 1;
        deform_groupbit = modify->fix[i]->groupbit;
      }
    }

  // region inits

  for (int i = 0; i < nregion; i++) regions[i]->init();
}

/* ----------------------------------------------------------------------
   set initial global box
   assumes boxlo/hi and triclinic tilts are already set
------------------------------------------------------------------------- */

void Domain::set_initial_box()
{
  // error checks for orthogonal and triclinic domains

  if (boxlo[0] >= boxhi[0] || boxlo[1] >= boxhi[1] || boxlo[2] >= boxhi[2])
    error->one(FLERR,"Box bounds are invalid");

  if (domain->dimension == 2 && (xz != 0.0 || yz != 0.0))
    error->all(FLERR,"Cannot skew triclinic box in z for 2d simulation");

  // error check or warning on triclinic tilt factors

  if (triclinic) {
    if ((fabs(xy/(boxhi[0]-boxlo[0])) > 0.5 && xperiodic) ||
        (fabs(xz/(boxhi[0]-boxlo[0])) > 0.5 && xperiodic) ||
        (fabs(yz/(boxhi[1]-boxlo[1])) > 0.5 && yperiodic)) {
      if (tiltsmall)
      error->all(FLERR,"Triclinic box skew is too large");
      else if (comm->me == 0)
        error->warning(FLERR,"Triclinic box skew is large");
    }
  }

  // set small based on box size and SMALL
  // this works for any unit system

  small[0] = SMALL * (boxhi[0] - boxlo[0]);
  small[1] = SMALL * (boxhi[1] - boxlo[1]);
  small[2] = SMALL * (boxhi[2] - boxlo[2]);

  // adjust box lo/hi for shrink-wrapped dims

  if (boundary[0][0] == 2) boxlo[0] -= small[0];
  else if (boundary[0][0] == 3) minxlo = boxlo[0];
  if (boundary[0][1] == 2) boxhi[0] += small[0];
  else if (boundary[0][1] == 3) minxhi = boxhi[0];

  if (boundary[1][0] == 2) boxlo[1] -= small[1];
  else if (boundary[1][0] == 3) minylo = boxlo[1];
  if (boundary[1][1] == 2) boxhi[1] += small[1];
  else if (boundary[1][1] == 3) minyhi = boxhi[1];

  if (boundary[2][0] == 2) boxlo[2] -= small[2];
  else if (boundary[2][0] == 3) minzlo = boxlo[2];
  if (boundary[2][1] == 2) boxhi[2] += small[2];
  else if (boundary[2][1] == 3) minzhi = boxhi[2];
}

/* ----------------------------------------------------------------------
   set global box params
   assumes boxlo/hi and triclinic tilts are already set
------------------------------------------------------------------------- */

void Domain::set_global_box()
{

  prd[0] = xprd = boxhi[0] - boxlo[0];
  prd[1] = yprd = boxhi[1] - boxlo[1];
  prd[2] = zprd = boxhi[2] - boxlo[2];

  h[0] = xprd;
  h[1] = yprd;
  h[2] = zprd;
  h_inv[0] = 1.0/h[0];
  h_inv[1] = 1.0/h[1];
  h_inv[2] = 1.0/h[2];

  prd_half[0] = xprd_half = 0.5*xprd;
  prd_half[1] = yprd_half = 0.5*yprd;
  prd_half[2] = zprd_half = 0.5*zprd;

  if (triclinic) {
    h[3] = yz;
    h[4] = xz;
    h[5] = xy;
    h_inv[3] = -h[3] / (h[1]*h[2]);
    h_inv[4] = (h[3]*h[5] - h[1]*h[4]) / (h[0]*h[1]*h[2]);
    h_inv[5] = -h[5] / (h[0]*h[1]);

    boxlo_bound[0] = MIN(boxlo[0],boxlo[0]+xy);
    boxlo_bound[0] = MIN(boxlo_bound[0],boxlo_bound[0]+xz);
    boxlo_bound[1] = MIN(boxlo[1],boxlo[1]+yz);
    boxlo_bound[2] = boxlo[2];

    boxhi_bound[0] = MAX(boxhi[0],boxhi[0]+xy);
    boxhi_bound[0] = MAX(boxhi_bound[0],boxhi_bound[0]+xz);
    boxhi_bound[1] = MAX(boxhi[1],boxhi[1]+yz);
    boxhi_bound[2] = boxhi[2];
  }
}

/* ----------------------------------------------------------------------
   set lamda box params
   assumes global box is defined and proc assignment has been made
   uses comm->xyz_split to define subbox boundaries in consistent manner
------------------------------------------------------------------------- */

void Domain::set_lamda_box()
{
  int *myloc = comm->myloc;
  double *xsplit = comm->xsplit;
  double *ysplit = comm->ysplit;
  double *zsplit = comm->zsplit;

  sublo_lamda[0] = xsplit[myloc[0]];
  subhi_lamda[0] = xsplit[myloc[0]+1];

  sublo_lamda[1] = ysplit[myloc[1]];
  subhi_lamda[1] = ysplit[myloc[1]+1];

  sublo_lamda[2] = zsplit[myloc[2]];
  subhi_lamda[2] = zsplit[myloc[2]+1];
}

/* ----------------------------------------------------------------------
   set local subbox params for orthogonal boxes
   assumes global box is defined and proc assignment has been made
   uses comm->xyz_split to define subbox boundaries in consistent manner
   insure subhi[max] = boxhi
------------------------------------------------------------------------- */

void Domain::set_local_box()
{
  int *myloc = comm->myloc;
  int *procgrid = comm->procgrid;
  double *xsplit = comm->xsplit;
  double *ysplit = comm->ysplit;
  double *zsplit = comm->zsplit;

  if (triclinic == 0) {
    sublo[0] = boxlo[0] + xprd*xsplit[myloc[0]];
    if (myloc[0] < procgrid[0]-1)
      subhi[0] = boxlo[0] + xprd*xsplit[myloc[0]+1];
    else subhi[0] = boxhi[0];

    sublo[1] = boxlo[1] + yprd*ysplit[myloc[1]];
    if (myloc[1] < procgrid[1]-1)
      subhi[1] = boxlo[1] + yprd*ysplit[myloc[1]+1];
    else subhi[1] = boxhi[1];

    sublo[2] = boxlo[2] + zprd*zsplit[myloc[2]];
    if (myloc[2] < procgrid[2]-1)
      subhi[2] = boxlo[2] + zprd*zsplit[myloc[2]+1];
    else subhi[2] = boxhi[2];

  }
}

/* ----------------------------------------------------------------------
   reset global & local boxes due to global box boundary changes
   if shrink-wrapped, determine atom extent and reset boxlo/hi
   for triclinic, atoms must be in lamda coords (0-1) before reset_box is called
------------------------------------------------------------------------- */

void Domain::reset_box()
{
  // perform shrink-wrapping
  // compute extent of atoms on this proc
  // for triclinic, this is done in lamda space

  if (nonperiodic == 2) {
    double extent[3][2],all[3][2];

    extent[2][0] = extent[1][0] = extent[0][0] = BIG;
    extent[2][1] = extent[1][1] = extent[0][1] = -BIG;

    double **x = atom->x;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      extent[0][0] = MIN(extent[0][0],x[i][0]);
      extent[0][1] = MAX(extent[0][1],x[i][0]);
      extent[1][0] = MIN(extent[1][0],x[i][1]);
      extent[1][1] = MAX(extent[1][1],x[i][1]);
      extent[2][0] = MIN(extent[2][0],x[i][2]);
      extent[2][1] = MAX(extent[2][1],x[i][2]);
    }

    modify->box_extent(extent[0][0],extent[0][1],
                       extent[1][0],extent[1][1],
                       extent[2][0],extent[2][1]);

    // compute extent across all procs
    // flip sign of MIN to do it in one Allreduce MAX

    extent[0][0] = -extent[0][0];
    extent[1][0] = -extent[1][0];
    extent[2][0] = -extent[2][0];

    MPI_Allreduce(extent,all,6,MPI_DOUBLE,MPI_MAX,world);

    // for triclinic, convert back to box coords before changing box

    if (triclinic) lamda2x(atom->nlocal);

    // in shrink-wrapped dims, set box by atom extent
    // if minimum set, enforce min box size settings
    // for triclinic, convert lamda extent to box coords, then set box lo/hi
    // decided NOT to do the next comment - don't want to sneakily change tilt
    // for triclinic, adjust tilt factors if 2nd dim is shrink-wrapped,
    //   so that displacement in 1st dim stays the same

    if (triclinic == 0) {
      if (xperiodic == 0) {
        if (boundary[0][0] == 2) boxlo[0] = -all[0][0] - small[0];
        else if (boundary[0][0] == 3)
          boxlo[0] = MIN(-all[0][0]-small[0],minxlo);
        if (boundary[0][1] == 2) boxhi[0] = all[0][1] + small[0];
        else if (boundary[0][1] == 3) boxhi[0] = MAX(all[0][1]+small[0],minxhi);
        if (boxlo[0] > boxhi[0]) error->all(FLERR,"Illegal simulation box");
      }
      if (yperiodic == 0) {
        if (boundary[1][0] == 2) boxlo[1] = -all[1][0] - small[1];
        else if (boundary[1][0] == 3)
          boxlo[1] = MIN(-all[1][0]-small[1],minylo);
        if (boundary[1][1] == 2) boxhi[1] = all[1][1] + small[1];
        else if (boundary[1][1] == 3) boxhi[1] = MAX(all[1][1]+small[1],minyhi);
        if (boxlo[1] > boxhi[1]) error->all(FLERR,"Illegal simulation box");
      }
      if (zperiodic == 0) {
        if (boundary[2][0] == 2) boxlo[2] = -all[2][0] - small[2];
        else if (boundary[2][0] == 3)
          boxlo[2] = MIN(-all[2][0]-small[2],minzlo);
        if (boundary[2][1] == 2) boxhi[2] = all[2][1] + small[2];
        else if (boundary[2][1] == 3) boxhi[2] = MAX(all[2][1]+small[2],minzhi);
        if (boxlo[2] > boxhi[2]) error->all(FLERR,"Illegal simulation box");
      }

    } else {
      double lo[3],hi[3];
      if (xperiodic == 0) {
        lo[0] = -all[0][0]; lo[1] = 0.0; lo[2] = 0.0;
        lamda2x(lo,lo);
        hi[0] = all[0][1]; hi[1] = 0.0; hi[2] = 0.0;
        lamda2x(hi,hi);
        if (boundary[0][0] == 2) boxlo[0] = lo[0] - small[0];
        else if (boundary[0][0] == 3) boxlo[0] = MIN(lo[0]-small[0],minxlo);
        if (boundary[0][1] == 2) boxhi[0] = hi[0] + small[0];
        else if (boundary[0][1] == 3) boxhi[0] = MAX(hi[0]+small[0],minxhi);
        if (boxlo[0] > boxhi[0]) error->all(FLERR,"Illegal simulation box");
      }
      if (yperiodic == 0) {
        lo[0] = 0.0; lo[1] = -all[1][0]; lo[2] = 0.0;
        lamda2x(lo,lo);
        hi[0] = 0.0; hi[1] = all[1][1]; hi[2] = 0.0;
        lamda2x(hi,hi);
        if (boundary[1][0] == 2) boxlo[1] = lo[1] - small[1];
        else if (boundary[1][0] == 3) boxlo[1] = MIN(lo[1]-small[1],minylo);
        if (boundary[1][1] == 2) boxhi[1] = hi[1] + small[1];
        else if (boundary[1][1] == 3) boxhi[1] = MAX(hi[1]+small[1],minyhi);
        if (boxlo[1] > boxhi[1]) error->all(FLERR,"Illegal simulation box");
        //xy *= (boxhi[1]-boxlo[1]) / yprd;
      }
      if (zperiodic == 0) {
        lo[0] = 0.0; lo[1] = 0.0; lo[2] = -all[2][0];
        lamda2x(lo,lo);
        hi[0] = 0.0; hi[1] = 0.0; hi[2] = all[2][1];
        lamda2x(hi,hi);
        if (boundary[2][0] == 2) boxlo[2] = lo[2] - small[2];
        else if (boundary[2][0] == 3) boxlo[2] = MIN(lo[2]-small[2],minzlo);
        if (boundary[2][1] == 2) boxhi[2] = hi[2] + small[2];
        else if (boundary[2][1] == 3) boxhi[2] = MAX(hi[2]+small[2],minzhi);
        if (boxlo[2] > boxhi[2]) error->all(FLERR,"Illegal simulation box");
        //xz *= (boxhi[2]-boxlo[2]) / xprd;
        //yz *= (boxhi[2]-boxlo[2]) / yprd;
      }
    }
  }

  // reset box whether shrink-wrapping or not

  set_global_box();
  if (decide_loadbalance())
    lbalance->loadbalance_local_boxes();
  else
    set_local_box();

  // if shrink-wrapped & kspace is defined (i.e. using MSM) call setup()

  if (nonperiodic == 2 && force->kspace) force->kspace->setup();

  // if shrink-wrapped & triclinic, re-convert to lamda coords for new box
  // re-invoke pbc() b/c x2lamda result can be outside [0,1] due to roundoff

  if (nonperiodic == 2 && triclinic) {
    x2lamda(atom->nlocal);
    pbc();
  }
}

int Domain::decide_loadbalance() // LB
{
   if (lbalance) return 1;
   return 0;
}

/* ----------------------------------------------------------------------
   enforce PBC and modify box image flags for each atom
   called every reneighboring and by other commands that change atoms
   resulting coord must satisfy lo <= coord < hi
   MAX is important since coord - prd < lo can happen when coord = hi
   if fix deform, remap velocity of fix group atoms by box edge velocities
   for triclinic, atoms must be in lamda coords (0-1) before pbc is called
   image = 10 bits for each dimension
   increment/decrement in wrap-around fashion
------------------------------------------------------------------------- */

void Domain::pbc()
{
  int i;
  tagint idim,otherdims;
  double *lo,*hi,*period;
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  tagint *image = atom->image;

  if (triclinic == 0) {
    lo = boxlo;
    hi = boxhi;
    period = prd;
  } else {
    lo = boxlo_lamda;
    hi = boxhi_lamda;
    period = prd_lamda;
  }

  for (i = 0; i < nlocal; i++) {
    if (xperiodic) {
      if (x[i][0] < lo[0]) {
        x[i][0] += period[0];
        if (deform_vremap && mask[i] & deform_groupbit) v[i][0] += h_rate[0];
        idim = image[i] & IMGMASK;
        otherdims = image[i] ^ idim;
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | idim;
      }
      if (x[i][0] >= hi[0]) {
        x[i][0] -= period[0];
        x[i][0] = MAX(x[i][0],lo[0]);
        if (deform_vremap && mask[i] & deform_groupbit) v[i][0] -= h_rate[0];
        idim = image[i] & IMGMASK;
        otherdims = image[i] ^ idim;
        idim++;
        idim &= IMGMASK;
        image[i] = otherdims | idim;
      }
    }

    if (yperiodic) {
      if (x[i][1] < lo[1]) {
        x[i][1] += period[1];
        if (deform_vremap && mask[i] & deform_groupbit) {
          v[i][0] += h_rate[5];
          v[i][1] += h_rate[1];
        }
        idim = (image[i] >> IMGBITS) & IMGMASK;
        otherdims = image[i] ^ (idim << IMGBITS);
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMGBITS);
      }
      if (x[i][1] >= hi[1]) {
        x[i][1] -= period[1];
        x[i][1] = MAX(x[i][1],lo[1]);
        if (deform_vremap && mask[i] & deform_groupbit) {
          v[i][0] -= h_rate[5];
          v[i][1] -= h_rate[1];
        }
        idim = (image[i] >> IMGBITS) & IMGMASK;
        otherdims = image[i] ^ (idim << IMGBITS);
        idim++;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMGBITS);
      }
    }

    if (zperiodic) {
      if (x[i][2] < lo[2]) {
        x[i][2] += period[2];
        if (deform_vremap && mask[i] & deform_groupbit) {
          v[i][0] += h_rate[4];
          v[i][1] += h_rate[3];
          v[i][2] += h_rate[2];
        }
        idim = image[i] >> IMG2BITS;
        otherdims = image[i] ^ (idim << IMG2BITS);
        idim--;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMG2BITS);
      }
      if (x[i][2] >= hi[2]) {
        x[i][2] -= period[2];
        x[i][2] = MAX(x[i][2],lo[2]);
        if (deform_vremap && mask[i] & deform_groupbit) {
          v[i][0] -= h_rate[4];
          v[i][1] -= h_rate[3];
          v[i][2] -= h_rate[2];
        }
        idim = image[i] >> IMG2BITS;
        otherdims = image[i] ^ (idim << IMG2BITS);
        idim++;
        idim &= IMGMASK;
        image[i] = otherdims | (idim << IMG2BITS);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   warn if image flags of any bonded atoms are inconsistent
   could be a problem when using replicate or fix rigid
------------------------------------------------------------------------- */

void Domain::image_check()
{
  int i,j,k;

  // only need to check if system is molecular and some dimension is periodic
  // if running verlet/split, don't check on KSpace partition since
  //    it has no ghost atoms and thus bond partners won't exist

  if (!atom->molecular) return;
  if (!xperiodic && !yperiodic && (dimension == 2 || !zperiodic)) return;
  if (strcmp(update->integrate_style,"verlet/split") == 0 &&
      universe->iworld != 0) return;

  // communicate unwrapped position of owned atoms to ghost atoms

  double **unwrap;
  memory->create(unwrap,atom->nmax,3,"domain:unwrap");

  double **x = atom->x;
  tagint *image = atom->image;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    unmap(x[i],image[i],unwrap[i]);

  comm->forward_comm_array(3,unwrap);

  // compute unwrapped extent of each bond
  // flag if any bond component is longer than 1/2 of periodic box length
  // flag if any bond component is longer than non-periodic box length
  //   which means image flags in that dimension were different

  int *num_bond = atom->num_bond;
  int **bond_atom = atom->bond_atom;

  double delx,dely,delz;

  int flag = 0;
  for (i = 0; i < nlocal; i++)
    for (j = 0; j < num_bond[i]; j++) {
      k = atom->map(bond_atom[i][j]);
      if (k == -1) error->one(FLERR,"Bond atom missing in image check");

      delx = unwrap[i][0] - unwrap[k][0];
      dely = unwrap[i][1] - unwrap[k][1];
      delz = unwrap[i][2] - unwrap[k][2];

      if (xperiodic && delx > xprd_half) flag = 1;
      if (xperiodic && dely > yprd_half) flag = 1;
      if (dimension == 3 && zperiodic && delz > zprd_half) flag = 1;
      if (!xperiodic && delx > xprd) flag = 1;
      if (!yperiodic && dely > yprd) flag = 1;
      if (dimension == 3 && !zperiodic && delz > zprd) flag = 1;
    }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);
  if (flagall && comm->me == 0)
    error->warning(FLERR,"Inconsistent image flags");

  memory->destroy(unwrap);
}

/* ----------------------------------------------------------------------
   warn if end atoms in any bonded interaction
     are further apart than half a periodic box length
   could cause problems when bonded neighbor list is built since
     closest_image() could return wrong image
------------------------------------------------------------------------- */

void Domain::box_too_small_check()
{
  int i,j,k;

  // only need to check if system is molecular and some dimension is periodic
  // if running verlet/split, don't check on KSpace partition since
  //    it has no ghost atoms and thus bond partners won't exist

  if (!atom->molecular) return;
  if (!xperiodic && !yperiodic && (dimension == 2 || !zperiodic)) return;
  if (strcmp(update->integrate_style,"verlet/split") == 0 &&
      universe->iworld != 0) return;

  // maxbondall = longest current bond length
  // if periodic box dim is tiny (less than 2 * bond-length),
  //   minimum_image() itself may compute bad bond lengths
  // in this case, image_check() should warn,
  //   assuming 2 atoms have consistent image flags

  int *num_bond = atom->num_bond;
  int **bond_atom = atom->bond_atom;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double delx,dely,delz,rsq;
  double maxbondme = 0.0;

  for (i = 0; i < nlocal; i++)
    for (j = 0; j < num_bond[i]; j++) {
      k = atom->map(bond_atom[i][j]);
      if (k == -1) error->one(FLERR,"Bond atom missing in box size check");
      delx = x[i][0] - x[k][0];
      dely = x[i][1] - x[k][1];
      delz = x[i][2] - x[k][2];
      minimum_image(delx,dely,delz);
      rsq = delx*delx + dely*dely + delz*delz;
      maxbondme = MAX(maxbondme,rsq);
    }

  double maxbondall;
  MPI_Allreduce(&maxbondme,&maxbondall,1,MPI_DOUBLE,MPI_MAX,world);
  maxbondall = sqrt(maxbondall);

  // maxdelta = furthest apart 2 atoms in a bonded interaction can be
  // include BONDSTRETCH factor to account for dynamics

  double maxdelta = maxbondall * BONDSTRETCH;
  if (atom->nangles) maxdelta = 2.0 * maxbondall * BONDSTRETCH;
  if (atom->ndihedrals) maxdelta = 3.0 * maxbondall * BONDSTRETCH;

  // warn if maxdelta > than half any periodic box length
  // since atoms in the interaction could rotate into that dimension

  int flag = 0;
  if (xperiodic && maxdelta > xprd_half) flag = 1;
  if (yperiodic && maxdelta > yprd_half) flag = 1;
  if (dimension == 3 && zperiodic && maxdelta > zprd_half) flag = 1;

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);
  if (flagall && comm->me == 0)
    error->warning(FLERR,
                   "Bond/angle/dihedral extent > half of periodic box length");
}

/* ----------------------------------------------------------------------
   minimum image convention
   use 1/2 of box size as test
   for triclinic, also add/subtract tilt factors in other dims as needed
------------------------------------------------------------------------- */

void Domain::minimum_image(double &dx, double &dy, double &dz)
{
  if (triclinic == 0) {
    if (xperiodic) {
      if (fabs(dx) > xprd_half) {
        if (dx < 0.0) dx += xprd;
        else dx -= xprd;
      }
    }
    if (yperiodic) {
      if (fabs(dy) > yprd_half) {
        if (dy < 0.0) dy += yprd;
        else dy -= yprd;
      }
    }
    if (zperiodic) {
      if (fabs(dz) > zprd_half) {
        if (dz < 0.0) dz += zprd;
        else dz -= zprd;
      }
    }

  } else {
    if (zperiodic) {
      if (fabs(dz) > zprd_half) {
        if (dz < 0.0) {
          dz += zprd;
          dy += yz;
          dx += xz;
        } else {
          dz -= zprd;
          dy -= yz;
          dx -= xz;
        }
      }
    }
    if (yperiodic) {
      if (fabs(dy) > yprd_half) {
        if (dy < 0.0) {
          dy += yprd;
          dx += xy;
        } else {
          dy -= yprd;
          dx -= xy;
        }
      }
    }
    if (xperiodic) {
      if (fabs(dx) > xprd_half) {
        if (dx < 0.0) dx += xprd;
        else dx -= xprd;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   minimum image convention
   use 1/2 of box size as test
   for triclinic, also add/subtract tilt factors in other dims as needed
------------------------------------------------------------------------- */

void Domain::minimum_image(double *delta)
{
  if (triclinic == 0) {
    if (xperiodic) {
      if (fabs(delta[0]) > xprd_half) {
        if (delta[0] < 0.0) delta[0] += xprd;
        else delta[0] -= xprd;
      }
    }
    if (yperiodic) {
      if (fabs(delta[1]) > yprd_half) {
        if (delta[1] < 0.0) delta[1] += yprd;
        else delta[1] -= yprd;
      }
    }
    if (zperiodic) {
      if (fabs(delta[2]) > zprd_half) {
        if (delta[2] < 0.0) delta[2] += zprd;
        else delta[2] -= zprd;
      }
    }

  } else {
    if (zperiodic) {
      if (fabs(delta[2]) > zprd_half) {
        if (delta[2] < 0.0) {
          delta[2] += zprd;
          delta[1] += yz;
          delta[0] += xz;
        } else {
          delta[2] -= zprd;
          delta[1] -= yz;
          delta[0] -= xz;
        }
      }
    }
    if (yperiodic) {
      if (fabs(delta[1]) > yprd_half) {
        if (delta[1] < 0.0) {
          delta[1] += yprd;
          delta[0] += xy;
        } else {
          delta[1] -= yprd;
          delta[0] -= xy;
        }
      }
    }
    if (xperiodic) {
      if (fabs(delta[0]) > xprd_half) {
        if (delta[0] < 0.0) delta[0] += xprd;
        else delta[0] -= xprd;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   return local index of atom J or any of its images that is closest to atom I
   if J is not a valid index like -1, just return it
------------------------------------------------------------------------- */

int Domain::closest_image(int i, int j)
{
  if (j < 0) return j;

  int *sametag = atom->sametag;
  double **x = atom->x;
  double *xi = x[i];

  int closest = j;
  double delx = xi[0] - x[j][0];
  double dely = xi[1] - x[j][1];
  double delz = xi[2] - x[j][2];
  double rsqmin = delx*delx + dely*dely + delz*delz;
  double rsq;

  while (sametag[j] >= 0) {
    j = sametag[j];
    delx = xi[0] - x[j][0];
    dely = xi[1] - x[j][1];
    delz = xi[2] - x[j][2];
    rsq = delx*delx + dely*dely + delz*delz;
    if (rsq < rsqmin) {
      rsqmin = rsq;
      closest = j;
    }
  }
  return closest;
}

/* ----------------------------------------------------------------------
   find and return Xj image = periodic image of Xj that is closest to Xi
   for triclinic, add/subtract tilt factors in other dims as needed
------------------------------------------------------------------------- */

void Domain::closest_image(const double * const xi, const double * const xj,
                           double * const xjimage)
{
  double dx = xj[0] - xi[0];
  double dy = xj[1] - xi[1];
  double dz = xj[2] - xi[2];

  if (triclinic == 0) {
    if (xperiodic) {
      if (dx < 0.0) {
        while (dx < 0.0) dx += xprd;
        if (dx > xprd_half) dx -= xprd;
      } else {
        while (dx > 0.0) dx -= xprd;
        if (dx < -xprd_half) dx += xprd;
      }
    }
    if (yperiodic) {
      if (dy < 0.0) {
        while (dy < 0.0) dy += yprd;
        if (dy > yprd_half) dy -= yprd;
      } else {
        while (dy > 0.0) dy -= yprd;
        if (dy < -yprd_half) dy += yprd;
      }
    }
    if (zperiodic) {
      if (dz < 0.0) {
        while (dz < 0.0) dz += zprd;
        if (dz > zprd_half) dz -= zprd;
      } else {
        while (dz > 0.0) dz -= zprd;
        if (dz < -zprd_half) dz += zprd;
      }
    }

  } else {
    if (zperiodic) {
      if (dz < 0.0) {
        while (dz < 0.0) {
          dz += zprd;
          dy += yz;
          dx += xz;
        }
        if (dz > zprd_half) {
          dz -= zprd;
          dy -= yz;
          dx -= xz;
        }
      } else {
        while (dz > 0.0) {
          dz -= zprd;
          dy -= yz;
          dx -= xz;
        }
        if (dz < -zprd_half) {
          dz += zprd;
          dy += yz;
          dx += xz;
        }
      }
    }
    if (yperiodic) {
      if (dy < 0.0) {
        while (dy < 0.0) {
          dy += yprd;
          dx += xy;
        }
        if (dy > yprd_half) {
          dy -= yprd;
          dx -= xy;
        }
      } else {
        while (dy > 0.0) {
          dy -= yprd;
          dx -= xy;
        }
        if (dy < -yprd_half) {
          dy += yprd;
          dx += xy;
        }
      }
    }
    if (xperiodic) {
      if (dx < 0.0) {
        while (dx < 0.0) dx += xprd;
        if (dx > xprd_half) dx -= xprd;
      } else {
        while (dx > 0.0) dx -= xprd;
        if (dx < -xprd_half) dx += xprd;
      }
    }
  }

  xjimage[0] = xi[0] + dx;
  xjimage[1] = xi[1] + dy;
  xjimage[2] = xi[2] + dz;
}

/* ----------------------------------------------------------------------
   remap the point into the periodic box no matter how far away
   adjust 3 image flags encoded in image accordingly
   resulting coord must satisfy lo <= coord < hi
   MAX is important since coord - prd < lo can happen when coord = hi
   for triclinic, point is converted to lamda coords (0-1) before doing remap
   image = 10 bits for each dimension
   increment/decrement in wrap-around fashion
------------------------------------------------------------------------- */

void Domain::remap(double *x, tagint &image)
{
  double *lo,*hi,*period,*coord;
  double lamda[3];
  tagint idim,otherdims;

  if (triclinic == 0) {
    lo = boxlo;
    hi = boxhi;
    period = prd;
    coord = x;
  } else {
    lo = boxlo_lamda;
    hi = boxhi_lamda;
    period = prd_lamda;
    x2lamda(x,lamda);
    coord = lamda;
  }

  if (xperiodic) {
    while (coord[0] < lo[0]) {
      coord[0] += period[0];
      idim = image & IMGMASK;
      otherdims = image ^ idim;
      idim--;
      idim &= IMGMASK;
      image = otherdims | idim;
    }
    while (coord[0] >= hi[0]) {
      coord[0] -= period[0];
      idim = image & IMGMASK;
      otherdims = image ^ idim;
      idim++;
      idim &= IMGMASK;
      image = otherdims | idim;
    }
    coord[0] = MAX(coord[0],lo[0]);
  }

  if (yperiodic) {
    while (coord[1] < lo[1]) {
      coord[1] += period[1];
      idim = (image >> IMGBITS) & IMGMASK;
      otherdims = image ^ (idim << IMGBITS);
      idim--;
      idim &= IMGMASK;
      image = otherdims | (idim << IMGBITS);
    }
    while (coord[1] >= hi[1]) {
      coord[1] -= period[1];
      idim = (image >> IMGBITS) & IMGMASK;
      otherdims = image ^ (idim << IMGBITS);
      idim++;
      idim &= IMGMASK;
      image = otherdims | (idim << IMGBITS);
    }
    coord[1] = MAX(coord[1],lo[1]);
  }

  if (zperiodic) {
    while (coord[2] < lo[2]) {
      coord[2] += period[2];
      idim = image >> IMG2BITS;
      otherdims = image ^ (idim << IMG2BITS);
      idim--;
      idim &= IMGMASK;
      image = otherdims | (idim << IMG2BITS);
    }
    while (coord[2] >= hi[2]) {
      coord[2] -= period[2];
      idim = image >> IMG2BITS;
      otherdims = image ^ (idim << IMG2BITS);
      idim++;
      idim &= IMGMASK;
      image = otherdims | (idim << IMG2BITS);
    }
    coord[2] = MAX(coord[2],lo[2]);
  }

  if (triclinic) lamda2x(coord,x);
}

/* ----------------------------------------------------------------------
   remap the point into the periodic box no matter how far away
   no image flag calculation
   resulting coord must satisfy lo <= coord < hi
   MAX is important since coord - prd < lo can happen when coord = hi
   for triclinic, point is converted to lamda coords (0-1) before remap
------------------------------------------------------------------------- */

void Domain::remap(double *x)
{
  double *lo,*hi,*period,*coord;
  double lamda[3];

  if (triclinic == 0) {
    lo = boxlo;
    hi = boxhi;
    period = prd;
    coord = x;
  } else {
    lo = boxlo_lamda;
    hi = boxhi_lamda;
    period = prd_lamda;
    x2lamda(x,lamda);
    coord = lamda;
  }

  if (xperiodic) {
    while (coord[0] < lo[0]) coord[0] += period[0];
    while (coord[0] >= hi[0]) coord[0] -= period[0];
    coord[0] = MAX(coord[0],lo[0]);
  }

  if (yperiodic) {
    while (coord[1] < lo[1]) coord[1] += period[1];
    while (coord[1] >= hi[1]) coord[1] -= period[1];
    coord[1] = MAX(coord[1],lo[1]);
  }

  if (zperiodic) {
    while (coord[2] < lo[2]) coord[2] += period[2];
    while (coord[2] >= hi[2]) coord[2] -= period[2];
    coord[2] = MAX(coord[2],lo[2]);
  }

  if (triclinic) lamda2x(coord,x);
}

/* ----------------------------------------------------------------------
   remap xnew to be within half box length of xold
   do it directly, not iteratively, in case is far away
   for triclinic, both points are converted to lamda coords (0-1) before remap
------------------------------------------------------------------------- */

void Domain::remap_near(double *xnew, double *xold)
{
  int n;
  double *coordnew,*coordold,*period,*half;
  double lamdanew[3],lamdaold[3];

  if (triclinic == 0) {
    period = prd;
    half = prd_half;
    coordnew = xnew;
    coordold = xold;
  } else {
    period = prd_lamda;
    half = prd_half_lamda;
    x2lamda(xnew,lamdanew);
    coordnew = lamdanew;
    x2lamda(xold,lamdaold);
    coordold = lamdaold;
  }

  // iterative form
  // if (xperiodic) {
  //   while (coordnew[0]-coordold[0] > half[0]) coordnew[0] -= period[0];
  //   while (coordold[0]-coordnew[0] > half[0]) coordnew[0] += period[0];
  // }

  if (xperiodic) {
    if (coordnew[0]-coordold[0] > period[0]) {
      n = static_cast<int> ((coordnew[0]-coordold[0])/period[0]);
      coordnew[0] -= n*period[0];
    }
    while (coordnew[0]-coordold[0] > half[0]) coordnew[0] -= period[0];
    if (coordold[0]-coordnew[0] > period[0]) {
      n = static_cast<int> ((coordold[0]-coordnew[0])/period[0]);
      coordnew[0] += n*period[0];
    }
    while (coordold[0]-coordnew[0] > half[0]) coordnew[0] += period[0];
  }

  if (yperiodic) {
    if (coordnew[1]-coordold[1] > period[1]) {
      n = static_cast<int> ((coordnew[1]-coordold[1])/period[1]);
      coordnew[1] -= n*period[1];
    }
    while (coordnew[1]-coordold[1] > half[1]) coordnew[1] -= period[1];
    if (coordold[1]-coordnew[1] > period[1]) {
      n = static_cast<int> ((coordold[1]-coordnew[1])/period[1]);
      coordnew[1] += n*period[1];
    }
    while (coordold[1]-coordnew[1] > half[1]) coordnew[1] += period[1];
  }

  if (zperiodic) {
    if (coordnew[2]-coordold[2] > period[2]) {
      n = static_cast<int> ((coordnew[2]-coordold[2])/period[2]);
      coordnew[2] -= n*period[2];
    }
    while (coordnew[2]-coordold[2] > half[2]) coordnew[2] -= period[2];
    if (coordold[2]-coordnew[2] > period[2]) {
      n = static_cast<int> ((coordold[2]-coordnew[2])/period[2]);
      coordnew[2] += n*period[2];
    }
    while (coordold[2]-coordnew[2] > half[2]) coordnew[2] += period[2];
  }

  if (triclinic) lamda2x(coordnew,xnew);
}

/* ----------------------------------------------------------------------
   unmap the point via image flags
   x overwritten with result, don't reset image flag
   for triclinic, use h[] to add in tilt factors in other dims as needed
------------------------------------------------------------------------- */

void Domain::unmap(double *x, tagint image)
{
  int xbox = (image & IMGMASK) - IMGMAX;
  int ybox = (image >> IMGBITS & IMGMASK) - IMGMAX;
  int zbox = (image >> IMG2BITS) - IMGMAX;

  if (triclinic == 0) {
    x[0] += xbox*xprd;
    x[1] += ybox*yprd;
    x[2] += zbox*zprd;
  } else {
    x[0] += h[0]*xbox + h[5]*ybox + h[4]*zbox;
    x[1] += h[1]*ybox + h[3]*zbox;
    x[2] += h[2]*zbox;
  }
}

/* ----------------------------------------------------------------------
   unmap the point via image flags
   result returned in y, don't reset image flag
   for triclinic, use h[] to add in tilt factors in other dims as needed
------------------------------------------------------------------------- */

void Domain::unmap(double *x, tagint image, double *y)
{
  int xbox = (image & IMGMASK) - IMGMAX;
  int ybox = (image >> IMGBITS & IMGMASK) - IMGMAX;
  int zbox = (image >> IMG2BITS) - IMGMAX;

  if (triclinic == 0) {
    y[0] = x[0] + xbox*xprd;
    y[1] = x[1] + ybox*yprd;
    y[2] = x[2] + zbox*zprd;
  } else {
    y[0] = x[0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
    y[1] = x[1] + h[1]*ybox + h[3]*zbox;
    y[2] = x[2] + h[2]*zbox;
  }
}

/* ----------------------------------------------------------------------
   adjust image flags due to triclinic box flip
   flip operation is changing box vectors A,B,C to new A',B',C'
     A' = A              (A does not change)
     B' = B + mA         (B shifted by A)
     C' = C + pB + nA    (C shifted by B and/or A)
   this requires the image flags change from (a,b,c) to (a',b',c')
   so that x_unwrap for each atom is same before/after
     x_unwrap_before = xlocal + aA + bB + cC
     x_unwrap_after = xlocal + a'A' + b'B' + c'C'
   this requires:
     c' = c
     b' = b - cp
     a' = a - (b-cp)m - cn = a - b'm - cn
   in other words, for xy flip, change in x flag depends on current y flag
   this is b/c the xy flip dramatically changes which tiled image of
     simulation box an unwrapped point maps to
------------------------------------------------------------------------- */

void Domain::image_flip(int m, int n, int p)
{
  tagint *image = atom->image;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    int xbox = (image[i] & IMGMASK) - IMGMAX;
    int ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
    int zbox = (image[i] >> IMG2BITS) - IMGMAX;

    ybox -= p*zbox;
    xbox -= m*ybox + n*zbox;

    image[i] = ((tagint) (xbox + IMGMAX) & IMGMASK) |
      (((tagint) (ybox + IMGMAX) & IMGMASK) << IMGBITS) |
      (((tagint) (zbox + IMGMAX) & IMGMASK) << IMG2BITS);
  }
}

/* ----------------------------------------------------------------------
   create a lattice
------------------------------------------------------------------------- */

void Domain::set_lattice(int narg, char **arg)
{
  if (lattice) delete lattice;
  lattice = new Lattice(lmp,narg,arg);
}

/* ----------------------------------------------------------------------
   create a new region
------------------------------------------------------------------------- */

void Domain::add_region(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal region command");

  if (strcmp(arg[1],"delete") == 0) {
    delete_region(narg,arg);
    return;
  }

  if (find_region(arg[0]) >= 0) error->all(FLERR,"Reuse of region ID");

  // extend Region list if necessary

  if (nregion == maxregion) {
    maxregion += DELTAREGION;
    regions = (Region **)
      memory->srealloc(regions,maxregion*sizeof(Region *),"domain:regions");
  }

  // create the Region

  if(!lmp->wb)
  {
    if (strcmp(arg[1],"none") == 0) error->all(FLERR,"Invalid region style");
    #define REGION_CLASS
    #define RegionStyle(key,Class) \
      else if (strcmp(arg[1],#key) == 0) \
        regions[nregion] = new Class(lmp,narg,arg);
    #include "style_region.h"
    #undef REGION_CLASS
    else error->all(FLERR,"Invalid region style");
  }

  else
  {
    if (strcmp(arg[1],"none") == 0) error->all(FLERR,"Invalid region style");
    #define REGION_CLASS
    #define RegionStyle(key,Class) \
      else if (strcmp(arg[1],#key) == 0) \
        regions[nregion] = new Class(lmp,narg,arg);
    #include "region_block.h"
    #include "region_cylinder.h"
    #include "region_intersect.h"
    #include "region_sphere.h"
    #include "region_union.h"
    #undef REGION_CLASS
    else error->all(FLERR,"Invalid region style");
  }

  // initialize any region variables via init()
  // in case region is used between runs, e.g. to print a variable

  regions[nregion]->init();
  nregion++;
}

/* ----------------------------------------------------------------------
   delete a region
------------------------------------------------------------------------- */

void Domain::delete_region(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal region command");

  int iregion = find_region(arg[0]);
  if (iregion == -1) error->all(FLERR,"Delete region ID does not exist");

  delete regions[iregion];
  regions[iregion] = regions[nregion-1];
  nregion--;
}

/* ----------------------------------------------------------------------
   return region index if name matches existing region ID
   return -1 if no such region
------------------------------------------------------------------------- */

int Domain::find_region(const char *name)
{
  for (int iregion = 0; iregion < nregion; iregion++)
    if (strcmp(name,regions[iregion]->id) == 0) return iregion;
  return -1;
}

/* ----------------------------------------------------------------------
   update regions upon (sub)domain change (e.g. due to load balancing)
------------------------------------------------------------------------- */

void Domain::update_all_regions()
{
  for (int i = 0; i < nregion; i++) regions[i]->rebuild();
}

/* ----------------------------------------------------------------------
   (re)set boundary settings
   flag = 0, called from the input script
   flag = 1, called from change box command
------------------------------------------------------------------------- */

void Domain::set_boundary(int narg, char **arg, int flag)
{
  if (narg != 3) error->all(FLERR,"Illegal boundary command");

  char c;
  for (int idim = 0; idim < 3; idim++)
    for (int iside = 0; iside < 2; iside++) {
      if (iside == 0) c = arg[idim][0];
      else if (iside == 1 && strlen(arg[idim]) == 1) c = arg[idim][0];
      else c = arg[idim][1];

      const int old = boundary[idim][iside];
      if (c == 'p') boundary[idim][iside] = 0;
      else if (c == 'f') boundary[idim][iside] = 1;
      else if (c == 's') boundary[idim][iside] = 2;
      else if (c == 'm') boundary[idim][iside] = 3;
      else {
        if (flag == 0) error->all(FLERR,"Illegal boundary command");
        if (flag == 1) error->all(FLERR,"Illegal change_box command");
      }
      // if the boundary was a periodic one, but no longer is one now, set all image flags to 0 in that dimension
      // and only do this in the change_box case
      if (iside == 1 && old == 0 && boundary[idim][iside] != 0 && flag == 1)
      {
        const int nlocal = atom->nlocal;
        tagint *image = atom->image;
        for (int i = 0; i < nlocal; i++) {
            tagint thisdim = 0;
            if (idim == 0)
            {
                thisdim = image[i] & IMGMASK;
                image[i] = (image[i] ^ thisdim) | IMGMAX;
            }
            else if (idim == 1)
            {
                thisdim = (image[i] >> IMGBITS) & IMGMASK;
                image[i] = (image[i] ^ (thisdim << IMGBITS)) | (IMGMAX << IMGBITS);
            }
            else
            {
                thisdim = image[i] >> IMG2BITS;
                image[i] = (image[i] ^ (thisdim << IMG2BITS)) | (IMGMAX << IMG2BITS);
            }
        }
      }
    }

  for (int idim = 0; idim < 3; idim++)
    if ((boundary[idim][0] == 0 && boundary[idim][1]) ||
        (boundary[idim][0] && boundary[idim][1] == 0))
      error->all(FLERR,"Both sides of boundary must be periodic");

  if (boundary[0][0] == 0) xperiodic = 1;
  else xperiodic = 0;
  if (boundary[1][0] == 0) yperiodic = 1;
  else yperiodic = 0;
  if (boundary[2][0] == 0) zperiodic = 1;
  else zperiodic = 0;

  periodicity[0] = xperiodic;
  periodicity[1] = yperiodic;
  periodicity[2] = zperiodic;

  nonperiodic = 0;
  if (xperiodic == 0 || yperiodic == 0 || zperiodic == 0) {
    nonperiodic = 1;
    if (boundary[0][0] >= 2 || boundary[0][1] >= 2 ||
        boundary[1][0] >= 2 || boundary[1][1] >= 2 ||
        boundary[2][0] >= 2 || boundary[2][1] >= 2) nonperiodic = 2;
  }
}

/* ----------------------------------------------------------------------
   set domain attributes
------------------------------------------------------------------------- */

void Domain::set_box(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal box command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"tilt") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal box command");
      if (strcmp(arg[iarg+1],"small") == 0) tiltsmall = 1;
      else if (strcmp(arg[iarg+1],"large") == 0) tiltsmall = 0;
      else error->all(FLERR,"Illegal box command");
      iarg += 2;
    } else error->all(FLERR,"Illegal box command");
  }
}

/* ----------------------------------------------------------------------
   print box info, orthogonal or triclinic
------------------------------------------------------------------------- */

void Domain::print_box(const char *str)
{
  if (comm->me == 0) {
    if (screen) {
      if (triclinic == 0)
        fprintf(screen,"%sorthogonal box = (%g %g %g) to (%g %g %g)\n",
                str,boxlo[0],boxlo[1],boxlo[2],boxhi[0],boxhi[1],boxhi[2]);
      else {
        char *format = (char *)
          "%striclinic box = (%g %g %g) to (%g %g %g) with tilt (%g %g %g)\n";
        fprintf(screen,format,
                str,boxlo[0],boxlo[1],boxlo[2],boxhi[0],boxhi[1],boxhi[2],
                xy,xz,yz);
      }
    }
    if (logfile) {
      if (triclinic == 0)
        fprintf(logfile,"%sorthogonal box = (%g %g %g) to (%g %g %g)\n",
                str,boxlo[0],boxlo[1],boxlo[2],boxhi[0],boxhi[1],boxhi[2]);
      else {
        char *format = (char *)
          "%striclinic box = (%g %g %g) to (%g %g %g) with tilt (%g %g %g)\n";
        fprintf(logfile,format,
                str,boxlo[0],boxlo[1],boxlo[2],boxhi[0],boxhi[1],boxhi[2],
                xy,xz,yz);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   format boundary string for output
   assume str is 9 chars or more in length
------------------------------------------------------------------------- */

void Domain::boundary_string(char *str)
{
  int m = 0;
  for (int idim = 0; idim < 3; idim++) {
    for (int iside = 0; iside < 2; iside++) {
      if (boundary[idim][iside] == 0) str[m++] = 'p';
      else if (boundary[idim][iside] == 1) str[m++] = 'f';
      else if (boundary[idim][iside] == 2) str[m++] = 's';
      else if (boundary[idim][iside] == 3) str[m++] = 'm';
    }
    str[m++] = ' ';
  }
  str[8] = '\0';
}

/* ----------------------------------------------------------------------
   convert triclinic 0-1 lamda coords to box coords for all N atoms
   x = H lamda + x0;
------------------------------------------------------------------------- */

void Domain::lamda2x(int n)
{
  double **x = atom->x;

  for (int i = 0; i < n; i++) {
    x[i][0] = h[0]*x[i][0] + h[5]*x[i][1] + h[4]*x[i][2] + boxlo[0];
    x[i][1] = h[1]*x[i][1] + h[3]*x[i][2] + boxlo[1];
    x[i][2] = h[2]*x[i][2] + boxlo[2];
  }
}

/* ----------------------------------------------------------------------
   convert box coords to triclinic 0-1 lamda coords for all N atoms
   lamda = H^-1 (x - x0)
------------------------------------------------------------------------- */

void Domain::x2lamda(int n)
{
  double delta[3];
  double **x = atom->x;

  for (int i = 0; i < n; i++) {
    delta[0] = x[i][0] - boxlo[0];
    delta[1] = x[i][1] - boxlo[1];
    delta[2] = x[i][2] - boxlo[2];

    x[i][0] = h_inv[0]*delta[0] + h_inv[5]*delta[1] + h_inv[4]*delta[2];
    x[i][1] = h_inv[1]*delta[1] + h_inv[3]*delta[2];
    x[i][2] = h_inv[2]*delta[2];
  }
}

/* ----------------------------------------------------------------------
   convert triclinic 0-1 lamda coords to box coords for one atom
   x = H lamda + x0;
   lamda and x can point to same 3-vector
------------------------------------------------------------------------- */

void Domain::lamda2x(double *lamda, double *x)
{
  x[0] = h[0]*lamda[0] + h[5]*lamda[1] + h[4]*lamda[2] + boxlo[0];
  x[1] = h[1]*lamda[1] + h[3]*lamda[2] + boxlo[1];
  x[2] = h[2]*lamda[2] + boxlo[2];
}

/* ----------------------------------------------------------------------
   convert box coords to triclinic 0-1 lamda coords for one atom
   lamda = H^-1 (x - x0)
   x and lamda can point to same 3-vector
------------------------------------------------------------------------- */

void Domain::x2lamda(double *x, double *lamda)
{
  double delta[3];
  delta[0] = x[0] - boxlo[0];
  delta[1] = x[1] - boxlo[1];
  delta[2] = x[2] - boxlo[2];

  lamda[0] = h_inv[0]*delta[0] + h_inv[5]*delta[1] + h_inv[4]*delta[2];
  lamda[1] = h_inv[1]*delta[1] + h_inv[3]*delta[2];
  lamda[2] = h_inv[2]*delta[2];
}

/* ----------------------------------------------------------------------
   convert box coords to triclinic 0-1 lamda coords for one atom
   use my_boxlo & my_h_inv stored by caller for previous state of box
   lamda = H^-1 (x - x0)
   x and lamda can point to same 3-vector
------------------------------------------------------------------------- */

void Domain::x2lamda(double *x, double *lamda,
                     double *my_boxlo, double *my_h_inv)
{
  double delta[3];
  delta[0] = x[0] - my_boxlo[0];
  delta[1] = x[1] - my_boxlo[1];
  delta[2] = x[2] - my_boxlo[2];

  lamda[0] = my_h_inv[0]*delta[0] + my_h_inv[5]*delta[1] + my_h_inv[4]*delta[2];
  lamda[1] = my_h_inv[1]*delta[1] + my_h_inv[3]*delta[2];
  lamda[2] = my_h_inv[2]*delta[2];
}

/* ----------------------------------------------------------------------
   convert 8 lamda corner pts of lo/hi box to box coords
   return bboxlo/hi = bounding box around 8 corner pts in box coords
------------------------------------------------------------------------- */

void Domain::bbox(double *lo, double *hi, double *bboxlo, double *bboxhi)
{
  double x[3];

  bboxlo[0] = bboxlo[1] = bboxlo[2] = BIG;
  bboxhi[0] = bboxhi[1] = bboxhi[2] = -BIG;

  x[0] = lo[0]; x[1] = lo[1]; x[2] = lo[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = hi[0]; x[1] = lo[1]; x[2] = lo[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = lo[0]; x[1] = hi[1]; x[2] = lo[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = hi[0]; x[1] = hi[1]; x[2] = lo[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = lo[0]; x[1] = lo[1]; x[2] = hi[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = hi[0]; x[1] = lo[1]; x[2] = hi[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = lo[0]; x[1] = hi[1]; x[2] = hi[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = hi[0]; x[1] = hi[1]; x[2] = hi[2];
  lamda2x(x,x);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);
}

/* ----------------------------------------------------------------------
   compute 8 corner pts of triclinic box
   8 are ordered with x changing fastest, then y, finally z
   could be more efficient if just coded with xy,yz,xz explicitly
------------------------------------------------------------------------- */

void Domain::box_corners()
{
  corners[0][0] = 0.0; corners[0][1] = 0.0; corners[0][2] = 0.0;
  lamda2x(corners[0],corners[0]);
  corners[1][0] = 1.0; corners[1][1] = 0.0; corners[1][2] = 0.0;
  lamda2x(corners[1],corners[1]);
  corners[2][0] = 0.0; corners[2][1] = 1.0; corners[2][2] = 0.0;
  lamda2x(corners[2],corners[2]);
  corners[3][0] = 1.0; corners[3][1] = 1.0; corners[3][2] = 0.0;
  lamda2x(corners[3],corners[3]);
  corners[4][0] = 0.0; corners[4][1] = 0.0; corners[4][2] = 1.0;
  lamda2x(corners[4],corners[4]);
  corners[5][0] = 1.0; corners[5][1] = 0.0; corners[5][2] = 1.0;
  lamda2x(corners[5],corners[5]);
  corners[6][0] = 0.0; corners[6][1] = 1.0; corners[6][2] = 1.0;
  lamda2x(corners[6],corners[6]);
  corners[7][0] = 1.0; corners[7][1] = 1.0; corners[7][2] = 1.0;
  lamda2x(corners[7],corners[7]);
}
