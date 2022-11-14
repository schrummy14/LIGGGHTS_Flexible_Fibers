/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-     DCS Computing GmbH, Linz

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "lmptype.h"
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "balance.h"
#include "atom.h"
#include "comm.h"
#include "irregular.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "error.h"
#include "pair.h" //NP modified C.K.
#include "neighbor.h" //NP modified C.K.
#include "vector_liggghts.h" //NP modified C.K.
#include "modify.h" //NP modified C.K.

using namespace LAMMPS_NS;

enum{NONE,UNIFORM,USER,DYNAMIC};
enum{X,Y,Z};

//#define BALANCE_DEBUG 1

/* ---------------------------------------------------------------------- */

Balance::Balance(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  memory->create(proccount,nprocs,"balance:proccount");
  memory->create(allproccount,nprocs,"balance:allproccount");

  user_xsplit = user_ysplit = user_zsplit = NULL;
  dflag = 0;

  fp = NULL;
  firststep = 1;
}

/* ---------------------------------------------------------------------- */

Balance::~Balance()
{
  memory->destroy(proccount);
  memory->destroy(allproccount);

  delete [] user_xsplit;
  delete [] user_ysplit;
  delete [] user_zsplit;

  if (dflag) {
    delete [] bdim;
    delete [] count;
    delete [] sum;
    delete [] target;
    delete [] onecount;
    delete [] lo;
    delete [] hi;
    delete [] losum;
    delete [] hisum;
  }

  if (fp) fclose(fp);
}

/* ----------------------------------------------------------------------
   called as balance command in input script
------------------------------------------------------------------------- */

void Balance::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Balance command before simulation box is defined");

  //NP modified C.K.
  //NP because would need irregular in any case
  if(disallow_irregular())
    error->all(FLERR,"Balance command cannot be used with meshes - please use fix balance instead");

  if (comm->me == 0 && screen) fprintf(screen,"Balancing ...\n");

  // parse arguments

  if (narg < 1) error->all(FLERR,"Illegal balance command");

  int dimension = domain->dimension;
  int *procgrid = comm->procgrid;
  xflag = yflag = zflag = NONE;
  dflag = 0;
  outflag = 0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"x") == 0) {
      if (xflag == DYNAMIC) error->all(FLERR,"Illegal balance command");
      if (strcmp(arg[iarg+1],"uniform") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal balance command");
        xflag = UNIFORM;
        iarg += 2;
      } else {
        if (1 + procgrid[0]-1 > narg)
          error->all(FLERR,"Illegal balance command");
        xflag = USER;
        delete [] user_xsplit;
        user_xsplit = new double[procgrid[0]+1];
        user_xsplit[0] = 0.0;
        iarg++;
        for (int i = 1; i < procgrid[0]; i++)
          user_xsplit[i] = force->numeric(FLERR,arg[iarg++]);
        user_xsplit[procgrid[0]] = 1.0;
      }
    } else if (strcmp(arg[iarg],"y") == 0) {
      if (yflag == DYNAMIC) error->all(FLERR,"Illegal balance command");
      if (strcmp(arg[iarg+1],"uniform") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal balance command");
        yflag = UNIFORM;
        iarg += 2;
      } else {
        if (1 + procgrid[1]-1 > narg)
          error->all(FLERR,"Illegal balance command");
        yflag = USER;
        delete [] user_ysplit;
        user_ysplit = new double[procgrid[1]+1];
        user_ysplit[0] = 0.0;
        iarg++;
        for (int i = 1; i < procgrid[1]; i++)
          user_ysplit[i] = force->numeric(FLERR,arg[iarg++]);
        user_ysplit[procgrid[1]] = 1.0;
      }
    } else if (strcmp(arg[iarg],"z") == 0) {
      if (zflag == DYNAMIC) error->all(FLERR,"Illegal balance command");
      if (strcmp(arg[iarg+1],"uniform") == 0) {
        if (iarg+2 > narg) error->all(FLERR,"Illegal balance command");
        zflag = UNIFORM;
        iarg += 2;
      } else {
        if (1 + procgrid[2]-1 > narg)
          error->all(FLERR,"Illegal balance command");
        zflag = USER;
        delete [] user_zsplit;
        user_zsplit = new double[procgrid[2]+1];
        user_zsplit[0] = 0.0;
        iarg++;
        for (int i = 1; i < procgrid[2]; i++)
          user_zsplit[i] = force->numeric(FLERR,arg[iarg++]);
        user_zsplit[procgrid[2]] = 1.0;
      }

    } else if (strcmp(arg[iarg],"dynamic") == 0) {
      if (xflag != NONE || yflag != NONE || zflag != NONE)
        error->all(FLERR,"Illegal balance command");
      if (iarg+4 > narg) error->all(FLERR,"Illegal balance command");
      dflag = 1;
      xflag = yflag = DYNAMIC;
      if (dimension == 3) zflag = DYNAMIC;
      if (strlen(arg[iarg+1]) > 3) error->all(FLERR,"Illegal balance command");
      strcpy(bstr,arg[iarg+1]);
      nitermax = force->inumeric(FLERR,arg[iarg+2]);
      if (nitermax <= 0) error->all(FLERR,"Illegal balance command");
      thresh = force->numeric(FLERR,arg[iarg+3]);
      if (thresh < 1.0) error->all(FLERR,"Illegal balance command");
      iarg += 4;

    } else if (strcmp(arg[iarg],"out") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal balance command");
      if (outflag) error->all(FLERR,"Illegal balance command");
      outflag = 1;
      if (me == 0) {
        fp = fopen(arg[iarg+1],"w");
        if (fp == NULL) error->one(FLERR,"Cannot open balance output file");
      }
      iarg += 2;
    } else error->all(FLERR,"Illegal balance command");
  }

  // error check

  if (zflag && dimension == 2)
    error->all(FLERR,"Cannot balance in z dimension for 2d simulation");

  if (xflag == USER)
    for (int i = 1; i <= procgrid[0]; i++)
      if (user_xsplit[i-1] >= user_xsplit[i])
        error->all(FLERR,"Illegal balance command");
  if (yflag == USER)
    for (int i = 1; i <= procgrid[1]; i++)
      if (user_ysplit[i-1] >= user_ysplit[i])
        error->all(FLERR,"Illegal balance command");
  if (zflag == USER)
    for (int i = 1; i <= procgrid[2]; i++)
      if (user_zsplit[i-1] >= user_zsplit[i])
        error->all(FLERR,"Illegal balance command");

  if (dflag) {
    for (size_t i = 0; i < strlen(bstr); i++) { //NP modified R.B.
      if (bstr[i] != 'x' && bstr[i] != 'y' && bstr[i] != 'z')
        error->all(FLERR,"Balance dynamic string is invalid");
      if (bstr[i] == 'z' && dimension == 2)
        error->all(FLERR,"Balance dynamic string is invalid");
      for (size_t j = i+1; j < strlen(bstr); j++) //NP modified R.B.
        if (bstr[i] == bstr[j])
          error->all(FLERR,"Balance dynamic string is invalid");
    }
  }

  // insure atoms are in current box & update box via shrink-wrap
  // no exchange() since doesn't matter if atoms are assigned to correct procs

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // imbinit = initial imbalance
  // use current splits instead of nlocal since atoms may not be in sub-box

  domain->x2lamda(atom->nlocal);
  int maxinit;
  double imbinit = imbalance_splits(maxinit);
  domain->lamda2x(atom->nlocal);

  // debug output of initial state

#ifdef BALANCE_DEBUG
  if (me == 0 && fp) dumpout(update->ntimestep,fp);
#endif

  int niter = 0;

  // explicit setting of sub-domain sizes

  if (xflag == UNIFORM) {
    for (int i = 0; i < procgrid[0]; i++)
      comm->xsplit[i] = i * 1.0/procgrid[0];
    comm->xsplit[procgrid[0]] = 1.0;
  }

  if (yflag == UNIFORM) {
    for (int i = 0; i < procgrid[1]; i++)
      comm->ysplit[i] = i * 1.0/procgrid[1];
    comm->ysplit[procgrid[1]] = 1.0;
  }

  if (zflag == UNIFORM) {
    for (int i = 0; i < procgrid[2]; i++)
      comm->zsplit[i] = i * 1.0/procgrid[2];
    comm->zsplit[procgrid[2]] = 1.0;
  }

  if (xflag == USER)
    for (int i = 0; i <= procgrid[0]; i++) comm->xsplit[i] = user_xsplit[i];

  if (yflag == USER)
    for (int i = 0; i <= procgrid[1]; i++) comm->ysplit[i] = user_ysplit[i];

  if (zflag == USER)
    for (int i = 0; i <= procgrid[2]; i++) comm->zsplit[i] = user_zsplit[i];

  // static load-balance of sub-domain sizes

  if (dflag) {
    static_setup(bstr);
    niter = dynamic();
  }

  // output of final result

  if (outflag && me == 0) dumpout(update->ntimestep,fp);

  // reset comm->uniform flag if necessary

  if (comm->uniform) {
    if (xflag == USER || xflag == DYNAMIC) comm->uniform = 0;
    if (yflag == USER || yflag == DYNAMIC) comm->uniform = 0;
    if (zflag == USER || zflag == DYNAMIC) comm->uniform = 0;
  } else {
    if (dimension == 3) {
      if (xflag == UNIFORM && yflag == UNIFORM && zflag == UNIFORM)
        comm->uniform = 1;
    } else {
      if (xflag == UNIFORM && yflag == UNIFORM) comm->uniform = 1;
    }
  }

  // reset proc sub-domains

  if (domain->triclinic) domain->set_lamda_box();
  domain->set_local_box();
  domain->update_all_regions();

  // move atoms to new processors via irregular()

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  Irregular *irregular = new Irregular(lmp);
  irregular->migrate_atoms();
  delete irregular;
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // check if any atoms were lost

  bigint natoms;
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (natoms != atom->natoms) {
    char str[128];
    sprintf(str,"Lost atoms via balance: original " BIGINT_FORMAT
            " current " BIGINT_FORMAT,atom->natoms,natoms);
    error->all(FLERR,str);
  }

  // imbfinal = final imbalance based on final nlocal

  int maxfinal;
  double imbfinal = imbalance_nlocal(maxfinal);

  if (me == 0) {
    if (screen) {
      fprintf(screen,"  iteration count = %d\n",niter);
      fprintf(screen,"  initial/final max atoms/proc = %d %d\n",
              maxinit,maxfinal);
      fprintf(screen,"  initial/final imbalance factor = %g %g\n",
              imbinit,imbfinal);
    }
    if (logfile) {
      fprintf(logfile,"  iteration count = %d\n",niter);
      fprintf(logfile,"  initial/final max atoms/proc = %d %d\n",
              maxinit,maxfinal);
      fprintf(logfile,"  initial/final imbalance factor = %g %g\n",
              imbinit,imbfinal);
    }
  }

  if (me == 0) {
    if (screen) {
      fprintf(screen,"  x cuts:");
      for (int i = 0; i <= comm->procgrid[0]; i++)
        fprintf(screen," %g",comm->xsplit[i]);
      fprintf(screen,"\n");
      fprintf(screen,"  y cuts:");
      for (int i = 0; i <= comm->procgrid[1]; i++)
        fprintf(screen," %g",comm->ysplit[i]);
      fprintf(screen,"\n");
      fprintf(screen,"  z cuts:");
      for (int i = 0; i <= comm->procgrid[2]; i++)
        fprintf(screen," %g",comm->zsplit[i]);
      fprintf(screen,"\n");
    }
    if (logfile) {
      fprintf(logfile,"  x cuts:");
      for (int i = 0; i <= comm->procgrid[0]; i++)
        fprintf(logfile," %g",comm->xsplit[i]);
      fprintf(logfile,"\n");
      fprintf(logfile,"  y cuts:");
      for (int i = 0; i <= comm->procgrid[1]; i++)
        fprintf(logfile," %g",comm->ysplit[i]);
      fprintf(logfile,"\n");
      fprintf(logfile,"  z cuts:");
      for (int i = 0; i <= comm->procgrid[2]; i++)
        fprintf(logfile," %g",comm->zsplit[i]);
      fprintf(logfile,"\n");
    }
  }
}

/* ----------------------------------------------------------------------
   calculate imbalance based on nlocal
   return max = max atom per proc
   return imbalance factor = max atom per proc / ave atom per proc
------------------------------------------------------------------------- */

double Balance::imbalance_nlocal(int &max)
{
  MPI_Allreduce(&atom->nlocal,&max,1,MPI_INT,MPI_MAX,world);
  double imbalance = 1.0;
  if (max) imbalance = max / (1.0 * atom->natoms / nprocs);
  return imbalance;
}

/* ----------------------------------------------------------------------
   calculate imbalance based on processor splits in 3 dims
   atoms must be in lamda coords (0-1) before called
   map atoms to 3d grid of procs
   return max = max atom per proc
   return imbalance factor = max atom per proc / ave atom per proc
------------------------------------------------------------------------- */

double Balance::imbalance_splits(int &max)
{
  double *xsplit = comm->xsplit;
  double *ysplit = comm->ysplit;
  double *zsplit = comm->zsplit;

  int nx = comm->procgrid[0];
  int ny = comm->procgrid[1];
  int nz = comm->procgrid[2];

  for (int i = 0; i < nprocs; i++) proccount[i] = 0;

  double **x = atom->x;
  int nlocal = atom->nlocal;
  int ix,iy,iz;

  for (int i = 0; i < nlocal; i++) {
    ix = binary(x[i][0],nx,xsplit);
    iy = binary(x[i][1],ny,ysplit);
    iz = binary(x[i][2],nz,zsplit);
    proccount[iz*nx*ny + iy*nx + ix]++;
  }

  MPI_Allreduce(proccount,allproccount,nprocs,MPI_INT,MPI_SUM,world);
  max = 0;
  for (int i = 0; i < nprocs; i++) max = MAX(max,allproccount[i]);
  double imbalance = 1.0;
  if (max) imbalance = max / (1.0 * atom->natoms / nprocs);
  return imbalance;
}

/* ----------------------------------------------------------------------
   setup static load balance operations
   called from command
   set rho = 0 for static balancing
------------------------------------------------------------------------- */

void Balance::static_setup(char *str)
{
  ndim = strlen(str);
  bdim = new int[ndim];

  for (size_t i = 0; i < strlen(str); i++) { //NP modified R.B.
    if (str[i] == 'x') bdim[i] = X;
    if (str[i] == 'y') bdim[i] = Y;
    if (str[i] == 'z') bdim[i] = Z;
  }

  int max = MAX(comm->procgrid[0],comm->procgrid[1]);
  max = MAX(max,comm->procgrid[2]);

  count = new bigint[max];
  onecount = new bigint[max];
  sum = new bigint[max+1];
  target = new bigint[max+1];
  lo = new double[max+1];
  hi = new double[max+1];
  losum = new bigint[max+1];
  hisum = new bigint[max+1];

  rho = 0;
}

/* ----------------------------------------------------------------------
   setup dynamic load balance operations
   called from fix balance
   set rho = 1 for dynamic balancing after call to dynamic_setup()
------------------------------------------------------------------------- */

void Balance::dynamic_setup(char *str, int nitermax_in, double thresh_in)
{
  dflag = 1;
  static_setup(str);
  nitermax = nitermax_in;
  thresh = thresh_in;

  rho = 1;
}

/* ----------------------------------------------------------------------
   load balance by changing xyz split proc boundaries in Comm
   called one time from input script command or many times from fix balance
   return niter = iteration count
------------------------------------------------------------------------- */

int Balance::dynamic()
{
  int i,j,k,m,np,max;
  double *split = NULL, *split_old = NULL; //NP modified C.K.

  // no balancing if no atoms

  bigint natoms = atom->natoms;
  if (natoms == 0) return 0;

  // set delta for 1d balancing = root of threshhold
  // root = # of dimensions being balanced on

  double delta = pow(thresh,1.0/ndim) - 1.0;
  int *procgrid = comm->procgrid;

  // all balancing done in lamda coords

  domain->x2lamda(atom->nlocal);

  // loop over dimensions in balance string

  int niter = 0;
  for (int idim = 0; idim < ndim; idim++) {

    // split = ptr to xyz split in Comm

    if (bdim[idim] == X) split = comm->xsplit;
    else if (bdim[idim] == Y) split = comm->ysplit;
    else if (bdim[idim] == Z) split = comm->zsplit;

    // intial count and sum

    np = procgrid[bdim[idim]];
    tally(bdim[idim],np,split);

    //NP modified C.K. begin
    memory->destroy(split_old);
    memory->create(split_old,np+1,"balance:split_old");
    vectorCopyN(split,split_old,np+1);
    //NP modified C.K. end

    // target[i] = desired sum at split I

    for (i = 0; i < np; i++)
      target[i] = static_cast<int> (1.0*natoms/np * i + 0.5);
    target[np] = natoms;

    // lo[i] = closest split <= split[i] with a sum <= target
    // hi[i] = closest split >= split[i] with a sum >= target

    lo[0] = hi[0] = 0.0;
    lo[np] = hi[np] = 1.0;
    losum[0] = hisum[0] = 0;
    losum[np] = hisum[np] = natoms;

    for (i = 1; i < np; i++) {
      for (j = i; j >= 0; j--)
        if (sum[j] <= target[i]) {
          lo[i] = split[j];
          losum[i] = sum[j];
          break;
        }
      for (j = i; j <= np; j++)
        if (sum[j] >= target[i]) {
          hi[i] = split[j];
          hisum[i] = sum[j];
          break;
        }
    }

    // iterate until balanced

#ifdef BALANCE_DEBUG
    if (me == 0) debug_output(idim,0,np,split);
#endif

    int doneflag;
    int change = 1;
    for (m = 0; m < nitermax; m++) {
      change = adjust(np,split);
      tally(bdim[idim],np,split);
      niter++;

#ifdef BALANCE_DEBUG
      if (me == 0) debug_output(idim,m+1,np,split);
      if (me == 0 && fp) dumpout(update->ntimestep,fp);
#endif

      // stop if no change in splits, b/c all targets are met exactly

      if (!change) break;

      // stop if all split sums are within delta of targets
      // this is a 1d test of particle count per slice
      // assumption is that this is sufficient accuracy
      //   for 3d imbalance factor to reach threshhold

      doneflag = 1;
      for (i = 1; i < np; i++)
        if (fabs(1.0*(sum[i]-target[i]))/target[i] > delta) doneflag = 0;
      if (doneflag) break;
    }

    // eliminate final adjacent splits that are duplicates
    // can happen if particle distribution is narrow and Nitermax is small
    // set lo = midpt between splits
    // spread duplicates out evenly between bounding midpts with non-duplicates
    // i,j = lo/hi indices of set of duplicate splits
    // delta = new spacing between duplicates
    // bounding midpts = lo[i-1] and lo[j]

    int duplicate = 0;
    for (i = 1; i < np-1; i++)
      if (split[i] == split[i+1]) duplicate = 1;
    if (duplicate) {
      for (i = 0; i < np; i++)
        lo[i] = 0.5 * (split[i] + split[i+1]);
      i = 1;
      while (i < np-1) {
        j = i+1;
        while (split[j] == split[i]) j++;
        j--;
        if (j > i) {
          delta = (lo[j] - lo[i-1]) / (j-i+2);
          for (k = i; k <= j; k++)
            split[k] = lo[i-1] + (k-i+1)*delta;
        }
        i = j+1;
      }
    }

    //NP modified C.K. begin
    //NP in case irregular is disallowed, take care borders are not shifted too far
    if(disallow_irregular())
    {
        //NP skin and cut in lamda coords
        double skin = neighbor->skin / (domain->boxhi[bdim[idim]]-domain->boxlo[bdim[idim]]);
        double cut = force->pair->cutforce / (domain->boxhi[bdim[idim]]-domain->boxlo[bdim[idim]]);
        for (i = 1; i < np; i++)
        {
            if(cut > 0. && split[i]-split[i-1] < 1.5*cut)
                split[i] = split[i-1]+1.5*cut;
            if(split[i] < split_old[i-1]+skin)
                split[i] = split_old[i-1]+skin;
            if(split[i] > split_old[i+1]-skin)
                split[i] = split_old[i+1]-skin;

            //split[i] = split_old[i];
            /*NL*/ //if (screen) fprintf(screen,"proc %d - cut %f\n",me,cut);
            /*NL*/ //if (screen) fprintf(screen,"proc %d - dim %d i %d: old %f, new %f\n",me,idim,i,split_old[i],split[i]);
        }
    }
    //NP modified C.K. end

    // sanity check on bad duplicate or inverted splits
    // zero or negative width sub-domains will break Comm class
    // should never happen if recursive multisection algorithm is correct

    int bad = 0;
    for (i = 0; i < np; i++)
      if (split[i] >= split[i+1]) bad = 1;
    if (bad) error->all(FLERR,"Balance produced bad splits");
    /*
      if (me == 0) {
      printf("BAD SPLITS %d %d %d\n",np+1,niter,delta);
      for (i = 0; i < np+1; i++)
      printf(" %g",split[i]);
      printf("\n");
      }
    */

    // stop at this point in bstr if imbalance factor < threshhold
    // this is a true 3d test of particle count per processor

    double imbfactor = imbalance_splits(max);
    if (imbfactor <= thresh) break;
  }

  memory->destroy(split_old);   //NP modified C.K.

  // restore real coords

  domain->lamda2x(atom->nlocal);

  return niter;
}

/* ----------------------------------------------------------------------
   count atoms in each slice, based on their dim coordinate
   N = # of slices
   split = N+1 cuts between N slices
   return updated count = particles per slice
   retrun updated sum = cummulative count below each of N+1 splits
   use binary search to find which slice each atom is in
------------------------------------------------------------------------- */

void Balance::tally(int dim, int n, double *split)
{
  for (int i = 0; i < n; i++) onecount[i] = 0;

  double **x = atom->x;
  int nlocal = atom->nlocal;
  int index;

  for (int i = 0; i < nlocal; i++) {
    index = binary(x[i][dim],n,split);
    onecount[index]++;
  }

  MPI_Allreduce(onecount,count,n,MPI_LMP_BIGINT,MPI_SUM,world);

  sum[0] = 0;
  for (int i = 1; i < n+1; i++)
    sum[i] = sum[i-1] + count[i-1];
}

/* ----------------------------------------------------------------------
   adjust cuts between N slices in a dim via recursive multisectioning method
   split = current N+1 cuts, with 0.0 and 1.0 at end points
   sum = cummulative count up to each split
   target = desired cummulative count up to each split
   lo/hi = split values that bound current split
   update lo/hi to reflect sums at current split values
   overwrite split with new cuts
     guaranteed that splits will remain in ascending order,
     though adjacent values may be identical
   recursive bisectioning zooms in on each cut by halving lo/hi
   return 0 if no changes in any splits, b/c they are all perfect
------------------------------------------------------------------------- */

int Balance::adjust(int n, double *split)
{
  int i;
  double fraction;

  // reset lo/hi based on current sum and splits
  // insure lo is monotonically increasing, ties are OK
  // insure hi is monotonically decreasing, ties are OK
  // this effectively uses info from nearby splits
  // to possibly tighten bounds on lo/hi

  for (i = 1; i < n; i++) {
    if (sum[i] <= target[i]) {
      lo[i] = split[i];
      losum[i] = sum[i];
    }
    if (sum[i] >= target[i]) {
      hi[i] = split[i];
      hisum[i] = sum[i];
    }
  }
  for (i = 1; i < n; i++)
    if (lo[i] < lo[i-1]) {
      lo[i] = lo[i-1];
      losum[i] = losum[i-1];
    }
  for (i = n-1; i > 0; i--)
    if (hi[i] > hi[i+1]) {
      hi[i] = hi[i+1];
      hisum[i] = hisum[i+1];
    }

  int change = 0;
  for (int i = 1; i < n; i++)
    if (sum[i] != target[i]) {
      change = 1;
      if (rho == 0) split[i] = 0.5 * (lo[i]+hi[i]);
      else {
        fraction = 1.0*(target[i]-losum[i]) / (hisum[i]-losum[i]);
        split[i] = lo[i] + fraction * (hi[i]-lo[i]);
      }
    }
  return change;
}

/* ----------------------------------------------------------------------
   OLD code: for local diffusion method that didn't work as well as RCB
   adjust cuts between N slices in a dim via diffusive method
   count = atoms per slice
   split = current N+1 cuts, with 0.0 and 1.0 at end points
   overwrite split with new cuts
   diffusion means slices with more atoms than their neighbors "send" atoms,
     by moving cut closer to sender, further from receiver
------------------------------------------------------------------------- */

void Balance::old_adjust(int iter, int n, bigint *count, double *split)
{
  // need to allocate this if start using it again

  double *cuts = NULL;

  // damping factor

  double damp = 0.5;

  // loop over slices
  // cut I is between 2 slices (I-1 and I) with counts
  // cut I+1 is between 2 slices (I and I+1) with counts
  // for a cut between 2 slices, only slice with larger count adjusts it
  // special treatment of end slices with only 1 neighbor

  bigint leftcount,mycount,rightcount;
  double rho,target; //NP modified R.B.

  for (int i = 0; i < n; i++) {
    if (i == 0) leftcount = MAXBIGINT;
    else leftcount = count[i-1];
    mycount = count[i];
    if (i == n-1) rightcount = MAXBIGINT;
    else rightcount = count[i+1];

    // middle slice is <= both left and right, so do nothing
    // special case if 2 slices both have count = 0 -> no change in cut

    if (mycount <= leftcount && mycount <= rightcount) {
      if (leftcount == 0) cuts[i] = split[i];
      if (rightcount == 0) cuts[i+1] = split[i+1];
      continue;
    }

    // rho = density of atoms in the slice

    rho = mycount / (split[i+1] - split[i]);

    // middle slice has more atoms than left or right slice
    // send atoms in that dir

    if (mycount > leftcount) {
      target = damp * 0.5*(mycount-leftcount);
      cuts[i] = split[i] + target/rho;
    }
    if (mycount > rightcount) {
      target = damp * 0.5*(mycount-rightcount);
      cuts[i+1] = split[i+1] - target/rho;
    }

    /*
    // middle slice has more atoms then left or right slice
    // if delta from middle to top slice > delta between top and bottom slice
    //   then send atoms both dirs to bring all 3 slices to same count
    // else bottom slice is very low, so send atoms only in that dir

    if (mycount > leftcount && mycount > rightcount) {
      if (mycount-MAX(leftcount,rightcount) >= fabs(leftcount-rightcount)) {
        if (leftcount <= rightcount) {
          targetleft = damp *
            (rightcount-leftcount + (mycount-rightcount)/3.0);
          targetright = damp * (mycount-rightcount)/3.0;
          cuts[i] = split[i] + targetleft/rho;
          cuts[i+1] = split[i+1] - targetright/rho;
        } else {
          targetleft = damp * (mycount-leftcount)/3.0;
          targetright = damp *
            (leftcount-rightcount + (mycount-leftcount)/3.0);
          cuts[i] = split[i] + targetleft/rho;
          cuts[i+1] = split[i+1] - targetright/rho;
        }
      } else if (leftcount < rightcount) {
        target = damp * 0.5*(mycount-leftcount);
        cuts[i] = split[i] + target/rho;
        cuts[i+1] = split[i+1];
      } else if (rightcount < leftcount) {
        target = damp * 0.5*(mycount-rightcount);
        cuts[i+1] = split[i+1] - target/rho;
        cuts[i] = split[i];
      }

    // middle slice has more atoms than only left or right slice
    // send atoms only in that dir

    } else if (mycount > leftcount) {
      target = damp * 0.5*(mycount-leftcount);
      cuts[i] = split[i] + target/rho;
    } else if (mycount > rightcount) {
      target = damp * 0.5*(mycount-rightcount);
      cuts[i+1] = split[i+1] - target/rho;
    }
    */
  }

  // overwrite adjustable splits with new cuts

  for (int i = 1; i < n; i++) split[i] = cuts[i];
}

/* ----------------------------------------------------------------------
   binary search for where value falls in N-length vec
   note that vec actually has N+1 values, but ignore last one
   values in vec are monotonically increasing, but adjacent values can be ties
   value may be outside range of vec limits
   always return index from 0 to N-1 inclusive
   return 0 if value < vec[0]
   reutrn N-1 if value >= vec[N-1]
   return index = 1 to N-2 inclusive if vec[index] <= value < vec[index+1]
   note that for adjacent tie values, index of lower tie is not returned
     since never satisfies 2nd condition that value < vec[index+1]
------------------------------------------------------------------------- */

int Balance::binary(double value, int n, double *vec)
{
  int lo = 0;
  int hi = n-1;

  if (value < vec[lo]) return lo;
  if (value >= vec[hi]) return hi;

  // insure vec[lo] <= value < vec[hi] at every iteration
  // done when lo,hi are adjacent

  int index = (lo+hi)/2;
  while (lo < hi-1) {
    if (value < vec[index]) hi = index;
    else if (value >= vec[index]) lo = index;
    index = (lo+hi)/2;
  }

  return index;
}

/* ----------------------------------------------------------------------
   write dump snapshot of line segments in Pizza.py mdump mesh format
   write xy lines around each proc's sub-domain for 2d
   write xyz cubes around each proc's sub-domain for 3d
   only called by proc 0
------------------------------------------------------------------------- */

void Balance::dumpout(bigint tstep, FILE *bfp)
{
  int dimension = domain->dimension;

  // write out one square/cube per processor for 2d/3d
  // only write once since topology is static

  if (firststep) {
    firststep = 0;
    fprintf(bfp,"ITEM: TIMESTEP\n");
    fprintf(bfp,BIGINT_FORMAT "\n",tstep);
    if (dimension == 2) fprintf(bfp,"ITEM: NUMBER OF SQUARES\n");
    else fprintf(bfp,"ITEM: NUMBER OF CUBES\n");
    fprintf(bfp,"%d\n",nprocs);
    if (dimension == 2) fprintf(bfp,"ITEM: SQUARES\n");
    else fprintf(bfp,"ITEM: CUBES\n");

    int nx = comm->procgrid[0] + 1;
    int ny = comm->procgrid[1] + 1;
    //int nz = comm->procgrid[2] + 1; //NP modified R.B.

    if (dimension == 2) {
      int m = 0;
      for (int j = 0; j < comm->procgrid[1]; j++)
        for (int i = 0; i < comm->procgrid[0]; i++) {
          int c1 = j*nx + i + 1;
          int c2 = c1 + 1;
          int c3 = c2 + nx;
          int c4 = c3 - 1;
          fprintf(bfp,"%d %d %d %d %d %d\n",m+1,m+1,c1,c2,c3,c4);
          m++;
        }

    } else {
      int m = 0;
      for (int k = 0; k < comm->procgrid[2]; k++)
        for (int j = 0; j < comm->procgrid[1]; j++)
          for (int i = 0; i < comm->procgrid[0]; i++) {
            int c1 = k*ny*nx + j*nx + i + 1;
            int c2 = c1 + 1;
            int c3 = c2 + nx;
            int c4 = c3 - 1;
            int c5 = c1 + ny*nx;
            int c6 = c2 + ny*nx;
            int c7 = c3 + ny*nx;
            int c8 = c4 + ny*nx;
            fprintf(bfp,"%d %d %d %d %d %d %d %d %d %d\n",
                    m+1,m+1,c1,c2,c3,c4,c5,c6,c7,c8);
            m++;
          }
    }
  }

  // write out nodal coords, can be different every call
  // scale xsplit,ysplit,zsplit values to full box
  // only implmented for orthogonal boxes, not triclinic

  int nx = comm->procgrid[0] + 1;
  int ny = comm->procgrid[1] + 1;
  int nz = comm->procgrid[2] + 1;

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  double *prd = domain->prd;

  fprintf(bfp,"ITEM: TIMESTEP\n");
  fprintf(bfp,BIGINT_FORMAT "\n",tstep);
  fprintf(bfp,"ITEM: NUMBER OF NODES\n");
  if (dimension == 2) fprintf(bfp,"%d\n",nx*ny);
  else fprintf(bfp,"%d\n",nx*ny*nz);
  fprintf(bfp,"ITEM: BOX BOUNDS\n");
  fprintf(bfp,"%g %g\n",boxlo[0],boxhi[0]);
  fprintf(bfp,"%g %g\n",boxlo[1],boxhi[1]);
  fprintf(bfp,"%g %g\n",boxlo[2],boxhi[2]);
  fprintf(bfp,"ITEM: NODES\n");

  if (dimension == 2) {
    int m = 0;
    for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++) {
        fprintf(bfp,"%d %d %g %g %g\n",m+1,1,
                boxlo[0] + prd[0]*comm->xsplit[i],
                boxlo[1] + prd[1]*comm->ysplit[j],
                0.0);
        m++;
      }
  } else {
    int m = 0;
    for (int k = 0; k < nz; k++)
      for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
          fprintf(bfp,"%d %d %g %g %g\n",m+1,1,
                  boxlo[0] + prd[0]*comm->xsplit[i],
                  boxlo[1] + prd[1]*comm->ysplit[j],
                  boxlo[2] + prd[2]*comm->zsplit[k]);
          m++;
      }
  }
}

/* ----------------------------------------------------------------------
   debug output for Idim and count
   only called by proc 0
------------------------------------------------------------------------- */

void Balance::debug_output(int idim, int m, int np, double *split)
{
  int i;
  const char *dim = NULL;

  double *boxlo = domain->boxlo;
  double *prd = domain->prd;

  if (bdim[idim] == X) dim = "X";
  else if (bdim[idim] == Y) dim = "Y";
  else if (bdim[idim] == Z) dim = "Z";
  printf("Dimension %s, Iteration %d\n",dim,m);

  printf("  Count:");
  for (i = 0; i < np; i++) printf(" " BIGINT_FORMAT,count[i]);
  printf("\n");
  printf("  Sum:");
  for (i = 0; i <= np; i++) printf(" " BIGINT_FORMAT,sum[i]);
  printf("\n");
  printf("  Target:");
  for (i = 0; i <= np; i++) printf(" " BIGINT_FORMAT,target[i]);
  printf("\n");
  printf("  Actual cut:");
  for (i = 0; i <= np; i++)
    printf(" %g",boxlo[bdim[idim]] + split[i]*prd[bdim[idim]]);
  printf("\n");
  printf("  Split:");
  for (i = 0; i <= np; i++) printf(" %g",split[i]);
  printf("\n");
  printf("  Low:");
  for (i = 0; i <= np; i++) printf(" %g",lo[i]);
  printf("\n");
  printf("  Low-sum:");
  for (i = 0; i <= np; i++) printf(" " BIGINT_FORMAT,losum[i]);
  printf("\n");
  printf("  Hi:");
  for (i = 0; i <= np; i++) printf(" %g",hi[i]);
  printf("\n");
  printf("  Hi-sum:");
  for (i = 0; i <= np; i++) printf(" " BIGINT_FORMAT,hisum[i]);
  printf("\n");
  printf("  Delta:");
  for (i = 0; i < np; i++) printf(" %g",split[i+1]-split[i]);
  printf("\n");

  bigint max = 0;
  for (i = 0; i < np; i++) max = MAX(max,count[i]);
  printf("  Imbalance factor: %g\n",1.0*max*np/target[np]);
}

//NP modified C.K.

/* ----------------------------------------------------------------------
   disallow irregular comm for some cases
------------------------------------------------------------------------- */

bool Balance::disallow_irregular()
{
    /*NL*///if (screen) fprintf(screen,"proc %d at step %d: result %d\n",me,update->ntimestep,modify->n_fixes_style("mesh"));
    //NP disallow if mesh is involved
    if(modify->n_fixes_style("mesh") > 0)
        return true;
    //NP disallow if multisphere is used
    //NP b/c body data is also communicated with exchange()-like
    if(modify->n_fixes_style("multisphere") > 0)
        return true;

    return false;
}
