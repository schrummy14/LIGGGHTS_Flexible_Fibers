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

#ifdef COMMAND_CLASS

CommandStyle(balance,Balance)

#else

#ifndef LMP_BALANCE_H
#define LMP_BALANCE_H

#include <stdio.h>
#include "pointers.h"

namespace LAMMPS_NS {

class Balance : protected Pointers {
 public:
  Balance(class LAMMPS *);
  ~Balance();
  void command(int, char **);
  void dynamic_setup(char *, int, double);
  int dynamic();
  double imbalance_nlocal(int &);
  void dumpout(bigint, FILE *);

  bool disallow_irregular();   //NP modified C.K.

 private:
  int me,nprocs;

  int xflag,yflag,zflag;                            // xyz LB flags
  double *user_xsplit,*user_ysplit,*user_zsplit;    // params for xyz LB

  int dflag;                 // dynamic LB flag
  int nitermax;              // params for dynamic LB
  double thresh;
  char bstr[4];

  int ndim;                  // length of balance string bstr
  int *bdim;                 // XYZ for each character in bstr
  bigint *count;             // counts for slices in one dim
  bigint *onecount;          // work vector of counts in one dim
  bigint *sum;               // cummulative count for slices in one dim
  bigint *target;            // target sum for slices in one dim
  double *lo,*hi;            // lo/hi split coords that bound each target
  bigint *losum,*hisum;      // cummulative counts at lo/hi coords
  int rho;                   // 0 for geometric recursion
                             // 1 for density weighted recursion

  int *proccount;            // particle count per processor
  int *allproccount;

  int outflag;               // for output of balance results to file
  FILE *fp;
  int firststep;

  void static_setup(char *);
  double imbalance_splits(int &);
  void tally(int, int, double *);
  int adjust(int, double *);
  void old_adjust(int, int, bigint *, double *);
  int binary(double, int, double *);
  void debug_output(int, int, int, double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Balance command before simulation box is defined

The balance command cannot be used before a read_data, read_restart,
or create_box command.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot open balance output file

Self-explanatory.

E: Cannot balance in z dimension for 2d simulation

Self-explanatory.

E: Balance dynamic string is invalid

The string can only contain the characters "x", "y", or "z".

E: Lost atoms via balance: original %ld current %ld

This should not occur.  Report the problem to the developers.

E: Balance produced bad splits

This should not occur.  It means two or more cutting plane locations
are on top of each other or out of order.  Report the problem to the
developers.

*/
