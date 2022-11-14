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

/* ----------------------------------------------------------------------
   Contributing authors:
   WangDi Wuhan University China
   E-mail: wangdi1010@whu.edu.cn
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(domain,FixDomain)

#else

#ifndef LMP_FIX_DOMAIN_H
#define LMP_FIX_DOMAIN_H

#include "fix.h"


namespace LAMMPS_NS
{
  class FixDomain : public Fix
  {
      public:
        FixDomain(class LAMMPS *, int, char **);
        ~FixDomain();

        int setmask();
        void init();
        virtual double compute_vector(int n);
        virtual void final_integrate();
      private:
        double xlo,xhi,ylo,yhi,zlo,zhi;
  };

} /* namespace LAMMPS_NS */

#endif
#endif
