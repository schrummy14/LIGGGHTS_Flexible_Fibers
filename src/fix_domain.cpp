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

#include "fix_domain.h"
#include <stdio.h>
#include "string.h"
#include "error.h"
#include "modify.h"
#include "comm.h"
#include "domain.h"
#include "fix_property_global.h"

using namespace LAMMPS_NS;
using namespace FixConst;
/* ---------------------------------------------------------------------- */

FixDomain::FixDomain(LAMMPS *lmp, int narg, char **arg)
: Fix(lmp, narg, arg)
{
    if (narg < 3) error->all(FLERR,"Illegal fix domain command");
    if (narg > 3) error->all(FLERR,"Illegal fix domain command");
    vector_flag = 1;
    size_vector = 6;
    global_freq = 1;
    extvector = 1;
    no_change_box = 1;

    xlo=0;
    xhi=0;
    ylo=0;
    yhi=0;
    zlo=0;
    zhi=0;

    if(domain->dimension ==3){
    xlo=domain->boxlo[0];
    xhi=domain->boxhi[0];
    ylo=domain->boxlo[1];
    yhi=domain->boxhi[1];
    zlo=domain->boxlo[2];
    zhi=domain->boxhi[2];
    }else error->all(FLERR,"Fix domain requests a 3 dimension box");

    int domain_check=0;
    if(xlo==xhi) domain_check=1;
    if(ylo==yhi) domain_check=1;
    if(zlo==zhi) domain_check=1;
    if(domain_check==1){
        if(xlo==0||ylo==0||zlo==0) error->all(FLERR,"Fail to get domain");
        else error->all(FLERR,"Illegal domain");
    }
}

/* ---------------------------------------------------------------------- */

FixDomain::~FixDomain()
{

}

/* ---------------------------------------------------------------------- */
void FixDomain::init()
{

}

/* ---------------------------------------------------------------------- */

int FixDomain::setmask()
{
    int mask = 0;
    mask |= FINAL_INTEGRATE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixDomain::final_integrate()
{
    int domain_check=0;
    if (domain->dimension ==3) {
        xlo=domain->boxlo[0];
        xhi=domain->boxhi[0];
        ylo=domain->boxlo[1];
        yhi=domain->boxhi[1];
        zlo=domain->boxlo[2];
        zhi=domain->boxhi[2];
    } else error->all(FLERR,"Fix domain requests a 3 dimension box");

    if(xlo==xhi) domain_check=1;
    if(ylo==yhi) domain_check=1;
    if(zlo==zhi) domain_check=1;
    if(domain_check==1){
        if(xlo==0||ylo==0||zlo==0) error->all(FLERR,"Fail to get domain");
        else error->all(FLERR,"Illegal domain");
    }
}

/* ----------------------------------------------------------------------
   return total force or torque component on body
------------------------------------------------------------------------- */

double FixDomain::compute_vector(int n)
{
  switch(n) {
      case 0 :
          return xlo;
      case 1 :
          return xhi;
      case 2 :
          return ylo;
      case 3 :
          return yhi;
      case 4 :
          return zlo;
      case 5 :
          return zhi;
  }
  error->all(FLERR,"More than 6 parameters request by fix domain");
}
