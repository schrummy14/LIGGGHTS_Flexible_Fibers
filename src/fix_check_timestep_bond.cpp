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

#include <string.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include "atom.h"
#include "update.h"
#include "bond.h"
#include "error.h"
#include "fix_check_timestep_bond.h"
#include "pair_gran.h"
#include "properties.h"
#include "fix_property_global.h"
#include "force.h"
#include "comm.h"
#include "modify.h"
#include "fix_wall_gran.h"
#include "fix_mesh_surface.h"
#include "force.h"
#include "neighbor.h"
#include "mpi_liggghts.h"
#include "property_registry.h"
#include "global_properties.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MODEL_PARAMS;

/* ---------------------------------------------------------------------- */

FixCheckTimestepBond::FixCheckTimestepBond(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg),
    warnflag(true),
    errorflag(false),
    fraction_bond(0.),
    fraction_skin(0.)
{
    if (narg < 5) error->all(FLERR,"Illegal fix check/timestep/bond command, not enough arguments");

    int iarg = 5;

    if(0 == strcmp("check_every_time", arg[3])) {
        nevery = atoi(arg[4]);
        fraction_bond_lim = atof(arg[5]);
        iarg = 6;
    } else {
        nevery = atoi(arg[3]);
        fraction_bond_lim = atof(arg[4]);
        iarg = 5;
    }

    while(iarg < narg) {
        if (strcmp(arg[iarg],"warn") == 0) {
            if (narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'warn'");
            if(0 == strcmp(arg[iarg+1],"no"))
                warnflag = false;
            else if(0 == strcmp(arg[iarg+1],"yes"))
                warnflag = true;
            else
                error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'warn'");
            iarg += 2;
        } else if (strcmp(arg[iarg],"error") == 0) {
            if (narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'error'");
            if(0 == strcmp(arg[iarg+1],"no"))
                errorflag = false;
            else if(0 == strcmp(arg[iarg+1],"yes"))
                errorflag = true;
            else
                error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'error'");
            iarg += 2;
        }
    }

    vector_flag = 1;
    size_vector = 3;
    global_freq = nevery;
    extvector = 1;
}

/* ---------------------------------------------------------------------- */

int FixCheckTimestepBond::setmask()
{
    int mask = 0;
    mask |= END_OF_STEP;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixCheckTimestepBond::init()
{
    //some error checks
    if(!atom->radius_flag || !atom->density_flag || !atom->molecule_flag)
        error->all(FLERR,"Fix check/timestep/bond can only be used together with atom style hybrid sphere gran/bond");
}

/* ---------------------------------------------------------------------- */

void FixCheckTimestepBond::end_of_step()
{
    calc_bond_estims();

    const double dt = update->dt;

    fraction_bond = dt/bond_time;

    if(errorflag || (warnflag && comm->me==0)) {
        char errstr[512];

        if(fraction_bond > fraction_bond_lim)
        {
            sprintf(errstr,"time-step is %f %% of bond time",fraction_bond*100.);
            if(errorflag)
                error->fix_error(FLERR,this,errstr);
            else
                error->warning(FLERR,errstr);
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixCheckTimestepBond::calc_bond_estims()
{
    bond_time = force->bond->getMinDt();
    MPI_Min_Scalar(bond_time,world);
}

/* ----------------------------------------------------------------------
   return fraction of bond time-step
------------------------------------------------------------------------- */

double FixCheckTimestepBond::compute_vector(int n)
{
    if(n == 0)      return fraction_bond;
    return 0.;
}
