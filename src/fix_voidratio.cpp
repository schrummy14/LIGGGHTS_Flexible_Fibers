
#include <stdio.h>
#include <climits>
#include <stdlib.h>
#include <string.h>
#include "domain.h"
#include "region.h"
#include "mpi_liggghts.h"
#include "comm.h"
#include "update.h"
#include "atom.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "fix_voidratio.h"
#include "respa.h"

using namespace LAMMPS_NS;
using namespace FixConst;

FixVoidRatio::FixVoidRatio(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
    fprintf(screen, "WARNING: Fix voidratio is not coded for superquadratics\n");
    if (narg < 9) error->all(FLERR,"Illegal fix voidratio command");
    iregion = -1;
    idregion = NULL;
    doEndOfRun = false;
    vector_flag = 1;
    size_vector = 4;
    global_freq = 1;
    extvector = 1;
    nextCheck = update->ntimestep+1;

    int iarg = 3;
    while (iarg < narg) {
        if (strcmp(arg[iarg],"nevery") == 0) {
            if (iarg+1 > narg) error->all(FLERR,"Illegal fix voidratio command");
            nevery = force->inumeric(FLERR,arg[iarg+1]);
            if (nevery == 0) doEndOfRun = true;
            if (nevery < 0) error->all(FLERR,"Illegal fix voidratio command: nevery cannot be negative");
            iarg += 2;
        } else if (strcmp(arg[iarg],"region") == 0) {
            if (iarg+1 > narg) error->all(FLERR,"Illegal fix voidratio command");
            iregion = domain->find_region(arg[iarg+1]);
            if (iregion == -1) error->all(FLERR, "Illegal fix voidration command: Could not find region");
            int n = strlen(arg[iarg+1]) + 1;
            idregion = new char[n];
            strcpy(idregion, arg[iarg+1]);
            iarg += 2;
        } else if (strcmp(arg[iarg],"ntry") == 0) {
            if (iarg+1 > narg) error->all(FLERR,"Illegal fix voidratio command");
            ntry = force->inumeric(FLERR,arg[iarg+1]);
            iarg += 2;
        } else error->all(FLERR,"Illegal fix voidratio command");
    }

    if (iregion == -1) error->all(FLERR,"Illegal fix voidratio command: Must provide a region");
    region_volume = -1.0;
}


FixVoidRatio::~FixVoidRatio()
{
    delete [] idregion;
}

int FixVoidRatio::setmask()
{
    int mask = 0;
    mask |= POST_FORCE;
    mask |= POST_FORCE_RESPA;
    mask |= MIN_POST_FORCE;
    return mask;
}

void FixVoidRatio::init()
{
  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixVoidRatio::setup(int vflag)
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

void FixVoidRatio::post_force(int vflag)
{
    if (doEndOfRun) {
        if (update->ntimestep == update->endstep) {
            if (comm->me == 0)
                fprintf(screen, "Updating voidratio at step %li\n", update->ntimestep);
            voidratio = getVoidRatio();
        }
    } else if (update->ntimestep == nextCheck /*|| update->ntimestep == update->endstep*/) {
        if (comm->me == 0)
            fprintf(screen, "Updating voidratio at step %li\n", update->ntimestep);
        nextCheck += static_cast<long int>(nevery);
        voidratio = getVoidRatio();
    }
}

void FixVoidRatio::post_run() //end_of_step()
{
    // if (!doEndOfRun) return;
    // if (comm->me == 0)
    //     fprintf(screen, "Updating voidratio\n");
    // voidratio = getVoidRatio();
}

double FixVoidRatio::getVoidRatio()
{
    int nlocal = atom->nlocal;
    double *r = atom->radius;
    double **x = atom->x;

    double randX[3];

    double dx, dy, dz, atom_volume_local;
    
    unsigned long pointInside = 0;

    // Get volume of region
    if (region_volume < 0 || domain->regions[iregion]->dynamic_check() == 1) {
        int volTry = INT_MAX < ntry ? INT_MAX : ntry;
        domain->regions[iregion]->volume_mc(volTry ,false, 0.0, region_volume, region_volume_local);
    }

    if (nlocal > 0 && region_volume_local > 0.0) {
        for (unsigned long i = 0; i < ntry; i++) {
            // Get random point in local region
            domain->regions[iregion]->generate_random(randX, true);
            for (int j = 0; j < nlocal; j++) {
                dx = randX[0] - x[j][0];
                dy = randX[1] - x[j][1];
                dz = randX[2] - x[j][2];

                // See if random point is inside an atom
                // For superquadrics
                // https://cse.buffalo.edu/~jryde/cse673/files/superquadrics.pdf Equations 2.16
                // https://link.springer.com/article/10.1007/s40571-016-0131-6 Equation 1
                if (dx*dx + dy*dy + dz*dz <= r[j]*r[j]) {
                    pointInside++;
                    break;
                }
            }
        }
    }

    // Calculate the local atom volume
    atom_volume_local = region_volume_local*static_cast<double>(pointInside)/static_cast<double>(ntry);

    // Sum up the local atom volumes to get total atom volume
    MPI_Sum_Scalar(atom_volume_local, atom_volume, world);

    // calculate the void volume from the region and atom volumes
    void_volume = region_volume - atom_volume;

    if (atom_volume == 0.0) return -2.0; // This should be something else...
    else return void_volume/atom_volume;

}

void FixVoidRatio::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

double FixVoidRatio::compute_vector(int index)
{
    if(index == 0) return voidratio;
    if(index == 1) return void_volume;
    if(index == 2) return atom_volume;
    if(index == 3) return voidratio/(1+voidratio);
    return 0.0;
}