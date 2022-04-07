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

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "fix_remove_molecule.h"
#include "update.h"
#include "respa.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "neighbor.h"
#include "error.h"
#include "modify.h"
#include "domain.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixRemoveMolecule::FixRemoveMolecule(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
    // error check
    if (atom->molecular == 0)
        error->all(FLERR,"Cannot use fix remove molecule with non-molecular systems");
}

/* ---------------------------------------------------------------------- */

FixRemoveMolecule::~FixRemoveMolecule()
{
}

/* ---------------------------------------------------------------------- */

int FixRemoveMolecule::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRemoveMolecule::init()
{
    if (strcmp(update->integrate_style,"respa") == 0)
        nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixRemoveMolecule::post_integrate()
{
    // Should try to use delete_atom class
    const int nbondlist = neighbor->nbondlist;
    int * const * const bondlist = neighbor->bondlist;
    int * const mol = atom->molecule;

    double * const * const x = atom->x;
    double * const radius = atom->radius;

    for (int k = 0; k < nbondlist; k++) {
        if (bondlist[k][3] != 1) continue;

        // now check if the atom is inside the boundary
        const int i1 = bondlist[k][0];
        const int i2 = bondlist[k][1];
        if (!((x[i1][0] - 0.25*radius[i1] < domain->boxlo[0] || x[i1][0] + 0.25*radius[i1] > domain->boxhi[0]) && domain->xperiodic == 0)) continue;
        if (!((x[i1][1] - 0.25*radius[i1] < domain->boxlo[1] || x[i1][1] + 0.25*radius[i1] > domain->boxhi[1]) && domain->yperiodic == 0)) continue;
        if (!((x[i1][2] - 0.25*radius[i1] < domain->boxlo[2] || x[i1][2] + 0.25*radius[i1] > domain->boxhi[2]) && domain->zperiodic == 0)) continue;

        if (!((x[i2][0] - 0.25*radius[i2] < domain->boxlo[0] || x[i2][0] + 0.25*radius[i2] > domain->boxhi[0]) && domain->xperiodic == 0)) continue;
        if (!((x[i2][1] - 0.25*radius[i2] < domain->boxlo[1] || x[i2][1] + 0.25*radius[i2] > domain->boxhi[1]) && domain->yperiodic == 0)) continue;
        if (!((x[i2][2] - 0.25*radius[i2] < domain->boxlo[2] || x[i2][2] + 0.25*radius[i2] > domain->boxhi[2]) && domain->zperiodic == 0)) continue;

        if (mol[i1] == mol[i2]) {
            remove_molecule(mol[i1]);
        } else {
            remove_molecule(mol[i1]);
            remove_molecule(mol[i2]);
        }

    }

}

/* ---------------------------------------------------------------------- */

void FixRemoveMolecule::remove_molecule(int mol_id)
{
    const int nlocal = atom->nlocal;

    int * const tag = atom->tag;
    int * const mol = atom->molecule;

    const int nprocs = comm->nprocs;

    int * mol_ids = new int[nprocs];

    // Need to get all molecule ids from all processors
    MPI_Allgather(&mol_id, 1, MPI_INT, mol_ids, 1, MPI_INT, world);

    int kk = 0;
    int * atom_tags = new int[512]; // remove upto 512 atoms each time...
    for (int m = 0; m < nprocs; m++) {
        for (int k = 0; k < nlocal; k++) {
            if (mol[k] == mol_ids[m] && kk < 512) {
                atom_tags[kk++] = tag[k];
            }
        }
    }

    remove_atoms(atom_tags, kk);

    // We have removed the atoms on our processor that had the above mol id...
    // Now we need to tell the other procs to do the same.

    delete [] mol_ids;
    delete [] atom_tags;
}

/* ---------------------------------------------------------------------- */

void FixRemoveMolecule::remove_atoms(int * atom_tag, int numAtoms)
{
    // Delete atom with the specific tag
    int nlocal = atom->nlocal;

    int * const tag = atom->tag;
    // delete local atoms flagged in dlist
    // reset nlocal
    bigint natoms_previous = atom->natoms;

    AtomVec *avec = atom->avec;

    int i = 0;
    while (i < nlocal) {
        bool isInAtomTag = false;
        for (int k = 0; k < numAtoms; k++) {
            if (tag[i] == atom_tag[k]) {
                isInAtomTag = true;
                break;
            }
        }
        if (isInAtomTag) {
            avec->copy(nlocal-1,i,1);
            nlocal--;
            break;
        } else i++;
    }

    atom->nlocal = nlocal;

    // reset atom->natoms
    // reset atom->map if it exists
    // set nghost to 0 so old ghosts of deleted atoms won't be mapped

    bigint nblocal = atom->nlocal;
    MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
    if (atom->map_style) {
        atom->nghost = 0;
        atom->map_init();
        atom->map_set();
    }

    // print before and after atom count

    bigint ndelete = natoms_previous - atom->natoms;

    if (comm->me == 0) {
        if (screen) fprintf(screen,"Deleted " BIGINT_FORMAT
                            " atoms, new total = " BIGINT_FORMAT "\n",
                            ndelete,atom->natoms);
        if (logfile) fprintf(logfile,"Deleted " BIGINT_FORMAT
                            " atoms, new total = " BIGINT_FORMAT "\n",
                            ndelete,atom->natoms);
    }
}

/* ---------------------------------------------------------------------- */

void FixRemoveMolecule::post_integrate_respa(int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_integrate();
}
