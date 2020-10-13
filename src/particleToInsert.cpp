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

#include "particleToInsert.h"
#include <cmath>
#include "error.h"
#include "update.h"
#include "domain.h"
#include "atom.h"
#include "atom_vec.h"
#include "fix_property_atom.h"
#include "vector_liggghts.h"
#include "math_extra_liggghts.h"
#include "modify.h"
#include "force.h"

using namespace LAMMPS_NS;

ParticleToInsert::ParticleToInsert(LAMMPS* lmp, int ns, FixPropertyAtom * const _fix_release) :
    Pointers(lmp),
    nspheres(ns),
    groupbit(0),
    atom_type(0),
    bond_type(0),
    density_ins(0.0),
    volume_ins(0.0),
    mass_ins(0.0),
    r_bound_ins(0.0),
    atom_type_vector_flag(false),
    id_ins(-1),
    fix_release(_fix_release),
    fix_property(NULL),
    n_fix_property(0),
    fix_property_nentry(NULL),
    fix_property_value(NULL),
    local_start(-1),
    needs_bonding(false),
    fix_template_(NULL)
{
    groupbit = 0;

    distorder = -1;

    nparticles = ns;

    memory->create(x_ins,nparticles,3,"x_ins");
    radius_ins = new double[nparticles];

    atom_type_vector = new int[nparticles];
    atom_type_vector_flag = false;

    if (fix_release && fix_release->get_nvalues() <= 14)
        error->all(FLERR, "Invalid fix_release, more than 14 entries required");
}

/* ---------------------------------------------------------------------- */

ParticleToInsert::~ParticleToInsert()
{
        memory->destroy(x_ins);
        delete []radius_ins;
        delete []atom_type_vector;
        if (fix_property_value)
        {
            for (int i = 0; i < n_fix_property; i++)
                delete [] fix_property_value[i];
            delete [] fix_property_value;
        }
        if (fix_property)
            delete [] fix_property;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::insert()
{
    // perform the actual insertion
    // add particles, set coordinate and radius
    // set group mask to "all" plus fix groups

    // Atom is the primative element (sphere)
    // Molecule is a mega-particle consisting of 1 or multiple atoms

    int inserted = 0;
    int nfix = modify->nfix;
    Fix **fix = modify->fix;

    // int molId;
    local_start = atom->nlocal;
    // if (atom->molecule_flag) {
    //     if (atom->nlocal > 0) molId = -abs(atom->molecule[atom->nlocal-1]) - 1;
    //     else molId = -3; // multisphere particles start with a value of -2, so we will start with a value of -3
    // }

    for(int i = 0; i < nparticles; i++)
    {
        /*NL*/ //if (screen) fprintf(screen,"proc %d tyring to insert particle at pos %f %f %f\n",comm->me,x_ins[i][0],x_ins[i][1],x_ins[i][2]);
        //NP do not need subdomain check any longer since have processor-local lists anyway
        //if (domain->is_in_extended_subdomain(x_ins[i]))
        //{
        /*NL*/ //if (screen) fprintf(screen,"   proc %d inserting particle at pos %f %f %f\n",comm->me,x_ins[i][0],x_ins[i][1],x_ins[i][2]);
        inserted++;
        if(atom_type_vector_flag)
            atom->avec->create_atom(atom_type_vector[i],x_ins[i]);
        else
            atom->avec->create_atom(atom_type,x_ins[i]);
        int m = atom->nlocal - 1;
        atom->mask[m] = 1 | groupbit;
        vectorCopy3D(v_ins,atom->v[m]);
        vectorCopy3D(omega_ins,atom->omega[m]);
        atom->radius[m] = radius_ins[i];
        atom->density[m] = density_ins;
        atom->rmass[m] = (nspheres==1)? (mass_ins) : (4.18879020479*radius_ins[i]*radius_ins[i]*radius_ins[i]*density_ins); // 4/3 * pi
        

        //pre_set_arrays() called via FixParticleDistribution
        for (int j = 0; j < nfix; j++)
            if (fix[j]->create_attribute) fix[j]->set_arrays(m);

        // apply fix property setting coming from fix insert
        // this overrides the set_arrays call above
        if(fix_property)
        {
            for (int j = 0; j < n_fix_property; j++)
            {
                if (fix_property_nentry[j] == 1)
                {
                    fix_property[j]->vector_atom[m] = fix_property_value[j][0];
                    if(strcmp(fix_property[j]->id,"bond_random_id") == 0)
                    {
                        if (atom->molecule_flag)
                        {
                            needs_bonding = true;
                            // use random part as base for dummy molecule ID
                            // see FixTemplateMultiplespheres::randomize_ptilist
                            double dmol = (fix_property_value[j][0] - static_cast<double>(update->ntimestep));
                            if(dmol > 1.0 || dmol < 0.0)
                                error->one(FLERR, "Internal error (particle to insert: mol id)");
                            atom->molecule[m] = -static_cast<int>(dmol * INT_MAX);
                            // actual molecule value needs to be created afterwards via atom->mol_extend()
                        }
                    }
                }
                else
                {
                    for (int k = 0; k < fix_property_nentry[j]; k++)
                        fix_property[j]->array_atom[m][k] = fix_property_value[j][k];
                }
            }
        }
        if (fix_template_)
            fix_template_->vector_atom[m] = (double)distorder;
        if (fix_release)
            fix_release->array_atom[m][14] = (double) id_ins;
    }
    
    return inserted;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::create_bond_partners(int*& npartner, int**&partner)
{
    npartner = new int[nspheres]();
    partner = new int*[nspheres];

    for(int i = 0; i < nspheres; ++i)
        partner[i] = new int[nspheres-1]();

    int create_bonds = 0;

    for(int i = 0; i < nspheres-1; ++i)
    {
        const double xtmp = x_ins[i][0];
        const double ytmp = x_ins[i][1];
        const double ztmp = x_ins[i][2];

        for(int j = i+1; j < nspheres; ++j)
        {
            const double max_bonding_dist = radius_ins[i] + radius_ins[j] + neighbor->skin;
            const double delx = xtmp - x_ins[j][0];
            const double dely = ytmp - x_ins[j][1];
            const double delz = ztmp - x_ins[j][2];
            const double rsq = delx * delx + dely * dely + delz * delz;

            if(rsq < max_bonding_dist*max_bonding_dist)
            {
                if(npartner[i] == nspheres-1 || npartner[j] == nspheres-1)
                {
                    // should not happen print warning
                    continue;
                }
                partner[i][npartner[i]] = j;
                partner[j][npartner[j]] = i;
                npartner[i]++;
                npartner[j]++;
                create_bonds = 2;
            }
        }
    }

    return create_bonds;
}

/* ---------------------------------------------------------------------- */

void ParticleToInsert::destroy_bond_partners(int* npartner, int** partner)
{
    for(int i = 0; i < nspheres; ++i)
        delete [] partner[i];
    delete [] partner;
    delete [] npartner;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::create_bonds(int *npartner, int **partner)
{
    // fprintf(screen, "\n--------------------------------------\nIn create bonds\n");
    if(nspheres == 1 || !needs_bonding || local_start < 0)
        return 0;

    needs_bonding = false; // reset in case pti gets reused

    int create_bonds = 1;

    // fprintf(screen, "Checking for npartner and partner\n");
    if(!npartner && !partner)
    {
        create_bonds = create_bond_partners(npartner, partner);
    }

    int ncreated = 0;

    // fprintf(screen, "Checking create_bonds\n");
    if(create_bonds)
    {
        // create bonds
        int **_bond_type = atom->bond_type;
        int **_bond_atom = atom->bond_atom;
        int *num_bond = atom->num_bond;
        int newton_bond = force->newton_bond;
        int n_bondhist = atom->n_bondhist;
        double ***bond_hist = atom->bond_hist;

        // fprintf(screen, "Loaded variables... Starting for loop\n");
        // fprintf(screen, "local_start = %i\n", local_start);

        // Note: atoms are created with dummy tag = 0, but
        //       actual tags must be available at this point,
        //       i.e. atom->tag_extend() must have been called
        for(int i = 0; i < nspheres; ++i)
        {
            // fprintf(screen, "\nChecking npartner == 0\n");
            if (npartner[i] == 0) continue;
            // fprintf(screen, "Looping through partners\n-----------------\n");
            for(int k = 0; k < npartner[i]; ++k)
            {
                // fprintf(screen, "Setting j and checking if i is less than j\n");
                const int j = partner[i][k];
                if (!newton_bond || i < j)
                {
                    const int ilocal = local_start + i;

                    if (num_bond[ilocal] == atom->bond_per_atom)
                    {
                        error->one(FLERR,"New bond exceeded bonds per atom in fix bond/create");
                    }

                    // fprintf(screen, "Setting bond_type\n");
                    // if (num_bond != NULL) fprintf(screen, "num_bond != NULL\n");
                    // if (_bond_type != NULL) fprintf(screen, "bond_type != NULL\n");
                    // if (_bond_atom != NULL) fprintf(screen, "bond_atom != NULL\n");
                    // fprintf(screen, "bond_type = %i\n", bond_type);
                    _bond_type[ilocal][num_bond[ilocal]] = bond_type;
                    // fprintf(screen, "Setting tag id\n");
                    // fprintf(screen, "local_start+j = %i\n", local_start+j);
                    // fprintf(screen, "atom->tag[%i] = %i\n",local_start+j, atom->tag[local_start+j]);
                    _bond_atom[ilocal][num_bond[ilocal]] = atom->tag[local_start+j];

                    // reset history
                    for (int ih = 0; ih < n_bondhist; ++ih)
                    {
                        bond_hist[ilocal][num_bond[ilocal]][ih] = 0.;
                    }
                    num_bond[ilocal]++;
                }

                if(i < j)
                    ++ncreated;
            }
        }
    }

    // fprintf(screen, "Destroy partners\n");
    if(create_bonds != 1) // create_bond_partners allocated memory
    {
        destroy_bond_partners(npartner, partner);
    }

    return ncreated;
}

/* ---------------------------------------------------------------------- */
//NP checks against xnear list
//NP returns # inerted spheres
//NP if >1, increase nbodies

/* ---------------------------------------------------------------------- */

/*
int ParticleToInsert::check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, double **xnear, int &nnear)
{
    if(nparticles > 1)
        return check_near_set_x_v_omega_ms(x,v, omega,quat,xnear,nnear);

    // check sphere against all others in xnear
    // if no overlap add to xnear
    double del[3], rsq, radsum;

    vectorCopy3D(x,x_ins[0]);

    for(int i = 0; i < nnear; i++)
    {
        vectorSubtract3D(x_ins[0],xnear[i],del);
        rsq = vectorMag3DSquared(del);
        
/*
        radsum = radius_ins[0] + xnear[i][3];

        // no success in overlap
        if (rsq <= radsum*radsum) return 0;
    }

    // no overlap with any other - success

    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    // add to xnear
    vectorCopy3D(x_ins[0],xnear[nnear]);
    xnear[nnear][3] = radius_ins[0];
    nnear++;

    return 1;
}*/

/* ---------------------------------------------------------------------- */

/*
int ParticleToInsert::check_near_set_x_v_omega_ms(double *x,double *v, double *omega, double *quat, double **xnear, int &nnear)
{
    // x is position where insertion should take place
    // v and omega are the velocity and omega for the newly inserted particles
    double rel[3],xins_j_try[3];
    double del[3], rsq, radsum;

    // check insertion position, take quat into account
    // relative position of spheres to each other already stored at this point
    // check sphere against all others in xnear
    for(int j = 0; j < nspheres; j++)
    {
        // take orientation into account; x_bound_ins is in the global coordinate system
        // calculate xins_j_try for every sphere and check if would work
        vectorSubtract3D(x_ins[j],x_bound_ins,rel);
        MathExtraLiggghts::vec_quat_rotate(rel,quat);
        vectorAdd3D(rel,x,xins_j_try);

        for(int i = 0; i < nnear; i++)
        {
           vectorSubtract3D(xins_j_try,xnear[i],del);
           rsq = vectorMag3DSquared(del);
           radsum = radius_ins[j] + xnear[i][3];

           // no success in overlap
           if (rsq <= radsum*radsum)
            return 0;
        }
    }

    // no overlap with any other - success
    // set x_ins, v_ins and omega_ins
    for(int j = 0; j < nspheres; j++)
    {
        vectorSubtract3D(x_ins[j],x_bound_ins,rel);
        MathExtraLiggghts::vec_quat_rotate(rel,quat);
        vectorAdd3D(rel,x,x_ins[j]);
    }
    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    // add to xnear for future checks
    for(int j = 0; j < nspheres; j++)
    {
        vectorCopy3D(x_ins[j],xnear[nnear]);
        xnear[nnear][3] = radius_ins[j];
        nnear++;
    }

    return nparticles;
}*/

/* ---------------------------------------------------------------------- */

int ParticleToInsert::check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, RegionNeighborList<interpolate_no> & neighList)
{
    if(nparticles > 1)
        return check_near_set_x_v_omega_ms(x,v, omega,quat,neighList);

    vectorCopy3D(x,x_ins[0]);

    if(neighList.hasOverlap(x_ins[0], radius_ins[0])) {
        return 0;
    }

    // no overlap with any other - success

    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    neighList.insert(x_ins[0], radius_ins[0]);

    return 1;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, RegionNeighborList<interpolate_no> & neighList, int reg_id) // Added by Matt Schramm
{
  if(nparticles > 1) {
    return check_near_set_x_v_omega_ms(x,v, omega,quat,neighList,reg_id);
  }

  vectorCopy3D(x,x_ins[0]);

  if(neighList.hasOverlap(x_ins[0], radius_ins[0])) {
    return 0;
  }

  // no overlap with any other - success

  vectorCopy3D(v,v_ins);
  vectorCopy3D(omega,omega_ins);

  neighList.insert(x_ins[0], radius_ins[0]);

  return 1;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::check_near_set_x_v_omega_ms(double *x,double *v, double *omega, double *quat, RegionNeighborList<interpolate_no> & neighList)
{
    // x is position where insertion should take place
    // v and omega are the velocity and omega for the newly inserted particles
    double rel[3],xins_j_try[3];
    //double del[3], rsq, radsum;

    // check insertion position, take quat into account
    // relative position of spheres to each other already stored at this point
    // check sphere against all others in xnear
    for(int j = 0; j < nparticles; j++)
    {
        // take orientation into account; x_bound_ins is in the global coordinate system
        // calculate xins_j_try for every sphere and check if would work
        vectorSubtract3D(x_ins[j],x_bound_ins,rel);
        MathExtraLiggghts::vec_quat_rotate(rel,quat);
        vectorAdd3D(rel,x,xins_j_try);

        if(neighList.hasOverlap(xins_j_try, radius_ins[j])) {
            return 0;
        }
    }

    // no overlap with any other - success
    // set x_ins, v_ins and omega_ins
    for(int j = 0; j < nparticles; j++)
    {
        vectorSubtract3D(x_ins[j],x_bound_ins,rel);
        MathExtraLiggghts::vec_quat_rotate(rel,quat);
        vectorAdd3D(rel,x,x_ins[j]);
    }
    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    // add to xnear for future checks
    for(int j = 0; j < nparticles; j++)
    {
        neighList.insert(x_ins[j], radius_ins[j]);
    }

    return nparticles;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::check_near_set_x_v_omega_ms(double *x,double *v, double *omega, double *quat, RegionNeighborList<interpolate_no> & neighList, int reg_id) // Added by Matt Schramm
{
  // x is position where insertion should take place
  // v and omega are the velocity and omega for the newly inserted particles
  double rel[3],xins_j_try[3];
  //double del[3], rsq, radsum;

  // check insertion position, take quat into account
  // relative position of spheres to each other already stored at this point
  // check sphere against all others in xnear
  for(int j = 0; j < nparticles; j++) {
    // take orientation into account; x_bound_ins is in the global coordinate system
    // calculate xins_j_try for every sphere and check if would work
    vectorSubtract3D(x_ins[j],x_bound_ins,rel);
    MathExtraLiggghts::vec_quat_rotate(rel,quat);
    vectorAdd3D(rel,x,xins_j_try);

    if(neighList.hasOverlap(xins_j_try, radius_ins[j])) {
//            fprintf(screen,"Failed to insert particle...overlaps\n");
        return 0;
    }
  }
  
  // Now check that all particles are inside the factory region???
  for(int j = 0; j < nparticles; j++) {
    vectorSubtract3D(x_ins[j],x_bound_ins,rel);
    MathExtraLiggghts::vec_quat_rotate(rel,quat);
    vectorAdd3D(rel,x,xins_j_try);
//        fprintf(screen,"x = %e, y = %e, z = %e\n",xins_j_try[0],xins_j_try[1],xins_j_try[2]);
    bool in_region = domain->regions[reg_id]->match(xins_j_try[0],xins_j_try[1],xins_j_try[2]);
//        fprintf(screen,"In region == %d\n",in_region);
    if(!in_region)
    {
        fprintf(screen,"Failed to insert particle...extends from region\n");
        return 0;
    }
  }



  // no overlap with any other - success
  // set x_ins, v_ins and omega_ins
  for(int j = 0; j < nparticles; j++) {
    vectorSubtract3D(x_ins[j],x_bound_ins,rel);
    MathExtraLiggghts::vec_quat_rotate(rel,quat);
    vectorAdd3D(rel,x,x_ins[j]);
  }
  vectorCopy3D(v,v_ins);
  vectorCopy3D(omega,omega_ins);

  // add to xnear for future checks
  for(int j = 0; j < nparticles; j++) {
    neighList.insert(x_ins[j], radius_ins[j]);
  }

  return nparticles;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::set_x_v_omega(double *x, double *v, double *omega, double *quat)
{
    double rel[3];

    // x is position where insertion should take place
    // v and omega are the velocity and omega for the newly inserted particles

    // add insertion position
    // relative position of spheres to each other already stored at this point
    // also take quat into account
    for(int j = 0; j < nparticles; j++)
    {
        // if only one sphere, then x_bound = x_ins and there is
        // no relevant orientation
        if(1 == nparticles)
            vectorAdd3D(x_ins[j],x,x_ins[j]);
        // if > 1 sphere, take orientation into account
        // x_bound_ins is in the global coordinate system
        else
        {
            vectorSubtract3D(x_ins[j],x_bound_ins,rel);
            MathExtraLiggghts::vec_quat_rotate(rel,quat);
            vectorAdd3D(rel,x,x_ins[j]);
        }
    }

    // set velocity and omega
    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    return nparticles;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::check_near_set_x_v_omega_dense(double *x,double *v, double *omega, double *quat, RegionNeighborList<interpolate_no> & neighList, int reg_id) // Added by Matt Schramm
{
    // x is position where insertion should take place
    // v and omega are the velocity and omega for the newly inserted particles
    double rel[3],xins_j_try[3];
    //double del[3], rsq, radsum;

    // check insertion position, take quat into account
    // relative position of spheres to each other already stored at this point
    // check sphere against all others in xnear
    //fprintf(screen,"Check if spheres touch\n");
    for(int j = 0; j < nparticles; j++)
    {
        // take orientation into account; x_bound_ins is in the global coordinate system
        // calculate xins_j_try for every sphere and check if would work
        vectorSubtract3D(x_ins[j],x_bound_ins,rel);
        MathExtraLiggghts::vec_quat_rotate(rel,quat);
        vectorAdd3D(rel,x,xins_j_try);

        //fprintf(screen,"xins_j_try == %f, radius_ins_j == %f\n",xins_j_try,radius_ins[j]);

        if(neighList.hasOverlap(xins_j_try, radius_ins[j])) {
            //fprintf(screen,"Failed to insert particle...overlaps\n");
            return 0;
        }
    }
    
    // Now check that all particles are inside the factory region???
    //fprintf(screen,"Check if fiber is in region\n");
    for(int j = 0; j < nparticles; j++)
    {
        vectorSubtract3D(x_ins[j],x_bound_ins,rel);
        MathExtraLiggghts::vec_quat_rotate(rel,quat);
        vectorAdd3D(rel,x,xins_j_try);
//        fprintf(screen,"x = %e, y = %e, z = %e\n",xins_j_try[0],xins_j_try[1],xins_j_try[2]);
        bool in_region = domain->regions[reg_id]->match(xins_j_try[0],xins_j_try[1],xins_j_try[2]);
//        fprintf(screen,"In region == %d\n",in_region);
        if(!in_region)
        {
            //fprintf(screen,"Failed to insert particle...extends from region\n");
            return 0;
        }
    }



    // no overlap with any other - success
    // set x_ins, v_ins and omega_ins
    //fprintf(screen,"So far so good, copy temp_vel to vel\n");
    for(int j = 0; j < nparticles; j++)
    {
        vectorSubtract3D(x_ins[j],x_bound_ins,rel);
        MathExtraLiggghts::vec_quat_rotate(rel,quat);
        vectorAdd3D(rel,x,x_ins[j]);
    }
    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    // add to xnear for future checks
    for(int j = 0; j < nparticles; j++)
    {
        //fprintf(screen,"j=%i, xins==%f, radius_ins==%f\n",j,x_ins[j],radius_ins[j]);
        neighList.insert(x_ins[j], radius_ins[j]);
    }

    return nparticles;
}

/* ---------------------------------------------------------------------- */

void ParticleToInsert::scale_pti(double r_scale)
{
    double r_scale3 = r_scale*r_scale*r_scale;

    for(int i = 0; i < nparticles; i++)
    {
        radius_ins[i] *= r_scale;
        vectorScalarMult3D(x_ins[i],r_scale);
    }

    volume_ins *= r_scale3;
    mass_ins *= r_scale3;

    r_bound_ins *= r_scale;

    vectorScalarMult3D(x_bound_ins,r_scale);
}
