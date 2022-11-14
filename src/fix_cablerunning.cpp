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

#include <sstream>
#include <iterator>
#include <algorithm>

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include "fix_cablerunning.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "force.h"
#include "domain.h"
#include "region.h"
#include "fix_property_global.h"
#include "modify.h"

#include "property_registry.h"
#include "global_properties.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MODEL_PARAMS;

/* ---------------------------------------------------------------------- */

FixCableRunning::FixCableRunning(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg),
    Y_cable(1.0e7),
    nu_cable(0.3),
    coeffRest_cable_conduit(0.03),
    coeffFriction_cable_conduit(0.03),
    Y_conduit(1.0e7),
    nu_conduit(0.3),
    conduitRadius(1.0),
    airStartTime(0.0),
    airVelocity(0.0),
    airDynamicViscosity(1.8e-7),
    airDensity(1.225),
    timeSinceStart(0.0),
    pushSpeed(0.0)
{
    double scaling = 1.0;
    if (narg < 7)
        error->all(FLERR,"Illegal fix cable running command\nMust provide file name and conduitRadius");

    bool hasFile = false;
    bool hasConduitRadius = false;
    int iarg = 3;
    while (iarg < narg) {
        if (0 == strcmp("file", arg[iarg])) {
            if (narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'file'");
            hasFile = true;
            const int n = strlen(arg[iarg+1]) + 1;
            this->pointsFileName = new char[n];
            strcpy(this->pointsFileName, arg[iarg+1]);
            iarg += 2;
        } else if (0 == strcmp("scaling", arg[iarg])) {
            if (narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'scaling'");
            scaling = force->numeric(FLERR, arg[iarg+1]);
            iarg += 2;
        } else if (0 == strcmp("pushSpeed", arg[iarg])) {
            if (narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'pushSpeed'");
            this->pushSpeed = force->numeric(FLERR, arg[iarg+1]);
            iarg += 2;
        } else if (0 == strcmp("conduitRadius", arg[iarg])) {
            if (narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'conduitRadius'");
            hasConduitRadius = true;
            this->conduitRadius = force->numeric(FLERR, arg[iarg+1]);
            iarg += 2;
        } else if (0 == strcmp("cableYoungs", arg[iarg])) {
            if (narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'cableYoungs'");
            this->Y_cable = force->numeric(FLERR, arg[iarg+1]);
            iarg += 2;
        } else if (0 == strcmp("conduitYoungs", arg[iarg])) {
            if (narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'conduitYoungs'");
            this->Y_conduit = force->numeric(FLERR, arg[iarg+1]);
            iarg += 2;
        } else if (0 == strcmp("cablePoissons", arg[iarg])) {
            if (narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'cablePoissons'");
            this->nu_cable = force->numeric(FLERR, arg[iarg+1]);
            iarg += 2;
        } else if (0 == strcmp("conduitPoissons", arg[iarg])) {
            if (narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'conduitPoissons'");
            this->nu_conduit = force->numeric(FLERR, arg[iarg+1]);
            iarg += 2;
        } else if (0 == strcmp("coeffFrctionCableConduit", arg[iarg])) {
            if (narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'coeffFrctionCableConduit'");
            this->coeffFriction_cable_conduit = force->numeric(FLERR, arg[iarg+1]);
            iarg += 2;
        } else if (0 == strcmp("coeffRestCableConduit", arg[iarg])) {
            if (narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'coeffRestCableConduit'");
            this->coeffRest_cable_conduit = force->numeric(FLERR, arg[iarg+1]);
            iarg += 2;
        } else if (0 == strcmp("airStartTime", arg[iarg])) {
            if (narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'airStartTime'");
            this->airStartTime = force->numeric(FLERR, arg[iarg+1]);
            iarg += 2;
        } else if (0 == strcmp("airVelocity", arg[iarg])) {
            if (narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'airVelocity'");
            this->airVelocity = force->numeric(FLERR, arg[iarg+1]);
            iarg += 2;
        } else if (0 == strcmp("airDynamicViscosity", arg[iarg])) {
            if (narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'airDynamicViscosity'");
            this->airDynamicViscosity = force->numeric(FLERR, arg[iarg+1]);
            iarg += 2;
        } else if (0 == strcmp("airDensity", arg[iarg])) {
            if (narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'airDensity'");
            this->airDensity = force->numeric(FLERR, arg[iarg+1]);
            iarg += 2;
        } else {
            error->fix_error(FLERR, this, "Incorrect options...");
        }
    }

    if (!hasFile) error->fix_error(FLERR,this,"Missing file name...");
    if (!hasConduitRadius) error->fix_error(FLERR,this,"Missing conduit radius...");

    // Need to load in the points and build the parameterization curves
    std::ifstream pointsFile (this->pointsFileName);
    if (!pointsFile.is_open())
        error->all(FLERR,"Failed to open file");

    std::string curLine;
    std::string space_delimiter = " ";
    while (pointsFile) {
        std::getline(pointsFile, curLine);

        if (curLine.length() < 3)
            continue;

        // Need the next three values
        size_t pos = 0;
        double * xyz = new double[3];
        std::vector<std::string> values{};
        pos = curLine.find(space_delimiter);
        xyz[0] = scaling*force->numeric(FLERR, curLine.substr(0, pos).c_str());
        curLine.erase(0, pos + space_delimiter.length());

        pos = curLine.find(space_delimiter);
        xyz[1] = scaling*force->numeric(FLERR, curLine.substr(0, pos).c_str());
        curLine.erase(0, pos + space_delimiter.length());

        xyz[2] = scaling*force->numeric(FLERR, curLine.substr(0, pos).c_str());

        // https://www.delftstack.com/howto/cpp/cpp-split-string-by-space/

        this->points.push_back(xyz);
    }

}

/* ---------------------------------------------------------------------- */

FixCableRunning::~FixCableRunning()
{
    delete [] this->pointsFileName;
    for (double * point : this->points) {
        delete [] point;
    }
}

/* ---------------------------------------------------------------------- */

int FixCableRunning::setmask()
{
    int mask = 0;
    mask |= POST_FORCE;
    mask |= POST_FORCE_RESPA;
    mask |= MIN_POST_FORCE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixCableRunning::init()
{
    if (strstr(update->integrate_style,"respa"))
        nlevels_respa = ((Respa *) update->integrate)->nlevels;

    // this->properties = atom->get_properties();
    // int max_type = this->properties->max_type();

    // this->Y = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0,style));
    // this->nu = static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0,style));
    // const double * ff = this->Y->get_values();

    // if(!this->Y || !this->nu)
    //     error->all(FLERR,"Fix cablerunning only works with a pair style that defines youngsModulus and poissonsRatio");
}

/* ---------------------------------------------------------------------- */

void FixCableRunning::setup(int vflag)
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

void FixCableRunning::min_setup(int vflag)
{
    post_force(vflag);
}

/* ---------------------------------------------------------------------- */
#define M_6_PI 18.8495559215387594307758602
#define M_8_PI 25.1327412287183459077011471
void FixCableRunning::post_force(int vflag)
{
    // apply drag force to atoms in group
    // direction is opposed to velocity vector
    // magnitude depends on atom type
    const double * const * const x = atom->x;
    const double * const * const omega = atom->omega;

    double **v = atom->v;
    double **f = atom->f;
    double **T = atom->torque;

    const double * const radius = atom->radius;
    const int * const mask = atom->mask;
    const int nlocal = atom->nlocal;

    const double sqrtFiveOverSix = 0.91287092917527685576161630466800355658790782499663875;
    const double Yi = this->Y_conduit;
    const double Yj = this->Y_cable;
    const double vi = this->nu_conduit;
    const double vj = this->nu_cable;
    const double Yeff = 1./((1.-pow(vi,2.))/Yi+(1.-pow(vj,2.))/Yj);
    const double Geff = 1./(2.*(2.-vi)*(1.+vi)/Yi+2.*(2.-vj)*(1.+vj)/Yj);
    const double coeffRestLog = log(this->coeffRest_cable_conduit);
    const double betaeff = coeffRestLog / sqrt(pow(coeffRestLog,2.)+pow(M_PI,2.));

    int oldTag = 0;
    int k_old = 1;
    for (int i = 0; i < nlocal; i++) {
        int k_start = 1;
        if (mask[i] & groupbit) {

            if (atom->tag[i] > oldTag) {
                oldTag = atom->tag[i];
                k_start = k_old;
            }

            // Find the section p1 and p2 that contains the current point
            double p0[3] = {0.0};
            double p1[3] = {0.0};
            double p2[3] = {0.0};
            bool foundPoints = false;
            double x1 = 0.0;
            double x2 = 0.0;
            double y1 = 0.0;
            double y2 = 0.0;
            double z1 = 0.0;
            double z2 = 0.0;
            for (int k = k_start; k < this->points.size(); k++) {
                p1[0] = this->points[k-1][0];
                p1[1] = this->points[k-1][1];
                p1[2] = this->points[k-1][2];

                p2[0] = this->points[k-0][0];
                p2[1] = this->points[k-0][1];
                p2[2] = this->points[k-0][2];

                x1 = MIN(p1[0], p2[0]);
                x2 = MAX(p1[0], p2[0]);

                y1 = MIN(p1[1], p2[1]);
                y2 = MAX(p1[1], p2[1]);

                z1 = MIN(p1[2], p2[2]);
                z2 = MAX(p1[2], p2[2]);

                if (x[i][0] < x1-this->conduitRadius) continue;
                if (x2+this->conduitRadius < x[i][0]) continue;
                if (x[i][1] < y1-this->conduitRadius) continue;
                if (y2+this->conduitRadius < x[i][1]) continue;
                if (x[i][2] < z1-this->conduitRadius) continue;
                if (z2+this->conduitRadius < x[i][2]) continue;

                foundPoints = true;
                k_old = k;
                break;
            }

            if (!foundPoints)
                error->all(FLERR,"Unable to find particle between two points...");

            // We have the two points
            // Need to calculate the distance from the line from p1 to p2, to the particle
            // d^2 = (|(p2-p1)x(p1-p0)|)^2 / (|p2-p1|)^2
            const double d21[3] = {
                p2[0]-p1[0],
                p2[1]-p1[1],
                p2[2]-p1[2]
            };

            const double d10[3] = {
                p1[0]-x[i][0],
                p1[1]-x[i][1],
                p1[2]-x[i][2]
            };

            const double d21_cross_d10[3] = {
                d21[1]*d10[2] - d21[2]*d10[1],
                d21[2]*d10[0] - d21[0]*d10[2],
                d21[0]*d10[1] - d21[1]*d10[0]
            };

            const double mag_d21_cross_d10_squared = (
                d21_cross_d10[0]*d21_cross_d10[0] +
                d21_cross_d10[1]*d21_cross_d10[1] +
                d21_cross_d10[2]*d21_cross_d10[2]
            );

            const double mag_d21 = (
                d21[0]*d21[0] +
                d21[1]*d21[1] +
                d21[2]*d21[2]
            );

            const double d2 = mag_d21_cross_d10_squared/mag_d21;
            const double d = std::sqrt(d2);

            // Check if this is a PUSHED particle (need to check if atom is in first section...)
            if (x[i][1] < 0.0) {
                // Zero out force and torque
                f[i][0] = 0.0;
                f[i][1] = 0.0;
                f[i][2] = 0.0;
                T[i][0] = 0.0;
                T[i][1] = 0.0;
                T[i][2] = 0.0;

                // Set atom velocity to pushed velocity
                const double pushedVelocity = this->pushSpeed;
                const double inv_sqrt_mag_d21 = 1.0/sqrt(mag_d21);
                const double pushed_ex = d21[0]*inv_sqrt_mag_d21;
                const double pushed_ey = d21[1]*inv_sqrt_mag_d21;
                const double pushed_ez = d21[2]*inv_sqrt_mag_d21;
                v[i][0] = pushed_ex*pushedVelocity;
                v[i][1] = pushed_ey*pushedVelocity;
                v[i][2] = pushed_ez*pushedVelocity;

                continue;
            }

            // Apply wind drag force now...
            // FORCES ----------------------------------------------------------
            if (this->timeSinceStart >= this->airStartTime) {
                // Calculate b and c coefficients for forces
                // Equations from "Classical Mechanics" by John R. Taylor
                // Linear term b = 3*pi*nu*D
                const double air_b = M_6_PI * this->airDynamicViscosity * radius[i]; // 6 * pi * nu * r
                // Quadratic term c = k*p*A
                const double air_c = M_PI_4 * this->airDensity *radius[i]*radius[i]; // 0.25 * rho * pi*r*r

                const double inv_sqrt_mag_d21 = 1.0/sqrt(mag_d21);
                const double air_ex = d21[0]*inv_sqrt_mag_d21;
                const double air_ey = d21[1]*inv_sqrt_mag_d21;
                const double air_ez = d21[2]*inv_sqrt_mag_d21;

                const double air_wx = this->airVelocity*air_ex;
                const double air_wy = this->airVelocity*air_ey;
                const double air_wz = this->airVelocity*air_ez;

                const double air_vrx = air_wx - v[i][0];
                const double air_vry = air_wy - v[i][1];
                const double air_vrz = air_wz - v[i][2];
                const double air_vel_norm = sqrt(air_vrx*air_vrx + air_vry*air_vry + air_vrz*air_vrz);

                const double temp = air_b + air_c*air_vel_norm;
                f[i][0] += air_vrx*temp;
                f[i][1] += air_vry*temp;
                f[i][2] += air_vrz*temp;

                // TORQUES ---------------------------------------------------------
                // Equations from "Viscous torque on a sphere under arbitrary rotation" by U. Lei et all
                const double airTorque_d = M_8_PI*this->airDynamicViscosity*radius[i]*radius[i]*radius[i]; // 8 * pi * nu * r*r*r
                T[i][0] -= airTorque_d*omega[i][0];
                T[i][1] -= airTorque_d*omega[i][1];
                T[i][2] -= airTorque_d*omega[i][2];
            }

            // Check if the sphere is touching the cable
            // const double atomCableOverlap = d+radius[i] - this->conduitRadius;
            if (d+radius[i] - this->conduitRadius < 0.0)
                continue;

            // (x1-x0) * (x2-x1) / mag_d21
            const double t = (d10[0]*d21[0] + d10[1]*d21[1] + d10[2]*d21[2]) / mag_d21;
            p0[0] = p1[0] + t*(p1[0]-p2[0]);
            p0[1] = p1[1] + t*(p1[1]-p2[1]);
            p0[2] = p1[2] + t*(p1[2]-p2[2]);

            double dx = x[i][0] - p0[0];
            double dy = x[i][1] - p0[1];
            double dz = x[i][2] - p0[2];
            double magd = sqrt(dx*dx + dy*dy + dz*dz);
            double invmagd = 1.0/magd;
            double ex = invmagd*dx;
            double ey = invmagd*dy;
            double ez = invmagd*dz;

            // Move p0 to the surface of the conduit
            p0[0] += this->conduitRadius*ex;
            p0[1] += this->conduitRadius*ey;
            p0[2] += this->conduitRadius*ez;

            // Calculate new dx
            dx = x[i][0] - p0[0];
            dy = x[i][1] - p0[1];
            dz = x[i][2] - p0[2];
            magd = sqrt(dx*dx + dy*dy + dz*dz);
            invmagd = 1.0/magd;
            const double atomCableOverlap = radius[i] > magd ? radius[i] - magd : 0.0;

            ex = invmagd*dx;
            ey = invmagd*dy;
            ez = invmagd*dz;

            const double vr1 = v[i][0]; // - this->conduit_velocity_x;
            const double vr2 = v[i][1]; // - this->conduit_velocity_x;
            const double vr3 = v[i][2]; // - this->conduit_velocity_x;
            const double vn = vr1*ex + vr2*ey + vr3*ez;
            const double vn1 = ex*vn;
            const double vn2 = ey*vn;
            const double vn3 = ez*vn;
            const double vt1 = vr1-vn1;
            const double vt2 = vr2-vn2;
            const double vt3 = vr3-vn3;

            // atomCableOverlap
            const double meff = 4.0 * M_PI * radius[i]*radius[i]*radius[i]*atom->density[i] / 3.0;
            const double sqrtval = sqrt(radius[i]*atomCableOverlap);
            const double Sn=2.*Yeff*sqrtval;
            const double St=8.*Geff*sqrtval;

            double kn=4./3.*Yeff*sqrtval;
            double kt=St;
            const double gamman=-2.*sqrtFiveOverSix*betaeff*sqrt(Sn*meff);
            const double gammat=-2.*sqrtFiveOverSix*betaeff*sqrt(St*meff);

            const double Fn_contact = kn*atomCableOverlap;
            const double Fn_damping = -gamman*vn;
            const double Fn = Fn_contact + Fn_damping;

            f[i][0] += (ex*Fn);
            f[i][1] += (ey*Fn);
            f[i][2] += (ez*Fn);


            // Apply tangent force update
            const double cr = radius[i] - 0.5*atomCableOverlap;
            const double rinv = 1.0/(radius[i] - atomCableOverlap);
            const double rvx = cr * omega[i][0] * rinv;
            const double rvy = cr * omega[i][1] * rinv;
            const double rvz = cr * omega[i][2] * rinv;

            const double vtr1 = vt1 - (dz*rvy - dy*rvz);
            const double vtr2 = vt2 - (dx*rvz - dz*rvx);
            const double vtr3 = vt3 - (dy*rvx - dx*rvy);

            const double oldShearX = 0.0;
            const double oldShearY = 0.0;
            const double oldShearZ = 0.0;

            double shearX = oldShearX + vtr1*update->dt; // Dont have access to oldshear unless we add it to history values...
            double shearY = oldShearY + vtr2*update->dt;
            double shearZ = oldShearZ + vtr3*update->dt;

            const double rsht = shearX*ex + shearY*ey + shearZ*ez;
            shearX -= rsht * ex;
            shearY -= rsht * ey;
            shearZ -= rsht * ez;

            const double shearmag = std::sqrt(shearX*shearX + shearY*shearY + shearZ*shearZ);

            double ftx = -kt*shearX;
            double fty = -kt*shearY;
            double ftz = -kt*shearZ;

            const double ft_shear = kt*shearmag;
            const double ft_friction = this->coeffFriction_cable_conduit*std::fabs(Fn);

            if (ft_shear > ft_friction) {
                if (shearmag != 0.0) {
                    const double ratio = ft_friction / ft_shear;

                    ftx *= ratio;
                    fty *= ratio;
                    ftz *= ratio;

                    shearX = -ftx/kt;
                    shearY = -fty/kt;
                    shearZ = -ftz/kt;
                } else {
                    ftx = 0.0;
                    fty = 0.0;
                    ftz = 0.0;
                }
            } else {
                ftx -= gammat*vtr1;
                fty -= gammat*vtr2;
                ftz -= gammat*vtr3;
            }

            const double torquex = ey*ftz - ez*fty;
            const double torquey = ez*ftx - ex*ftz;
            const double torquez = ex*fty - ey*ftx;

            const double torque1x = -cr*torquex;
            const double torque1y = -cr*torquey;
            const double torque1z = -cr*torquez;

            f[i][0] += ftx;
            f[i][1] += fty;
            f[i][2] += ftz;

            T[i][0] += torque1x;
            T[i][1] += torque1y;
            T[i][2] += torque1z;
        }
    }

    this->timeSinceStart += update->dt;
}

/* ---------------------------------------------------------------------- */

void FixCableRunning::post_force_respa(int vflag, int ilevel, int iloop)
{
    if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixCableRunning::min_post_force(int vflag)
{
    post_force(vflag);
}
