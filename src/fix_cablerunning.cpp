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

FixCableRunning::FixCableRunning(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
    double scaling = 0.3048;
    this->Y_cable = 1.0e7;
    this->nu_cable = 0.3;
    this->coeffRest_cable_conduit = 0.3;
    this->Y_conduit = 1.0e7;
    this->nu_conduit = 0.3;
    if (narg < 4)
        error->all(FLERR,"Illegal fix cable running command\nMust provide file name");
    const int n = strlen(arg[3]) + 1;
    this->pointsFileName = new char[n];
    strcpy(this->pointsFileName, arg[3]);
    this->cableRadius = force->numeric(FLERR, arg[4]);


    // if (narg < 5)
    //     error->all(FLERR,"Illegal fix air drag command\nMust have air_viscosity and air_density");

    // air_viscosity = force->numeric(FLERR,arg[3]);
    // air_density   = force->numeric(FLERR,arg[4]);

    // iregion = -1;
    // idregion = NULL;
    // wx = 0.0;
    // wy = 0.0;
    // wz = 0.0;

    // if (narg > 5) {
    //     if (narg < 8)
    //         error->all(FLERR,"Illegal fix air drag command\nMissing wx, wy, and wz"); // region
    //     wx =  force->numeric(FLERR,arg[5]);
    //     wy =  force->numeric(FLERR,arg[6]);
    //     wz =  force->numeric(FLERR,arg[7]);
    // }

    // if (narg > 8) {
    //     if (narg > 9) error->all(FLERR,"Illegal fix air drag command\nRegion variable name error");

    //     iregion = domain->find_region(arg[8]);
    //     if (iregion == -1) error->all(FLERR,"Region ID for fix setforce does not exist");
    //     const int n = strlen(arg[8]) + 1;
    //     idregion = new char[n];
    //     strcpy(idregion,arg[8]);
    // }

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
// Needed constants
void FixCableRunning::post_force(int vflag)
{
    // apply drag force to atoms in group
    // direction is opposed to velocity vector
    // magnitude depends on atom type
    double * const * const x = atom->x;
    double * const * const v = atom->v;
    double * const * const omega = atom->omega;

    double **f = atom->f;
    double **T = atom->torque;

    double * const radius = atom->radius;
    int * const mask = atom->mask;
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

    for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {

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
            for (int k = 1; k < this->points.size(); k++) {
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

                if (x[i][0] < x1-this->cableRadius) continue;
                if (x2+this->cableRadius < x[i][0]) continue;
                if (x[i][1] < y1-this->cableRadius) continue;
                if (y2+this->cableRadius < x[i][1]) continue;
                if (x[i][2] < z1-this->cableRadius) continue;
                if (z2+this->cableRadius < x[i][2]) continue;

                foundPoints = true;
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
            const double d = sqrt(d2);

            // Apply wind drag force now...

            // Check if the sphere is touching the cable
            // const double atomCableOverlap = d+radius[i] - this->cableRadius;
            if (d+radius[i] - this->cableRadius < 0.0)
                continue;

            // (x1-x0) * (x2-x1) / mag_d21
            const double t = (d10[0]*d21[0] + d10[1]*d21[1] + d10[2]*d21[2]) / mag_d21;
            p0[0] = p1[0] + t*(p1[0]-p2[0]);
            p0[1] = p1[1] + t*(p1[1]-p2[1]);
            p0[2] = p1[2] + t*(p1[2]-p2[2]);

            // For HM, Lets just try harmonic...
            // const double ddxx = x[i][0] - p0[0];
            // const double ddyy = x[i][1] - p0[1];
            // const double ddzz = x[i][2] - p0[2];
            // const double ddmm = sqrt(ddxx*ddxx + ddyy*ddyy + ddzz*ddzz);
            // const double eexx = ddxx/ddmm;
            // const double eeyy = ddyy/ddmm;
            // const double eezz = ddzz/ddmm;

            // p0[0] = p0[0] + eexx*this->cableRadius;
            // p0[1] = p0[1] + eeyy*this->cableRadius;
            // p0[2] = p0[2] + eezz*this->cableRadius;

            // const double dx = x[i][0] - p0[0];
            // const double dy = x[i][1] - p0[1];
            // const double dz = x[i][2] - p0[2];
            // const double dm2 = ddxx*ddxx + ddyy*ddyy + ddzz*ddzz;
            // const double invdm2 = 1./dm2;
            // const double dm = sqrt(dm2);
            // const double invdm = 1./dm;

            // // normal component
            // // double vnnr = vr1 * delx + vr2 * dely + vr3 * delz;
            // // vnnr *= rsqinv;
            // // const double vn1 = delx * vnnr;
            // // const double vn2 = dely * vnnr;
            // // const double vn3 = delz * vnnr;

            // // // tangential component at the center of the sphere
            // // const double vt1 = vr1 - vn1;
            // // const double vt2 = vr2 - vn2;
            // // const double vt3 = vr3 - vn3;

            // const double vr1 = v[i][0];
            // const double vr2 = v[i][1];
            // const double vr3 = v[i][2]; // - this->conduit_velocity_z
            // const double vnnr = (vr1*dx + vr2*dy + vr3*dz)*invdm2;
            // const double vn1 = dx*vnnr;
            // const double vn2 = dy*vnnr;
            // const double vn3 = dz*vnnr;
            // const double vt1 = vr1-vn1;
            // const double vt2 = vr2-vn2;
            // const double vt3 = vr3-vn3;

            // const double Yi = this->Y_conduit;
            // const double Yj = this->Y_cable;
            // const double vi = this->nu_conduit;
            // const double vj = this->nu_cable;
            // const double Yeff = 1./((1.-pow(vi,2.))/Yi+(1.-pow(vj,2.))/Yj);
            // const double Geff = 1./(2.*(2.-vi)*(1.+vi)/Yi+2.*(2.-vj)*(1.+vj)/Yj);
            // const double coeffRestLog = log(this->coeffRest_cable_conduit);
            // const double betaeff = coeffRestLog / sqrt(pow(coeffRestLog,2.)+pow(M_PI,2.));

            // const double meff = 4.0 * M_PI * radius[i]*radius[i]*radius[i]*atom->density[i] / 3.0;
            // const double sqrtval = sqrt(radius[i]*d);
            // const double Sn=2.*Yeff*sqrtval;
            // const double St=8.*Geff*sqrtval;

            // double kn=4./3.*Yeff*sqrtval;
            // double kt=St;
            // const double gamman=-2.*sqrtFiveOverSix*betaeff*sqrt(Sn*meff);
            // const double gammat=-2.*sqrtFiveOverSix*betaeff*sqrt(St*meff);

            // // convert Kn and Kt from pressure units to force/distance^2
            // kn /= force->nktv2p;
            // kt /= force->nktv2p;

            // const double Fn_damping = -gamman*sqrt(vn1*vn1 + vn2*vn2 + vn3*vn3);
            // const double Fn_contact = kn*d;
            // double Fn = Fn_contact + Fn_damping;

            // f[i][0] += (invdm*Fn*dx);
            // f[i][1] += (invdm*Fn*dy);
            // f[i][2] += (invdm*Fn*dz);

            // Harmonic update
            // fprintf(screen,"-----------------\n");
            // fprintf(screen, "P0 = {%e, %e, %e}\n", p0[0], p0[1], p0[2]);
            // fprintf(screen, "P1 = {%e, %e, %e}\n", p1[0], p1[1], p1[2]);
            // fprintf(screen, "P2 = {%e, %e, %e}\n", p2[0], p2[1], p2[2]);

            double dx = x[i][0] - p0[0];
            double dy = x[i][1] - p0[1];
            double dz = x[i][2] - p0[2];
            double magd = sqrt(dx*dx + dy*dy + dz*dz);
            double invmagd = 1.0/magd;
            double ex = invmagd*dx;
            double ey = invmagd*dy;
            double ez = invmagd*dz;

            // Move p0 to the surface of the conduit
            p0[0] += this->cableRadius*ex;
            p0[1] += this->cableRadius*ey;
            p0[2] += this->cableRadius*ez;

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
            // fprintf(
            //     screen,
            //     "------------------------\n" \
            //     "tag = %d, Fn = %e, kn = %e, Yeff = %e, sqrtval = %e\n" \
            //     "radius = %e, atomCableOverlap = %e\n" \
            //     "ex = %e, ey = %e, ez = %e\n" \
            //     "dx = %e, dy = %e, dz = %e\n" \
            //     "p0[0] = %e, p0[1] = %e, p0[2] = %e\n" \
            //     "p1[0] = %e, p1[1] = %e, p1[2] = %e\n" \
            //     "p2[0] = %e, p2[1] = %e, p2[2] = %e\n" \
            //     "x[0] = %e, x[1] = %e, x[2] = %e\n",
            //     atom->tag[i], Fn, kn, Yeff, sqrtval, \
            //     radius[i], atomCableOverlap, \
            //     ex, ey, ez, \
            //     dx, dy, dz, \
            //     p0[0], p0[1], p0[2], \
            //     p1[0], p1[1], p1[2],\
            //     p2[0], p2[1], p2[2], \
            //     x[i][0], x[i][1], x[i][2]
            // );
            f[i][0] += (ex*Fn);
            f[i][1] += (ey*Fn);
            f[i][2] += (ez*Fn);



            //limit force to avoid the artefact of negative repulsion force

            // Apply hertz mindlin contact model...
            // const double poi = Sn[type]/(2.0*St[type]) - 1.0; // 0.25; --> Sn/(2*St) - 1.0 = poi
            // const double one_minus_p2_inv = 1. / (1.0 - poi * poi);
            // const double deltan = r - bondLength;
            // const double reff = 0.5 * rout;
            // const double Yeff = 0.5 * Sn[type] * one_minus_p2_inv; // Need Poisson's ratio here also... set to 0.25...
            // const double sqrtval = sqrt(-reff * deltan);
            // const double kn = 4. / 3. * Yeff * sqrtval;
            // fn = MAX(-kn * deltan * rinv, 0.0);
            // f[i1][0] += fn * delx;
            // f[i1][1] += fn * dely;
            // f[i1][2] += fn * delz;

            // if (iregion >= 0 && !domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
            //     continue;

            // FORCES ----------------------------------------------------------
            // Calculate b and c coefficients for forces
            // Equations from "Classical Mechanics" by John R. Taylor
            // Linear term b = 3*pi*nu*D
            // b = M_6_PI * air_viscosity * radius[i]; // 6 * pi * nu * r
            // Quadratic term c = k*p*A
            // c = M_PI_4 * air_density *radius[i]*radius[i]; // 0.25 * rho * pi*r*r

            // vrx = wx - v[i][0];
            // vry = wy - v[i][1];
            // vrz = wz - v[i][2];
            // vel_norm = sqrt(vrx*vrx + vry*vry + vrz*vrz);

            // temp = b + c*vel_norm;
            // f[i][0] += vrx*temp;
            // f[i][1] += vry*temp;
            // f[i][2] += vrz*temp;

            // TORQUES ---------------------------------------------------------
            // Equations from "Viscous torque on a sphere under arbitrary rotation" by U. Lei et all
            // d = M_8_PI*air_viscosity*radius[i]*radius[i]*radius[i]; // 8 * pi * nu * r*r*r
            // T[i][0] -= d*omega[i][0];
            // T[i][1] -= d*omega[i][1];
            // T[i][2] -= d*omega[i][2];
        }
    }
    // delete [] Y_values;
    // delete [] nu_values;
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
