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

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Richard Berger (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef NORMAL_MODEL
NORMAL_MODEL(JKR,jkr,16)
#else
#ifndef NORMAL_MODEL_JKR_H_
#define NORMAL_MODEL_JKR_H_

#include <iostream>
#include <fstream>
#include <iomanip>

#include "global_properties.h"
#include "fix_property_atom.h"
#include <cmath>
#include "normal_model_base.h"

#include "fix_mesh_surface.h"
#include "jkr_lookup_table.h"
#include "neighbor.h"


namespace MODEL_PARAMS
{
  static ScalarProperty* createResolutionJKR(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
      return createScalarProperty(registry, "resolutionJKR", caller);
  }
  static MatrixProperty* createWorkOfAdhesion(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    return createPerTypePairProperty(registry, "workOfAdhesion", caller);
  }

}

namespace LIGGGHTS {

namespace ContactModels
{
  class ContactModelBase;

  template<>
  class NormalModel<JKR> : public NormalModelBase
  {
  public:
    NormalModel(LAMMPS * lmp, IContactHistorySetup* hsetup, class ContactModelBase *c) :
      NormalModelBase(lmp, hsetup, c),
      Yeff(NULL),
      Geff(NULL),
      betaeff(NULL),
      limitForce(false),
      displayedSettings(false),
      cmb(c),
      history_offset(0)
    {
      history_offset = hsetup->add_history_value("contflag", "0");
      debug = false;
    }

    void registerSettings(Settings & settings)
    {
      settings.registerOnOff("tangential_damping", tangential_damping, true);
      settings.registerOnOff("limitForce", limitForce);
    }

    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb)
    {
    }

    inline double get_max(double ** matrix,int max_type)
    {
      double maximum = 0.;
      for(int x=1; x<max_type+1; ++x)
        for(int y=1; y<max_type+1; ++y)
          maximum = fmax(matrix[x][y], maximum);
      return maximum;
    }
    inline double get_min(double ** matrix,int max_type)
    {
      double minimum = 1e50;
      for(int x=1; x<max_type+1; ++x)
        for(int y=1; y<max_type+1; ++y)
          minimum = fmin(matrix[x][y], minimum);

      return minimum;
    }
    void connectToProperties(PropertyRegistry & registry)
    {
      registry.registerProperty("Yeff", &MODEL_PARAMS::createYeff,"model jkr");
      registry.registerProperty("Geff", &MODEL_PARAMS::createGeff,"model jkr");
      registry.registerProperty("betaeff", &MODEL_PARAMS::createBetaEff,"model jkr");
      registry.registerProperty("workOfAdhesion", &MODEL_PARAMS::createWorkOfAdhesion,  "model jkr");
      registry.registerProperty("resolutionJKR", &MODEL_PARAMS::createResolutionJKR,  "model jkr");

      registry.connect("Yeff", Yeff,"model jkr");
      registry.connect("Geff", Geff,"model jkr");
      registry.connect("betaeff", betaeff,"model jkr");
      registry.connect("workOfAdhesion", workOfAdhesion,"model jkr");
      registry.connect("resolutionJKR", resolutionJKR,"model jkr");
      
      
      //Initalize the Look Up Table (LUT) 
      lut.init(resolutionJKR);

      //Display LUT properties to screen and logfile
      if (comm->me == 0) {
        if (screen)
          fprintf(screen,"JKR table created with:\nResolution = %g \nMaximum(delta_delta_c) = %g\nPoints = %i\n",
                  resolutionJKR,lut.get_max_delta_delta_c(),lut.get_nr_of_points());
        if (logfile)
          fprintf(logfile,"JKR table created with:\nResolution = %g \nMaximum(delta_delta_c) = %g\nPoints = %i\n",
                  resolutionJKR,lut.get_max_delta_delta_c(),lut.get_nr_of_points());
      }
      
      //Find max values for negative overlap threshold
      const double max_rad = registry.max_radius();
      const double min_rad = registry.min_radius();
      const int max_type = registry.max_type();
      double workOfAdhesion_max = get_max(workOfAdhesion,max_type);
      double Yeff_min = get_min(Yeff,max_type);
      double a_0_max = lut.calculate_a0(workOfAdhesion_max,max_rad,Yeff_min);
      double delta_c_max = lut.calculate_delta_c(a_0_max,max_rad);

      double contact_distance_factor = 1 + delta_c_max/min_rad;

      //If no work of adhesion cancel surfaceClose calculation
      if (workOfAdhesion_max == 0.0)
      {
        contact_distance_factor = 1.0;
      }
      neighbor->register_contact_dist_factor(contact_distance_factor);
    }
    
    inline double stressStrainExponent()
    {
      return 1.5;
    }

    inline void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      // store that surfaces intersect for necking later
      double * const contflag = &sidata.contact_history[history_offset];

      contflag[0] = 1.0;

      const int itype = sidata.itype;
      const int jtype = sidata.jtype;
      const double radi = sidata.radi;
      const double radj = sidata.radj;
      double reff = sidata.is_wall ? radi : (radi*radj/(radi+radj));


      #ifdef SUPERQUADRIC_ACTIVE_FLAG
      if(sidata.is_non_spherical && atom->superquadric_flag)
        reff = sidata.reff;
      #endif
      const double meff=sidata.meff;

      if(!displayedSettings)
      {
        displayedSettings = true;
      }

      //// convert Kn and Kt from pressure units to force/distance^2
      // kn /= force->nktv2p;
      // kt /= force->nktv2p;

      double Fn_contact,Fc,delta_c,a_0,a ;

      //If gamma_surf is set to zero then no adhesion is calculated (Hertz solution)
      double a_a0;
      if (workOfAdhesion[sidata.itype][sidata.jtype] == 0.0)
      {
        Fc = 0.0;
        delta_c = 0.0;
        a_0 = 1;
        a = sqrt(reff*sidata.deltan);
        Fn_contact = 4./3.*Yeff[itype][jtype]*a*sidata.deltan;
      }
      else
      {
        //Calculate conact force for JKR model
        a_0 = lut.calculate_a0(workOfAdhesion[sidata.itype][sidata.jtype],reff,Yeff[itype][jtype]);
        delta_c = lut.calculate_delta_c(a_0,reff);
        //Fc = 3.*M_PI*gamma_surf[sidata.itype][sidata.jtype]*reff;
        Fc = lut.calculate_fc(workOfAdhesion[sidata.itype][sidata.jtype],reff);
        Fn_contact = lut.get_fn_fc(sidata.deltan/delta_c)*Fc;

        a_a0 = lut.get_a_a0(sidata.deltan/delta_c);
        a = a_a0*a_0;
        //Updating kt to be correct for JKR surface as in Marshall 2009 eq 28 and kn as eq 21
      }

      const double kt = 8.*Geff[itype][jtype]*a;
      const double kn = 4./3.*Yeff[itype][jtype]*a;

      const double Sn=2.*Yeff[itype][jtype]*sqrt(reff*sidata.deltan);
      const double St=8.*Geff[itype][jtype]*sqrt(reff*sidata.deltan);
      const double sqrtFiveOverSix = 0.91287092917527685576161630466800355658790782499663875;
      const double gamman=-2.*sqrtFiveOverSix*betaeff[itype][jtype]*sqrt(Sn*meff);
      const double gammat= tangential_damping ? -2.*sqrtFiveOverSix*betaeff[itype][jtype]*sqrt(St*meff) : 0.0;
      double Fn_damping = -gamman*sidata.vn;

      double Fn = Fn_contact + Fn_damping;

      //limit force to avoid the artefact of negative repulsion force
      if(limitForce && (Fn<0.0) )
      {
          Fn = 0.0;
      }

      sidata.Fn = Fn;
      sidata.kn = kn;
      sidata.kt = kt;
      sidata.gamman = gamman;
      sidata.gammat = gammat;

      #ifdef NONSPHERICAL_ACTIVE_FLAG
          double torque_i[3] = {0., 0., 0.};
          double Fn_i[3] = { Fn * sidata.en[0], Fn * sidata.en[1], Fn * sidata.en[2]};
          if(sidata.is_non_spherical) {
            double xci[3];
            vectorSubtract3D(sidata.contact_point, atom->x[sidata.i], xci);
            vectorCross3D(xci, Fn_i, torque_i);
          }
          
      #endif

      // apply normal force
      if(sidata.is_wall) 
      {
        const double Fn_ = Fn * sidata.area_ratio;
        i_forces.delta_F[0] += Fn_ * sidata.en[0];
        i_forces.delta_F[1] += Fn_ * sidata.en[1];
        i_forces.delta_F[2] += Fn_ * sidata.en[2];

        #ifdef NONSPHERICAL_ACTIVE_FLAG
                if(sidata.is_non_spherical) {
                  //for non-spherical particles normal force can produce torque!
                  i_forces.delta_torque[0] += torque_i[0];
                  i_forces.delta_torque[1] += torque_i[1];
                  i_forces.delta_torque[2] += torque_i[2];
                }
        #endif
      } 
      else {
        i_forces.delta_F[0] += sidata.Fn * sidata.en[0];
        i_forces.delta_F[1] += sidata.Fn * sidata.en[1];
        i_forces.delta_F[2] += sidata.Fn * sidata.en[2];

        j_forces.delta_F[0] += -i_forces.delta_F[0];
        j_forces.delta_F[1] += -i_forces.delta_F[1];
        j_forces.delta_F[2] += -i_forces.delta_F[2];
    
        #ifdef NONSPHERICAL_ACTIVE_FLAG
          if(sidata.is_non_spherical) 
          {
            //for non-spherical particles normal force can produce torque!
            double xcj[3], torque_j[3];
            double Fn_j[3] = { -Fn_i[0], -Fn_i[1], -Fn_i[2]};
            vectorSubtract3D(sidata.contact_point, atom->x[sidata.j], xcj);
            vectorCross3D(xcj, Fn_j, torque_j);

            i_forces.delta_torque[0] += torque_i[0];
            i_forces.delta_torque[1] += torque_i[1];
            i_forces.delta_torque[2] += torque_i[2];

            j_forces.delta_torque[0] += torque_j[0];
            j_forces.delta_torque[1] += torque_j[1];
            j_forces.delta_torque[2] += torque_j[2];
          }
        #endif
      }
    }

    void surfacesClose(SurfacesCloseData &scdata, ForceData& i_forces, ForceData& j_forces)
    {

        if (scdata.contact_flags) *scdata.contact_flags |= CONTACT_NORMAL_MODEL;


        double * const contflag = &scdata.contact_history[history_offset];
        const double dx = scdata.delta[0];
        const double dy = scdata.delta[1];
        const double dz = scdata.delta[2];

        //If no surface energy between particles then contact breaks at delta=0
        if (workOfAdhesion[scdata.itype][scdata.jtype] == 0.0)
        {
          contflag[0] = 0.0;
        }

        //Checking if surface is close after contact
        if(MathExtraLiggghts::compDouble(contflag[0],1.0,1e-6))
        {
          const int itype = scdata.itype;
          const int jtype = scdata.jtype;

          const double radi = scdata.radi;
          const double radj = scdata.is_wall ? radi : scdata.radj;
          double reff = scdata.is_wall ? radi : (radi*radj/(radi+radj));
          double a_0 = lut.calculate_a0(workOfAdhesion[itype][jtype],reff,Yeff[itype][jtype]);
          double delta_c = lut.calculate_delta_c(a_0,reff);
          const double r = sqrt(scdata.rsq);
          
          //In this case deltan == - particle distance
          const double deltan =  scdata.is_wall ?  radi - r : (radi + radj) - r;

          //Necking stage of contact
          if (deltan > -delta_c )
          {

            double Fc = lut.calculate_fc(workOfAdhesion[itype][jtype],reff);
            double Fn = lut.get_fn_fc(deltan/delta_c)*Fc;

            const double rinv = 1.0 / r;

            const double enx = dx * rinv;
            const double eny = dy * rinv;
            const double enz = dz * rinv;

            // apply normal force
            scdata.has_force_update = true;
            if(scdata.is_wall)
            {
              const double Fn_ = Fn * scdata.area_ratio;
              i_forces.delta_F[0] += Fn_ * enx;
              i_forces.delta_F[1] += Fn_ * eny;
              i_forces.delta_F[2] += Fn_ * enz;
            }
            else {
              i_forces.delta_F[0] += Fn * enx;
              i_forces.delta_F[1] += Fn * eny;
              i_forces.delta_F[2] += Fn * enz;

              j_forces.delta_F[0] += -i_forces.delta_F[0];
              j_forces.delta_F[1] += -i_forces.delta_F[1];
              j_forces.delta_F[2] += -i_forces.delta_F[2];
            }
          }
          //Contact breaks
          else
          {
            contflag[0] = 0.0;
          }

        }
    }

    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

  protected:
    double ** Yeff;
    double ** Geff;
    double ** betaeff;
    bool limitForce;
    bool displayedSettings;
    class ContactModelBase *cmb;
    int history_offset;
    double ** workOfAdhesion;
    bool tangential_damping;
    LUT lut;
    double resolutionJKR;
    bool debug;
  };

}

}
#endif
#endif
