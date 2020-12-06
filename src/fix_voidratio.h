
#ifdef FIX_CLASS

FixStyle(voidratio,FixVoidRatio)

#else

#ifndef LMP_FIX_VOIDRATIO_H
#define LMP_FIX_VOIDRATIO_H

#include "fix.h"

namespace LAMMPS_NS {
class FixVoidRatio : public Fix  {
 public:
  FixVoidRatio(class LAMMPS *, int, char **);
  ~FixVoidRatio();
  int setmask();
  void init();
  void setup(int);
  // void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void post_run(); // end_of_step();
  double compute_vector(int index);
  double getVoidRatio();

 protected:
  bool doEndOfRun;
  double voidratio;
  double atom_volume, atom_volume_local;
  double void_volume, void_volume_local;
  double region_volume,region_volume_local;
  int nlevels_respa, nevery;
  long int nextCheck;
  unsigned long ntry;
  char *idregion;
  int varflag,iregion;
};
}
#endif
#endif