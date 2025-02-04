/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifdef SURF_COLLIDE_CLASS

SurfCollideStyle(evapref,SurfCollideEvapRef)

#else

#ifndef SPARTA_SURF_COLLIDE_EVAP_REF_H
#define SPARTA_SURF_COLLIDE_EVAP_REF_H

#include "surf_collide.h"
#include "surf.h"
#include "particle.h"

namespace SPARTA_NS {

class SurfCollideEvapRef : public SurfCollide {
 public:
  SurfCollideEvapRef(class SPARTA *, int, char **);
  SurfCollideEvapRef(class SPARTA *sparta) : SurfCollide(sparta) {} // needed Kokkos
  virtual ~SurfCollideEvapRef();
  virtual void init();
  Particle::OnePart *collide(Particle::OnePart *&, double &,
                             int, double *, int, int &);

 protected:
  double acc;                // surface accomodation coeff
  class RanKnuth *random;     // RNG for particle reflection

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Surf_collide diffuse rotation invalid for 2d

Specified rotation vector must be in z-direction.

E: Surf_collide diffuse variable name does not exist

Self-explanatory.

E: Surf_collide diffuse variable is invalid style

It must be an equal-style variable.

*/
