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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "surf_collide_evap_ref_part.h"
#include "surf.h"
#include "surf_react.h"
#include "input.h"
#include "variable.h"
#include "particle.h"
#include "domain.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "math_const.h"
#include "math_extra.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{NUMERIC,CUSTOM,VARIABLE,VAREQUAL,VARSURF};   // surf_collide classes

/* ---------------------------------------------------------------------- */

SurfCollideEvapRefPart::SurfCollideEvapRefPart(SPARTA *sparta, int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal surf_collide evapref command");

  parse_tsurf(arg[2]);

  acc = input->numeric(FLERR,arg[3]);
  if (acc < 0.0 || acc > 1.0)
    error->all(FLERR,"Illegal surf_collide diffuse command");
  if (narg == 4){
    direction = -1;
    liqmin = 0.0;
    liqmax = 0.0;
  } else if (narg == 7){
    direction = atoi(arg[4]);
    liqmin = atof(arg[5]);
    liqmax = atof(arg[6]);
    if (direction < 0 || direction > 2) {
      error->all(FLERR,"Illegal direction in surf_collide evapref command");
    }
    if (liqmin > liqmax) {
      error->all(FLERR,"liqmin cannot be greater than liqmax in surf_collide evapref command");
    }
    if (liqmin == liqmax) {
      error->all(FLERR,"liqmin cannot be equal to liqmax in surf_collide evapref command");
    }
  } else {
    error->all(FLERR,"Illegal surf_collide evapref command");
  }

  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);
  allowreact = 0;
}

/* ---------------------------------------------------------------------- */

SurfCollideEvapRefPart::~SurfCollideEvapRefPart()
{
  if (copy) return;

  delete random;
}

/* ---------------------------------------------------------------------- */

void SurfCollideEvapRefPart::init()
{
  SurfCollide::init();
  check_tsurf();
}

/* ----------------------------------------------------------------------
   particle collision with surface with optional chemistry
   ip = particle with current x = collision pt, current v = incident v
   isurf = index of surface element
   norm = surface normal unit vector
   isr = index of reaction model if >= 0, -1 for no chemistry
   ip = reset to NULL if destroyed by chemistry
   return jp = new particle if created by chemistry
   return reaction = index of reaction (1 to N) that took place, 0 = no reaction
   resets particle(s) to post-collision outward velocity
------------------------------------------------------------------------- */

Particle::OnePart *SurfCollideEvapRefPart::
collide(Particle::OnePart *&ip, double &,
        int isurf, double *norm, int isr, int &reaction)
{
  nsingle++;
  double *x = ip->x; // particle position
  
  if ( direction >= 0 && (x[direction] < liqmin || x[direction] > liqmax )){
    MathExtra::reflect3(ip->v,norm);
  } else if (random->uniform() < acc) {
    ip = NULL;
  } else {
    double *v = ip->v; // particle velocity
    double tangent1[3],tangent2[3];
    Particle::Species *species = particle->species;
    int ispecies = ip->ispecies;

    double vrm = sqrt(2.0*update->boltz * tsurf / species[ispecies].mass);
    double vperp = vrm * sqrt(-log(random->uniform()));

    double theta = MY_2PI * random->uniform();
    double vtangent = vrm * sqrt(-log(random->uniform()));
    double vtan1 = vtangent * sin(theta);
    double vtan2 = vtangent * cos(theta);

    double dot = MathExtra::dot3(v,norm);
    double beta_un,normalized_distbn_fn;

    tangent1[0] = v[0] - dot*norm[0];
    tangent1[1] = v[1] - dot*norm[1];
    tangent1[2] = v[2] - dot*norm[2];

    if (MathExtra::lensq3(tangent1) == 0.0) {
      tangent2[0] = random->uniform();
      tangent2[1] = random->uniform();
      tangent2[2] = random->uniform();
      MathExtra::cross3(norm,tangent2,tangent1);
    }

    MathExtra::norm3(tangent1);
    MathExtra::cross3(norm,tangent1,tangent2);
    v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0];
    v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1];
    v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2];

    // initialize rot/vib energy
    ip->erot = particle->erot(ispecies,tsurf,random);
    ip->evib = particle->evib(ispecies,tsurf,random);
  }
  return NULL;
}