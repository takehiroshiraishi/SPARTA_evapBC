# SPARTA_inletBC

This repository contains the implementation of various surface collision models for the SPARTA (Stochastic PArallel Rarefied-gas Time-accurate Analyzer) simulator. These models are used to simulate the interaction of particles with surfaces under different conditions.

## File Structure

- `surf_collide_evap_ref_part.h` and `surf_collide_evap_ref_part.cpp`: Implementation of the `SurfCollideEvapRefPart` class, which handles particle collisions with surfaces that have evaporation and reflection properties.
- `surf_collide_evap_ref.h` and `surf_collide_evap_ref.cpp`: Implementation of the `SurfCollideEvapRef` class, which handles particle collisions with surfaces that have evaporation properties.
  
## Usage

To use these collision models in SPARTA, include the appropriate header files and instantiate the classes as needed in your simulation setup.

## License

This software is distributed under the GNU General Public License. See the LICENSE file for more details.

## Contact

For questions or issues, please contact Steve Plimpton (sjplimp@gmail.com) or Michael Gallis (magalli@sandia.gov) at Sandia National Laboratories.
