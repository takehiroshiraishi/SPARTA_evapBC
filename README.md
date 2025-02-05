# SPARTA Boundary Condition for Evaporating Surface

This repository contains the implementation of various surface collision models for the SPARTA (Stochastic PArallel Rarefied-gas Time-accurate Analyzer) simulator. These models are used to simulate the interaction of particles with surfaces under different conditions.

## File Structure

- `surf_collide_evap_ref_part.h` and `surf_collide_evap_ref_part.cpp`: Implementation of the `SurfCollideEvapRefPart` class, which handles particle collisions with surfaces that have evaporation and reflection properties.
- `surf_collide_evap_ref.h` and `surf_collide_evap_ref.cpp`: Implementation of the `SurfCollideEvapRef` class, which handles particle collisions with surfaces that have evaporation properties.
  
## Usage
### evaprefpart

When the whole wall is liquid surface

```cpp
double tsurf; // Surface temperature
double acc; // Condensation coefficient
```
Example: All of ylo is liquid surface at 353.15 K with condensation coefficient of 0.25.
```cpp
surf_collide        inlet evaprefpart 353.15 0.25
bound_modify        ylo collide inlet
```
When the wall is partly liquid surface
```cpp
double tsurf; // Surface temperature
double acc; // Condensation coefficient
int direction; // 0: x, 1: y, 2: 3
double liqmin; // Min value of liquid surface [length]
double liqmax; // Max value of liquid surface [length]
```
Example: 0 < x < 5 um part of ylo is liquid surface at 353.15 K with condensation coefficient of 0.25.
```cpp
variable            coeff equal 0.25 // evaporation, condensation coefficient
variable            nw equal ${coeff}*9.725e24 // effective number density of the liquid surface

species             ../air.species H2O // read water properties
mixture             water H2O vstream 0.0 0 0 temp 353.15 nrho ${nw} // half maxwellian at the liquid surface
region              liquid block 0 5e-6 -1.5e-8 1.5e-8 INF INF // define liquid surface (ylo at 0 < x < 5 um)
fix                 in emit/face water ylo region liquid
surf_collide        inlet evaprefpart 353.15 ${coeff} 0 0.0 5e-6
bound_modify        ylo collide inlet
```
### evaprefpart


## License

This software is distributed under the GNU General Public License. See the LICENSE file for more details.

## Contact

For questions or issues, please contact Steve Plimpton (sjplimp@gmail.com) or Michael Gallis (magalli@sandia.gov) at Sandia National Laboratories.
