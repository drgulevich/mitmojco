# mitmojco

MiTMoJCo (Microscopic Tunneling Model for Josephson Contacts) is C code designed for modeling superconducting
Josephson junctions within the formalism of microscopic tunneling theory. The purpose of the code is to offer 
implementation of a computationally demanding part of this calculation which is evaluation of the superconducting
pair and quasiparticle tunnel currents. 

The source code contains 5 examples of modeling common cases of Josephson contacts:

1. Current-biased SIS junction.
2. Voltage-biased SIS junction under ac drive.
3. Sine-Gordon breather in long Josephson junction.
4. Fluxon in an Annular Josephson junction.
5. Flux Flow Oscillator.

To compile a particular example, type

$ make example-#

substituting # for the example number, or,

$ make all

to compile all the provided examples.

See the doc folder for the User Guide.

References:
D. R. Gulevich, V. P. Koshelets, and F. V. Kusmartsev, Phys. Rev. B 96, 024215 (2017). https://arxiv.org/abs/1704.03045
