# MiTMoJCo version 1.1

MiTMoJCo (Microscopic Tunneling Model for Josephson Contacts) is C library designed for modeling superconducting
Josephson junctions within the formalism of microscopic tunneling theory. The purpose of the code is to offer 
implementation of a computationally demanding part of this calculation which is evaluation of the superconducting
pair and quasiparticle tunnel currents. 

The source code contains 5 examples of modeling common cases of Josephson contacts:

1. Current-biased SIS junction.
2. Voltage-biased SIS junction under ac drive.
3. Sine-Gordon breather in long Josephson junction.
4. Fluxon in an Annular Josephson junction.
5. Flux Flow Oscillator.

**Installation**

$ git clone https://github.com/drgulevich/mitmojco

Modify CMakeLists.txt if needed (default installation paths are
    /usr/local/lib for the library and /usr/local/include/mitmojco for headers)

From mitmojco directory enter the (empty) build folder. From there, execute cmake, make and install the library files: 

    $ cd build
    $ cmake ..
    $ make
    $ sudo make install

To check the installation by running one of the provided examples. To compile the first example
type

    $ make example-1

Use similar command for the others.

See the doc folder for the User Guide.

**References**

D. R. Gulevich, V. P. Koshelets, and F. V. Kusmartsev, Phys. Rev. B 96, 024215 (2017). https://arxiv.org/abs/1704.03045
