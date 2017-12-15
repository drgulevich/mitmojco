# MiTMoJCo version 1.1

MiTMoJCo (Microscopic Tunneling Model for Josephson Contacts) is C library designed for modeling superconducting
Josephson junctions within the formalism of microscopic tunneling theory. The purpose of the code is to offer 
implementation of a computationally demanding part of this calculation which is evaluation of the superconducting
pair and quasiparticle tunnel currents. 

The source code contains 5 examples of modeling common cases of Josephson contacts:

1. Current-biased SIS junction.

2. Voltage-biased SIS junction under ac drive.
![Alt text](/examples/figures/example-2.png?raw=true "Example 2")

3. Sine-Gordon breather in long Josephson junction.
![Alt text](/examples/figures/example-3.png?raw=true "Example 3")

4. Fluxon in an Annular Josephson junction.
![Alt text](/examples/figures/example-4.png?raw=true "Example 4")

5. Flux Flow Oscillator.
![Alt text](/examples/figures/example-5.png?raw=true "Example 5")

**Installation**

    $ git clone https://github.com/drgulevich/mitmojco

Modify CMakeLists.txt to suit your needs. Default installation paths are set in CMakeLists.txt to
    /usr/local/lib for the shared library and /usr/local/include/mitmojco for the headers.

In the mitmojco directory create and enter the build/ folder. Then execute cmake, make and the installation procedure: 

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make
    $ sudo make install

Check the installation by running one of the provided examples. To compile the first example type

    $ make example-1

Use similar command for the others.

See the doc folder for the User Guide.

**References**

1. MiTMoJCo 1.1 User Guide: https://www.researchgate.net/publication/318380494_User_Guide_for_MiTMoJCo_Microscopic_Tunneling_Model_for_Josephson_Contacts.

2. Original paper: D. R. Gulevich, V. P. Koshelets, and F. V. Kusmartsev, Phys. Rev. B 96, 024215 (2017); https://arxiv.org/abs/1704.03045.
