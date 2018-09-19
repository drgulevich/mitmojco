# MiTMoJCo version 1.1

MiTMoJCo (Microscopic Tunneling Model for Josephson Contacts) is C library designed for modeling superconducting
Josephson junctions within the formalism of microscopic tunneling theory. The purpose of the code is to offer 
implementation of a computationally demanding part of this calculation which is evaluation of the superconducting
pair and quasiparticle tunnel currents. 

**Description of source files**

``amplitudes/``: folder containing a library of pre-calculated fits of tunnel current amplitudes for common types of Nb and NbN junctions

``doc/``: folder containing documentation

``examples/``: folder containing examples of using MiTMoJCo C library

``include/``: folder with C headers ``mitmojco/mitmojco.h`` and ``mitmojco/opt_filter.h``

``src/``: folder with the C source code ``mitmojco.c`` and supplementary optimum filtration routine ``opt_filter.c``

``CMakeLists.txt``: input to CMake at installation

``LICENSE``: GNU General Public License information

``README.md``: general information about the code

``amplitudes.ipynb``: Jupyter notebook demonstrating the use of ``mitmojco`` Python module

``module.py``: supplementary Python module for working with tunnel current amplitudes

**Examples**

The source code contains 6 examples of modeling common cases of Josephson contacts:

1. Current-biased SIS junction.

2. Voltage-biased SIS junction under ac drive.

3. Sine-Gordon breather in long Josephson junction.

4. Fluxon in an Annular Josephson junction.

5. Flux Flow Oscillator.

6. deal.II+MiTMoJCo: 2D model of Josephson junction

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

1. D. R. Gulevich, V. P. Koshelets, and F. V. Kusmartsev, *Josephson Flux Flow Oscillator: the Microscopic Tunneling Approach*, Phys. Rev. B 96, 024215 (2017); https://arxiv.org/abs/1704.03045.

2. D. R. Gulevich, *MiTMoJCo: Microscopic Tunneling Model for Josephson Contacts*, https://arxiv.org/abs/1809.04706.

3. MiTMoJCo 1.1 User Guide: https://www.researchgate.net/publication/318380494_User_Guide_for_MiTMoJCo_Microscopic_Tunneling_Model_for_Josephson_Contacts.

