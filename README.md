# MiTMoJCo version 1.2

MiTMoJCo (Microscopic Tunneling Model for Josephson Contacts) 
represents a C library and collection of Python tools designed to assist modeling superconducting Josephson junctions within the formalism of microscopic tunneling theory. The purpose of the C code is to offer implementation of a computationally demanding part of this calculation which is evaluation of the superconducting pair and quasiparticle tunnel currents. 
The tunnel currents calculated by MiTMoJCo are offered to user's disposal to be employed in a specialized ODE/PDE solver or within a finite difference or finite element scheme in a custom C code.
The source code contains a collection of Python tools to work with tunnel current amplitudes which characterize specific superconducting materials constituting a tunnel junction and 
examples of modeling some common cases of Josephson contacts.

**Description of the source files**

``amplitudes/``: folder containing a library of pre-calculated fits of tunnel current amplitudes for common types of Nb and NbN junctions.

``doc/``: folder containing documentation.

``examples/``: folder containing examples of using MiTMoJCo C library.

``figures/``: illustrations to features of the Python ``mitmojco`` module.

``include/``: folder with C headers ``mitmojco/mitmojco.h`` and ``mitmojco/opt_filter.h``.

``src/``: folder with the C source code ``mitmojco.c`` and supplementary optimum filtration routine ``opt_filter.c``.

``CMakeLists.txt``: input to CMake at installation.

``LICENSE``: GNU General Public License information.

``README.md``: general information about the code.

``amplitudes.ipynb``: Jupyter notebook demonstrating the use of ``mitmojco`` Python module.

``module.py``: supplementary Python module for calculation of tunnel current amplitudes (TCAs) from the BCS expressions, smoothing Riedel peaks and creation of custom fits of TCAs by a sum of complex exponentials used in the C library.  

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

You may need to run 

    $ sudo ldconfig
	
to reload the list of system-wide library paths. If the library path is loaded correctly, the following output is expected:

    $ ldconfig -p | grep mitmojco
      libmitmojco.so (libc6,x86-64) => /usr/local/lib/libmitmojco.so

Check the installation by running one of the provided examples. 

**Examples**

The source code contains 6 examples of modeling common cases of Josephson contacts:

1. Current-biased SIS junction.

2. Voltage-biased SIS junction under ac drive.

3. Sine-Gordon breather in long Josephson junction.

4. Fluxon in an Annular Josephson junction.

5. Flux Flow Oscillator.

6. Gap resonance (D. V. Averin, *Gap resonance in the classical dynamics of the current-biased Josephson tunnel junctions*, arXiv:2111.07206; https://arxiv.org/abs/2111.07206).

7. deal.II+MiTMoJCo: 2D model of Josephson junction.

To compile the first example type from the `examples/` folder

    $ make example-1

Use similar command to compile examples 1-5. The last example 6 requires `deal.II` finite element library (https://www.dealii.org/).

**Features of the Python ``mitmojco`` module**

1. Calculation of BCS tunnel current amplitudes.

![Alt text](/figures/NbNbN_4K2.png?raw=true "BCS Tunnel Current Amplitudes")

2. Smoothing Riedel peaks in BCS tunnel current amplitudes.

![Alt text](/figures/NbNbN_4K2_smoothed.png?raw=true "Smoothed BCS Tunnel Current Amplitudes")

3. Fit by a sum of complex exponentials.

![Alt text](/figures/NbNbN_4K2_smoothed_fit.png?raw=true "Smoothed BCS Tunnel Current Amplitudes")

4. Saving, loading and comparison of tunnel current amplitude fits.

![Alt text](/figures/TCA_comparison.png?raw=true "Comparison of different tunnel current amplitude fits")

**References**

1. D. R. Gulevich, *MiTMoJCo: Microscopic Tunneling Model for Josephson Contacts*, Comput. Phys. Commun. 251, 107091 (2020); https://doi.org/10.1016/j.cpc.2019.107091.

2. MiTMoJCo 1.2 User Guide: https://www.researchgate.net/publication/318380494_User_Guide_for_MiTMoJCo_Microscopic_Tunneling_Model_for_Josephson_Contacts.

**Research using MiTMoJCo**

1. D. R. Gulevich, V. P. Koshelets, and F. V. Kusmartsev, *Josephson Flux Flow Oscillator: the Microscopic Tunneling Approach*, Phys. Rev. B 96, 024215 (2017); https://doi.org/10.1103/PhysRevB.96.024515.

2. D. R. Gulevich, V. P. Koshelets, F. V. Kusmartsev, *Bridging the Terahertz gap for chaotic sources with superconducting junctions*, Phys. Rev. B 99, 060501(R) (2019); https://doi.org/10.1103/PhysRevB.99.060501.

3. D. R. Gulevich, L. V. Filippenko, V. P. Koshelets, *Microscopic tunneling model of Nb-AlN-NbN Josephson flux-flow oscillator*, J. Low Temp. Phys. 194, 312 (2019); https://doi.org/10.1007/s10909-018-2106-x.

4. MiTMoJCo tunnel current amplitude ``NbNb_4K2_001.fit`` enters *PSCAN2 Superconductor Circuit Simulator* under the TJM approximation name ``tjm1``; http://www.pscan2sim.org.

5. J. A. Delport, *Simulation and verification software for superconducting electronic circuits*, PhD thesis, Stellenbosch University (2019); https://scholar.sun.ac.za/handle/10019.1/106048.

6. D. R. Gulevich, *MiTMoJCo: Microscopic Tunneling Model for Josephson Contacts*, Comput. Phys. Commun. 251, 107091 (2020); https://doi.org/10.1016/j.cpc.2019.107091.

7. H. G. Ahmad, L. Di Palma, R. Caruso, A. Pal, G. P. Pepe, M. G. Blamire, F. Tafuri, and D. Massarotti, *Critical Current Suppression in Spin-Filter Josephson Junctions*, J. Supercond. Nov. Magn. 33, 3043 (2020); https://doi.org/10.1007/s10948-020-05577-0.

8. Д. Р. Гулевич, *Динамика нелинейных систем с памятью*, Университет ИТМО, Санкт-Петербург (2020) [Translation from Russian: D. R. Gulevich, Dynamics of nonlinear systems with memory, the ITMO University, Saint-Petersburg (2020)]. ISBN 978-5-7577-0642-9.
