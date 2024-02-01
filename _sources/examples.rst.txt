Examples
========

The best way to learn the code is to work through the examples in the directory ``examples``. The current set of examples illustrates the range of the functionality of UppASD, and gives a feel for the data analysis that is necessary in order to interpret the results of the simulations. The examples are grouped into the subdirectories ``Mappings``, ``PhaseDiagrams``, ``SimpleSystems``, ``SpecialFeatures``, ``SpinLattice``, and ``SpinWaves``.

bcc Fe
------

The ``SimpleSystems/bccFe`` directory contains the necessary input in order to simulate bcc Fe. There are a couple of Matlab scripts enclosed that plot quantities such as the average magnetization as a function of time. By varying input parameters such as the temperature, one may get a feel for how these affect these observable quantities.

The ``PhaseDiagrams/bccFe-Tsweep`` directory contains very similar input to the ``SimpleSystems/bccFe`` directory. However, a small wrapper is included in order to run SD simulations as a function of temperature. This is often useful for checking the transition temperature of a given system.

Heisenberg Spin Chain
---------------------

The ``SimpleSystems/HeisChain`` directory contains possibly the simplest system that can be set up with UppASD: a Heisenberg ferromagnetic chain. This example is set up to allow the spin wave spectrum of the system to be probed. Since the dispersion of this system can be derived analytically (as discussed any many textbooks concerning solid state physics), it is a handy example to get acquainted with spin correlations and dynamics in UppASD.

By copying over the ``sqw.simid.out`` and ``inpsd.dat`` files into the ``Sqw`` directory, and running the ``sqw_map.sh`` script, :math:`S(\mathbf{q},\omega)` as a function of :math:`\mathbf{q}` is plotted, giving rise to the spin wave dispersion. The dynamics of spin waves have been studied in this way in the paper by Bergman *et al* [Bergman2010]_ .

Two-dimensional Systems
-----------------------

The ``SimpleSystems/fcc001`` directory contains input files for simulating the surface of the Co/Cu(001) monolayer system.

.. For both the ``fcc001`` and ``sc`` examples, a directory ``Snapshots`` is included, that allows snapshots of the lattices to be taken. The scripts are written in VTK, and require ``coord.simid.out`` and ``moments.simid.out`` files to work.

fcc Co
------

Another very common structure in solid state physics is the fcc structure. The ``SimpleSystems/fccCo`` example allows the measurement of the dyamical structure factor, :math:`S(\mathbf{q},\omega)`, for fcc Co, as in the Heisenberg spin chain example.

GaMnAs
------

The ``GaMnAs`` directory contains input data that sets up a GaMnAs dilute magnetic semiconductor [Hellsvik2008]_ .

.. Random Alloy
.. ------------
..
.. The ``Ralloy`` directory contains input data necessary to set up a dilute random alloy.

SKKR Input
----------

The ``xctensortest`` directory contains the same data as in the bccFe directory, but set up in the tensorial format that arises from the Vienna-Budapest SKKR code [Udvardi2003]_ .
