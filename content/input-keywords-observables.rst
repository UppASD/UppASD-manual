inpsd.dat keywords: observables
===============================

Parameters for measuring of observables
---------------------------------------

Typically the measurement of each observable is controlled by two parameters in a combination as follows; ``do_observable`` that enables the measurement and ``observable_step`` that determines the frequency of the measurements. Here the ``observable`` should be replaced by the internal name of the wanted quantity i.e. ``do_avrg`` and ``avrg_step`` for the average magnetization.

.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------------------------------------------------+
|  plotenergy   |    Flag to enable the calculation of the energy of the system projected to the different components of |
|               |    the Hamiltonian. ``{0=off}/1=on``)                                                            .     |
+---------------+--------------------------------------------------------------------------------------------------------+
|  do_avrg      |    Sample and print average magnetization, and its higher order moments. ``Y/N``                       |
+---------------+--------------------------------------------------------------------------------------------------------+
|  avrg_step    |    Number of time steps between sampling of averages. Set to 100 by default.                           |
+---------------+--------------------------------------------------------------------------------------------------------+
|  avrg_buff    |    Number of samplings of averages to buffer between printing to file. Set to 10 by default.           |
+---------------+--------------------------------------------------------------------------------------------------------+
|  do_proj_avrg |    Sample and print type (*i.e*. sublattice) projected average moments. (``Y/N/A``).                   |
+---------------+--------------------------------------------------------------------------------------------------------+
| do_projch_avrg|    Sample and print chemical (*i.e.*} sublattice) projected average moments (``Y/N/A``).               |
+---------------+--------------------------------------------------------------------------------------------------------+
|  do_cumu      |    Sample cumulants (Y/N). Automatically enabled for Monte Carlo simulations.                          |
+---------------+--------------------------------------------------------------------------------------------------------+
|  cumu_step    |    Number of time steps between sampling of cumulants. Set to 25 by default.                           |
+---------------+--------------------------------------------------------------------------------------------------------+
|  do_tottraj   |    Sample and print all trajectories (moments) in the system. (Y/N). Generates the (rather large)      |
|               |    ``moments.simid.out`` file.                                                                         |
+---------------+--------------------------------------------------------------------------------------------------------+
|  tottraj_step |    Number of time steps between samplings of moments. Set to 1000 by default.                          |
+---------------+--------------------------------------------------------------------------------------------------------+
|  tottraj_buff |    Number of samplings of moments to buffer between printing to file. Set to 10 by default.            |
+---------------+--------------------------------------------------------------------------------------------------------+
|  ntraj        |    Number of individual trajectories to sample and print. Followed by ``ntraj`` lines describing atoms |
|               |    to sample, time step between samples and steps to buffer. Set to 0 by default.                      |
+---------------+--------------------------------------------------------------------------------------------------------+
|  do_pol       |    Sample and print average ferroelectric polarization (Y/N) according to the expression               |
|               |    :math:`P\propto \gamma\sum_{i,j}\hat{\mathbf{e}}_{ij}\times(\mathbf{m}_i\times\mathbf{m}_j)`.       |
|               |    Uses the neighbour lists set up for exchange but here the sum is performed up to a threshold        |
|               |    ``max_pol_nn``. For this construction to work, it is important to set the flag ``do_sortcoup N``.   |
+---------------+--------------------------------------------------------------------------------------------------------+
|  max_pol_nn   |    Number of neighbours to use when evaluating the polarization.                                       |
+---------------+--------------------------------------------------------------------------------------------------------+
|  pol_step     |    Number of time steps between sampling of polarization averages. Set to 100 by default.              |
+---------------+--------------------------------------------------------------------------------------------------------+
|  pol_buff     |  Number of samplings of polarization averages to buffer between printing to file. Set to 10 by default.|
+---------------+--------------------------------------------------------------------------------------------------------+
|  do_stiffness |    Calculation of spin-wave stiffness (and tensor) and micromagnetic exchange constant (Y/N).          |
+---------------+--------------------------------------------------------------------------------------------------------+
|  eta_min      |    Lowest value of auxiliary convergence parameter in stiffness calculation (recommended around 6-8)   |
+---------------+--------------------------------------------------------------------------------------------------------+
|  eta_max      |    Largest value of auxiliary convergence parameter in stiffness calculation (recommended around 10-12)|
+---------------+--------------------------------------------------------------------------------------------------------+
|  alat         |    Lattice constant (in m) for calculation of exchange stiffness                                       |
+---------------+--------------------------------------------------------------------------------------------------------+
