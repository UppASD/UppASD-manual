inpsd.dat keywords: simulation and phases
=========================================

General simulation parameters
-----------------------------

.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------------------------------------------------+
|  do_ralloy    |    Flag to set if a random alloy is being simulated (*0=off*/1=on).                                    |
+---------------+--------------------------------------------------------------------------------------------------------+
|  aunits       |    Implement atomic units, *i.e.* set :math:`k_B`, :math:`\hbar`, ... :math:`=1` (Y/*N*). If this      |
|               |    is switched on, the ``timestep`` in SD mode should be roughly :math:`0.1J_{ij}`.                    |
+---------------+--------------------------------------------------------------------------------------------------------+
|  sdealgh      |    Switch for choosing SDE solver (*1=Midpoint*, 4=Heun, 5=Depondt-Mertens). The default option        |
|               |    runs the semi-implicit midpoint solver developed by Mentink *et al.* [Mentink2010]_.                |
|               |    In this case, as when using the Depondt-Mertens solver [Depondt2009]_, the ``timestep``             |
|               |    can be as large as :math:`10^{-16}` seconds, but this should *always* be checked carefully          |
+---------------+--------------------------------------------------------------------------------------------------------+
|  mensemble    |    Number of ensembles to simulate. The default value is 1, but this may be increased to improve       |
|               |    statistics, especially if investigating laterally confined systems, such as finite                  |
|               |    clusters or other low-dimensional systems.                                                          |
+---------------+--------------------------------------------------------------------------------------------------------+
|  tseed        |    Random number seed for the stochastic field simulating the fluctuations due to temperature.         |
|               |    Default value is 1.                                                                                 |
+---------------+--------------------------------------------------------------------------------------------------------+
|  do_sortcoup  |    Flag to specify if the arrays of couplings should be sorted or not (*Y=yes*, N=no). Of              |
|               |    importance for sampling of polarization. Could be very slow if long range interactions.             |
+---------------+--------------------------------------------------------------------------------------------------------+


Initialization parameters
-------------------------

.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------------------------------------------------+
|  initmag      |    Switch for setting up the initial configuration of the magnetic moments (1=Random distribution,     |
|               |    2=Cone, 3=aligned along direction defined in momfile, *4=Read from restartfile*).                   |
+---------------+--------------------------------------------------------------------------------------------------------+
|  restartfile  |    External file containing stored snapshot from previous simulation (used when initmag=4).            |
|               |    The format coincides with the format of the output file ``restart.simid.out``.                      |
+---------------+--------------------------------------------------------------------------------------------------------+
|  mseed        |    Random number seed for magnetic moments if initmag=1. Set to 1 by default.                          |
+---------------+--------------------------------------------------------------------------------------------------------+
|  theta0       |    If ``initmag`` 2, the magnetic moments are randomly distributed in a cone                           |
|               |    prescribed by this angle, and ``phi0``. Set to 0 by default.                                        |
+---------------+--------------------------------------------------------------------------------------------------------+
|  phi0         |    Cone angle for ``initmag`` 2. Set to 0 by default.                                                  |
+---------------+--------------------------------------------------------------------------------------------------------+
|  roteul       |    Perform global rotation of magnetization. Set to 0 by default.                                      |
+---------------+--------------------------------------------------------------------------------------------------------+
|  rotang       |    Euler angles describing the rotation if ``roteul`` 1.                                               |
+---------------+--------------------------------------------------------------------------------------------------------+
|  initexc      |    Perform initial excitation of the spin system (*N=none*}, I=Vacancies,                              |
|               |    R=Two magnon Raman scattering).                                                                     |
+---------------+--------------------------------------------------------------------------------------------------------+
|  initconc     |    Concentration of vacancies or two magnon spin scattering.                                           |
+---------------+--------------------------------------------------------------------------------------------------------+
|  initneigh    |    eighbour index referring to the list of neighbours for Heisenberg exchange. Determines which spins  |
|               |    to swap in two magnon spin scattering.                                                              |
+---------------+--------------------------------------------------------------------------------------------------------+


Initial phase parameters
------------------------

.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------------------------------------------------+
|  ip_mode      |    Mode for initial phase run (S=SD, M=Monte Carlo, H=Heat bath Monte Carlo, *N=none*).                |
+---------------+--------------------------------------------------------------------------------------------------------+
|  ip_temp      |    Temperature for initial phase run if Monte Carlo (``ip_mode`` M or H).                              |
+---------------+--------------------------------------------------------------------------------------------------------+
|  ip_hfield    |    External applied field (in units of Tesla) for initial phase run, given in Cartesian coordinates,   |
|               |    *e.g.* ``hfield   1.0   0.0   0.0``.                                                                |
+---------------+--------------------------------------------------------------------------------------------------------+
|  ip_mcnstep   |    Number of Monte Carlo sweeps (MCS) over the system if ip_mode=M or H.                               |
+---------------+--------------------------------------------------------------------------------------------------------+
|  ip_damping   |    Damping parameter :math:`\alpha` for SD initial phase. Default value is 0.05.                       |
+---------------+--------------------------------------------------------------------------------------------------------+
|  ip_nphase    |    Number of initial phases to be done with SD.                                                        |
+---------------+--------------------------------------------------------------------------------------------------------+

This must be followed by ``ip_nphase`` lines containing number of steps, temperature, timestep and damping for each phase. An example (for an initialization with the temperature decreasing from 300 K to 10 K) can look like::

  ip_nphase 3
  20000 300.0 1.0d-16  0.1
  20000 100.0 1.0d-16  0.1
  30000 010.0 1.0d-16  0.1

.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------------------------------------------------+
|  ip_mcanneal  |     Number of initial phases to be done with MC.                                                       |
+---------------+--------------------------------------------------------------------------------------------------------+

This must be followed by ``ip_mcanneal`` lines containing number of steps and temperature for each phase. An example (for an initialization with the temperature decreasing from 300 K to 10 K) can look like::

  ip_mcanneal 3
  20000 300.0
  20000 100.0 
  30000 010.0 


Measurement phase parameters
----------------------------

.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------------------------------------------------+
|  mode         |    Mode for measurement phase run (*S=SD*, M=Monte Carlo, H=Heat bath Monte Carlo).                    |
+---------------+--------------------------------------------------------------------------------------------------------+
|  temp         |    Temperature for measurement phase.                                                                  |
+---------------+--------------------------------------------------------------------------------------------------------+
|  hfield       |    External applied field (in units of Tesla) for measurement phase.                                   |
+---------------+--------------------------------------------------------------------------------------------------------+
|  mcnstep      |    Number of Monte Carlo sweeps (MCS) over the system if mode=M or H.                                  |
+---------------+--------------------------------------------------------------------------------------------------------+
|  damping      |    Damping parameter :math:`\alpha` for SD measurement phase. Default value is 0.05.                   |
+---------------+--------------------------------------------------------------------------------------------------------+
|  timestep     |    Time step between SD iterations. Unless ``aunits Y``, this should typically be set to a value       |
|               |    between :math:`10^{-17}` and :math:`10^{-15}` seconds, depending on the system and SDE solver.      |
+---------------+--------------------------------------------------------------------------------------------------------+
|  relaxtime    |    Relaxation time in LLG+I equation (if sdealgh=11).                                                  |
+---------------+--------------------------------------------------------------------------------------------------------+
|  set_bpulse   |    Add magnetic field pulse ``0=no``, :math:`1-4` for different shapes)                                |
+---------------+--------------------------------------------------------------------------------------------------------+
