Detailed background and inputs
==============================

This chapter consolidates simulation methods, external control modules, and observable analysis into one place. Use it as the deep-dive reference when selecting algorithms, driving simulations with stimuli, and interpreting outputs.

Simulation Methods and Algorithms
---------------------------------

Fundamental and specialized approaches for running magnetic dynamics simulations.

Stochastic Integration
~~~~~~~~~~~~~~~~~~~~~~

Numerical integration schemes (solvers) for solving the stochastic Landau-Lifshitz-Gilbert equation:

- :doc:`methods/solvers`

Monte Carlo simulations and minimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Equilibrium sampling techniques for thermodynamic properties:

- :doc:`methods/monte-carlo`

Advanced Sampling and Equilibration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Specialized methods for exploring phase space and computing thermodynamic properties:

- :doc:`methods/monte-carlo`
- :doc:`methods/wang-landau`
- :doc:`methods/replica-exchange`

Path-Finding and Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Methods for finding transitions between magnetic configurations and optimizing spin textures:

- :doc:`methods/gneb`
- :doc:`methods/spin-spiral`

Rare Event and Dynamics Methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Techniques for studying infrequent events and non-equilibrium dynamics:

- :doc:`methods/kmc`

Coupled Spin-Lattice Dynamics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Coupled spin and lattice evolution including magnetoelastic effects:

- :doc:`methods/spin-lattice`

Multiscale Simulations
~~~~~~~~~~~~~~~~~~~~~~

Hybrid atomistic–continuum simulations for large-scale magnetic systems:

- :doc:`methods/multiscale`


External Stimuli and Control
----------------------------

Apply external driving forces and thermal control modules.

Spatially and temporally varying magnetic fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- :doc:`stimuli/fields`
- :doc:`stimuli/microwave_fields`

Spin-polarized currents and torques
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- :doc:`stimuli/currents`

Temperature control and three-temperature model (3TM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- :doc:`stimuli/temperature-3tm`

Temperature gradients:
~~~~~~~~~~~~~~~~~~~~~~

- :doc:`stimuli/temperature-gradients`

Observable Measurements and Analysis
------------------------------------

Theoretical foundations for measuring and analyzing magnetic observables and configurations.

Magnetization Observables
~~~~~~~~~~~~~~~~~~~~~~~~~

Measurements of magnetization and thermodynamic moments:

- :doc:`observables/averages`

Spin-wave Spectra (AMS/LSWT)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Adiabatic magnon spectra and linear spin-wave theory:

- :doc:`observables/ams`

Correlations and Structure Factors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Space and time-displaced correlation functions:

- :doc:`observables/correlations`
- :doc:`observables/autocorrelation`

Polarization and Chirality
~~~~~~~~~~~~~~~~~~~~~~~~~~

Ferroelectric polarization and local chirality:

- :doc:`observables/polarization`

Configuration Analysis and Topology
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Characterization and topological analysis of magnetic configurations:

- :doc:`methods/clusters`
- :doc:`observables/topology`

Micromagnetic Stiffness and Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Extraction of continuum parameters from atomistic simulations:

- :doc:`observables/stiffness`

.. toctree::
   :maxdepth: 2
   :hidden:

   methods/index
   stimuli/index
   observables/index
