Detailed background and inputs
==============================

This chapter consolidates simulation methods, external control modules, and observable analysis into one place. Use it as the deep-dive reference when selecting algorithms, driving simulations with stimuli, and interpreting outputs.

Simulation Methods and Algorithms
---------------------------------

Fundamental and specialized approaches for running magnetic dynamics simulations.

Stochastic Integration
~~~~~~~~~~~~~~~~~~~~~~

Numerical integration schemes (solvers) for solving the stochastic Landau-Lifshitz-Gilbert equation:

- :doc:`input-keywords-solvers`

Monte Carlo simulations and minimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Equilibrium sampling techniques for thermodynamic properties:

- :doc:`input-keywords-montecarlo`

Advanced Sampling and Equilibration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Specialized methods for exploring phase space and computing thermodynamic properties:

- :doc:`input-keywords-montecarlo`
- :doc:`input-keywords-wanglandau`
- :doc:`input-keywords-replicaexchange`

Path-Finding and Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Methods for finding transitions between magnetic configurations and optimizing spin textures:

- :doc:`input-keywords-gneb`
- :doc:`input-keywords-spinspiral`

Rare Event and Dynamics Methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Techniques for studying infrequent events and non-equilibrium dynamics:

- :doc:`input-keywords-kmc`

Spin-Lattice Coupling
~~~~~~~~~~~~~~~~~~~~~

Coupled spin and lattice evolution including magnetoelastic effects:

- :doc:`input-keywords-sld`

Multiscale Coupling
~~~~~~~~~~~~~~~~~~~

Hybrid atomistic–continuum simulations for large-scale magnetic systems:

- :doc:`input-keywords-multiscale`

External Stimuli and Control
----------------------------

Apply external driving forces and thermal control modules.

Spatially and temporally varying magnetic fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- :doc:`input-keywords-fields`

Spin-polarized currents and torques
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- :doc:`input-keywords-currents`

Temperature control and three-temperature model (3TM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- :doc:`input-keywords-temperature-3tm`

Temperature gradients:
~~~~~~~~~~~~~~~~~~~~~~

- :doc:`input-keywords-temperature-gradients`

Observable Measurements and Analysis
------------------------------------

Theoretical foundations for measuring and analyzing magnetic observables and configurations.

Magnetization Observables
~~~~~~~~~~~~~~~~~~~~~~~~~

Measurements of magnetization and thermodynamic moments:

- :doc:`input-keywords-averages`

Spin-wave Spectra (AMS/LSWT)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Adiabatic magnon spectra and linear spin-wave theory:

- :doc:`input-keywords-ams`

Correlations and Structure Factors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Space and time-displaced correlation functions:

- :doc:`input-keywords-correlations`
- :doc:`input-keywords-autocorrelation`

Polarization and Chirality
~~~~~~~~~~~~~~~~~~~~~~~~~~

Ferroelectric polarization and local chirality:

- :doc:`input-keywords-polarization`

Configuration Analysis and Topology
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Characterization and topological analysis of magnetic configurations:

- :doc:`input-keywords-clusters`
- :doc:`input-keywords-topology`

Micromagnetic Stiffness and Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Extraction of continuum parameters from atomistic simulations:

- :doc:`input-keywords-stiffness`

.. toctree::
   :maxdepth: 2
   :hidden:

   dib-simulation-methods
   dib-external-stimuli
   dib-observables-analysis
