Introduction
============


Theoretical overview
--------------------

This user guide describes the essential features of the Uppsala Atomistic Spin Dynamics program (UppASD). The emphasis is on the input files necessary to run calculations using UppASD and the output files it generates. Some information concerning the analysis of data generated by the code is also given. The ASD method and our implementation of it is described in the article by Skubic *et al* [Skubic2008]_. The underlying theory is described in the articles by Antropov *et al.* [Antropov1996]_ and García-Palacios and Lázaro [Garcia-Palacios1998]_. An old, yet remarkably lucid overview is also given by Watson *et al.* [Watson1969]_. A thorough review of the method and relevant applications is given in the book by Eriksson *et al.* [Eriksson2017]_.

The program evolves the equations of motion for atomic magnetic moments in a solid. These take the form of the Landau-Lifshitz-Gilbert (LLG) equation

.. math::
   
   \frac{d\mathbf{m}_i}{dt}=-\frac{\gamma}{1+\alpha^2} \mathbf{m}_i \times [\mathbf{B}_{i}+\mathbf{b}_{i}(t)]-\frac{\gamma}{m_i} \frac{\alpha}{1+\alpha^2} \mathbf{m}_i \times (\mathbf{m}_i \times [\mathbf{B}_{i}+\mathbf{b}_{i}(t)]).

In this expression, :math:`\gamma` is the gyromagnetic ratio and :math:`\mathbf{b}_{i}(t)` is a stochastic magnetic field with a Gaussian distribution, the magnitude of which is related to the damping parameter :math:`\alpha`, which eventually brings the system into thermal equilibrium. The typical time step for solving the differential equations is :math:`\Delta t=0.1` femtoseconds, *i.e.* :math:`10^{-16}` seconds.

The effective field, :math:`\mathbf{B}_i`, experienced by each atom :math:`\textit{i}` is given by the partial derivative of the Hamiltonian :math:`\mathcal{H}` with respect to the local magnetic moment :math:`\mathbf{m}_i`

.. .. _effectivefield:

.. math::
  \mathbf{B}_i=-\frac{ \partial \mathcal{H} }{ \partial \mathbf{m}_i },

where :math:`\mathcal{H}` is the spin Hamiltonian that takes all relevant interactions into account. Note that in the input files to UppASD, the convention is that it is the unit vectors :math:`\mathbf{e}_i=\frac{\mathbf{m}_i}{m_i}` that enter the spin Hamiltonian

.. .. _spinHamiltonian:
.. Eqn. :ref:`test <effectivefield>`: .

.. math::   
   \mathcal{H}=-\sum_{i,j} J_{ij}\mathbf{e}_i \cdot \mathbf{e}_j - \sum_{i,j} \mathbf{D}_{ij}\cdot(\mathbf{e}_i \times \mathbf{e}_j)-\sum_i K_i (\hat{\mathbf{e}}_i \cdot \mathbf{e}_i^K)^2-\sum_i \mathbf{B}^{ext}\cdot\mathbf{e}_i  + \ldots

This is consistent with most electronic structure codes that output :math:`J_{ij}` and other exchange interactions but care should be taken since in many other models found in the litterature, the Hamiltonian depends explicitly of :math:`\mathbf{m}_i` instead. The most important contribution to the Hamiltonian is typically given by the Heisenberg exchange Hamiltonian, given by the first term in the spin Hamiltonian. There, :math:`i` and :math:`j` are atomic indices, and :math:`J_{ij}` is the strength of the exchange interaction. These exchange interactions can be obtained from first-principles calculations, or alternatively be inferred from experiments. It is also possible (and in many cases even essential) to include other terms to the Hamiltonian, including Dzyaloshinskii-Moriya exchange, magnetic anisotropies and external magnetic fields, as are also exemplified in the spin Hamiltonian. There are currently several other interactions available in the UppASD code and additional interactions can be implemented quite straightforwardly. Please note that the format of the Hamiltonian can be defined differently regarding prefactors, inclusion of moment magnitudes, summation convention and more. The input format in UppASD is conformal with most commonly used electronic structure codes that have the capability of calculating :math:`J_{ij}` (and sometimes :math:`\mathbf{D}_{ij}`), *i.e.* following the same convention of the Heisenberg Hamiltonian as Liechtenstein *et al.* [Lichtenstein1987]_.


License
-------

The UppASD code is developed by the Division of Materials Theory, in the Department of Physics and Astronomy at Uppsala University, Sweden. The copyright of the code is held by the developers but the program is open for use and distribution according to the GPLv3 license.

Further information concerning the license and contact information of the developers may be found on the UppASD webpage https://www.physics.uu.se/UppASD.

.. The current version of the code (5.0) is still under active development.


Installation
------------

The source code is distributed on https://github.com/UppASD/UppASD along with documentation and a growing set of examples. To install, perform the following actions

  - Obtain the code, by downloading and unpacking a release::

      wget https://github.com/UppASD/UppASD/archive/refs/tags/v6.0.2.tar.gz
      tar xvzf v6.0.2.tar.gz
      cd UppASD-6.0.2

    or by cloning the git repository::

      git clone https://github.com/UppASD/UppASD.git
      cd UppASD

  - Generate the dependencies needed for compiling the code::

      make deps

    (Optional) Perform a system check for available compiler profiles::

      make probe

    Compile the code with the selected compiler profile::

      make <profile>

    where ``<profile>`` is the name of the profile, i.e. ``ifort``, ``ifort-cuda``, ``gfortran``,
    ``gfortran-osx``, and so on, e.g. ``make ifort``.
    
  - Test the compiled program against a selection of realistic runs::

      make asd-tests

In addition to the source files, the UppASD distribution also contains several examples (in the directory ``examples/``), documentation,

.. including this file (in  ``docs/``)

and routines and reference data (``tests/``) for validating the installation of the UppASD program.


Principles of the Code
----------------------

When run, UppASD essentially goes through three stages:

- Initialization: the system is set up.
- Initial phase: an optional stage in which the system is brought into thermal equilibrium, with limited data sampling.
- Measurement phase: the system is evolved in time, with complete data sampling being made.

During the initialization phase, all the parameters necessary to describe the system of interest, such as its geometry, dimensions, exchange couplings and boundary conditions, are set up. In addition, the initial phase also sets the simulation parameters, such as the number of simulation steps to record data over, which stochastic differential equation (SDE) integrator to use, and the temperature at which the simulation should be run.

The initial phase, which is optional, is typically performed in order to bring the system into thermal equilibrium, so that the data recorded in the measurement phase is for a thermalized system. Obviously, if one is interested in out-of-equilibrium dynamics, then there is no need to perform this phase. The initial phase can either be performed using Spin Dynamics (SD), or the Metropolis or Heatbath Monte Carlo (MC) algorithms [Binder2009]_. The latter is convenient for ground state searches, provided the system is not too complex, for instance a system with a spin glass phase.

During the measurement phase, the data sampling is performed. Simulations can be run in either MC or SD mode. In MC mode only magnetization averages and static correlation functions may be measured. In SD mode, a much richer set of observables are measured, including the dynamical structure factor.
