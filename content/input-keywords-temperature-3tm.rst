Three-temperature models
========================

This section documents the implementation of temperature–bath models used to
describe energy flow between electrons, spins, and lattice degrees of freedom.
Two related models are supported:

* Conventional three-temperature model (3TM)
* Heat-conserving three-temperature model (HC-3TM)

Both are used primarily in ultrafast laser-excitation and non-equilibrium
thermalization studies.

Overview
--------

In multireservoir models, the system is divided into three coupled subsystems:

* Electron bath (e)
* Spin bath (s)
* Lattice bath (l)

Each bath is assigned a temperature :math:`T_e, T_s, T_l` and exchanges energy
with the others. These temperatures then control the stochastic terms in the
equations of motion.

Two philosophies exist:

* **3TM**: temperatures are governed by coupled differential equations with
  phenomenological coupling constants.
* **HC-3TM**: spin and lattice temperatures are measured dynamically during
  simulation, and the electron temperature is adjusted to conserve total energy.

Supported modes
----------------

* Conventional 3TM with fixed coupling constants
* Conventional 3TM with temperature-dependent heat capacities
* Optional inclusion of spin–lattice coupling through explicit SLD

Keyword summary
----------------

``do_3tm``
  Activate 3TM: `Y` = full, `E` = electron-only, `N` = off.

``temp_elec_init``
  Initial electron temperature (K).

``temp_spin_init``
  Initial spin temperature (K).

``temp_latt_init``
  Initial lattice temperature (K).

``Gel``
  Electron–lattice coupling (W/m³/K).

``Ges``
  Electron–spin coupling (W/m³/K).

``Gsl``
  Spin–lattice coupling (W/m³/K).

``P_pulse``
  Pulse fluence (J/m²).

``t0_pulse``
  Pulse center time (s).

``sigma_pulse``
  Pulse width (s).

``gamma_Ce``
  Electron Sommerfeld constant (J/m³/K²).

``Cs``
  Spin heat capacity (J/m³/K).

``Cl``
  Lattice heat capacity (J/m³/K).

``G_cool``
  Cooling rate to ambient (W/m³/K).

``print_3tm_step``
  Print interval (number of simulation steps).

Conventional three-temperature model (3TM)
------------------------------------------

In the conventional 3TM, temperatures are governed by coupled differential
equations [Beaurepaire1996]

.. math::

  \begin{aligned}
  C_e(T_e)\frac{dT_e}{dt} &= -G_{el}(T_e-T_l) - G_{es}(T_e-T_s) + P(t) \\
  C_s(T_s)\frac{dT_s}{dt} &= -G_{es}(T_s-T_e) - G_{sl}(T_s-T_l) \\
  C_l(T_l)\frac{dT_l}{dt} &= -G_{el}(T_l-T_e) - G_{sl}(T_l-T_s)
  \end{aligned}

Here:

* :math:`C_e, C_s, C_l` are heat capacities.
* :math:`G_{el}, G_{es}, G_{sl}` are phenomenological coupling constants.
* :math:`P(t)` is the laser source term.

In this approach:

* Spin and lattice temperatures are not measured from dynamics.
* Energy flow is imposed through fitted coupling parameters.
* Results depend sensitively on the chosen :math:`G` values.

This model is useful for:

* Simple ultrafast demagnetization studies
* Benchmarking against earlier literature
* Systems where couplings are well constrained

Heat-conserving three-temperature model (HC-3TM)
-------------------------------------------------

The HC-3TM removes the phenomenological inter-bath coupling constants and
instead enforces energy conservation dynamically [Pankratova2022]_.

The HC-3TM removes the phenomenological inter-bath coupling constants and
instead enforces energy conservation dynamically.

Philosophy
^^^^^^^^^^

* Spin and lattice temperatures are measured directly from simulation:
  
  .. math::

     T_l = \langle E_{kin}^{lattice} \rangle / k_B

     T_s = \langle E_{ex}^{spin} \rangle / k_B

* These temperatures emerge from Langevin dynamics driven by the electron bath.
* The electron temperature is then adjusted to conserve total energy.

Electron temperature evolution:

.. math::

   \Delta T_e(t) =
   -\frac{C_l(T_l)}{C_e(T_e)}T_l(t)
   -\frac{C_s(T_s)}{C_e(T_e)}T_s(t)
   +\frac{W(t)}{C_e(T_e)}

where :math:`W(t)` is the laser source term.

Algorithm
^^^^^^^^^

At each time step:

1. Propagate spin and lattice with stochastic terms set by current :math:`T_e`.
2. Measure :math:`T_s` and :math:`T_l` from instantaneous energies.
3. Update :math:`T_e` so that total energy is conserved.

No explicit :math:`G_{el}, G_{es}, G_{sl}` are used.

Key consequences:

* Heat flow is emergent, not imposed.
* Fewer poorly known parameters.
* Faster initial demagnetization due to direct electron-controlled noise.

Comparison: 3TM vs HC-3TM
-------------------------

+----------------------------+--------------------------+-----------------------------+
| Feature                    | 3TM                      | HC-3TM                      |
+============================+==========================+=============================+
| Energy exchange            | Phenomenological G’s     | Emergent from dynamics      |
+----------------------------+--------------------------+-----------------------------+
| Spin temperature           | Solved from ODEs         | Measured from spin energy   |
+----------------------------+--------------------------+-----------------------------+
| Lattice temperature        | Solved from ODEs         | Measured from kinetic energy|
+----------------------------+--------------------------+-----------------------------+
| Electron temperature       | ODE with couplings       | Adjusted for energy balance |
+----------------------------+--------------------------+-----------------------------+
| Number of fit parameters   | Several (Gel, Ges, Gsl)  | None (besides damping)      |
+----------------------------+--------------------------+-----------------------------+
| Ultrafast demagnetization  | Often too slow           | Correct sub-ps timescale    |
+----------------------------+--------------------------+-----------------------------+

The faster demagnetization in HC-3TM originates from the fact that the stochastic
fields acting on spins are governed directly by the electron temperature, which
rises rapidly after laser excitation.

Laser source term
-----------------

The laser pulse enters through:

.. math::

   W(t) = A \exp\left[-\frac{(t-t_0)^2}{2\sigma^2}\right]

where:

* :math:`A` is related to the absorbed fluence.
* :math:`t_0` is the pulse center.
* :math:`\sigma` is the temporal width.

This affects only the electron bath, in both 3TM and HC-3TM.

Heat capacities
----------------

Available models:

* Electron:
  :math:`C_e = \gamma_e T_e`
* Lattice:
  Debye model
* Spin:
  * Classical
  * Quantum (Bose-Einstein)
  * Mixed

Temperature-dependent capacities are strongly recommended, especially near
:math:`T_C`.

Spin–lattice coupling
---------------------

If explicit spin–lattice dynamics (SLD) is enabled, exchange parameters depend on
atomic displacements. This allows direct microscopic energy exchange between
spin and lattice subsystems.

In HC-3TM:

* This provides an additional physical channel.
* For many systems (e.g. fcc Ni) the effect is small compared to electron-spin
  energy transfer.

Example: Conventional 3TM
--------------------------

::

   do_3tm              Y
   temp_elec_init      300
   temp_spin_init      300
   temp_latt_init      300

   Gel                 8.0e17
   Ges                 6.0e17
   Gsl                 0.3e17

   P_pulse             5.0
   t0_pulse            1.0e-12
   sigma_pulse         1.0e-13

Example: 3TM with temperature-dependent heat capacities
----------------------------------------------------------

::

   do_3tm              Y
   temp_elec_init      300
   temp_spin_init      300
   temp_latt_init      300

   ! Use temperature-dependent heat capacities
   do_cs_temp          Y
   cv_spinfile         spin_cv.dat
   do_cl_temp          Y
   cv_lattfile         latt_cv.dat

   Gel                 8.0e17
   Ges                 6.0e17
   Gsl                 0.3e17

   P_pulse             5.0
   t0_pulse            1.0e-12
   sigma_pulse         1.0e-13

Notes and recommendations
--------------------------

* Use 3TM for ultrafast demagnetization studies with fitted coupling constants.
* Temperature-dependent heat capacities are recommended, especially near :math:`T_C`.
* Treat coupling constants (Gel, Ges, Gsl) as fit parameters to experimental data.
* For simple studies, use constant heat capacities (gamma_Ce, Cs, Cl).
* Monitor output in ``temperature_3tm.<simid>.out`` to verify multi-timescale dynamics.

Related keywords and cross-references
-------------------------------------

- ``temp``: Equilibrium temperature for non-3TM simulations (see :doc:`input-keywords-simulation`)
- ``damping``: Gilbert damping :math:`\alpha` influences spin-lattice relaxation
  timescales (see :doc:`input-keywords-simulation`)
- ``alat``: Lattice parameter, required for volume/density calculations
- ``Mensemble``: Ensemble averaging for thermal fluctuations

References
==========

See the centralized :doc:`references` for full bibliographic entries:

- [Beaurepaire1996]_ - Ultrafast spin dynamics in ferromagnetic nickel
- [Koopmans2010]_ - Diversity of ultrafast laser-induced demagnetization
- [Mentink2012]_ - Ultrafast spin dynamics in multisublattice magnets
- [Rethfeld2002]_ - Ultrafast dynamics of nonequilibrium electrons in metals
- [Atxitia2010]_ - Multiscale modeling of ultrafast element-specific magnetization dynamics
- [Pankratova2022]_ - Heat-conserving three-temperature model for ultrafast demagnetization
