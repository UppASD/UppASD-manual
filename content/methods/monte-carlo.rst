Monte Carlo simulations
=======================

Parameters for Monte Carlo simulations and thermal equilibration
----------------------------------------------------------------

The UppASD code supports multiple Monte Carlo (MC) algorithms for thermal
equilibration and equilibrium phase studies of atomistic magnetic systems. The
implementation includes [Landau2014]_, [Binney1992]_:

- :ref:`Metropolis algorithm for general Heisenberg and Ising models <mc-metropolis>`
- :ref:`Heat bath (Glauert) algorithm for improved rejection-free sampling <mc-heatbath>`
- :ref:`Glauber dynamics for non-equilibrium relaxation studies <mc-glauber>`
- :ref:`Spin-ice and loop algorithms for constrained systems <mc-spinice>`
- :ref:`Local spin fluctuation (LSF) and induced moment coupling <mc-lsf-induced>`

All algorithms operate through single-spin or multi-spin updates evaluated with
the full interatomic Hamiltonian, including exchange, anisotropy, dipolar, and
relativistic interactions. Monte Carlo is essential for accurate thermodynamic
measurements and phase transition characterization [Binder1981]_.

.. important::

   The code automatically enables cumulant measurements (``do_cumu Y``) for all
   Monte Carlo simulations, as equilibrium statistical mechanics requires
   comprehensive sampling of magnetization fluctuations for accurate thermodynamic
   properties and critical behavior analysis.

-------------------------------------------------
Canonical Monte Carlo acceptance criterion
-------------------------------------------------

For standard single-spin flip Monte Carlo, the **Metropolis** acceptance
probability is [Metropolis1953]_:

.. math::

   P_{\mathrm{acc}} = \min(1, \exp(-\beta \Delta E)),

where:

- :math:`\beta = 1/(k_B T)` is the inverse temperature
- :math:`\Delta E = E_{\mathrm{trial}} - E_{\mathrm{current}}` is the energy difference
- :math:`k_B` is Boltzmann's constant

For the **Heat bath** (or Glauert) algorithm, the acceptance probability depends
on the total effective field [Glauert1974]_:

.. math::

   P_{\mathrm{acc}} = \frac{1}{1 + \exp(\beta \Delta E)},

which eliminates rejections and provides better thermal statistics, particularly
near critical temperatures.

The **Glauber dynamics** algorithm uses:

.. math::

   P_{\mathrm{acc}} = \frac{1}{1 + \exp(\beta \Delta E)},

identical to heat bath but historically derived from Brownian dynamics [Glauber1963]_.

.. _mc-metropolis:

-------------------------------------------------
Metropolis algorithm
-------------------------------------------------

The Metropolis algorithm generates trial moment configurations by random rotation
according to the Hinzke-Nowak update scheme [Hinzke2000]_:

Mode M
   Select this mode in ``mode`` parameter to enable Metropolis algorithm. One random
   spin per atom is attempted each MC sweep. Trial configurations are generated as a
   mixture of uniform random rotations (probability 1/3), Gaussian rotations around
   current direction (1/3), and spin flips (1/3). After evaluation via Metropolis
   criterion, the spin is updated or rejected.

The trial moment direction is chosen uniformly on the unit sphere with
:math:`\approx 33\%` probability, or via a Gaussian perturbation (thermal width
:math:`\delta \propto T^{0.2}`) with same probability, or via antipodal flip
(remaining 33%). This mixed scheme provides efficient sampling across all
temperature regimes without manual tuning of step sizes.

A full MC **sweep** consists of one trial flip attempt per atom averaged over
all atoms in the system.

.. _mc-heatbath:

-------------------------------------------------
Heat bath algorithm
-------------------------------------------------

The heat bath (also called microcanonical ensemble or "mcheatbath") algorithm
performs rejection-free Monte Carlo by drawing the new moment direction directly
from the Boltzmann distribution of the effective magnetic field [Glauert1974]_.

Mode H
   Select this mode to enable heat bath algorithm. Computes the total effective field
   acting on each spin and samples the new moment orientation directly from the
   canonical distribution. Provides higher acceptance rates than Metropolis, typically
   leading to faster equilibration. Particularly effective near phase transitions where
   energy barriers are significant.

**Advantages of heat bath over Metropolis:**

- Zero rejection rate (all spins are updated)
- Better equilibration near critical temperatures
- Reduced autocorrelation times in equilibrium
- More accurate critical exponents in finite-size scaling

**Computational cost:** Slightly higher per-spin due to Boltzmann sampling but
offset by lack of rejections.

.. _mc-glauber:

-------------------------------------------------
Glauber dynamics (non-equilibrium relaxation)
-------------------------------------------------

The Glauber algorithm (mode D) uses the Glauber acceptance criterion to study
non-equilibrium relaxation and aging in magnetic systems [Glauber1963]_.

Mode D
   Select this mode for Glauber dynamics, implementing relaxation toward equilibrium.
   Useful for non-equilibrium studies, aging, and transient dynamics. Related to
   Langevin spin dynamics in the heavily overdamped limit with time rescaled by
   microscopic rates.

Glauber dynamics provides access to:

- Aging exponents and effective temperature in supercooled states
- Coarsening dynamics during phase separation
- Relaxation times in glassy systems

.. _mc-spinice:

-------------------------------------------------
Spin-ice and loop algorithms
-------------------------------------------------

For constrained model systems (spin ice with ice-rule constraints, frustrated
dipolar lattices), the loop algorithm provides efficient global updates that
respect conservation laws [Evertz1993]_.

Mode L
   Select for loop algorithm or spin-ice updates. Generates large-scale correlated
   moves respecting the constraint structure. Essential for systems with frustration
   where single-spin flips are inefficient. Automatically invokes spin-ice
   (``mc_update_spinice``) routines.

The loop algorithm is particularly effective for:

- Spin-ice models with ice-rule constraints
- Frustrated lattices (triangular, kagomé)
- Artificial spin ice structures
- Systems with topological defects (monopoles, vortices)

.. _mc-ising:

-------------------------------------------------
Ising model mode
-------------------------------------------------

Mode I
   Select for Ising model dynamics. Restricts spin orientations to discrete values
   (up/down along easy axis) and uses specialized flip routines (``Ising_random_flip``).
   More efficient than full Heisenberg sampling for systems where anisotropy freezes
   out transverse fluctuations.

.. _mc-lsf-induced:

-------------------------------------------------
Local spin fluctuation (LSF) and induced moments
-------------------------------------------------

The code supports two advanced MC modes for treating electronic and magnetic
disorder:

1. **Local Spin Fluctuation (LSF)** [Ruban2004]_: Models the electronic spin
   fluctuations by allowing the magnetic moment magnitude at each site to vary
   according to a specified distribution of configurations. Implemented in
   ``mc_update_LSF``.

2. **Induced Moments** [Ebert2010]_: Treats systems with induced magnetic moments
   on non-magnetic sites, coupled self-consistently to the fixed moments.
   Implemented in ``mc_update_ind_mom``.

do_lsf
   Enable local spin fluctuation calculations (Y=yes, *N=no*). When enabled, each
   atom can sample multiple magnetic moment magnitudes from a configuration list. The
   sampling weights are determined by internal energy. Requires specification of
   ``conf_num`` (number of configurations per site). Output includes averaged energies
   for each configuration.

conf_num
   Number of magnetic moment magnitude configurations per atom in LSF mode (integer,
   ≥ 1). Typical values: 1–10. Each configuration represents an electronic state with
   different magnetic moment. Higher values improve description of disorder but
   increase computational cost. Default: 1 (pure Heisenberg model, no fluctuations).

lsf_metric
   Metric for phase space integration in LSF (1=Murata-Doniach, 2=Jacobian, *default
   1*). Determines the mathematical weighting when transforming between moment
   magnitude and energy probability distributions.

lsf_window
   Range of moment variation in LSF as fraction of average moment (0 < value ≤ 1,
   default 0.2). Sets the width of allowed moment fluctuations around the nominal
   value. Larger windows allow more magnetic disorder but may reduce thermodynamic
   stability.

lsf_interpolate
   Interpolate LSF energy between configurations (Y=yes, *N=no*, default N). When
   enabled, energies are smoothly interpolated to improve accuracy of derived
   quantities (specific heat, susceptibility).

lsf_field
   Choose field for LSF energy calculation (L=local field on atom, T=total system
   field, default T). Local field improves accuracy for dilute alloys; total field
   for concentrated systems.

exc_inter
   Interpolate exchange coupling between ferromagnetic (FM) and disordered local
   moment (DLM) limits (Y=yes, *N=no*, default N). Used in combination with LSF to
   smoothly transition exchange strength based on local moment orientation.


-------------------------------------------------
Trial spin selection and update scheme
-------------------------------------------------

Each Monte Carlo sweep consists of the following stages:

**1. Random atom sequence:** Generate a random permutation of atoms using
``choose_random_atom_x``. This ensures each atom is attempted exactly once per
sweep, preventing artificial correlations from sequential updating.

**2. Trial move generation:** Generate trial moment directions using
``choose_random_flip`` (Metropolis/Glauber) or geometry-dependent sampling (heat bath).

**3. Energy calculation:** Compute energy difference :math:`\Delta E` for the
trial configuration using ``calculate_energy``. This includes:

   - Heisenberg exchange interactions
   - Dzyaloshinskii-Moriya interactions (DMI)
   - Pseudo-dipolar interactions
   - Biquadratic exchange
   - Scalar chirality
   - Symmetric anisotropic exchange
   - Magnetic anisotropy (uniaxial, cubic, combined)
   - Dipolar interactions (brute-force or macrocell-accelerated)
   - External magnetic fields

**4. Acceptance decision:** Apply acceptance criterion (Metropolis, heat bath, or
Glauber) via ``flip_a``, ``flip_g``, or ``flip_h`` routines.

**5. Update:** If accepted, update moment vectors ``emom`` and magnitudes ``emomM``.

-------------------------------------------------
Energy calculation and Hamiltonian interactions
-------------------------------------------------

The total energy for a configuration is [Eriksson2017]_:

.. math::

   E = -\frac{1}{2}\sum_{i,j} J_{ij} \mathbf{m}_i \cdot \mathbf{m}_j
   - \sum_i \mathbf{B}_{\mathrm{ext}} \cdot \mathbf{m}_i
   + \sum_i E_{\mathrm{ani}}(\mathbf{m}_i)
   + E_{\mathrm{DM}} + E_{\mathrm{DIP}} + \ldots

The energy difference only requires summing over neighbors of the flipped atom
(``nlistsize``), reducing computational cost from O(N²) to O(1) per update.

**Supported interactions:**

- **Exchange:** Isotropic Heisenberg or exchange tensor (SKKR style)
- **Anisotropy:** Uniaxial (first and second order K₁, K₂), cubic, or combined
- **Relativistic:** Dzyaloshinskii-Moriya, pseudo-dipolar, scalar chirality
- **Higher order:** Biquadratic exchange, four-spin ring exchange
- **Dipolar:** Brute-force or FFT-accelerated macrocell method
- **Induced moments:** Self-consistent coupling to non-magnetic sites

-------------------------------------------------
Temperature control and thermodynamic sampling
-------------------------------------------------

-------------------------------------------------
Temperature control and thermodynamic sampling
-------------------------------------------------

Temp
   Simulation temperature in Kelvin (real number). Controls the Metropolis acceptance
   probability via :math:`\beta = 1/(k_B T)`. Lower temperatures require larger energy
   barriers to be accepted, leading to slower equilibration but more accurate ground
   state sampling. Temperature can be swept externally by running multiple simulations
   at different ``Temp`` values.

compensate_drift
   Enable quantum heat bath (QHB) thermostat for improved temperature control
   (Y=yes, *N=no*). When enabled, introduces a Nosé-Hoover-type thermostat that
   rescales effective temperature to account for drift in energy-based thermostats.
   Requires specification of ``temprescale`` and ``temprescalegrad``.

temprescale
   Temperature rescaling factor for QHB (real, typically 0.1–1.0). Multiplicative
   factor in effective temperature :math:`T_{\mathrm{eff}} = T_{\mathrm{rescale}}\times
   T + T_{\mathrm{rescalegrad}} \times t`. Only used if ``compensate_drift Y``.

temprescalegrad
   Temperature rescaling gradient for QHB thermostat (real, typically ≤ 0.01 K/step).
   Rate of change of effective temperature during simulation. Allows temperature
   scheduling without external modification of input files.

.. _mc-annealing:

-------------------------------------------------
Simulated annealing
-------------------------------------------------

Simulated annealing gradually reduces temperature to drive the system toward its
ground state or metastable configurations [Kirkpatrick1983]_. The UppASD code
supports multi-phase annealing schedules where temperature and iteration count
are specified for each phase independently.

ip_mcanneal
   Number of temperature phases for simulated annealing schedule (integer, ≥ 0).
   Set to 0 to disable annealing (default). When set to :math:`N > 0`, the code
   expects :math:`N` pairs of input lines following this keyword, each containing
   ``(mcnstep_i, Temp_i)``. The simulation runs for ``mcnstep_i`` Monte Carlo
   sweeps at temperature ``Temp_i`` in phase :math:`i`.

**Annealing schedule format:**

If ``ip_mcanneal`` is set to 3, the input file should contain::

   ip_mcanneal         3
   5000  500.0          ! Phase 1: 5000 sweeps at 500 K
   5000  250.0          ! Phase 2: 5000 sweeps at 250 K
   10000  50.0          ! Phase 3: 10000 sweeps at 50 K

The temperature is held constant throughout each phase, then switched discretely
to the next phase temperature. For continuous temperature reduction, divide the
annealing into more phases with smaller temperature steps.

**Annealing strategy guidance:**

- **Initial temperature:** Should be high enough to overcome energy barriers
  (typically 2–5 times estimated critical temperature)
- **Cooling rate:** Slower cooling (more phases, smaller :math:`\Delta T`) generally
  gives better results but increases wall-clock time. A logarithmic schedule
  :math:`T_i \propto \log(i)` is often effective
- **Final temperature:** Set to target temperature for ground state or equilibrium
  measurement. For zero-temperature structures, use a small finite temperature
  (e.g., 1 K) for numerical stability
- **Number of sweeps per phase:** Should be sufficient for system to relax at each
  temperature; typical values are :math:`\geq N_{\text{atom}}` sweeps per phase

**Example 10-phase logarithmic annealing schedule (1000 K → 10 K):**

.. code-block:: text

   ip_mcanneal         10
   2000  1000.0
   2000  562.0
   2000  316.0
   2000  178.0
   2000  100.0
   2000  56.2
   2000  31.6
   2000  17.8
   2000  10.0
   5000  10.0          ! Extra sweeps at final temperature for convergence

**Important notes on annealing:**

- Annealing schedules are typically **system and Hamiltonian specific**; optimal
  cooling rates should be determined empirically for your problem
- Cumulant measurements are automatically enabled during annealing but may not be
  meaningful during temperature ramps; analyze final phase measurements separately
- Multiple independent annealing runs with different random seeds should be
  performed to assess solution quality and variance
- For complex landscapes with many metastable states, consider escaping local
  minima using parallel tempering (``Mensemble`` > 1) instead of single-run annealing

-------------------------------------------------
Output and statistical measurements
-------------------------------------------------

All Monte Carlo simulations automatically enable cumulant measurements (equivalent
to ``do_cumu Y``) regardless of input setting. This is essential because:

- Equilibrium MC generates uncorrelated samples from the Boltzmann distribution
- Binder cumulants provide unbiased phase transition characterization
- Magnetic susceptibility from fluctuations is automatically available
- Specific heat is computed from energy fluctuations

Additional output is controlled by standard measurement keywords:

- ``do_avrg Y`` enables time-averaged magnetization output
- ``avrg_step``, ``avrg_buff`` control sampling frequency and buffering
- ``plotenergy 1`` computes and outputs energy components for specific heat

See :doc:`../observables/averages` for detailed measurement configuration.

-------------------------------------------------
Notes and best practices
-------------------------------------------------

**Choosing an algorithm:**

- **Metropolis (M):** Simplest, most flexible; good for prototyping and low temperatures
- **Heat bath (H):** Better rejection-free sampling; recommended for critical phenomena
- **Glauber (D):** For non-equilibrium studies and relaxation analysis
- **Spin-ice (L):** Essential for constrained systems; uses different acceptance rules
- **Ising (I):** For pure Ising models with strong uniaxial anisotropy

**Equilibration:**

- Run for :math:`\gtrsim N_{\text{atom}} \times 100` sweeps before measuring
- Monitor cumulant output for convergence of :math:`\langle m^2 \rangle` and :math:`\langle m^4 \rangle`
- For systems near :math:`T_c`, increase equilibration time proportionally to :math:`L^z` where :math:`L` is system size

**Statistical accuracy:**

- Perform multiple independent runs with different random seeds
- Use ensemble averaging (``Mensemble`` > 1) for reduced variance
- For finite-size scaling, run multiple system sizes at same temperature

**Computational efficiency:**

- Heat bath typically 2–3× more expensive per sweep than Metropolis
- Use macrocell dipolar acceleration (``do_dip 2``) for large systems
- Reduce measurement frequency (larger ``avrg_step``) to minimize I/O overhead

-------------------------------------------------
Example inpsd.dat snippet
-------------------------------------------------

**Basic Metropolis equilibration:**

.. code-block:: text

   ! Monte Carlo equilibration
   mode          M                ! Metropolis algorithm
   Temp          300.0            ! Temperature 300 K
   mcnstep       10000            ! 10,000 MC sweeps

   ! Enable measurements
   do_avrg       Y
   avrg_step     10               ! Sample every 10 sweeps
   cumu_buff     10               ! Cumulants output every 100 sweeps

**Heat bath with fast equilibration:**

.. code-block:: text

   mode          H                ! Heat bath (rejection-free)
   Temp          1.0              ! Reduced units
   mcnstep       5000             ! Fewer sweeps needed due to better sampling

   ! Cumulant measurements (automatic)
   plotenergy    1                ! Compute energy for specific heat
   do_avrg       Y
   avrg_step     20

**Temperature-dependent phase transition study:**

.. code-block:: text

   ! Run at multiple temperatures for phase diagram
   mode          H
   Temp          2.5              ! Near critical temperature
   mcnstep       50000            ! Long runs for accuracy

   ! High-resolution cumulant measurement
   do_avrg       Y
   avrg_step     100
   cumu_buff     50               ! Output every 5000 sweeps

   ! Finite-size scaling
   plotenergy    1
   do_cumu       Y                ! Redundant but explicit

**Local spin fluctuation with disordered moments:**

.. code-block:: text

   mode          M
   Temp          200.0
   mcnstep       20000

   ! Enable LSF disorder
   do_lsf        Y
   conf_num      4                ! 4 moment configurations per site
   lsf_window    0.25             ! ±25% moment variation
   lsf_field     L                ! Use local field

   plotenergy    1                ! Monitor disorder effects
   do_avrg       Y

**Induced moments in partially filled systems:**

.. code-block:: text

   mode          M
   Temp          100.0
   mcnstep       15000

   ! Induced moments on non-magnetic sites
   ind_mom_flag  Y                ! Enable induced moment coupling
   ind_mom_type  2                ! Ebert-type treatment

   do_avrg       Y
   avrg_step     50

-------------------------------------------------
Related keywords and cross-references
-------------------------------------------------

- ``Temp``: Temperature in Kelvin (see :doc:`../input/simulation`)
- ``mcnstep``: Number of Monte Carlo sweeps (see :doc:`../input/simulation`)
- ``Mensemble``: Number of ensemble replicas for parallel tempering (see :doc:`../input/simulation`)
- ``do_avrg``, ``avrg_step``: Average magnetization measurements (see :doc:`../observables/averages`)
- ``do_cumu``: Cumulant measurements (automatic for MC; see :doc:`../observables/averages`)
- ``plotenergy``: Energy component output for specific heat (see :doc:`../input/observables`)
- Hamiltonian interactions: exchange, anisotropy, dipolar (see :doc:`../input/hamiltonian`)

Output files use the simulation identifier ``simid`` (see :doc:`../input/system`).


References
----------

See the centralized :doc:`../references` for full bibliographic entries:

- [Landau2014]_ - Guide to Monte Carlo simulation in statistical physics
- [Binney1992]_ - Theory of critical phenomena and renormalization group
- [Metropolis1953]_ - Equation of state calculations by fast computing machines
- [Glauert1974]_ - Nonuniversal critical dynamics in Monte Carlo simulations
- [Glauber1963]_ - Time-dependent statistics of the Ising model
- [Hinzke2000]_ - Stochastic dynamics of magnetic nanoparticles
- [Evertz1993]_ - Cluster algorithm for vertex models
- [Ruban2004]_ - Surface segregation energies in transition-metal alloys
- [Ebert2010]_ - Calculating condensed matter properties using KKR-Green's function method
- [Kirkpatrick1983]_ - Optimization by simulated annealing
- [Eriksson2017]_ - Atomistic spin dynamics foundations and applications
- [Binder1981]_ - Finite size scaling and simulation of first order phase transitions

