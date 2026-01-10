Averages and cumulants
======================

The UppASD code supports comprehensive statistical measurements of magnetization
and thermodynamic observables during spin dynamics simulations. The averaging and
cumulant measurement subsystem provides [Binder1981]_, [Landau2014]_:

- :ref:`Average magnetization with standard deviation <avrg-basic>`
- :ref:`Type-projected (sublattice) average magnetizations <avrg-projected>`
- :ref:`Chemically-projected average magnetizations for random alloys <avrg-chemical>`
- :ref:`Binder cumulants for phase transition characterization <cumulants-binder>`
- :ref:`Magnetic susceptibility and specific heat <cumulants-thermal>`
- :ref:`Antiferromagnetic order parameter measurements <cumulants-afm>`

All measurements support ensemble averaging for improved statistical accuracy in
Monte Carlo and Langevin dynamics simulations. The subsystem employs a buffered
output strategy to minimize I/O overhead during long simulations.

.. tip::

   **Input keywords:** For a complete list of all average and cumulant input
   parameters (``do_avrg``, ``avrg_step``, ``do_cumu``, etc.), see the
   comprehensive reference in :doc:`input-keywords-observables`.

-------------------------------------------------
Canonical measurement definitions
-------------------------------------------------

The average magnetization is defined as the ensemble-averaged total magnetic moment:

.. math::

   \langle \mathbf{M} \rangle = \frac{1}{N_{\text{ens}}} \sum_{k=1}^{N_{\text{ens}}}
   \sum_{i=1}^{N_{\text{atom}}} \mathbf{m}_i^{(k)},

where :math:`N_{\text{ens}}` is the number of ensembles (replicas) and
:math:`\mathbf{m}_i^{(k)}` is the magnetic moment vector of atom :math:`i` in
ensemble :math:`k`. The magnitude of the average magnetization per atom is:

.. math::

   m_{\text{avg}} = \frac{|\langle \mathbf{M} \rangle|}{N_{\text{atom}}}.

The **Binder cumulant** (also called Binder parameter) is a fourth-order
statistical moment used to identify phase transitions [Binder1981]_:

.. math::

   U_L = 1 - \frac{\langle m^4 \rangle}{3 \langle m^2 \rangle^2},

where :math:`m = |\mathbf{M}|/N_{\text{atom}}` is the magnetization per atom.
At a second-order phase transition, :math:`U_L` exhibits universal behavior
independent of system size for sufficiently large systems. Typical values:

- :math:`U_L \approx 0` for disordered (paramagnetic) phase
- :math:`U_L \approx 2/3` for ordered (ferromagnetic) phase in 3D Heisenberg model
- :math:`U_L` crosses at :math:`T_c` for different system sizes (finite-size scaling)

The **magnetic susceptibility** is calculated from magnetization fluctuations:

.. math::

   \chi = \frac{\langle m^2 \rangle - \langle m \rangle^2}{k_B T} \cdot N_{\text{atom}} \mu_B^2,

with units of :math:`k_B` (Boltzmann constant). The **specific heat** is
obtained from energy fluctuations:

.. math::

   C_V = \frac{\langle E^2 \rangle - \langle E \rangle^2}{k_B T^2}
   \cdot N_{\text{atom}} \left(\frac{\partial T_{\text{eff}}}{\partial T}\right),

where the derivative term accounts for temperature rescaling in quantum heat bath
(QHB) thermostats [Evans1985]_.

For **antiferromagnetic systems**, the staggered magnetization (Néel order parameter)
is measured using a sublattice structure vector :math:`\boldsymbol{\ell}`:

.. math::

   \mathbf{L} = \sum_{i=1}^{N_A} \ell_i \sum_{\mathbf{R}} \mathbf{m}_{i,\mathbf{R}},

where :math:`N_A` is the number of atoms per unit cell, :math:`\mathbf{R}` runs
over all unit cells, and :math:`\ell_i = \pm 1` defines the sublattice structure.

.. _avrg-basic:

-------------------------------------------------
Basic average magnetization measurements
-------------------------------------------------

The fundamental measurement mode samples the total magnetization vector and its
magnitude at regular intervals during the simulation. Statistics are accumulated
across all ensembles and output is buffered to reduce file I/O.

``do_avrg``
   Enable sampling and printing of average magnetization (Y=yes, *N=no*). When enabled, the code
   calculates :math:`\langle \mathbf{M} \rangle`, its magnitude, and standard deviation across
   ensembles. Output is written to ``averages.simid.out``.

``avrg_step``
   Sampling interval in simulation steps (integer). Determines how frequently the magnetization is
   measured. Default: 100. Smaller values provide finer temporal resolution but increase file size.
   For dynamics simulations, set this to capture relevant timescales (e.g., precession periods).

``avrg_buff``
   Buffer size for average magnetization (integer). Number of samples to accumulate in memory before
   writing to file. Default: 10. Larger values reduce I/O overhead but increase memory usage.
   Total samples in output file = ``avrg_buff`` :math:`\times` (number of flush operations).

The output file ``averages.simid.out`` contains columns:

- Column 1: Simulation step (or time in seconds if ``real_time_measure Y``)
- Columns 2–4: :math:`\langle M_x \rangle`, :math:`\langle M_y \rangle`, :math:`\langle M_z \rangle`
- Column 5: :math:`\sqrt{\langle M_x \rangle^2 + \langle M_y \rangle^2 + \langle M_z \rangle^2}` (total magnitude)
- Column 6: Standard deviation of magnetization magnitude across ensembles

All magnetization values are in units of total Bohr magnetons (:math:`\mu_B`).

.. _avrg-projected:

-------------------------------------------------
Type-projected (sublattice) averages
-------------------------------------------------

For systems with multiple magnetic atom types or sublattices, type-projected
averages provide separate statistics for each magnetic species. This is essential
for analyzing ferrimagnetic, antiferromagnetic, or multi-component systems.

``do_proj_avrg``
   Enable type-projected average measurements (Y=by type, A=by site, *N=no*). When set to Y, averages
   are calculated separately for each atom type (chemical species). When set to A, averages are
   calculated for each inequivalent site in the unit cell. Automatically enables ``do_avrg``.
   Output is written to ``projavgs.simid.out``.

**Mode Y (by type)**: Groups atoms by their type index (as specified in the
``momfile`` or ``posfile``). All atoms of type :math:`k` are summed together,
accounting for periodic repetitions via :math:`N_1 \times N_2 \times N_3`.

**Mode A (by site)**: Groups atoms by their position within the unit cell (site
index). Useful when different sites have the same chemical type but distinct
magnetic environments.

The output file ``projavgs.simid.out`` contains:

- Column 1: Simulation step (or time)
- Column 2: Type or site index
- Column 3: Average magnetization magnitude for this type/site
- Column 4: Standard deviation
- Columns 5–7: :math:`\langle M_x \rangle`, :math:`\langle M_y \rangle`, :math:`\langle M_z \rangle` for this type/site

All values are normalized per unit cell (divided by :math:`N_1 N_2 N_3`).

.. _avrg-chemical:

-------------------------------------------------
Chemically-projected averages for random alloys
-------------------------------------------------

For random alloy simulations (``do_ralloy 1``), the chemical composition varies
from site to site. Chemically-projected averages track magnetization separately
for each chemical component.

``do_projch_avrg``
   Enable chemically-projected average measurements for random alloys (Y=yes, *N=no*). Requires
   ``do_ralloy 1``. Calculates averages for each chemical species (element) present in the alloy,
   accounting for concentration fluctuations. Output is written to ``projchavgs.simid.out``.

The chemical species are indexed by ``achem_ch`` values read from the chemical
configuration file (``chemfile``). The measurement accounts for the actual spatial
distribution of chemical species.

The output file ``projchavgs.simid.out`` has the same format as ``projavgs.simid.out``
but with an additional column:

- Column 8: Sum of magnetization magnitudes :math:`\sum_k \langle |M_k| \rangle`

This total can differ from the overall system average due to non-collinear magnetic
structures or concentration correlations in the alloy.

.. _cumulants-binder:

-------------------------------------------------
Binder cumulant and phase transition detection
-------------------------------------------------

The Binder cumulant is a powerful tool for locating phase transitions and
characterizing critical behavior. It is size-independent at the critical point,
making it ideal for finite-size scaling analyses.

``do_cumu``
   Enable cumulant measurements (Y=standard, A=antiferromagnetic, *N=no*). Calculates Binder cumulant,
   magnetic susceptibility, and specific heat from magnetization and energy fluctuations. For Monte Carlo
   simulations, this is automatically set to Y regardless of input. Output written to
   ``cumulants.simid.out`` and ``cumulants.simid.json``.

``cumu_step``
   Sampling interval for cumulant measurements (integer). Default: 50. Should be chosen to ensure
   statistical independence between samples. For equilibrium MC, set to :math:`\sim N_{\text{atom}}` or larger.
   For dynamics, consider correlation times.

``cumu_buff``
   Buffer size for cumulant output (integer). Default: 10. Number of cumulant samples to accumulate before
   writing to file. Larger values smooth statistical noise but delay output updates.

``do_cumu_proj``
   Enable type-projected cumulant measurements (Y=yes, *N=no*). Calculates Binder cumulant and
   susceptibility separately for each atom type. Useful for analyzing sublattice ordering in complex
   magnetic structures. Output written to ``projcumulants.simid.out``.

**Mode Y (standard)**: Measures total magnetization cumulants using
:math:`m = |\mathbf{M}|/N_{\text{atom}}`. Appropriate for ferromagnetic and
paramagnetic phases.

.. _cumulants-afm:

**Mode A (antiferromagnetic)**: Measures staggered magnetization cumulants using
the Néel order parameter :math:`\mathbf{L}`. Currently uses a fixed sublattice
structure suitable for bipartite antiferromagnets (future versions will support
custom AFM vectors).

The output file ``cumulants.simid.out`` contains:

- Column 1: Cumulative sample count (normalized by ``Mensemble``)
- Column 2: :math:`\langle m \rangle`
- Column 3: :math:`\langle m^2 \rangle`
- Column 4: :math:`\langle m^4 \rangle`
- Column 5: Binder cumulant :math:`U_L`
- Column 6: Magnetic susceptibility :math:`\chi` (units of :math:`k_B`)
- Column 7: Specific heat :math:`C_V` (units of :math:`k_B`)
- Column 8: Average total energy :math:`\langle E \rangle`
- Column 9: Average exchange energy :math:`\langle E_{\text{exc}} \rangle`
- Column 10: Average local spin field energy :math:`\langle E_{\text{lsf}} \rangle`

Energy columns (8–10) are only populated if ``plotenergy 1``.

The JSON output file ``cumulants.simid.json`` provides a machine-readable summary
with keys: ``temperature``, ``magnetization``, ``binder_cumulant``, ``energy``,
``susceptibility``, ``specific_heat``, and optionally ``skyrmion_num`` if
skyrmion counting is enabled (see :doc:`input-keywords-topology`).

.. _cumulants-thermal:

-------------------------------------------------
Thermodynamic observables and finite-size scaling
-------------------------------------------------

The susceptibility and specific heat are intensive thermodynamic quantities that
exhibit characteristic behavior near phase transitions:

- **Susceptibility** :math:`\chi` diverges at :math:`T_c` as
  :math:`\chi \sim |T - T_c|^{-\gamma}` in the critical region.
- **Specific heat** :math:`C_V` shows a peak (continuous transition) or jump
  (first-order transition) at :math:`T_c`.

For **quantum heat bath (QHB)** thermostats (``compensate_drift`` enabled), the
effective temperature depends on the Nosé-Hoover chain parameters. The code
automatically accounts for :math:`\partial T_{\text{eff}}/\partial T` using
``temprescale`` and ``temprescalegrad``.

For **finite-size scaling**, the Binder cumulant crossing method is recommended:

1. Run simulations at multiple temperatures near the suspected :math:`T_c`
2. Use several system sizes :math:`L`
3. Plot :math:`U_L(T)` for each size
4. The crossing point approximates :math:`T_c(L \to \infty)`
5. Extrapolate to infinite size using :math:`T_c(L) = T_c(\infty) + a L^{-1/\nu}`

where :math:`\nu` is the correlation length critical exponent.

-------------------------------------------------
Weight and ensemble averaging details
-------------------------------------------------

The cumulant calculation uses a cumulative weighted average scheme to handle
variable timesteps and non-uniform sampling in adaptive integrators:

.. math::

   \langle O \rangle_{\text{cum}} = \frac{\sum_i w_i O_i}{\sum_i w_i},

where :math:`w_i` is the weight of sample :math:`i` (currently set to unity for
constant timestep integrators). This infrastructure supports future extensions to
variable-timestep adaptive schemes where weights would depend on :math:`\Delta t`.

All ensemble averages are computed by summing over replicas:

.. math::

   \langle O \rangle = \frac{1}{N_{\text{ens}}} \sum_{k=1}^{N_{\text{ens}}} O^{(k)},

ensuring proper statistical error reduction in parallel tempering and replica
exchange simulations.

-------------------------------------------------
Interaction with other measurement modes
-------------------------------------------------

The averaging subsystem is independent of but compatible with:

- **Trajectory output** (``do_tottraj``): Full moment configurations can be saved
  separately at different sampling rates.
- **Correlation functions** (``do_sc``): Space-time correlations are measured
  independently but can share the same timestep if ``avrg_step`` equals ``sc_step``.
- **Real-time measurements**: Set ``real_time_measure Y`` to output physical time
  (seconds) instead of simulation steps in all averaging output files. This requires
  specifying the timestep ``delta_t`` in the simulation parameters.

-------------------------------------------------
Notes and best practices
-------------------------------------------------

**For Monte Carlo simulations:**

- Perform sufficiently long runs to reduce statistical noise in :math:`\langle m^4 \rangle`
- Check convergence by monitoring running averages in ``cumulants.simid.out``

**For spin dynamics simulations:**

- Set ``avrg_step`` to sample faster than the slowest relevant timescale (e.g., domain wall motion)
- Equilibrate thoroughly before enabling cumulant measurements

**For parallel tempering:**

- Ensure ``Mensemble`` is sufficient for accurate ensemble averaging (typically ≥ 16)

-------------------------------------------------
Example inpsd.dat snippet
-------------------------------------------------

**Minimal configuration for average magnetization:**

.. code-block:: text

   do_avrg       Y
   avrg_step     100
   avrg_buff     10

**Complete configuration with cumulants and projections:**

.. code-block:: text

   # Enable basic averages
   do_avrg       Y
   avrg_step     50         # Sample every 50 steps
   avrg_buff     20         # Buffer 20 samples before writing

   # Enable type-projected averages (ferrimagnetic system)
   do_proj_avrg  Y

   # Enable cumulants for phase transition study
   do_cumu       Y
   cumu_step     100        # Sample every 100 steps for statistics
   cumu_buff     10         # Output every 1000 steps

   # Plot energy components
   plotenergy    1

   # Use real time output (requires delta_t specification)
   real_time_measure Y

**Finite-size scaling study (Monte Carlo):**

.. code-block:: text

   # High-precision cumulant measurement
   do_avrg       Y
   avrg_step     1000       # Sample once per 1000 MC sweeps
   avrg_buff     100        # Accumulate 100 samples

   do_cumu       Y          # Automatically enabled for MC
   cumu_step     1000       # Match avrg_step for consistency
   cumu_buff     50         # Output every 50000 sweeps

   plotenergy    1          # Calculate specific heat

   # Temperature is set in separate section
   Temp          1.5        # Near suspected T_c

**Antiferromagnetic order parameter tracking:**

.. code-block:: text

   do_avrg       Y
   avrg_step     100

   # Measure staggered magnetization cumulant
   do_cumu       A          # AFM mode
   cumu_step     100
   cumu_buff     10

   # Also measure sublattice-resolved cumulants
   do_cumu_proj  Y

-------------------------------------------------
Related keywords and cross-references
-------------------------------------------------

The averages and cumulants subsystem interacts with several global parameters:

- ``Temp``: Temperature for susceptibility and specific heat calculations (see :doc:`input-keywords-simulation`)
- ``delta_t``: Timestep for real-time measurements (see :doc:`input-keywords-simulation`)
- ``Mensemble``: Number of ensembles for statistical averaging (see :doc:`input-keywords-simulation`)
- ``plotenergy``: Enables energy measurements for specific heat (see :doc:`input-keywords-observables`)
- ``do_ralloy``: Enables random alloy simulations for chemical projections (see :doc:`input-random-alloys`)
- ``compensate_drift``: Activates QHB thermostat with temperature rescaling (see :doc:`input-keywords-simulation`)

Output files are always named with the simulation identifier ``simid`` (see :doc:`input-keywords-system`).
