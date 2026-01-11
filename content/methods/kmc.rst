.. _input-keywords-kmc:

Kinetic Monte Carlo (KMC)
======================================

Overview
--------

The **Kinetic Monte Carlo (KMC)** module in UppASD implements a **rare event dynamics** approach
for modeling slow particle motion coupled to magnetic dynamics. This is particularly relevant for
systems where certain degrees of freedom (e.g., charge carriers, vacancies, interstitials) evolve
on much longer timescales than spin dynamics.

**Primary applications:**

- **Magnetic polarons**: Charge carriers (electrons/holes) that locally distort magnetic structure
- **Vacancy diffusion**: Point defects hopping between lattice sites
- **Interstitial migration**: Diffusion of interstitial atoms in magnetic materials
- **Spin-lattice coupling**: Systems where atomic motion affects magnetic interactions

**Key concept:** Particles overcome energy barriers between lattice sites with transition rates
governed by the Arrhenius law. The magnetic Hamiltonian is dynamically updated as particles move,
allowing self-consistent treatment of coupled spin-charge or spin-lattice dynamics.

.. note::
   **Inherent symmetry breaking:** Due to the nature of particle localization, the KMC module
   is **not compatible** with ``do_reduced='Y'`` simulations.


Physical Background
-------------------

Arrhenius Law and Transition Rates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Thermally activated hopping processes follow the Arrhenius law for transition rates:

.. math::

   \Gamma_{i \to j} = \nu_0 \exp\left(-\frac{\Delta E_{ij}}{k_B T}\right)

where:

- :math:`\Gamma_{i \to j}`: Transition rate from site :math:`i` to site :math:`j` (Hz)
- :math:`\nu_0`: Attempt frequency (``rate0`` parameter), typically :math:`10^{12}` Hz for phonon-assisted processes
- :math:`\Delta E_{ij}`: Energy barrier between sites (from ``barrfile``)
- :math:`k_B`: Boltzmann constant
- :math:`T`: Temperature

**Physical interpretation:** :math:`\nu_0` represents the fundamental vibrational frequency of the
particle in its local potential well (harmonic approximation). For atomic processes, this is typically
the Debye frequency; for electronic processes (polarons), it may be related to phonon modes.

Electric Field Effects
^^^^^^^^^^^^^^^^^^^^^^

For charged particles (polarons), an external electric field modifies the barrier:

.. math::

   \Delta E_{ij}^{\text{eff}} = \Delta E_{ij} + q \mathbf{E} \cdot \mathbf{r}_{ij}

where:

- :math:`q \mathbf{E} \cdot \mathbf{r}_{ij}`: Electrostatic energy change (``do_efield='Y'``)
- :math:`\mathbf{E}`: Applied electric field vector (``efield``)
- :math:`\mathbf{r}_{ij}`: Displacement vector between sites

This creates a bias in hopping rates, driving **directed migration** (drift) in addition to diffusion.

Kinetic Monte Carlo Algorithm
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Standard KMC (``kmc_method=1``) proceeds as follows:

**1. Calculate transition rates:**
   
   For each KMC particle at site :math:`i`, compute rates to all neighboring sites :math:`j`:

   .. math::

      \Gamma_{i \to j} = \nu_0 \exp\left(-\frac{\Delta E_{ij} + \mathbf{E} \cdot \mathbf{r}_{ij}}{k_B T}\right)

   **Site blocking:** If target site :math:`j` is already occupied, :math:`\Gamma_{i \to j} = 0`.

**2. Select event stochastically:**
   
   Total rate: :math:`\Gamma_{\text{tot}} = \sum_j \Gamma_{i \to j}`

   Randomly select jump :math:`i \to k` with probability:

   .. math::

      P(i \to k) = \frac{\Gamma_{i \to k}}{\Gamma_{\text{tot}}}

**3. Calculate waiting time:**
   
   The time until the next jump is exponentially distributed:

   .. math::

      \tau_{\text{KMC}} = -\frac{1}{\Gamma_{\text{tot}}} \ln(\text{rand})

   Convert to simulation time steps:

   .. math::

      \Delta t_{\text{steps}} = \text{nint}\left(\frac{\tau_{\text{KMC}}}{\delta t}\right) + m_{\text{step}}

   where :math:`\delta t` is the integration timestep and :math:`m_{\text{step}}` is the current step.

**4. Update Hamiltonian:**
   
   When the particle moves:
   
   - Swap atomic types, chemical identities, and moment magnitudes between sites :math:`i` and :math:`k`
   - Update exchange couplings :math:`J_{ij}` for neighboring bonds affected by the particle
   - Recalculate effective fields and transition rates

Dynamic Hamiltonian Coupling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The KMC module implements **full self-consistency** between particle positions and magnetic interactions:

**Central site update:**

- Atomic type ``atype(i)`` ↔ ``atype(k)``
- Chemical type ``achem_ch(i)`` ↔ ``achem_ch(k)``
- Moment magnitude ``mmom(i)`` ↔ ``mmom(k)``

**Neighbor coupling update:**

For all neighbors :math:`n` of both sites :math:`i` and :math:`k`, exchange interactions
:math:`J_{in}` and :math:`J_{kn}` are swapped if they connect to equivalent bonding vectors.

This ensures that the **magnetic environment follows the particle**, enabling accurate treatment
of magnetic polarons where local exchange is modified by the charge carrier.

**Tolerance:** Bonding vectors are compared with tolerance :math:`\epsilon = 0.005` to identify equivalent neighbors.


Input Parameters
----------------

The following keywords control KMC behavior in the ``inpsd.dat`` file:

.. list-table:: KMC Keywords
   :widths: 20 15 50 15
   :header-rows: 1

   * - Keyword
     - Type
     - Description
     - Default
   * - ``do_kmc``
     - Character(1)
     - Enable Kinetic Monte Carlo (``'Y'``/``'N'``)
     - ``'N'``
   * - ``kmc_method``
     - Integer
     - Algorithm: ``1`` = Standard KMC
     - ``1``
   * - ``rate0``
     - Real
     - Attempt frequency :math:`\nu_0` (Hz)
     - ``0.0``
   * - ``barrfile``
     - Character(35)
     - Filename with energy barriers
     - (required)
   * - ``kmc_posfile``
     - Character(35)
     - Filename with initial particle positions
     - (required)
   * - ``do_efield``
     - Character(1)
     - Enable external electric field (``'Y'``/``'N'``)
     - ``'N'``
   * - ``efield``
     - Real(3)
     - Electric field vector (eV/Å or dimensionless)
     - ``0.0 0.0 0.0``
   * - ``kmc_step``
     - Integer
     - Measurement interval (time steps)
     - ``100``
   * - ``kmc_buff``
     - Integer
     - Buffer size for trajectory output
     - ``10``
   * - ``do_prn_kmc``
     - Character(1)
     - Output mode: ``'Y'`` = regular, ``'D'`` = adaptive
     - ``'N'``

Parameter Details
^^^^^^^^^^^^^^^^^

**``rate0`` - Attempt Frequency**

Typical values depend on the physical process:

- **Phonon-assisted atomic diffusion:** :math:`\nu_0 \sim 10^{12} - 10^{13}` Hz (Debye frequency)
- **Magnetic polaron hopping:** :math:`\nu_0 \sim 10^{11} - 10^{12}` Hz (magnon-electron coupling)
- **Vacancy migration in metals:** :math:`\nu_0 \sim 10^{12}` Hz

.. note::
   ``rate0`` must be specified in **Hz**. It fundamentally sets the timescale for rare event dynamics.

**``barrfile`` - Energy Barrier File**

Similar format to ``jfile`` (exchange file):

.. code-block:: none

   # Format (without random alloy):
   isite  jsite  rx  ry  rz  barrier(eV)
   
   # Format (with random alloy):
   isite  jsite  ichem  jchem  rx  ry  rz  barrier(eV)

- ``isite``, ``jsite``: Atom indices
- ``rx, ry, rz``: Neighbor vector (Cartesian or Direct coordinates matching ``posfiletype``)
- ``barrier``: Energy barrier :math:`\Delta E_{ij}` in **eV**

**Example:**

.. code-block:: none

   1   2   1.0  0.0  0.0   0.15
   1   3   0.0  1.0  0.0   0.15
   1   5  -1.0  0.0  0.0   0.20

**Physical meaning:** Barriers represent saddle point energies for particle migration. For polarons,
these include both elastic distortion energy and magnetic exchange differences.

**``kmc_posfile`` - Initial Particle Positions**

Specifies starting positions of KMC particles:

.. code-block:: none

   # Format (Cartesian coordinates):
   site  type  Rx  Ry  Rz
   
   # Format (Direct coordinates):
   site  type  rx  ry  rz
   
   # Format (random alloy):
   site  type  chem  conc  Rx  Ry  Rz

- ``site``: Atom index where particle is initially located
- ``type``: Particle type (1=polaron site, others=host)
- ``chem``, ``conc``: Chemical type and concentration (for random alloys)
- Coordinates: Match ``posfiletype`` (``'C'`` or ``'D'``)

**Example (Cartesian):**

.. code-block:: none

   512  2  5.0  5.0  5.0
   1024 2  10.0 10.0 10.0

This places two KMC particles at atoms 512 and 1024.

**``efield`` - Electric Field**

Vector format: ``efield(1) efield(2) efield(3)``

.. code-block:: none

   efield   0.01   0.0   0.0

Applies field along :math:`x`-axis. **Units:** The energy term :math:`\mathbf{E} \cdot \mathbf{r}`
should give eV when :math:`\mathbf{r}` is in Ångströms.

**Typical values:**

- Lab-scale fields: :math:`E \sim 0.001 - 0.1` eV/Å (:math:`\sim 10^8 - 10^{10}` V/m)
- Effective polaron fields: :math:`E \sim 0.01` eV/Å


Output Files
------------

``kmc_info.<simid>.out``
^^^^^^^^^^^^^^^^^^^^^^^^

Main trajectory output file with format:

.. code-block:: none

   # Columns:
   # 1: Particle index
   # 2: Time step (or real time if real_time_measure='Y')
   # 3: Current site index
   # 4-6: Cartesian coordinates (x, y, z)
   # 7: KMC time step counter (next jump at this step)
   # 8-10: Average magnetic moment in polaron cloud (mx, my, mz)

**Example:**

.. code-block:: none

   1     1000    512    5.0000    5.0000    5.0000     1250    0.9850   0.0000   0.0000
   2     1000   1024   10.0000   10.0000   10.0000     1300   -0.9750   0.0000   0.0000
   1     1250    513    6.0000    5.0000    5.0000     1480    0.9700   0.1200   0.0000

**Interpretation:**

- Particle 1 jumped from site 512 → 513 at step 1250
- Average moment in polaron 1's neighborhood: :math:`\mathbf{m} = (0.97, 0.12, 0)` at step 1250
- Next jump for particle 1 scheduled at step 1480

**Printing modes:**

- ``do_prn_kmc='Y'``: Output every ``kmc_step`` time steps
- ``do_prn_kmc='D'``: Adaptive printing only when particle moves (saves disk space)

``struct_kmc.<simid>.out``
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Structural information** (written once at initialization):

.. code-block:: none

   Sorted data from KMC
   ----------------------------------
   Atom=    512    No neigh=     12
                     513   514   ...
         1.5000E-01   1.5000E-01   1.5000E-01   ...
   ----------------------------------
   Atom=    513    No neigh=     12
   ...

- Lists neighbor indices and corresponding energy barriers for each site
- Useful for **debugging** barrier connectivity


Physical Interpretation
-----------------------

Diffusion Coefficients
^^^^^^^^^^^^^^^^^^^^^^

For isotropic barriers, the diffusion coefficient can be estimated from KMC trajectories:

.. math::

   D = \frac{\langle |\mathbf{r}(t) - \mathbf{r}(0)|^2 \rangle}{6t}

Track particle positions from ``kmc_info`` output and calculate mean-square displacement vs. time.

**Arrhenius plot:**

.. math::

   D(T) = D_0 \exp\left(-\frac{E_a}{k_B T}\right)

Plot :math:`\ln D` vs. :math:`1/T` to extract activation energy :math:`E_a` (should match average barrier).

Polaron Cloud Size
^^^^^^^^^^^^^^^^^^

The **average moment** output (columns 8-10 in ``kmc_info``) represents the local magnetization
around the particle. For a magnetic polaron:

- **FM host:** :math:`|\mathbf{m}| \approx 1` → particle weakly perturbs spins
- **AFM host:** :math:`|\mathbf{m}| \ll 1` → strong local spin canting (polaron formation)

Monitor :math:`|\mathbf{m}(t)|` to track polaron stability and distortion dynamics.

Residence Time Distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The distribution of waiting times :math:`\tau_{\text{KMC}}` before jumps follows:

.. math::

   P(\tau) = \Gamma_{\text{tot}} \exp(-\Gamma_{\text{tot}} \tau)

Extract from differences in column 7 (``kmc_time_steps``). Deviations from exponential indicate:

- **Trapped particles:** :math:`\Gamma_{\text{tot}} \to 0` at some sites (high barriers or blocked neighbors)
- **Correlated jumps:** Multiple particles interacting via excluded volume


Examples
--------

Example 1: Single Vacancy Diffusion in FM Ni
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**System:** FCC Ni (``a = 3.52`` Å), single vacancy, nearest-neighbor barriers.

**``inpsd.dat``:**

.. code-block:: none

   do_kmc          Y
   kmc_method      1
   rate0           1.0e12           # Debye frequency (Hz)
   barrfile        barriers.dat
   kmc_posfile     vacancy.dat
   kmc_step        100
   kmc_buff        50
   do_prn_kmc      Y

**``barriers.dat``:**

.. code-block:: none

   # Nearest-neighbor hopping (12 neighbors in FCC)
   512  513   1.0  0.0  0.0   1.05
   512  514  -1.0  0.0  0.0   1.05
   512  515   0.0  1.0  0.0   1.05
   512  516   0.0 -1.0  0.0   1.05
   512  517   0.0  0.0  1.0   1.05
   512  518   0.0  0.0 -1.0   1.05
   ...
   # (Repeat for all sites with consistent 1.05 eV barrier)

**``vacancy.dat``:**

.. code-block:: none

   512  1  17.6  17.6  17.6

**Expected results:**

- Vacancy performs random walk with average rate :math:`\Gamma \approx 10^{12} \exp(-1.05/0.026) \sim 10^{-5}` Hz at 300 K
- Mean-square displacement: :math:`\langle r^2(t) \rangle \sim 6Dt` with :math:`D \sim 10^{-16}` m²/s
- Trajectory output shows isotropic diffusion

**Physical validation:**

Compare extracted :math:`D(T)` with experimental vacancy diffusion in Ni (literature: :math:`E_a \sim 1.0 - 1.1` eV).


Example 2: Magnetic Polaron with Electric Field
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**System:** Skyrmion-hosting chiral magnet, single polaron, applied electric field.

**``inpsd.dat``:**

.. code-block:: none

   do_kmc          Y
   kmc_method      1
   rate0           5.0e11
   barrfile        polaron_barriers.dat
   kmc_posfile     polaron_init.dat
   do_efield       Y
   efield          0.02  0.0  0.0    # Field along x
   kmc_step        50
   kmc_buff        100
   do_prn_kmc      D                 # Adaptive printing

**``polaron_barriers.dat``:**

.. code-block:: none

   # Lower barriers in FM regions, higher in skyrmion cores
   256  257   1.0  0.0  0.0   0.08   # FM region
   512  513   1.0  0.0  0.0   0.25   # Skyrmion boundary
   768  769   1.0  0.0  0.0   0.50   # Skyrmion core

**``polaron_init.dat``:**

.. code-block:: none

   256  2  10.0  10.0  10.0

**Expected results:**

- **Drift velocity:** :math:`v_d \propto E` along :math:`x`-axis (extract from mean displacement :math:`\langle x(t) \rangle`)
- **Skyrmion trapping:** Polaron slows near skyrmion cores (high barriers :math:`\sim 0.5` eV)
- **Local moment distortion:** Columns 8-10 show reduced :math:`|\mathbf{m}|` when polaron enters skyrmion

**Physical insight:** Electric field drives polaron motion, but topological barriers from skyrmions
create **pinning sites**. Useful for modeling skyrmionics applications.


Example 3: Two-Polaron System with Coulomb Repulsion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**System:** Two like-charged polarons with excluded volume interactions.

**``inpsd.dat``:**

.. code-block:: none

   do_kmc          Y
   kmc_method      1
   rate0           1.0e12
   barrfile        two_polaron_barr.dat
   kmc_posfile     two_polarons.dat
   kmc_step        100
   do_prn_kmc      Y

**``two_polarons.dat``:**

.. code-block:: none

   128  2   5.0   5.0   5.0
   384  2  15.0  15.0  15.0

**Algorithm behavior:**

- When polaron 1 attempts to jump to an occupied site → :math:`\Gamma = 0` (blocked)
- **Correlation effects:** Particles avoid each other, leading to reduced effective diffusion
- Output: Track separation :math:`|\mathbf{r}_1(t) - \mathbf{r}_2(t)|` to study pair correlation

**Physical relevance:** Models **polaron-polaron interactions** in manganites or cuprates, where
Coulomb repulsion prevents site double occupancy.


Example 4: Random Alloy with Site-Dependent Barriers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**System:** Fe₅₀Co₅₀ random alloy, barriers depend on local chemical environment.

**``inpsd.dat``:**

.. code-block:: none

   do_kmc          Y
   do_ralloy       1
   kmc_method      1
   rate0           8.0e11
   barrfile        alloy_barriers.dat
   kmc_posfile     alloy_vacancy.dat
   kmc_step        100

**``alloy_barriers.dat``:**

.. code-block:: none

   # Site  Neigh  iChem  jChem  rx  ry  rz  Barrier
   512    513      1      1    1.0  0.0  0.0   0.95   # Fe-Fe
   512    514      1      2    1.0  0.0  0.0   1.10   # Fe-Co
   512    515      2      1    1.0  0.0  0.0   1.10   # Co-Fe
   512    516      2      2    1.0  0.0  0.0   1.20   # Co-Co

**``alloy_vacancy.dat``:**

.. code-block:: none

   512  1  1  0.5  10.0  10.0  10.0   # At Fe site with 50% Fe concentration

**Expected results:**

- **Heterogeneous diffusion:** Vacancy hops faster in Fe-rich regions (lower barriers)
- Extract site-resolved diffusion coefficient from trajectory
- Compare with compositional dependence of diffusion in Fe-Co alloys


Troubleshooting
---------------

Issue: Very Long or Very Short KMC Time Steps
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Symptom:** ``kmc_time_steps`` in output are extremely large (>10⁶) or very small (<10).

**Causes:**

1. **Incorrect ``rate0``:** Check units (must be Hz). Too small → unrealistic long waits.
2. **Unphysical barriers:** Very high :math:`\Delta E` (>2 eV) → exponentially suppressed rates.
3. **Timestep mismatch:** Ensure ``delta_t`` is appropriate for KMC timescale.

**Solution:**

- Verify :math:`\Gamma_{\text{tot}} \sim \nu_0 \exp(-\Delta E / k_B T)` is in reasonable range (10² - 10⁶ Hz)
- Adjust ``rate0`` or barriers to match physical system
- Check that :math:`\delta t \ll \tau_{\text{KMC}}` (otherwise KMC steps skip MD dynamics)

Issue: Particle Appears Trapped
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Symptom:** Particle does not move for entire simulation (``site`` column in output stays constant).

**Causes:**

1. **All neighbors blocked:** In multi-particle system, surrounding sites may be occupied.
2. **No neighbors defined:** ``barrfile`` missing entries for current particle site.
3. **Infinitely high barriers:** All :math:`\Gamma_{i \to j} \to 0`.

**Solution:**

- Check ``struct_kmc`` output to verify neighbor list is correct
- Ensure ``barrfile`` covers all reachable states
- Inspect barrier magnitudes: typical activated processes have :math:`\Delta E = 0.5 - 2.0` eV

Issue: Hamiltonian Update Errors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Symptom:** Simulation crashes or produces unphysical moment configurations after KMC jump.

**Causes:**

1. **Inconsistent neighbor maps:** ``barrfile`` and ``jfile`` use different coordinate conventions.
2. **Tolerance too strict:** Bonding vector comparison fails (``tol = 0.005`` in code).
3. **Missing random alloy data:** Chemical types not properly swapped for ``do_ralloy=1``.

**Solution:**

- Ensure ``barrfile`` and ``jfile`` use same ``posfiletype`` (``'C'`` or ``'D'``)
- Verify coordinate transformations (Cartesian ↔ Direct) are consistent
- For random alloys, check ``kmc_posfile`` includes chemical types

Issue: No Output in ``kmc_info`` File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Symptom:** File ``kmc_info.<simid>.out`` is empty or not created.

**Causes:**

1. **``do_prn_kmc='N'``:** Printing is disabled.
2. **Adaptive mode with no jumps:** ``do_prn_kmc='D'`` only writes when particles move.
3. **Buffer not flushed:** Simulation ended before ``kmc_buff`` entries accumulated.

**Solution:**

- Set ``do_prn_kmc='Y'`` for regular output
- Ensure simulation runs long enough for particles to jump
- Check ``kmc_step`` is not too large (reduce for more frequent measurements)


Performance Considerations
--------------------------

Computational Cost
^^^^^^^^^^^^^^^^^^

KMC adds **minimal overhead** to UppASD simulations:

- **Rate calculation:** :math:`O(N_{\text{part}} \times N_{\text{neigh}})` per time step
- **Hamiltonian update:** :math:`O(N_{\text{neigh}}^2)` when particle jumps

For typical systems (:math:`N_{\text{part}} \sim 1 - 10`, :math:`N_{\text{neigh}} \sim 12`),
KMC overhead is < 1% of total simulation time.

Timescale Separation
^^^^^^^^^^^^^^^^^^^^

KMC is efficient when:

.. math::

   \tau_{\text{KMC}} \gg \delta t_{\text{MD}}

Typical: :math:`\tau_{\text{KMC}} \sim 10^{-8}` s vs. :math:`\delta t_{\text{MD}} \sim 10^{-15}` s
→ :math:`10^7` MD steps between KMC events.

**Caution:** If barriers are very low (:math:`\Delta E < 0.1` eV), KMC becomes inefficient
(frequent jumps every few MD steps). In this limit, use full molecular dynamics instead.

Parallelization
^^^^^^^^^^^^^^^

Current implementation:

- **Setup phases** (``setup_barriers``, rate calculation) are OpenMP parallelized
- **Hamiltonian updates** are serial (due to data dependencies)

For large ``NA_KMC`` (>100 particles), consider subdividing system if particles are spatially separated.


References
----------

**KMC Theory:**

1. A. F. Voter, "Introduction to the Kinetic Monte Carlo Method," *Radiation Effects in Solids* (Springer, 2007)
2. A. P. J. Jansen, *An Introduction to Kinetic Monte Carlo Simulations of Surface Reactions* (Springer, 2012)

**Magnetic Polarons:**

3. E. Dagotto, "Nanoscale Phase Separation and Colossal Magnetoresistance," *Springer Series in Solid-State Sciences* (2003)
4. M. B. Salamon and M. Jaime, "The physics of manganites: Structure and transport," *Rev. Mod. Phys.* **73**, 583 (2001)

**Vacancy Diffusion:**

5. H. Mehrer, *Diffusion in Solids* (Springer, 2007)
6. G. H. Vineyard, "Frequency factors and isotope effects in solid state rate processes," *J. Phys. Chem. Solids* **3**, 121 (1957)

**UppASD Implementation:**

7. J. Chico et al., "KMC module in UppASD for coupled spin-charge dynamics" (in code documentation)


Summary of Workflow
-------------------

1. **Prepare input files:**
   
   - ``inpsd.dat``: Set ``do_kmc='Y'``, specify ``rate0``, file names
   - ``barrfile``: Energy barriers for all hopping connections
   - ``kmc_posfile``: Initial particle positions

2. **Run simulation:**
   
   - KMC particles jump stochastically based on Arrhenius rates
   - Magnetic Hamiltonian updates after each jump
   - Spins evolve via Landau-Lifshitz-Gilbert between KMC events

3. **Analyze output:**
   
   - ``kmc_info``: Extract trajectories, calculate :math:`D(T)`, plot mean-square displacement
   - Monitor polaron cloud magnetization (columns 8-10)
   - Check ``struct_kmc`` for barrier connectivity

4. **Physical validation:**
   
   - Compare :math:`E_a` from Arrhenius plot with known activation energies
   - Verify diffusion coefficients match experimental data
   - Check that polaron formation energy is consistent with model parameters

.. seealso::

   - :doc:`../stimuli/temperature-gradients` for coupling KMC with spatially varying temperatures
   - :doc:`spin-lattice` for full spin-lattice dynamics (alternative to KMC for faster motion)
   - :doc:`../input/simulation` for temperature control and ensemble settings in equilibrium simulations

.. note::
   **Future extensions:** Planned features include ``kmc_method=2`` (dynamic barriers calculated
   from current spin configuration) and ``kmc_method=3`` (adiabatic exchange interpolation). Currently
   only ``kmc_method=1`` (fixed barriers from input) is implemented.

