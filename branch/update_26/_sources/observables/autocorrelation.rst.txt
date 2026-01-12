Autocorrelation functions
=========================

Parameters for time-resolved autocorrelation measurements
---------------------------------------------------------

The UppASD code supports measurement of autocorrelation functions to characterize
the relaxation dynamics and equilibration timescales of magnetic systems.
Autocorrelation measurements provide [Reif1965]_, [Evans2000]_:

- :ref:`Spin autocorrelation functions at multiple waiting times <autocorr-definition>`
- :ref:`Site-resolved and macrocell-averaged spatial correlations <autocorr-spatial>`
- :ref:`Flexible waiting time sampling with external specification <autocorr-waiting-times>`

The autocorrelation measurement scheme allows extraction of dynamical properties
such as relaxation times, correlation decay rates, and aging effects in
non-equilibrium simulations. Results are written to ``autocorr.simid.out`` with
optional spatial binning into macrocells for coarse-grained analysis.

.. tip::

   **Input keywords:** For a complete list of all autocorrelation measurement
   parameters (``do_autocorr``, ``acfile``, ``ac_step``, etc.), see the
   comprehensive reference in :doc:`../input/observables`.

.. _autocorr-definition:

-------------------------------------------------
Autocorrelation function definition
-------------------------------------------------

The spin autocorrelation function measures the persistence of magnetic moment
orientations as the system evolves. For a system with :math:`N_{\text{atom}}` atoms,
the global autocorrelation is defined as [Reif1965]_:

.. math::

   C(t_w, t) = \frac{1}{N_{\text{atom}}} \sum_{i=1}^{N_{\text{atom}}}
   \hat{\mathbf{m}}_i(t_w) \cdot \hat{\mathbf{m}}_i(t_w + t),

where:

- :math:`t_w` is the **waiting time** (initial measurement time)
- :math:`t` is the **time separation** (evolution time)
- :math:`\hat{\mathbf{m}}_i` is the unit magnetic moment vector at atom :math:`i`

The autocorrelation probes how well the system "remembers" its earlier state.
Typical behavior:

- :math:`C(t_w, 0) = 1` at :math:`t = 0` (perfect initial correlation)
- :math:`C(t_w, t) \to 0` as :math:`t \to \infty` (decorrelation in equilibrium)
- :math:`C(t_w, t) > 0` during thermalization (slow relaxation)
- :math:`C(t_w, t)` exhibits subdiffusive or glassy decay in disordered systems

The implementation in UppASD computes correlations only over ensemble 1
(first replica) and averages over all atoms to improve statistics.

.. _autocorr-spatial:

-------------------------------------------------
Spatial binning and macrocell averaging
-------------------------------------------------

For spatially-resolved analysis, autocorrelation can be computed separately within
each **macrocell** (computational domain subdivision). 

When ``do_macro_cells Y`` is set in the input, the code:

1. Partitions the system into equal-sized cubic macrocells
2. Computes autocorrelation separately within each macrocell
3. Outputs per-macrocell correlations to files ``autocorr_loc.XXXX.simid.out``

This approach is useful for:

- Studying domain wall motion and propagation
- Detecting spatial phase separation
- Analyzing non-uniform relaxation in thin films or defect-containing systems
- Validating assumptions of spatial homogeneity

.. _autocorr-waiting-times:

-------------------------------------------------
Waiting time specification
-------------------------------------------------

The autocorrelation function is computed at :math:`N_{\text{wait}}` different
waiting times :math:`t_{w,i}`. These times can be:

- Equally spaced (simple case)
- Logarithmically spaced (to cover multiple decades of timescale)
- Irregularly spaced (to focus on regions of interest)

All waiting times are specified externally in a dedicated text file (``acfile``)
containing one waiting time per line (in simulation steps). This design allows
reuse across different simulations and avoids hardcoding timescales.

-------------------------------------------------
Basic autocorrelation keywords
-------------------------------------------------

``do_autocorr``
   Enable autocorrelation function sampling (Y=yes, *N=no*). When enabled, the code measures spin
   autocorrelation at specified waiting times. Waiting times must be provided in an external file
   specified by ``acfile``. Output is written to ``autocorr.simid.out``.

``acfile``
   Name of external file containing waiting times for autocorrelation measurements. Each line should
   contain one integer (simulation step number). At least 2 waiting times should be provided to
   establish a meaningful correlation decay. The file is read once at initialization.

``ac_step``
   Sampling interval in simulation steps (integer). Determines how frequently autocorrelation samples
   are accumulated in the measurement buffer. Default: 100. Must be set if ``do_autocorr Y``. Smaller
   values increase temporal resolution but increase memory and I/O overhead.

``ac_buff``
   Buffer size for autocorrelation samples (integer). Number of samples to accumulate before writing
   to file. Default: 10. Larger buffers reduce I/O overhead but increase memory. Total output samples
   = ``ac_buff`` :math:`\times` (number of flush cycles).

.. important::

   **Waiting time file format**: The external file specified by ``acfile`` should
   contain one or more waiting times (simulation step numbers), one per line.
   Example file ``twfile`` with waiting times at steps 0, 100, 1000, 10000::

       0
       100
       1000
       10000

   The number of waiting times (``nspinwait``) is inferred from the file length
   automatically.

-------------------------------------------------
Workflow for autocorrelation measurements
-------------------------------------------------

The measurement process follows these stages:

**1. Initialization phase** (``initmag`` = new simulation)

   - At step :math:`mstep = 0`, the code reads the current moment configuration
   - For each atom :math:`i` and waiting time :math:`w_i`, stores
     :math:`\hat{\mathbf{m}}_i(t = t_{w,i})`
   - Writes initial moments to auxiliary files ``spin.simid.J.out`` (one per
     waiting time :math:`J = 1, \ldots, N_{\text{wait}}`)

**2. Sampling phase** (during simulation)

   - At steps matching the waiting times (``spinwaitt(j)``), the current moment
     configuration is stored in ``spin.simid.J.out``
   - At regular intervals (``ac_step``), autocorrelations are computed:

     .. math::

        \text{corr}(j) = \frac{1}{N_{\text{atom}}} \sum_{i=1}^{N_{\text{atom}}}
        \hat{\mathbf{m}}_i(t_{w,j}) \cdot \hat{\mathbf{m}}_i(\text{current step})

   - Computed correlations are buffered in memory

**3. Output phase** (every ``ac_buff`` samples)

   - Buffered autocorrelations are written to ``autocorr.simid.out``
   - One row per sampling point
   - Columns: iteration (or time) + correlations for each waiting time
   - Optional spatial data written to ``autocorr_loc.XXXX.simid.out`` if
     ``do_macro_cells Y``

-------------------------------------------------
Output file formats
-------------------------------------------------

**Global autocorrelation file: ``autocorr.simid.out``**

Columns:

- Column 1: Simulation step or real time (seconds if ``real_time_measure Y``)
- Columns 2–(1+N_wait): :math:`C(t_{w,j}, t)` for each waiting time :math:`j = 1, \ldots, N_{\text{wait}}`

Example output (3 waiting times: 0, 100, 1000 steps)::

       0        1.000000e+00  1.000000e+00  1.000000e+00
     100        9.834562e-01  1.000000e+00  1.000000e+00
    1000        8.245671e-01  9.921342e-01  1.000000e+00
   10000        6.123456e-01  8.756234e-01  9.876543e-01

**Spatial macrocell files: ``autocorr_loc.XXXX.simid.out``**

One file per macrocell (XXXX = 4-digit macrocell index). Same format as global
file but autocorrelation normalized by ``max_num_atom_macro_cell`` (maximum
number of atoms in any macrocell).

-------------------------------------------------
Practical considerations and timescale selection
-------------------------------------------------

**Choosing waiting times:**

The choice of waiting times depends on the timescale of interest [Evans2000]_:

- **Fast dynamics**: Use waiting times spaced by :math:`\sim 10`–:math:`100` steps
- **Aging/slow dynamics**: Use logarithmically spaced times: :math:`t_w \propto 10^k`
- **Equilibrium studies**: Use early waiting time (:math:`t_w \approx 0`) and track decay
- **Non-equilibrium**: Use multiple waiting times to identify aging exponents

**Memory and I/O overhead:**

- Auxiliary storage for moment snapshots: :math:`N_{\text{wait}} \times N_{\text{atom}} \times 8` bytes per snapshot
- Files ``spin.simid.J.out`` are overwritten at each waiting time (not accumulated)
- Autocorrelation buffer: :math:`N_{\text{wait}} \times \text{ac\_buff} \times 8` bytes

For large :math:`N_{\text{wait}}`, ensure adequate disk space for auxiliary files.

**Performance notes:**

- Autocorrelation sampling only uses ensemble 1 (single replica)
- Spatial binning (macrocells) incurs modest computational overhead
- For very long simulations, consider reducing ``ac_buff`` to lower memory pressure

-------------------------------------------------
Interactions with other measurement modes
-------------------------------------------------

Autocorrelation measurements are independent but compatible with:

- **Averages** (``do_avrg``): Can use overlapping ``avrg_step`` and ``ac_step``
  for simultaneous measurements at same frequencies
- **Cumulants** (``do_cumu``): Run separately; no resource conflicts
- **Trajectory output** (``do_tottraj``): Complementary; trajectory files
  provide per-atom moment evolution if full detail needed

.. note::

   Unlike averages and cumulants which average over ensembles, autocorrelation
   currently samples only ensemble 1. Future versions may support ensemble
   averaging of autocorrelations for improved statistics.

-------------------------------------------------
Restart and continuation
-------------------------------------------------

When restarting a simulation (``initmag`` :math:`\neq 1`):

- The code reads previously saved ``spin.simid.J.out`` files
- Determines the latest completed waiting time index ``n0spinwait``
- Continues autocorrelation computation from that point forward
- Existing autocorrelation output is appended to ``autocorr.simid.out``

This allows extending autocorrelation measurements across multiple simulation
segments without loss of continuity.

-------------------------------------------------
Example inpsd.dat snippet
-------------------------------------------------

**Basic autocorrelation measurement:**

.. code-block:: text

   do_autocorr   Y
   acfile        twfile          # External file with waiting times
   ac_step       100             # Sample autocorr every 100 steps
   ac_buff       10              # Buffer 10 samples before writing

**Fine-grained temporal sampling:**

.. code-block:: text

   do_autocorr   Y
   acfile        twfile_fine     # File with densely-spaced waiting times
   ac_step       50              # High-frequency sampling
   ac_buff       20

**With spatial macrocell analysis:**

.. code-block:: text

   do_autocorr       Y
   acfile            twfile
   ac_step           100
   ac_buff           10

   # Enable spatial binning
   do_macro_cells    Y
   macro_cell_size   10.0         # Macrocell edge length (Angstrom)

**Aging study (logarithmic waiting times):**

.. code-block:: text

   do_autocorr   Y
   acfile        twfile_aging    # Logarithmically spaced: 1, 10, 100, 1000, ...
   ac_step       10              # Fine temporal resolution
   ac_buff       50              # Large buffer for stability

**Example ``twfile`` for logarithmic sampling (10 decades):**

.. code-block:: text

   0
   1
   10
   100
   1000
   10000
   100000

-------------------------------------------------
Relationship to other dynamics measures
-------------------------------------------------

Autocorrelation complements other dynamical quantities measured in UppASD:

- **Mean squared displacement** (if implemented): Tracks position changes; autocorrelation tracks orientation
- **Energy fluctuations**: Measure thermodynamic response; autocorrelation measures relaxation
- **Correlation functions** (``do_sc``): Measure spatial structure at fixed times; autocorrelation measures temporal decay

For a complete dynamical characterization, consider enabling multiple measurement
modes simultaneously.

-------------------------------------------------
Related keywords and cross-references
-------------------------------------------------

- ``delta_t``: Timestep for real-time measurements (see :doc:`../input/simulation`)
- ``real_time_measure``: Output physical time instead of steps (see :doc:`../input/simulation`)
- ``do_macro_cells``: Enable spatial macrocell subdivision (see :doc:`../input/system`)
- ``macro_cell_size``: Size of macrocells in physical units (see :doc:`../input/system`)
- ``do_avrg``, ``avrg_step``: Companion average magnetization measurements (see :doc:`averages`)

Output files use the simulation identifier ``simid`` (see :doc:`../input/system`).


References
----------

See the centralized :doc:`../references` for full bibliographic entries:

- [Reif1965]_ - Statistical and thermal physics fundamentals
- [Evans2000]_ - Statistical mechanics of nonequilibrium systems

