Wang-Landau sampling
=====================

Parameters for Wang-Landau Monte Carlo sampling
-------------------------------------------------

The UppASD code supports the Wang-Landau Monte Carlo algorithm [WangLandau2001]_
for calculating the density of states :math:`g(E)` of spin Hamiltonians.
The Wang-Landau method is a powerful technique that allows for:

- Direct calculation of the **density of states** :math:`g(E)` over a wide energy range
- Extraction of **thermodynamic properties** (free energy, entropy, heat capacity) from the DOS via integration
- **Efficient sampling** of low-probability energy regions that are poorly accessed by canonical ensemble simulations
- **Flat energy histograms** through adaptive biasing of energy visitation

The method is particularly useful for systems with:

- First-order phase transitions
- Rugged energy landscapes
- Strong metastability
- Low-temperature phases inaccessible by standard Monte Carlo

The Wang-Landau algorithm is invoked by setting the simulation mode to ``W``
in the measurement phase of the simulation.

--------------------------------------------------
Wang-Landau algorithm and entropy representation
--------------------------------------------------

In the canonical Wang-Landau approach, a random walk in energy space is performed
with a visitation probability proportional to :math:`1/g(E)`. The density of
states :math:`g(E)` is iteratively refined by multiplying it by a modification
factor :math:`f > 1` each time an energy :math:`E` is visited:

.. math::

   g(E) \rightarrow g(E) \times f

Initially, :math:`f = f_0 = e`, and a histogram :math:`H(E)` records the number
of visits to each energy. Once the histogram becomes approximately **flat**
(within a tolerance, typically 70–80%), the histogram is reset and :math:`f`
is reduced:

.. math::

   f \rightarrow \sqrt{f}

The procedure is repeated until :math:`f` approaches unity (typically 
:math:`f_{\text{final}} \approx 1.000001`), at which point the density of
states is converged.

For numerical stability, UppASD works with the **logarithm** of the density
of states:

.. math::

   S(E) = \ln g(E)

During the random walk, acceptance of a Monte Carlo trial move from energy
:math:`E_0` to :math:`E_1` is determined by the modified acceptance probability:

.. math::

   P_{\text{accept}} = \min\left(1, \exp\left[S(E_0) - S(E_1)\right]\right)

This ensures that all energies are visited with equal probability when the
simulation converges.

--------------------------------------------------
Energy histogram and discretization
--------------------------------------------------

The energy range :math:`[E_{\text{min}}, E_{\text{max}}]` is divided into
``wl_nhist`` discrete bins of width :math:`\Delta E`. The bin centers are:

.. math::

   E_i = E_{\text{min}} + i \cdot \Delta E, \quad i = 1, 2, \ldots, N_{\text{hist}}

Each bin accumulates a histogram count :math:`H(E_i)` and an associated
DOS estimate :math:`g(E_i)`.

--------------------------------------------------
Broadening kernel
--------------------------------------------------

To reduce correlations and improve convergence, a **Gaussian broadening kernel**
is applied when updating the DOS. The width of the kernel is controlled by
the parameter ``wl_sigma``, which specifies the broadening as a fraction of
the total energy window:

.. math::

   \sigma_E = \text{wl_sigma} \times N_{\text{hist}} \times \Delta E

This broadening helps smooth the DOS and ensures that nearby energy bins are
updated simultaneously, reducing statistical noise.

--------------------------------------------------
Automatic energy window determination
--------------------------------------------------

By default, UppASD performs an **initial phase** (``ip_mode W``) to automatically
determine the accessible energy range :math:`[E_{\text{min}}, E_{\text{max}}]`.
This phase involves:

1. Starting from a random spin configuration
2. Performing a one-sided minimization in the positive energy direction to find :math:`E_{\text{max}}`
3. Repeating the minimization in the negative energy direction to find :math:`E_{\text{min}}`

The final energy window is then adjusted by cutoff factors ``wl_lcut`` and ``wl_hcut``:

.. math::

   E_{\text{min}}^{\text{final}} &= \text{wl_lcut} \times E_{\text{min}} \\
   E_{\text{max}}^{\text{final}} &= \text{wl_hcut} \times E_{\text{max}}

If ``wl_emin`` and ``wl_emax`` are explicitly provided in the input file, the
automatic energy window determination is skipped.

--------------------------------------------------
Thermodynamic integration
--------------------------------------------------

Once the density of states :math:`g(E)` is obtained, thermodynamic quantities
can be computed via integration. The partition function at temperature :math:`T` is:

.. math::

   Z(T) = \sum_E g(E) e^{-E/(k_B T)}

From this, the free energy, internal energy, entropy, and heat capacity are
derived as:

.. math::

   F(T) &= -k_B T \ln Z(T) \\
   U(T) &= \frac{\sum_E E \, g(E) e^{-E/(k_B T)}}{Z(T)} \\
   S(T) &= \frac{U(T) - F(T)}{T} \\
   C_V(T) &= \frac{\partial U}{\partial T}

UppASD automatically performs this integration over the temperature range
:math:`[T_{\text{min}}, T_{\text{max}}]` and outputs the results to file.

--------------------------------------------------
Wang-Landau keywords
--------------------------------------------------

wl_emin
   Minimum energy (in mRy) for the Wang–Landau histogram. If unset, this value is determined automatically during the initial phase. If both ``wl_emin`` and ``wl_emax`` are provided, the initial phase is skipped. Default: *automatic determination*.

wl_emax
   Maximum energy (in mRy) for the Wang–Landau histogram. If unset, this value is determined automatically during the initial phase. Default: *automatic determination*.

wl_nhist
   Number of bins in the Wang–Landau energy histogram. If unset, the number is chosen automatically such that the bin width is approximately 1 mRy. The automatic choice is :math:`N_{\\text{hist}} = N_{\\text{threads}} \\times |E_{\\text{max}} - E_{\\text{min}}|`, where :math:`N_{\\text{threads}}` is the number of OpenMP threads. Default: *automatic* (**0**).

wl_nloop
   Number of Wang–Landau histogram loops (refinement iterations). For each loop, the modification factor :math:`f` is reduced to :math:`\\sqrt{f}`. Starting from :math:`f_0 = e`, after 10 loops :math:`f \\approx 1.003`, after 20 loops :math:`f \\approx 1.000001`. The error in the DOS scales as :math:`\\log f`, so more loops give higher accuracy. Default: **20**.

wl_lcut
   Lower cutoff factor for the automatic energy window. The automatically determined :math:`E_{\\text{min}}` is multiplied by this factor to set the final lower bound. Values less than 1.0 reduce the energy window. Default: **0.975**.

wl_hcut
   Higher cutoff factor for the automatic energy window. The automatically determined :math:`E_{\\text{max}}` is multiplied by this factor to set the final upper bound. Values less than 1.0 reduce the energy window. Default: **0.950**.

wl_sigma
   Broadening width (in units of the energy window) for the Gaussian kernel applied during DOS updates. The broadening in energy units is :math:`\\sigma_E = \\text{wl_sigma} \\times N_{\\text{hist}} \\times \\Delta E`. Larger values smooth the DOS more aggressively but may reduce fine structure. Default: **0.005**.

wl_gfac
   Global prefactor applied to the final density of states. Used for rescaling or normalization. Default: **1.0**.

wl_minT
   Minimum temperature (in K) for thermodynamic integration. The free energy, entropy, internal energy, and heat capacity are computed over the range :math:`[T_{\\text{min}}, T_{\\text{max}}]`. Default: **1.0** K.

wl_maxT
   Maximum temperature (in K) for thermodynamic integration. Default: **1500.0** K.

--------------------------------------------------
Simulation mode and phase structure
--------------------------------------------------

Wang-Landau simulations are typically split into two phases:

1. **Initial phase** (``ip_mode W``): Automatically determines the accessible energy range by performing directional minimizations. This phase outputs :math:`E_{\text{min}}` and :math:`E_{\text{max}}` and sets up the histogram binning.

2. **Measurement phase** (``mode W``): Performs the Wang-Landau sampling to build the density of states. The histogram is reset and the modification factor :math:`f` is reduced once flatness criteria are met.

The number of Monte Carlo steps in each phase is controlled by the standard
``ipmcnstep`` (initial phase) and ``mcnstep`` (measurement phase) keywords.

--------------------------------------------------
Output files
--------------------------------------------------

The Wang-Landau simulation generates several output files:

- ``wlhistogram.<simid>.out``: Histogram data for each refinement loop, showing energy bins, DOS, histogram counts, and average magnetization per bin. This file is appended after each histogram reset.

- ``wlfinal.<simid>.out``: Final converged density of states :math:`g(E)` vs energy :math:`E`, along with histogram counts and average magnetization.

- ``wl_integration.<simid>.out``: Thermodynamic properties (free energy, internal energy, entropy, heat capacity) as a function of temperature, obtained by integrating the DOS.

All energy values are in **milli-Rydberg (mRy)** unless atomic units are enabled
(``aunits Y``), in which case energies are in **Joules**.

--------------------------------------------------
Notes on parallelization and performance
--------------------------------------------------

The Wang-Landau implementation in UppASD is **OpenMP-parallelized**. Multiple
threads simultaneously perform independent random walks in energy space, each
updating a shared DOS array. This parallelization significantly reduces
computational time but requires careful synchronization to avoid race conditions.

The automatic histogram binning (when ``wl_nhist = 0``) scales with the number
of threads to ensure sufficient bins for accurate parallel sampling:

.. math::

   N_{\text{hist}} = N_{\text{threads}} \times |E_{\text{max}} - E_{\text{min}}|

For large systems or tight convergence requirements, it may be beneficial to
manually set ``wl_nhist`` to a larger value.

--------------------------------------------------
Example: Wang-Landau simulation of a spin system
--------------------------------------------------

Minimal ``inpsd.dat`` snippet for a Wang-Landau simulation with automatic energy
window determination:

.. code-block:: text

   simid    WL_test

   !! Initial phase: determine energy window
   ip_mode    W
   ip_nphase  1
   10000   0.00001   1.0e-15   0.1

   !! Measurement phase: Wang-Landau sampling
   mode       W
   mcnstep    500000
   temp       100.0

   !! Wang-Landau parameters
   wl_nloop   20
   wl_lcut    0.975
   wl_hcut    0.950
   wl_sigma   0.005
   wl_nhist   0
   wl_minT    1.0
   wl_maxT    1500.0

Example with **explicit energy window** (skips initial phase):

.. code-block:: text

   simid    WL_explicit

   !! Measurement phase: Wang-Landau sampling
   mode       W
   mcnstep    500000
   temp       100.0

   !! Explicit energy range
   wl_emin    -500.0
   wl_emax     100.0
   wl_nhist    600
   wl_nloop    25
   wl_sigma    0.003

In this case, the energy window :math:`[-500, 100]` mRy is divided into 600 bins,
and 25 refinement loops are performed.

--------------------------------------------------
Details: Flatness criterion
--------------------------------------------------

The flatness of the histogram :math:`H(E)` is checked every 10,000 Monte Carlo
steps. The flatness is defined as:

.. math::

   \text{flatness} = \frac{\min H(E)}{\langle H(E) \rangle} \times 100\%

where the minimum and average are taken over the **interior bins** (excluding
the first and last bins, which may be undersampled). Once the flatness exceeds
a threshold (typically 70% plus a small increment that grows with the loop count),
the histogram is reset and :math:`f` is reduced.

The acceptance rate is also monitored and used to adaptively adjust the
trial move step size, improving sampling efficiency.

--------------------------------------------------
Details: Magnetization histogram
--------------------------------------------------

In addition to the energy histogram, UppASD also accumulates a **magnetization
histogram** :math:`\mathbf{M}(E)`, which records the average magnetization
vector in each energy bin. This allows for the calculation of magnetic
susceptibilities and order parameters as a function of energy and temperature.

The average magnetization per energy bin is written to the output files and
can be used to identify magnetic phases and transitions.

--------------------------------------------------
References and further reading
--------------------------------------------------


The Wang-Landau algorithm was introduced by Wang and Landau [WangLandau2001]_.

See also
--------

- :doc:`input-keywords-montecarlo` (general Monte Carlo keywords)
- :doc:`input-keywords-system` (system setup and lattice parameters)


References
----------

See the centralized :doc:`references` for full bibliographic entries:

- [WangLandau2001]_ - Efficient multiple-range random walk algorithm to calculate density of states

