Replica exchange (parallel tempering)
======================================

Parameters for Replica Exchange / Parallel Tempering Monte Carlo
------------------------------------------------------------------

The UppASD code implements the **Replica Exchange** (also known as **Parallel Tempering**)
Monte Carlo method [Swendsen1986]_, [Geyer1995]_, [Hukushima1996]_ for enhanced sampling
of complex energy landscapes in magnetic systems. This method is particularly effective for:

- **Overcoming energy barriers** that trap conventional Monte Carlo in metastable states
- **Accelerating equilibration** in systems with slow relaxation dynamics
- **Sampling rare configurations** in systems with rough free energy landscapes
- **Determining phase transition temperatures** and critical behavior
- **Finding low-energy magnetic structures** in frustrated or competing-interaction systems

The replica exchange method runs multiple **replicas** (independent copies) of the system
simultaneously at different temperatures. Periodically, configurations at neighboring
temperatures are exchanged according to a Metropolis-like criterion, allowing high-temperature
replicas to escape local minima and low-temperature replicas to sample rare states.

In UppASD, replica exchange is enabled by setting the initial phase mode to ``P``
(``ip_mode P``), with one temperature per replica specified via the ``ip_nphase`` keyword.

--------------------------------------------------
Replica exchange algorithm
--------------------------------------------------

The replica exchange method maintains :math:`N_{\text{rep}}` independent replicas of the
system at temperatures :math:`T_1 < T_2 < \cdots < T_{N_{\text{rep}}}`. Each replica
evolves independently according to the Metropolis or heat-bath Monte Carlo algorithm
at its assigned temperature.

At regular intervals (every ``pt_step`` Monte Carlo steps), neighboring replicas at
temperatures :math:`T_i` and :math:`T_{i+1}` attempt to **exchange configurations**.
The exchange is accepted with probability:

.. math::

   P_{\text{exchange}}(i, i+1) = \min\left(1, \exp\left[(\beta_{i+1} - \beta_i)(E_i - E_{i+1})\right]\right)

where :math:`\beta_i = 1/(k_B T_i)` is the inverse temperature and :math:`E_i` is the
total energy of replica :math:`i`. This acceptance criterion ensures detailed balance
and preserves the correct canonical distribution at each temperature.

The exchange attempts are performed in two sweeps:

1. **Even sweep**: Attempt exchanges between replicas :math:`(1,2), (3,4), \ldots`
2. **Odd sweep**: Attempt exchanges between replicas :math:`(2,3), (4,5), \ldots`

This two-sweep strategy ensures that all neighboring pairs have an opportunity to
exchange within one cycle.

--------------------------------------------------
Temperature ladder selection
--------------------------------------------------

The choice of temperature ladder :math:`\{T_1, T_2, \ldots, T_{N_{\text{rep}}}\}` is
critical for efficient sampling. The **acceptance rate** for exchanges between
neighboring temperatures depends on the energy fluctuations and the temperature spacing.

A commonly used heuristic is to aim for an acceptance rate of 20–40% between neighboring
replicas. For systems with smooth energy landscapes, a **geometric temperature ladder**
often works well:

.. math::

   T_{i+1} = T_i \times r, \quad r = \left(\frac{T_{\text{max}}}{T_{\text{min}}}\right)^{1/(N_{\text{rep}}-1)}

For systems with rough energy landscapes or first-order transitions, **adaptive schemes**
or empirical tuning may be required.

In UppASD, the temperature ladder is specified via the ``ip_nphase`` keyword, with each
phase corresponding to one replica temperature.

--------------------------------------------------
Equilibration and measurement
--------------------------------------------------

The replica exchange simulation in UppASD is structured as an **initial phase only**
(``ip_mode P``). This phase serves to:

- Equilibrate all replicas at their respective temperatures
- Allow sufficient exchange attempts to decorrelate the configurations
- Identify the replica with the **lowest energy** for subsequent measurements

At the end of the initial phase, the replica with the minimum energy is automatically
selected and its configuration is copied to the main simulation arrays (``emom``, ``mmom``).
This configuration can then be used as the starting point for a subsequent measurement
phase (e.g., spin dynamics, further Monte Carlo sampling, or correlation measurements).

The exchange attempts are **disabled during the final 10%** of the initial phase to
ensure that each replica reaches thermal equilibrium at its assigned temperature before
the final selection.

--------------------------------------------------
Replica exchange keywords
--------------------------------------------------

Replica exchange keywords
--------------------------------------------------

ip_mode
   Initial phase mode. Set to **P** to enable Replica Exchange / Parallel Tempering. This runs multiple replicas at different temperatures with periodic configuration exchanges. Default: **S** (spin dynamics). See also: M=Metropolis MC, H=Heat bath MC, W=Wang-Landau.

ip_nphase
   Number of replicas (temperature points) for replica exchange. Each replica corresponds to one temperature. The temperatures are specified in subsequent lines, one per phase. For PT, set this to the desired number of replicas :math:`N_{\text{rep}}`. Default: **1** (single phase, no PT).

ip_nstep
   Number of Monte Carlo steps per replica in each phase line. This is the number of MC sweeps performed before checking for replica exchanges. The format is: ``ip_nstep  temp  timestep  damping``, repeated for each phase. For PT, all replicas should use the same number of steps in the single phase block. Default: **0** (no steps).

pt_step
   Interval (in Monte Carlo steps) between replica exchange attempts. Every ``pt_step`` steps, the code calculates energies for all replicas and attempts configuration exchanges between neighboring temperatures. Smaller values increase exchange frequency but add computational overhead. Typical range: 50–500 steps. Default: **100**.

**Note:** The temperature for each replica is specified implicitly via the phase structure.
Each line in the ``ip_nphase`` block corresponds to one replica and includes the temperature
for that replica. The first column after ``ip_nphase`` is the number of MC steps, followed
by the temperature.

--------------------------------------------------
Simulation mode and MC algorithm
--------------------------------------------------

Replica exchange in UppASD currently supports only the **Heat Bath Monte Carlo** algorithm
(``mode H``). The initial phase automatically uses heat-bath updates for all replicas,
regardless of the global ``mode`` setting.

The general structure for a replica exchange simulation is:

.. code-block:: text

   simid    PT_simulation

   !! Initial phase: Replica Exchange / Parallel Tempering
   ip_mode    P
   ip_nphase  N_rep
   nstep_1   temp_1   timestep   damping
   nstep_2   temp_2   timestep   damping
   ...
   nstep_N   temp_N   timestep   damping

   !! Optional: subsequent measurement phase with selected replica
   mode       M     ! or S, H, etc.
   nstep      10000
   temp       temp_measurement

The ``ip_nphase`` keyword specifies the number of replicas, and each subsequent line
defines the number of MC steps, temperature, timestep, and damping for that replica.

--------------------------------------------------
Output files
--------------------------------------------------

The replica exchange simulation produces the following output files:

- ``ptinitial.<simid>.out``: Summary statistics for the initial phase, including average magnetization, Binder cumulant, and susceptibility for each replica (currently minimal output).

- ``pt_restart.<simid>.out``: Detailed restart file containing the final magnetic configuration for all replicas, with temperatures and energies. This file can be used to inspect the state of each replica at the end of the initial phase.

- Standard restart file (``restart.<simid>.out``): Contains the magnetic configuration of the **lowest-energy replica** selected at the end of the initial phase. This configuration is automatically copied to the main simulation arrays for subsequent use.

The code also prints progress information to stdout, including:

- Exchange hit rate (percentage of successful exchanges)
- Number of exchange attempts
- Final energy and magnetization for each replica

--------------------------------------------------
Performance and parallelization
--------------------------------------------------

The replica exchange implementation in UppASD is designed for **multi-replica parallelism**
but currently runs replicas **sequentially** within a single process. Each replica performs
its Monte Carlo evolution independently, and energies are computed for all replicas before
exchange attempts.

For large systems or many replicas, the dominant computational cost is:

1. **Monte Carlo sweeps**: :math:`\mathcal{O}(N_{\text{atoms}} \times N_{\text{rep}} \times N_{\text{steps}})`
2. **Energy evaluations**: :math:`\mathcal{O}(N_{\text{atoms}} \times N_{\text{rep}})` per exchange attempt
3. **Exchange attempts**: :math:`\mathcal{O}(N_{\text{rep}})` per ``pt_step`` interval

The exchange interval ``pt_step`` should be large enough to allow decorrelation between
exchanges but small enough to ensure efficient mixing. A typical rule of thumb is to
set ``pt_step`` to 50–200 MC sweeps.

--------------------------------------------------------------------------------
Example: Replica exchange for spin glass equilibration
--------------------------------------------------------------------------------

Minimal ``inpsd.dat`` snippet for replica exchange with 8 replicas spanning temperatures
from 10 K to 500 K:

.. code-block:: text

   simid    SpinGlass_PT

   !! Replica Exchange initial phase
   ip_mode    P
   ip_nphase  8
   50000   10.0    1.0e-15   0.1
   50000   20.0    1.0e-15   0.1
   50000   40.0    1.0e-15   0.1
   50000   80.0    1.0e-15   0.1
   50000   160.0   1.0e-15   0.1
   50000   250.0   1.0e-15   0.1
   50000   350.0   1.0e-15   0.1
   50000   500.0   1.0e-15   0.1

   !! Exchange interval
   pt_step    100

   !! Optional: measurement phase with selected configuration
   mode       M
   mcnstep    100000
   temp       10.0

In this example:

- 8 replicas run at geometrically spaced temperatures
- Each replica performs 50,000 MC sweeps
- Exchanges are attempted every 100 steps
- The lowest-energy replica is selected for a subsequent MC measurement phase at 10 K

**Geometric temperature ladder** with :math:`r = (500/10)^{1/7} \approx 1.93`.

--------------------------------------------------
Example: Replica exchange for skyrmion nucleation
--------------------------------------------------

Example for finding low-energy skyrmion configurations in a frustrated magnet:

.. code-block:: text

   simid    Skyrmion_PT

   !! Replica Exchange initial phase
   ip_mode    P
   ip_nphase  6
   100000   1.0     1.0e-16   0.05
   100000   5.0     1.0e-16   0.05
   100000   20.0    1.0e-16   0.05
   100000   50.0    1.0e-16   0.05
   100000   100.0   1.0e-16   0.05
   100000   200.0   1.0e-16   0.05

   pt_step    50

   !! Measurement: spin dynamics with selected configuration
   mode       S
   SDEalgh    1
   nstep      500000
   timestep   1.0e-16
   temp       1.0
   damping    0.05

Here, the replica exchange runs at 6 temperatures from 1 K to 200 K. The low-temperature
replicas explore metastable skyrmion states, while high-temperature replicas allow rapid
exploration of configuration space. After equilibration, the lowest-energy configuration
(likely a skyrmion) is used as the initial state for spin dynamics at 1 K.

--------------------------------------------------
Details: Exchange acceptance rate
--------------------------------------------------

The acceptance rate for replica exchanges is printed during the simulation. A typical
target is **20–40%** for efficient sampling. If the acceptance rate is:

- **Too high (>60%)**: Temperatures are too close; reduce the number of replicas or increase the temperature spacing.
- **Too low (<10%)**: Temperatures are too far apart; add more intermediate replicas or adjust the temperature ladder.

The code prints the hit rate every 10% of the simulation:

.. code-block:: text

   IP PT  10% done.  No. trials:   4500  Hitrate:  32.4 %
   IP PT  20% done.  No. trials:   9000  Hitrate:  31.8 %
   ...

Monitor these statistics to assess the efficiency of the temperature ladder.

--------------------------------------------------
Details: Replica selection criterion
--------------------------------------------------

At the end of the initial phase, UppASD selects the replica with the **minimum energy**:

.. math::

   i_{\text{selected}} = \arg\min_i E_i

The magnetic configuration (``emom``, ``mmom``) of this replica is copied to the main
simulation arrays. This ensures that subsequent measurement phases start from a
low-energy state, which is particularly useful for:

- Finding ground-state configurations
- Equilibrating systems with slow dynamics
- Initializing spin dynamics simulations near energy minima

If a specific temperature is desired for measurements (rather than the lowest energy),
the restart file ``pt_restart.<simid>.out`` contains all replica configurations and
can be manually inspected.


--------------------------------------------------
References and further reading
--------------------------------------------------

The replica exchange method (parallel tempering) is discussed in classic and review literature [Swendsen1986]_, [Geyer1995]_, [Hukushima1996]_. For sampling strategies and practical guidance see the overview by Earl and Deem [Earl2005]_ and the feedback-optimization study by Katzgraber *et al.* [Katzgraber2006]_.

See also
--------

- :doc:`monte-carlo` (general Monte Carlo keywords)
- :doc:`wang-landau` (Wang-Landau sampling)
- :doc:`../input/system` (system setup and temperature control)


References
----------

See the centralized :doc:`../references` for full bibliographic entries:

- [Swendsen1986]_ - Replica Monte Carlo simulation of spin-glasses
- [Geyer1995]_ - Markov chain Monte Carlo maximum likelihood
- [Hukushima1996]_ - Exchange Monte Carlo method and application to spin glass simulations
- [Earl2005]_ - Parallel tempering theory applications and perspectives
- [Katzgraber2006]_ - Feedback-optimized parallel tempering Monte Carlo

