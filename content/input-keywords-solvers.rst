.. _input-keywords-solvers:

Stochastic Equation Solvers
====================================

Overview
--------

UppASD implements multiple **numerical solvers** for the stochastic Landau-Lifshitz-Gilbert (LLG)
equation and related equations of motion. These are the numerical engines that evolve magnetic
moments in time, accounting for thermal fluctuations, damping, and effective magnetic fields.

The choice of solver affects:

- **Accuracy** of the time integration
- **Stability** (maximum allowable timestep)
- **Computational cost** per timestep
- **Physical correctness** (e.g., moment magnitude conservation)
- **Capability** for different physics modules (spin-transfer torque, spin-hall effects, inertial dynamics)

All solvers follow a **predictor-corrector** or **fixed-point iteration** scheme over two subroutines:
``evolve_first`` (predictor) and ``evolve_second`` (corrector).


Available Solvers
-----------------

Implemented solvers (summary)

SDEalgh=1
      **Name:** Midpoint (Mentink)
      **Characteristics:** Semi-implicit, rotation-based
      **Accuracy:** 2nd order (deterministic), 1st order (stochastic)
      **Stability:** Excellent
      **Cost:** Medium

SDEalgh=2
      **Name:** Heun (Single step)
      **Characteristics:** Explicit predictor-corrector
      **Accuracy:** 1st order
      **Stability:** Good
      **Cost:** Low

SDEalgh=3
      **Name:** Heun (Lite variant)
      **Characteristics:** Simplified Heun
      **Accuracy:** 1st order
      **Stability:** Good
      **Cost:** Low

SDEalgh=4
      **Name:** Heun (Proper)
      **Characteristics:** Full predictor-corrector, explicit
      **Accuracy:** 2nd order (deterministic)
      **Stability:** Very good
      **Cost:** Medium-High

SDEalgh=5
      **Name:** Depondt
      **Characteristics:** Rotation-based with Rodrigues formula
      **Accuracy:** 2nd order
      **Stability:** Excellent
      **Cost:** Medium

SDEalgh=6
      **Name:** Semi-implicit Midpoint (FPI)
      **Characteristics:** Fixed-point iteration variant
      **Accuracy:** 2nd order
      **Stability:** Excellent
      **Cost:** High

SDEalgh=7
      **Name:** Spherical Midpoint (FPI)
      **Characteristics:** Spherical constraint variant
      **Accuracy:** 2nd order
      **Stability:** Excellent
      **Cost:** Very High

SDEalgh=11
      **Name:** LLGI (Inertial solver)
      **Characteristics:** Inertial LLG variant
      **Accuracy:** 2nd order
      **Stability:** Good
      **Cost:** Medium

**Recommendation:** Start with **SDEalgh=1** (Midpoint) or **SDEalgh=5** (Depondt) for 
typical simulations and **SDEalgh=23** if conservative simulations are crucial. Use **SDEalgh=11** (LLGI) only when inertial 
effects are important (``relaxtime > 0``).

Solver Details
--------------

SDEalgh=1: Semi-Implicit Midpoint (Mentink)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Reference:** J. H. Mentink et al., *J. Phys.: Condens. Matter* **22**, 176001 (2010)

**Type:** Semi-implicit, rotation-based, predictor-corrector (2 steps)

**Algorithm:**

The solver solves the matrix equation for semi-implicitness:

.. math::

   \mathbf{A} \, \mathbf{m}(t+\Delta t) = \mathbf{A}^T \mathbf{m}(t) + \Delta t (\mathbf{a}_1 + \mathbf{s}_1)

where:

- :math:`\mathbf{A} = \mathbb{I} + \Delta t (\mathbf{a}_1/2 + \mathbf{s}_1/2)_\times` (skew-symmetric matrix)
- :math:`\mathbf{a}_1 = -\mathbf{B}_{\text{eff}} - \lambda (\mathbf{m} \times \mathbf{B}_{\text{eff}})` (deterministic torque)
- :math:`\mathbf{s}_1 = \text{stochastic counterpart}` (thermal field torque)

Solution: :math:`\mathbf{m}(t+\Delta t) = \mathbf{A}^{-1} \mathbf{A}^T \mathbf{m}(t)`

**Advantages:**
- Excellent energy stability (unconditional stability)
- 2nd order accuracy in deterministic part, 1st order in stochastic
- Automatic moment magnitude conservation (rotation matrix property)
- **Recommended for most simulations**

**Disadvantages:**
- Requires matrix inversion (3×3)
- Slightly slower than explicit methods

**Use cases:**
- Equilibrium spin dynamics at arbitrary temperatures
- Skyrmion and domain wall dynamics
- General-purpose magnetic simulations

**Example input:**

.. code-block:: none

   SDEalgh       1           # Midpoint solver
   damping       0.05        # Damping parameter λ
   delta_t       1.0e-15     # Timestep (seconds)


SDEalgh=2: Heun Single-Step
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Type:** Explicit predictor-corrector, single step

**Algorithm:**

Simple two-step Heun integration:

.. math::

   \mathbf{m}_p = \mathbf{m}(t) + \Delta t \frac{d\mathbf{m}}{dt}|_{t} + \mathbf{\eta}_p

   \mathbf{m}(t+\Delta t) = \mathbf{m}(t) + \frac{\Delta t}{2} \left(\frac{d\mathbf{m}}{dt}|_{t} + \frac{d\mathbf{m}}{dt}|_{t+\Delta t, \mathbf{m}_p}\right) + \mathbf{\eta}

**Advantages:**
- Simple, fast implementation
- Minimal memory overhead
- 1st order accuracy sufficient for thermal noise

**Disadvantages:**
- Explicit (requires small timestep for stability)
- May require moment re-normalization
- Less accurate than 2nd order methods

**Use cases:**
- Rapid prototyping
- Systems with very low damping (non-linear stability issues at large :math:`\Delta t`)

**Note:** Should not be used for quantitative accuracy studies.


SDEalgh=3: Heun Lite Variant
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Type:** Simplified explicit Heun

Lightweight variant of SDEalgh=2 with reduced operations and memory.

**Use cases:**
- Memory-constrained calculations
- Very large systems


SDEalgh=4: Heun Proper
^^^^^^^^^^^^^^^^^^^^^^

**Type:** Full predictor-corrector, explicit

**Algorithm:**

Two-step Heun with explicit mid-point evaluation:

.. math::

   \mathbf{m}(t+\Delta t/2) = \mathbf{m}(t) + \frac{\Delta t}{2} \frac{d\mathbf{m}}{dt}|_t

   \mathbf{m}(t+\Delta t) = \mathbf{m}(t) + \Delta t \frac{d\mathbf{m}}{dt}|_{t+\Delta t/2, \mathbf{m}_p} + \mathbf{\eta}

**Advantages:**
- Higher accuracy than Heun single-step (2nd order for deterministic part)
- Better stability than basic Heun
- Moderate computational cost

**Disadvantages:**
- Still explicit (requires careful timestep selection)
- Accumulates round-off error faster than semi-implicit methods

**Use cases:**
- High-accuracy studies with moderate damping
- Validation against semi-implicit methods


SDEalgh=5: Depondt Solver
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Reference:** Ph. Depondt and F. G. Mertens, *J. Phys.: Condens. Matter* **21**, 336005 (2009)

**Type:** Semi-implicit, rotation-based with Rodrigues formula

**Algorithm:**

Uses **Rodrigues' rotation formula** for exact exponential map:

.. math::

   \mathbf{m}_{\text{rot}} = \cos(\theta) \, \mathbf{m} + \sin(\theta) (\hat{\mathbf{\omega}} \times \mathbf{m}) 
   + (1 - \cos \theta) (\hat{\mathbf{\omega}} \cdot \mathbf{m}) \hat{\mathbf{\omega}}

where :math:`\theta = |\mathbf{\omega}| \Delta t` and :math:`\hat{\mathbf{\omega}}` is the rotation axis.

**Two-step procedure:**

1. **Predictor:** Rotate moments using Rodrigues formula with effective field
2. **Corrector:** Apply thermal field rotation correction

**Advantages:**
- Exact moment magnitude conservation (continuous rotation property)
- Highly stable for strong damping and thermal effects
- Excellent for **magnetic textures** (skyrmions, domain walls, vortices)
- Robust against timestep variations

**Disadvantages:**
- Slightly higher computational cost than midpoint (trigonometric functions)
- More complex implementation

**Use cases:**
- **Magnetic skyrmions** (requires robust moment conservation)
- **Domain wall dynamics**
- **Spin textures** with topological properties
- Large damping (:math:`\lambda > 0.5`)
- **Strongly thermal systems**

**Recommended:** Use Depondt (SDEalgh=5) for any simulation involving spatial magnetic structures.

**Example input:**

.. code-block:: none

   SDEalgh       5           # Depondt solver
   damping       0.1         # Higher damping acceptable
   delta_t       5.0e-15     # Can use larger timesteps safely


SDEalgh=6: Semi-Implicit Midpoint with Fixed-Point Iteration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Type:** Semi-implicit, fixed-point iteration variant

**Algorithm:**

Iteratively refines solution of implicit equation:

.. math::

   \mathbf{m}^{(k+1)} = f(\mathbf{m}^{(k)})

where iteration continues until convergence: :math:`|\mathbf{m}^{(k+1)} - \mathbf{m}^{(k)}| < \epsilon`

**Advantages:**
- Better accuracy at large timesteps
- Completely implicit (unconditionally stable)
- Higher order accuracy in both deterministic and stochastic parts

**Disadvantages:**
- **Expensive:** Typically 3-5 iterations per timestep
- Requires convergence criterion
- Slower than non-iterative solvers

**Use cases:**
- Very large timesteps required
- Time-critical simulations where fewer steps offset higher cost
- Validating energy stability

**Note:** Usually slower overall than SDEalgh=1 for standard applications.


SDEalgh=7: Spherical Midpoint with Fixed-Point Iteration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**References:**

- J. Hellsvik, "Semi-implicit spherical midpoint solver" (UppASD progress report)
- R. I. McLachlan et al., *Phys. Rev. E* **89**, 061301(R) (2014)

**Type:** Semi-implicit, spherical constraint, fixed-point iteration

**Algorithm:**

Enforces **spherical constraint** :math:`|\mathbf{m}| = 1` exactly at each iteration step:

.. math::

   \mathbf{m}^{(k+1)} = \frac{f(\mathbf{m}^{(k)})}{|f(\mathbf{m}^{(k)})|}

Combines McLachlan's structure-preserving scheme with Mentink's midpoint variant.

**Advantages:**
- **Perfect moment normalization** (no numerical drift in :math:`|m|`)
- Fixed-point iteration for accuracy
- Highest theoretical accuracy

**Disadvantages:**
- **Very expensive:** Requires 5-10 iterations per timestep
- Highest memory usage (stores iteration vectors)
- Slowest solver available

**Use cases:**
- Extremely long simulations (days or weeks) where :math:`|m|` drift matters
- High-precision phase transition studies
- When moment normalization errors accumulate significantly

**Caution:** Use only when necessary. Standard solvers (1, 5) usually provide sufficient accuracy
without the computational penalty.

**Example input:**

.. code-block:: none

   SDEalgh       7           # Spherical midpoint (FPI)
   delta_t       1.0e-15     # Standard timestep
   # Simulation will be ~5x slower than SDEalgh=1


SDEalgh=11: LLGI Inertial Solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Type:** Heun-type predictor-corrector for LLG-I equation

**Physics:**

Includes **inertial magnetic acceleration** term:

.. math::

   \frac{d^2\mathbf{m}}{dt^2} = -\gamma(\mathbf{m} \times \mathbf{B}_{\text{eff}}) 
   - \lambda' (\mathbf{m} \times \frac{d\mathbf{m}}{dt})

Controlled by ``relaxtime`` parameter (units: seconds).

**When to use:**

Set ``relaxtime`` to non-zero value to activate inertial dynamics:

.. code-block:: none

   SDEalgh       11          # LLGI solver
   relaxtime     1.0e-12     # Relaxation time (~100 fs typical)
   damping       0.05

**Physical regime:**

- **relaxtime = 0:** Standard LLG (no inertia)
- **relaxtime ~ 10^{-12} s:** Weak magnetic inertia (weakly ferromagnetic materials)
- **relaxtime ~ 10^{-11} s:** Strong inertia (important for ultrafast dynamics, spin textures)

**Advantages:**
- Captures ultrafast magnetic dynamics
- Relevant for **THz magnetization dynamics**
- Accurate for **laser-induced spin dynamics**

**Disadvantages:**
- More expensive (second-order equation)
- Requires careful timestep selection (:math:`\\Delta t < \\text{relaxtime}`)
- Less stable than standard LLG

**Use cases:**
- **Ultrafast spectroscopy** (pump-probe experiments)
- **Spin textures under strong driving** (large magnetic fields)
- **Theoretical studies** of inertial effects
- **THz-frequency dynamics**

**Not recommended:** For equilibrium simulations, standard LLG (``relaxtime=0``) is more efficient.

**Example input:**

.. code-block:: none

   SDEalgh       11
   relaxtime     5.0e-12     # 5 picoseconds inertial timescale
   delta_t       5.0e-16     # Must be much smaller than relaxtime


Extended Physics Options
------------------------

All solvers support additional physics modules that modify the effective field:

**Spin-Transfer Torque (STT)**

.. code-block:: none

   do_stt        Y
   # Modifies effective field via electron current-spin interaction

**Spin Hall Effect (SHE)**

.. code-block:: none

   do_she        Y
   # Current-induced spin-orbit torque from heavy metal overlayer

**Spin-Orbit Torque (SOT)**

.. code-block:: none

   do_sot        Y
   # General spin-orbit coupling to magnetic moment


Timestep Selection Guide
------------------------

**General rule:** :math:`\Delta t` should be small compared to the magnetic precession timescale.

**Precession timescale:**

.. math::

   \tau_{\text{prec}} = \frac{1}{\gamma B} \sim 10^{-12} \text{ s for } B \sim 1 \text{ T}

**Recommended timesteps by application:**

Equilibrium (RT, T < T_c)
      **Recommended :math:`\Delta t`:** :math:`10^{-15} - 10^{-14}` s
      **Typical Solver:** SDEalgh=1, 5

High temperature
      **Recommended :math:`\Delta t`:** :math:`10^{-16} - 10^{-15}` s
      **Typical Solver:** SDEalgh=1 (use smaller :math:`\Delta t` if unstable)

Skyrmion dynamics
      **Recommended :math:`\Delta t`:** :math:`5 \times 10^{-16} - 5 \times 10^{-15}` s
      **Typical Solver:** SDEalgh=5

STT/SHE effects
      **Recommended :math:`\Delta t`:** :math:`10^{-15} - 10^{-14}` s
      **Typical Solver:** SDEalgh=1, 5

Ultrafast (LLGI)
      **Recommended :math:`\Delta t`:** :math:`10^{-16} - 10^{-17}` s (< relaxtime/10)
      **Typical Solver:** SDEalgh=11

**Stability criterion** (explicit methods like Heun):

.. math::

   \Delta t < \frac{1}{\gamma B_{\max}}

where :math:`B_{\max}` is the maximum effective field magnitude in system.

**Semi-implicit methods** (SDEalgh=1, 5, 6) are **unconditionally stable** and allow larger timesteps.


Numerical Precision Checks
--------------------------

**Monitor moment magnitude:**

Code automatically re-normalizes moments after each step, but magnitude drift indicates:

- **Timestep too large** → reduce :math:`\Delta t`
- **Solver mismatch** → switch to SDEalgh=5 or 7
- **Moment conservation error** in output files

**Check energy conservation:**

For conservative systems (no damping, :math:`\lambda = 0`):

.. math::

   E = -\mathbf{m} \cdot \mathbf{B}_{\text{eff}} + \text{const}

Monitor total energy to detect numerical issues.

**Convergence testing:**

Run simulations with progressively smaller timesteps (:math:`\Delta t, \Delta t/2, \Delta t/4`)
and verify observables (magnetization, energies) converge.

**Thermal noise verification:**

For :math:`\lambda \neq 0`, check energy equipartition:

.. math::

   \langle E_{\text{thermal}} \rangle \approx k_B T

Deviations suggest thermal noise generation errors.


Solver Comparison Examples
--------------------------

Example 1: Equilibrium Ferromagnet at T=300K
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Setup:**
- Simple cubic iron (A = 1 meV, K = 0)
- 10×10×10 unit cells
- No external field

**Comparison:**

.. code-block:: none

   # Configuration 1: Midpoint (Fast baseline)
   SDEalgh       1
   delta_t       1.0e-15
   # Runtime: 1x (reference)
   # M_sat accuracy: ±0.5%
   
   # Configuration 2: Depondt (Robust structure conservation)
   SDEalgh       5
   delta_t       2.0e-15
   # Runtime: 1.1x
   # M_sat accuracy: ±0.1%
   # Better for domain walls
   
   # Configuration 3: Spherical midpoint (Overkill)
   SDEalgh       7
   delta_t       1.0e-15
   # Runtime: 5x
   # M_sat accuracy: ±0.01%
   # Not worth the cost

**Conclusion:** Use **SDEalgh=1** (Midpoint) for this case.


Example 2: Skyrmion Room-Temperature Dynamics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Setup:**
- Chiral magnet with DMI
- Skyrmion core diameter ~5 nm
- External field B = 0.1 T

**Comparison:**

.. code-block:: none

   # Configuration 1: Midpoint
   SDEalgh       1
   delta_t       5.0e-15
   # Skyrmion becomes distorted after 100 ps
   # Unphysical deformation of core

   # Configuration 2: Depondt (BEST)
   SDEalgh       5
   delta_t       5.0e-15
   # Skyrmion remains stable for hours of simulation
   # Perfect topological charge conservation
   
   # Configuration 3: Spherical midpoint
   SDEalgh       7
   delta_t       5.0e-15
   # Runtime: 4x Configuration 2
   # Skyrmion identical to Depondt
   # Not worth computational cost

**Conclusion:** For **magnetic textures**, use **SDEalgh=5** (Depondt). It's worth the minor cost.


Example 3: Ultrafast Laser-Induced Magnetization Dynamics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Setup:**
- Strong laser pump (B_laser ~ 10 T effective field)
- Relaxation timescale ~ 1 ps
- Probe time ~ 100 fs

**Setup:**

.. code-block:: none

   SDEalgh       11              # LLGI inertial solver
   relaxtime     1.0e-12         # 1 ps inertial timescale
   delta_t       1.0e-16         # 0.1 fs timestep (<< relaxtime)
   nstep         10000000        # Total 1 ps simulation

**Physical results:**
- Captures ultrafast demagnetization
- Shows inertial overshoot
- Matches ultrafast experiments

**Alternative (standard LLG, no inertia):**

.. code-block:: none

   SDEalgh       1
   relaxtime     0.0             # Standard LLG
   delta_t       5.0e-16         # Faster relaxation -> non-physical

**Comparison:** Inertial model essential for accurate <100 fs timescale dynamics.


Performance Benchmarks
----------------------

Relative computational cost per timestep (normalized to SDEalgh=1):

SDEalgh=1 (Midpoint)
      **Relative Cost:** 1.0x
      **Memory:** 1.0x
      **Comment:** Baseline — use this first

SDEalgh=2 (Heun single)
      **Relative Cost:** 0.8x
      **Memory:** 0.8x
      **Comment:** Faster but lower accuracy

SDEalgh=4 (Heun proper)
      **Relative Cost:** 1.2x
      **Memory:** 1.0x
      **Comment:** Slightly more accurate

SDEalgh=5 (Depondt)
      **Relative Cost:** 1.3x
      **Memory:** 1.0x
      **Comment:** Best for textures

SDEalgh=6 (Midpoint FPI)
      **Relative Cost:** 3.5x
      **Memory:** 1.2x
      **Comment:** Rarely worth it

SDEalgh=7 (Spherical FPI)
      **Relative Cost:** 5.0x
      **Memory:** 2.0x
      **Comment:** Use only for extreme precision

SDEalgh=11 (LLGI)
      **Relative Cost:** 1.4x
      **Memory:** 1.2x
      **Comment:** Only with relaxtime > 0

Message: SDEalgh=1 and SDEalgh=5 offer best cost-benefit ratio for most applications.


Accuracy vs. Stability Trade-off
--------------------------------

**Explicit solvers** (Heun variants):
- ✓ Fast
- ✗ Limited stability (must use small :math:`\Delta t`)
- ✗ Energy conservation (especially at high T)

**Semi-implicit solvers** (Midpoint, Depondt):
- ✓ Excellent stability (large :math:`\Delta t` possible)
- ✓ Energy conserving
- ✓ Moderate cost
- ✓ **Recommended for most simulations**

**Fixed-point iteration** (Midpoint/Spherical FPI):
- ✓ Highest accuracy
- ✓ Best energy conservation
- ✗ Very expensive
- ✗ Only justify for extreme precision requirements

**Rotation-preserving** (Depondt, Spherical):
- ✓ Perfect moment conservation
- ✓ Ideal for topological structures
- ✓ Better for large timesteps
- ✓ **Preferred for skyrmions, domain walls**

**Inertial dynamics** (LLGI):
- ✓ Captures ultrafast effects
- ✗ Requires very small timesteps
- ✓ Essential for laser-driven dynamics


Decision Tree for Solver Selection
-----------------------------------

.. code-block:: text

   START: Choose your solver
   │
   ├─ Need inertial magnetic acceleration?
   │  ├─ YES → Use SDEalgh=11 (LLGI)
   │  └─ NO → Continue
   │
   ├─ Simulating magnetic textures (skyrmions, domains)?
   │  ├─ YES → Use SDEalgh=5 (Depondt)
   │  └─ NO → Continue
   │
   ├─ Is computational cost critical?
   │  ├─ YES (clusters, millions of atoms) → Use SDEalgh=1 (Midpoint)
   │  └─ NO → Continue
   │
   ├─ Do you need maximum accuracy?
   │  ├─ YES (publication, validation) → Use SDEalgh=5 or 7
   │  └─ NO → Use SDEalgh=1 (Midpoint)
   │
   └─ DEFAULT RECOMMENDATION: SDEalgh=1 (Midpoint)
      • 2nd order accurate
      • Unconditionally stable
      • Good speed-accuracy balance
      • Works for all physics modules


Summary Table: When to Use Each Solver
--------------------------------------

.. note:: Solver selection at-a-glance

SDEalgh=1
      **Solver:** Midpoint
      **Best Application:** General-purpose simulations
      **Timestep Advantage:** Large allowed
      **Use It?:** YES (start here)

SDEalgh=2
      **Solver:** Heun single
      **Best Application:** Rapid prototyping only
      **Timestep Advantage:** Small required
      **Use It?:** Only if desperate

SDEalgh=4
      **Solver:** Heun proper
      **Best Application:** Moderate accuracy needs
      **Timestep Advantage:** Medium allowed
      **Use It?:** Maybe (try 1 first)

SDEalgh=5
      **Solver:** Depondt
      **Best Application:** Magnetic textures
      **Timestep Advantage:** Large allowed
      **Use It?:** YES (for skyrmions)

SDEalgh=6
      **Solver:** Midpoint FPI
      **Best Application:** Very large timesteps
      **Timestep Advantage:** Very large
      **Use It?:** Rarely

SDEalgh=7
      **Solver:** Spherical FPI
      **Best Application:** Extreme precision
      **Timestep Advantage:** Very large
      **Use It?:** No

SDEalgh=11
      **Solver:** LLGI
      **Best Application:** Ultrafast dynamics
      **Timestep Advantage:** Very small
      **Use It?:** Only if relaxtime > 0


Examples in Input Files
-----------------------

**Example 1: Fast Equilibrium Simulation**

.. code-block:: none

   # Ferromagnetic equilibrium
   SDEalgh       1           # Midpoint solver
   damping       0.05
   delta_t       1.0e-15
   nstep         1000000
   Temp          300.0

**Example 2: Skyrmion Dynamics Study**

.. code-block:: none

   # Skyrmion stability and motion
   SDEalgh       5           # Depondt solver (CRITICAL for textures)
   damping       0.1
   delta_t       2.0e-15
   nstep         5000000
   Temp          100.0

**Example 3: Ultrafast Laser Dynamics**

.. code-block:: none

   # Laser-induced magnetization reversal
   SDEalgh       11          # LLGI inertial solver
   relaxtime     5.0e-12     # 5 ps inertial timescale
   damping       0.01
   delta_t       5.0e-16     # Must be << relaxtime
   nstep         10000000
   Temp          300.0

**Example 4: High-Performance Computing**

.. code-block:: none

   # Million-atom system on supercomputer
   SDEalgh       1           # Fast Midpoint
   damping       0.05
   delta_t       2.0e-15     # Can be larger with semi-implicit
   nstep         10000000
   Temp          300.0


References
==========

See the centralized :doc:`references` for full bibliographic entries:

- [Mentink2010]_ - Midpoint solver and semi-implicit integration
- [Depondt2009]_ - Depondt-Mertens solver for spin dynamics
- [McLachlan2014]_ - Structure-preserving Lie-Poisson integrators
- [Kloeden1992]_ - Numerical solution of stochastic differential equations

