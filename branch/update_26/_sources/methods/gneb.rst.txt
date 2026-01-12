Geodesic nudged elastic band (GNEB)
====================================

Parameters for minimum energy path (MEP) calculations using GNEB
------------------------------------------------------------------

The UppASD code implements the **Geodesic Nudged Elastic Band (GNEB)** method
[Bessarab2012]_, [Bessarab2015]_ for calculating the **minimum energy path (MEP)**
between two magnetic configurations. This enables computation of:

- **Energy barriers** for magnetic transitions
- **Transition rates** via Harmonic Transition State Theory (HTST)
- **Saddle point (transition state) configurations**
- **Prefactors** for transition rates and lifetimes

The GNEB method is particularly useful for:

- Finding energy barriers in systems with competing interactions (frustrated magnets)
- Studying skyrmion creation and annihilation pathways
- Analyzing magnetic domain wall motion
- Computing transition rates in magnetic switching processes
- Identifying transition states for rare events in spin systems

GNEB combines the **nudged elastic band (NEB)** method, which connects images
along a path between two endpoints, with the **geodesic constraint** to ensure
calculations remain on the energy landscape surface. This is complemented by
For the original NEB method, see [Henkelman2000]_. This is complemented by
the **Velocity Projection Optimization (VPO)** algorithm used in the initial
phase for relaxing endpoint configurations.

---------------------------------------------------
Overview of the GNEB method
---------------------------------------------------

The GNEB method constructs a discrete path connecting two magnetic configurations
(initial and final states) using a series of :math:`N_{\text{images}}` intermediate
configurations called **images**. The path is evolved to the **minimum energy path (MEP)**
by simultaneously:

1. **Relaxing energy minima**: Images at local minima relax to nearby saddle points
2. **Constraining the path**: Spring forces keep images connected along the path
3. **Geodesic projection**: The constraint is applied only tangent to the energy landscape

The MEP connects the initial state (image 1) to the final state (image :math:`N_{\text{images}}`).
Intermediate images (2 to :math:`N_{\text{images}}-1`) are evolved toward the MEP.

---------------------------------------------------
The VPO minimizer
---------------------------------------------------

The **Velocity Projection Optimization (VPO)** algorithm [Bessarab2012]_ is used in
the initial phase (``ip_mode G`` or ``ip_mode E``) to relax individual magnetic configurations
to nearby energy minima. VPO is based on a velocity-Verlet integration scheme applied
to the equation of motion on the energy landscape.

The algorithm treats the magnetic moments as point masses on the energy surface and
integrates a damped equation of motion:

.. math::

   m \frac{d\mathbf{v}}{dt} = \mathbf{F} - \zeta \mathbf{v}

where :math:`\mathbf{v}` is the velocity, :math:`\mathbf{F}` is the effective force
(derived from the effective magnetic field), :math:`m` is a fictitious mass, and
:math:`\zeta` is a damping coefficient.

The effective force on each magnetic moment :math:`\mathbf{m}_i` is computed as the
component of the effective field perpendicular to the moment:

.. math::

   \mathbf{F}_i = \mathbf{m}_i \times (\mathbf{m}_i \times \mathbf{B}_{\text{eff},i})

This ensures that moments remain on the unit sphere during optimization.

The VPO algorithm includes:

1. **Velocity updates** based on the effective force
2. **Moment updates** via rotation around a computed axis
3. **Velocity projection** to maintain the geodesic constraint
4. **Convergence criterion**: Maximum effective field force falls below ``min_ftol``

VPO is controlled by two key parameters:

- ``vpo_mass``: Fictitious mass :math:`m` for the equation of motion. Larger values lead to slower updates but more stable convergence.
- ``vpo_dt``: Integration timestep :math:`\Delta t`. Should be small enough for stability (typically 0.001–0.01).

---------------------------------------------------
Initial phase: VPO relaxation of endpoints
---------------------------------------------------

The initial phase (``ip_mode G`` or ``ip_mode E``) uses VPO to relax the initial
and final configurations to nearby **energy minima**. This is important because:

- GNEB is most efficient when endpoints are at local minima
- It ensures the MEP connects true local minima
- It reduces the number of GNEB iterations needed

The VPO initial phase:

1. Reads initial and final configurations from the configuration files
2. Performs independent VPO minimization on each endpoint
3. Saves the minimized configurations for use in the GNEB measurement phase

If ``ip_mode`` is set to a Monte Carlo mode (M, H) instead, GNEB will still work
but the endpoints are not optimized and may not be at true local minima.

---------------------------------------------------
GNEB measurement phase: Finding the MEP
---------------------------------------------------

In the measurement phase (``mode G`` for standard GNEB, or ``mode E`` for climbing
image GNEB), the code:

1. **Generates an initial path** connecting the two endpoint configurations (interpolation)
2. **Evolves the path** toward the MEP using VPO optimization on all images
3. **Applies spring forces** to keep images evenly distributed
4. **Computes the saddle point** (highest energy along the MEP)
5. **Calculates the energy barrier** as the difference between saddle point and endpoints
6. **Optionally computes the Hessian matrix** at key points for HTST calculations

Spring force and elastic band
--------------------------------------------------

The initial path connecting the two endpoints can be generated in two ways:

1. **Geodesic interpolation** (``initpath = 1``): A geodesic curve on the energy landscape connecting the endpoints
2. **External path file** (``initpath = 2``): Read a path from the restart file specified in ``restartfile``

The **geodesic path** is computed as the shortest path on the unit sphere connecting
the two configurations. For each atom, the interpolation is:

.. math::

   \mathbf{m}_i(\lambda) = \frac{\sin((1-\lambda)\theta_i)}{\sin(\theta_i)} \mathbf{m}_i^{\text{init}} + \frac{\sin(\lambda \theta_i)}{\sin(\theta_i)} \mathbf{m}_i^{\text{final}}

where :math:`\lambda \in [0, 1]` is the path parameter and :math:`\theta_i` is the
angle between the initial and final moment directions for atom :math:`i`.

This ensures that the initial path is a **geodesic** connecting the endpoints on the
configuration space.

Spring force and elastic band
--------------------------------------------------

To keep images distributed along the path, **spring forces** are applied between
neighboring images. These springs have force constant ``spring`` (or ``kappa`` in code).
The spring force between images :math:`i` and :math:`i+1` is:

.. math::

   \mathbf{F}_{\text{spring},i} = \text{spring} \times (r_{i+1} - r_i) \times \frac{\mathbf{r}_{i+1} - \mathbf{r}_i}{|\mathbf{r}_{i+1} - \mathbf{r}_i|}

where :math:`r_i` is the path position (arc length) of image :math:`i`. Larger spring
constants prevent images from clustering.

Climbing Image GNEB (CI-GNEB)
--------------------------------------------------

the image with the highest energy into a **climbing image**. This image:

1. Is pushed uphill along the band direction
2. Converges directly to the **saddle point** (transition state)
3. Provides exact location and energy of the saddle point

The acceptance criterion for the climbing image is inverted: instead of being pushed
down, it is pushed up along the MEP direction, climbing toward the saddle point energy.

CI-GNEB is activated by setting ``do_gneb_ci = Y`` and ``do_gneb = N``.

Tangent vector and path constraints
------------------------------------

GNEB and MEP keywords
---------------------

ip_mode
   Initial phase mode. Set to **G** for VPO-based energy minimization of endpoint configurations. Default: **S** (spin dynamics). See also: M=Metropolis MC, H=Heat bath MC, W=Wang-Landau, P=Parallel Tempering, E=Endpoints-only minimization.

mode
   Measurement phase mode. Set to **G** for GNEB MEP calculation (standard GNEB). For climbing image variant, use ``do_gneb_ci = Y`` and ``do_gneb = N``. Default: **S** (spin dynamics).

minalgo
   Minimization algorithm. Currently, only **1** (VPO algorithm) is supported. Default: **1**.

min_ftol
   Convergence tolerance for the force (effective field magnitude in mRy/μ_B). The minimizer stops when the maximum force drops below this value. Smaller values ensure more accurate relaxation but require more iterations. Default: **1e-5**.

vpo_mass
   Fictitious mass for VPO equation of motion. Controls the "inertia" of the minimization. Larger values smooth the path but slow convergence. Typical range: 0.1–1.0. Default: **1.0**.

vpo_dt
   Integration timestep for VPO. Must be small enough for stability. Typical range: 0.001–0.05. Default: **0.01**.

do_gneb
   Enable standard GNEB MEP calculation (Y=yes, **N=no**). If **Y**, the code finds the MEP connecting the two endpoints. Default: **Y** for ``mode G``.

do_gneb_ci
   Enable Climbing Image GNEB variant (Y=yes, **N=no**). If **Y**, one image (the highest-energy image) climbs to the saddle point. Default: **N**.

mepitrmax
   Maximum number of iterations for GNEB MEP relaxation. Default: **1000**.

mepftol
   Convergence tolerance for GNEB (force criterion). GNEB stops when maximum force < this value. Default: **1e-4**.

spring
   Spring constant (kappa) for elastic band. Controls the strength of spring forces keeping images connected. Typical range: 0.1–1.0. Default: **0.5**.

initpath
   Initial path generation method: **1**: geodesic interpolation (recommended), **2**: read from external file. Default: **1**.

fixed_if
   Fix the endpoint configurations (Y=yes, **N=no**). If **Y**, endpoints remain at their initial positions during GNEB. If **N**, endpoints can move along energy isocontours. Default: **Y** (endpoints fixed).

Mensemble
   Number of images in the MEP path. This is equivalent to ``Num_images`` in some references. Must be > 2. Typical values: 10–50 depending on path complexity. Default: **11**.

Additional GNEB-specific keywords:

- ``en_zero``: Energy reference level (default: **E**, initial state)
- ``meptraj_step``: Save MEP trajectory every ``meptraj_step`` iterations
- ``sample_num``: Number of points for energy interpolation along the MEP
- ``prn_gneb_fields``: Print effective fields at each site along the path (Y/N)
- ``do_hess_ini``, ``do_hess_fin``, ``do_hess_sp``: Compute Hessian matrix at initial, final, and saddle point configurations (Y/N)

which defines the path direction. The tangent is computed as the normalized direction
connecting neighboring images:

.. math::

   \hat{\tau}_i = \frac{\mathbf{r}_{i+1} - \mathbf{r}_{i-1}}{|\mathbf{r}_{i+1} - \mathbf{r}_{i-1}|}

for interior images, with special handling at the endpoints.

All forces (effective field forces and spring forces) are projected perpendicular
to the tangent vector:

.. math::

   \mathbf{F}_{\perp,i} = \mathbf{F}_i - (\mathbf{F}_i \cdot \hat{\tau}_i) \hat{\tau}_i

This ensures that the path evolves along the landscape without being pushed along
the tangent direction, which would just slide images without finding the MEP.

--------------------------------------------------
Output files
--------------------------------------------------

The GNEB simulation produces several output files:

- ``mep_<simid>.out``: Final path configuration (all images in sequence)
- ``en_path.<simid>.out``: Energy at each image along the MEP
- ``enfit_path.<simid>.out``: Interpolated energy along the path using Hermite fitting
- ``en_sp.<simid>.out``: Energy of the identified saddle point
- ``sp_<simid>.out``: Magnetic configuration of the saddle point
- ``force_min.<simid>.out``: Force history during minimization (initial phase)
- ``beff_<image>.<simid>.out``: Effective field components at each site (optional, if ``prn_gneb_fields = Y``)

---------------------------------------------------
Harmonic Transition State Theory (HTST)
---------------------------------------------------

Once the saddle point and MEP are known, the **Harmonic Transition State Theory (HTST)**
can be used to compute transition rates and lifetimes (see [Eyring1935]_). The transition rate is:

.. math::

   k = \nu_0 \frac{\det H_{\text{init}}}{\det H_{\text{sp}}} \exp\left(-\frac{E_a}{k_B T}\right)

where:

- :math:`E_a` is the activation energy (barrier height)
- :math:`\nu_0 \approx k_B T / h` is the attempt frequency
- :math:`H_{\text{init}}` is the Hessian matrix at the initial state
- :math:`H_{\text{sp}}` is the Hessian matrix at the saddle point
- The ratio of determinants is the **prefactor**

To compute the Hessian, set:

.. code-block:: text

   do_hess_ini  Y
   do_hess_sp   Y

The code will automatically compute the Hessian matrices and save eigenvalues/eigenvectors
for subsequent HTST analysis.

---------------------------------------------------
Example: MEP for skyrmion nucleation
---------------------------------------------------

Example ``inpsd.dat`` for finding the MEP for skyrmion creation in a magnetic film:

.. code-block:: text

   simid    Skyrmion_MEP

   !! Initial phase: VPO relaxation of endpoints
   ip_mode    G
   ip_nstep   50000

   !! GNEB measurement phase
   mode       G
   mcnstep    500000

   !! GNEB parameters
   do_gneb         Y
   do_gneb_ci      N
   mepitrmax       500
   mepftol         1e-4
   spring          0.5

   !! Path setup
   initpath        1
   fixed_if        Y
   Mensemble       21

   !! Minimization parameters
   minalgo         1
   min_ftol        1e-5
   vpo_mass        1.0
   vpo_dt          0.01

   !! Output options
   meptraj_step    50
   sample_num      200
   do_hess_sp      Y

In this example:

- 21 images are used to represent the path
- The VPO minimizer relaxes initial and final configurations
- GNEB evolves the path with spring constant 0.5
- The saddle point is identified and the Hessian is computed for HTST

The energy barriers, MEP, and Hessian eigenvalues are output to files for analysis.

---------------------------------------------------
Example: CI-GNEB for transition state refinement
---------------------------------------------------

Climbing Image GNEB directly converges to the saddle point:

.. code-block:: text

   simid    Skyrmion_CIGNEB

   !! Initial phase: VPO relaxation
   ip_mode    G
   ip_nstep   50000

   !! CI-GNEB measurement phase
   mode       G
   mcnstep    300000

   !! CI-GNEB parameters
   do_gneb         N
   do_gneb_ci      Y
   mepitrmax       300
   mepftol_ci      1e-5
   spring          0.7

   !! Path setup
   initpath        1
   fixed_if        Y
   Mensemble       11

   !! Minimization
   minalgo         1
   min_ftol        1e-5
   vpo_mass        1.0
   vpo_dt          0.01

Here, ``mepftol_ci`` is tighter than the standard GNEB tolerance, ensuring the
climbing image converges precisely to the saddle point.

---------------------------------------------------
Performance and convergence
---------------------------------------------------

GNEB convergence can be improved by:

1. **Starting from relaxed endpoints**: Use VPO initial phase (``ip_mode G``) to ensure endpoints are at local minima
2. **Geodesic path initialization**: Start with geodesic interpolation (``initpath = 1``)
3. **Appropriate spring constant**: Too large → image bunching; too small → divergence
4. **Sufficient images**: More images (larger ``Mensemble``) resolve sharp barriers
5. **Climbing Image variant**: Use ``do_gneb_ci = Y`` if the saddle point is difficult to locate

For very steep or rough energy landscapes, adaptive spring constants or temperature-based sampling may be needed.

See :doc:`../references` for bibliography entries used in this page.

See also
--------

- :doc:`monte-carlo` (Monte Carlo sampling)
- :doc:`../input/system` (system setup)
- :doc:`../input/hamiltonian` (magnetic interactions)


References
----------

See the centralized :doc:`../references` for full bibliographic entries:

- [Bessarab2012]_ - Method for finding magnetic transition mechanisms and activation energies
- [Bessarab2015]_ - Harmonic transition-state theory for thermal spin transitions
- [Henkelman2000]_ - Climbing image nudged elastic band method
- [Eyring1935]_ - Activated complex theory for chemical reactions


