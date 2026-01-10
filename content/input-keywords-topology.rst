Topological measurements and observables
=========================================

.. _input-topology:

Overview
--------

UppASD provides comprehensive tools for characterizing topological magnetic
textures, primarily through the calculation of the skyrmion number (Pontryagin
index) and related local density measures. Two complementary methods are
implemented:

1. **Gradient-based (finite-difference)** approach using derivatives of the
   magnetization field
2. **Triangulation-based (volume integral)** approach using Delaunay simplices

Both methods compute the topological charge density integrated over the system,
but differ in numerical accuracy, computational cost, and sensitivity to lattice
discretization.

Skyrmion number: definition and physical meaning
-------------------------------------------------

The skyrmion number :math:`Q` (or Pontryagin index) quantifies the topological
winding of a magnetization texture in two-dimensional (or quasi-2D) systems. It
counts how many times the magnetization vector wraps around the unit sphere when
mapped over real space:

.. math::

   Q = \frac{1}{4\pi} \int\!\int \mathbf{m}\cdot\Bigl(\frac{\partial\mathbf{m}}{\partial x}
   \times \frac{\partial\mathbf{m}}{\partial y}\Bigr)\,dx\,dy,

where :math:`\mathbf{m}(\mathbf{r})` is the unit magnetization field. For a
discrete lattice, this becomes a sum over sites or simplices. A **skyrmion** has
:math:`Q=-1` (or :math:`+1` depending on convention); an **antiskyrmion** has
opposite sign; ferromagnetic domains have :math:`Q=0`.

Gradient-based method (finite-difference stencils)
--------------------------------------------------

The gradient-based approach evaluates the topological charge density at each
lattice site using finite-difference approximations of :math:`\nabla\mathbf{m}`:

.. math::

   Q = \frac{1}{4\pi} \sum_i \mathbf{m}_i \cdot \Bigl(\frac{\partial\mathbf{m}}{\partial x}\Big|_i
   \times \frac{\partial\mathbf{m}}{\partial y}\Big|_i\Bigr) \Delta A,

where :math:`\Delta A` is the effective area per site. The code computes
:math:`\nabla\mathbf{m}` via central differences over the neighbor list.

Implemented in ``pontryagin_no(Natom, Mensemble, emomM, grad_mom)`` and
``pontryagin_no_density(iatom, Natom, Mensemble, emomM, grad_mom)``.

**Advantages:**

- Straightforward implementation using existing neighbor lists and gradient routines
- Naturally provides site-resolved topological density :math:`q_i` for visualization
- Works for arbitrary lattices (not restricted to triangular grids)
- Accurate when the magnetization field is smooth on the lattice scale

**Disadvantages:**

- Requires pre-computed gradients (``grad_mom`` array), adding memory overhead
- Accuracy depends on stencil order (typically 2nd-order central differences)
- Sensitive to lattice discretization; coarse grids introduce numerical errors
- May be less accurate near domain walls or defects where :math:`\mathbf{m}` varies rapidly

Triangulation-based method (volume integrals)
---------------------------------------------

The triangulation method divides the 2D lattice into Delaunay triangles (simplices)
and computes the solid angle subtended by the magnetization triple
:math:`(\mathbf{m}_1,\mathbf{m}_2,\mathbf{m}_3)` at each triangle vertex on the
unit sphere:

.. math::

   Q = \frac{1}{4\pi}\sum_{\text{triangles}} \Omega(\mathbf{m}_1,\mathbf{m}_2,\mathbf{m}_3),

where the solid angle is

.. math::

   \Omega = 2\arctan\Bigl(\frac{\mathbf{m}_1\cdot(\mathbf{m}_2\times\mathbf{m}_3)}
   {1 + \mathbf{m}_1\cdot\mathbf{m}_2 + \mathbf{m}_1\cdot\mathbf{m}_3
   + \mathbf{m}_2\cdot\mathbf{m}_3}\Bigr).

This formula (the Berg-Lüscher winding number) is exact for piecewise-linear
interpolation of :math:`\mathbf{m}` over triangles.

Implemented in ``pontryagin_tri(Natom, Mensemble, emom)`` and related functions.
The Delaunay triangulation is precomputed via ``delaunay_tri_tri(nx, ny, nz, NT)``
for triangular lattices with periodic boundary conditions.

**Advantages:**

- Numerically robust; no explicit gradients needed
- Exact for smooth textures when the mesh resolves the structure
- Less sensitive to local noise in :math:`\mathbf{m}` compared to finite differences
- Natural for triangular or hexagonal lattices

**Disadvantages:**

- Requires a triangulation (currently hard-coded for triangular lattices)
- More complex bookkeeping (simplex indices, neighbor connectivity)
- Not easily extended to non-regular or 3D grids without external triangulation libraries
- Slight computational overhead from trigonometric functions (``atan``)

Comparison and recommendations
------------------------------

.. list-table:: Gradient vs Triangulation methods
   :widths: 25 35 40
   :header-rows: 1

   * - Aspect
     - Gradient (finite-difference)
     - Triangulation (volume integral)
   * - **Accuracy**
     - 2nd-order (central diff); sensitive to grid spacing
     - Exact for piecewise-linear :math:`\mathbf{m}`
   * - **Lattice restriction**
     - Any lattice with neighbors
     - Currently triangular only
   * - **Memory**
     - Requires ``grad_mom(3,3,Natom,Mensemble)``
     - Requires ``simp(3,nsimp)`` index array
   * - **Local density**
     - Direct per-site :math:`q_i`
     - Indirect (sum over simplices sharing site :math:`i`)
   * - **Speed**
     - Fast (simple dot/cross products)
     - Slightly slower (``atan`` calls)
   * - **Robustness**
     - Sensitive to noise and sharp features
     - More stable for smooth textures

**Guidelines:**

- Use **gradient method** for quick screening, arbitrary lattices, or when local
  density maps are the primary goal.
- Use **triangulation method** for high-accuracy skyrmion counting on triangular/hexagonal
  lattices, especially for publication-quality topological charge.
- For spin dynamics (SD), gradient-based is more common; for Monte Carlo (MC),
  both are viable (triangulation may be preferred for better convergence).

Keyword reference
-----------------

skyno
   Enable skyrmion number measurement (Y/*N*). Set to ``Y`` to activate
   triangulation-based calculation via ``pontryagin_tri`` or gradient-based via
   ``pontryagin_no``, depending on internal flags.

skyno_step
   Number of time steps between skyrmion number samples (integer, default 100).
   Controls sampling frequency for time-resolved topological dynamics.

skyno_buff
   Number of samples to buffer before writing to file (integer, default 10).
   Reduces I/O overhead for large simulations.

do_proj_skyno
   Type-projected (sublattice) skyrmion number (Y/*N*). Computes
   :math:`Q_{\alpha}` for each atom type :math:`\alpha` using
   ``proj_pontryagin_no``. Useful for multi-component systems (e.g., Néel vs
   Bloch skyrmions on different sublattices).

do_skyno_den
   Site-resolved skyrmion number density (Y/*N*). Outputs local topological
   charge :math:`q_i` at each site via ``pontryagin_tri_dens`` or
   ``pontryagin_no_density``. Essential for visualizing skyrmion cores, domain
   walls, and merons.

do_skyno_cmass
   Center-of-mass skyrmion tracking (Y/*N*). Experimental feature for tracking
   skyrmion positions and velocities. Requires spatial integration of
   :math:`q_i`.

Gradient calculation and stencils
---------------------------------

When using gradient-based methods, UppASD computes :math:`\nabla\mathbf{m}` via
central finite differences:

.. math::

   \frac{\partial\mathbf{m}}{\partial x}\Big|_i \approx
   \frac{\mathbf{m}_{i+\hat{x}} - \mathbf{m}_{i-\hat{x}}}{2\Delta x},

where neighbors are identified from the exchange neighbor list. The accuracy is
:math:`O(\Delta x^2)`. For anisotropic lattices, coordinate transformations
ensure consistent units.

**Improving gradient accuracy:**

- Use denser lattices (smaller :math:`\Delta x`) for sharper textures
- Ensure exchange cutoff ``rcutoff`` includes sufficient neighbors for stencil
- For higher-order schemes, custom gradient routines can replace the default
  2-point stencil (not currently automated in UppASD)

Delaunay triangulation details
-------------------------------

The ``delaunay_tri_tri`` routine constructs a triangulation for a regular
:math:`N_x\times N_y\times N_z` lattice (with :math:`N_{\text{type}}` atoms per
cell) assuming a triangular base:

- Each unit cell contributes two triangles
- Periodic boundary conditions wrap indices via ``wrap_idx``
- The simplex array ``simp(3,nsimp)`` stores vertex indices

For non-triangular lattices, external triangulation (e.g., Qhull, Triangle) is
required and must be interfaced manually.

Output files and analysis
-------------------------

skyrmion.<simid>.out
   Time-series of total skyrmion number :math:`Q(t)` and ensemble statistics.
   Format: ``step``, :math:`Q_{\text{mean}}`, :math:`\sigma_Q`.

skyrmion_proj.<simid>.out
   Type-projected skyrmion numbers :math:`Q_{\alpha}(t)` for each sublattice.

skyrmion_dens.<simid>.out
   Site-resolved topological density :math:`q_i` at each atom. Used for
   real-space visualization of skyrmion cores.

Example setup
-------------

**Gradient-based skyrmion counting:**

.. code-block:: text

   ! Enable gradient calculation in Hamiltonian setup
   do_jtensor   1          ! Ensure gradients are computed

   ! Skyrmion measurements
   skyno        Y          ! Enable skyrmion number
   skyno_step   50         ! Sample every 50 steps
   skyno_buff   20         ! Buffer 20 samples
   do_skyno_den Y          ! Output local density

**Triangulation-based skyrmion counting (triangular lattice):**

.. code-block:: text

   skyno        Y
   skyno_step   100
   skyno_buff   10
   do_proj_skyno Y         ! Type-projected Q

   ! Ensure lattice is triangular/hexagonal for automatic triangulation

Related keywords
----------------

- ``do_jtensor``: Enable gradient calculations (see :doc:`input-keywords-hamiltonian`)
- ``max_no_neigh``: Sufficient neighbors for gradient stencils
- Visualization tools: use ``skyrmion_dens.<simid>.out`` with ParaView or
  matplotlib for real-space topological maps

References
----------

Key references for topological charge calculations:

- Berg, B. and Lüscher, M., *Definition and statistical distributions of a
  topological number in the lattice O(3) sigma-model*, Nucl. Phys. B **190**, 412 (1981).
- Nagaosa, N. and Tokura, Y., *Topological properties and dynamics of magnetic
  skyrmions*, Nat. Nanotech. **8**, 899 (2013).
- Rohart, S. and Thiaville, A., *Skyrmion confinement in ultrathin film
  nanostructures in the presence of Dzyaloshinskii-Moriya interaction*, Phys.
  Rev. B **88**, 184422 (2013).
