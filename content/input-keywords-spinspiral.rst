Spin-spiral minimization (qminimizer / sx_driver)
=================================================

.. _input-spinspiral:

Overview
--------

Spin-spiral minimization uses the ``qminimizer`` routines (``qminimizer.f90``) and,
for parallel tempering seeding, ``sx_driver.f90``. The goal is to find the lowest
energy spiral over a user-supplied :math:`\mathbf{q}` mesh and spiral basis
vectors. Good results require a physically motivated :math:`\mathbf{q}` mesh and
orthonormal spiral directions.

Spin-spiral definition
----------------------

A single-:math:`\mathbf{q}` spiral is parameterized by an in-plane unit vector
:math:`\mathbf{s}` and a normal (rotation axis) unit vector :math:`\mathbf{n}`
with :math:`\mathbf{s}\cdot\mathbf{n}=0` and :math:`\lVert\mathbf{s}\rVert=
\lVert\mathbf{n}\rVert=1`.

.. math::

   \mathbf{m}(\mathbf{r}) = \mathbf{n} \cos\bigl(2\pi\,\mathbf{q}\cdot\mathbf{r} + \phi\bigr)
   + \mathbf{s} \sin\bigl(2\pi\,\mathbf{q}\cdot\mathbf{r} + \phi\bigr)

Here :math:`\mathbf{q}` is the propagation vector (in reciprocal lattice units)
and :math:`\phi` is a phase offset (often zero). For cycloidal spirals
(:math:`\mathbf{n}\parallel\mathbf{q}`) set ``qm_type C``; for helical spirals
(:math:`\mathbf{n}\perp\mathbf{q}`) set ``qm_type H``. If you supply explicit
``qm_svec`` and ``qm_nvec``, the code normalizes them and uses them directly.

Modes and drivers
-----------------

mode / ip_mode Q
   Line-search minimization over the provided :math:`\mathbf{q}` list (``qpoints``
   grid or explicit ``qm_qvec``). Uses ``sweep_q2`` and writes
   ``qm_sweep.<simid>.out`` and ``qm_minima.<simid>.out``.

ip_mode Y
   Alternative line search (``sweep_q3``) for multi-:math:`\mathbf{q}` scans or
   when seeding a 3Q texture via ``plot_q3``.

ip_mode Z
   Cube search (``sweep_cube``) that samples a Cartesian :math:`q_x,q_y,q_z` mesh
   when the propagation direction is unknown.

ip_mode SX
   Parallel tempering setup: builds 1Q and 3Q minima via ``qminimizer``, plus
   random and ferromagnetic replicas, before entering ``sx_iphase``.

Keyword reference
-----------------

qm_qvec
   Explicit :math:`\mathbf{q}`-vectors (three components per line) that override
   the default ``qpoints`` mesh. Provide one line per :math:`\mathbf{q}`; the
   minimizer scans them in order.

qm_svec
   Spiral in-plane unit vector(s). Defines the spiral plane together with
   ``qm_nvec``. If omitted, vectors are constructed from ``qm_type``.

qm_nvec
   Spiral normal (rotation axis). Must be orthogonal to ``qm_svec`` for a pure
   helix or cycloid; normalized on read.

qm_type
   Orientation helper: ``C`` = cycloidal (:math:`\mathbf{n}\parallel\mathbf{q}`),
   ``H`` = helical (:math:`\mathbf{n}\perp\mathbf{q}`). Default ``N`` uses your
   provided ``qm_svec``/``qm_nvec`` unchanged.

qm_rot
   If ``Y``, rotate the existing texture atom-by-atom instead of constructing a
   perfect spiral. Each atom :math:`i` at position :math:`\mathbf{r}_i` is rotated
   by an angle :math:`2\pi\,\mathbf{q}\cdot\mathbf{r}_i` around the axis
   :math:`\mathbf{n}` from its initial moment. Useful to embed a spiral modulation
   in a disordered or non-collinear background.

qm_cellrot
   If ``Y``, rotate the texture cell-by-cell rather than atom-by-atom. For each
   atom :math:`i` belonging to unit cell :math:`I`, the rotation angle is
   :math:`2\pi\,\mathbf{q}\cdot\mathbf{R}_I` where :math:`\mathbf{R}_I` is the
   position of the first atom in cell :math:`I`. All atoms in the same unit cell
   rotate by the same angle. This preserves the internal structure of complex
   basis sets while applying a long-wavelength spiral modulation.

qm_oaxis
   Enforce a rotation axis perpendicular to the average magnetization; helpful
   when starting from a finite-moment state.

qm_relax, qm_relax_mode, qm_relax_steps, qm_relax_temp
   Optional local MC relaxation after each trial spiral. Mode ``M`` = Metropolis,
   ``H`` = heat bath. Use small ``qm_relax_temp`` for near-zero-temperature
   minimization.

Selective minimization with qm_exclude
--------------------------------------

qm_exclude
   Number of atom types to exclude from spiral rotations, followed by that many
   type indices on separate lines. Excluded atoms keep their original moments
   during ``sweep_q2`` / ``sweep_q3`` / ``sweep_cube``; only non-excluded atoms
   are rotated when constructing the trial spiral and during optional relaxation.

- Use this to pin sublattices, non-magnetic species, or boundary atoms while
  minimizing over the remaining sites.
- The exclusion mask is built from atom types read in the structure; double-check
  type numbering in your geometry file.
- Exclusions apply to both pure spiral construction and rotation-based modes
  (``qm_rot`` / ``qm_cellrot``).

Guidance on q-mesh and directions
---------------------------------

- Provide a dense, physically motivated :math:`\mathbf{q}` mesh (via ``qpoints``
  presets or explicit ``qm_qvec``) around expected ordering vectors.
- Ensure ``qm_svec`` and ``qm_nvec`` are orthonormal; non-orthogonal inputs yield
  distorted spirals and spurious minima.
- For unknown directions, start with ``ip_mode Z`` (cube search), then refine with
  a focused ``qm_qvec`` list.

Rotation modes: atomwise vs cellwise
------------------------------------

When ``qm_rot`` or ``qm_cellrot`` is ``Y``, the code rotates an existing magnetic
configuration instead of building a pure spiral from the formula
:math:`\mathbf{m}(\mathbf{r}) = \mathbf{n}\cos(2\pi\mathbf{q}\cdot\mathbf{r}+\phi)
+ \mathbf{s}\sin(2\pi\mathbf{q}\cdot\mathbf{r}+\phi)`. This is useful to apply
spiral modulations to disordered, non-collinear, or multi-sublattice backgrounds.

- **qm_rot = Y (atomwise):** Each atom :math:`i` at position :math:`\mathbf{r}_i`
  is rotated from its starting moment :math:`\mathbf{m}_i^{(0)}` by angle
  :math:`\theta_i=2\pi\,\mathbf{q}\cdot\mathbf{r}_i` around axis :math:`\mathbf{n}`:

  .. math::

     \mathbf{m}_i = \mathbf{R}(\mathbf{n},\theta_i)\,\mathbf{m}_i^{(0)},

  where :math:`\mathbf{R}(\mathbf{n},\theta)` is the Rodrigues rotation matrix.
  Each site rotates independently, so intra-cell non-collinearity is preserved but
  modulated by the spiral wavevector.

- **qm_cellrot = Y (cellwise):** For atom :math:`i` in unit cell :math:`I`, the
  rotation angle is determined by the position :math:`\mathbf{R}_I` of the first
  atom in that cell (index :math:`((i-1)/N_{\text{basis}}) \times N_{\text{basis}}
  + 1`):

  .. math::

     \mathbf{m}_i = \mathbf{R}(\mathbf{n},\,2\pi\,\mathbf{q}\cdot\mathbf{R}_I)\,
     \mathbf{m}_i^{(0)}.

  All atoms in the same unit cell rotate rigidly by the same angle. This preserves
  complex spin textures within the basis while imposing a long-wavelength envelope.
  It is especially useful for multi-atom unit cells (e.g., skyrmion or hedgehog
  lattices) where the internal texture must remain intact.

**Use cases:**

- **qm_rot:** Spin-glass-like backgrounds, random initial states, embedding spirals
  in disordered alloys.
- **qm_cellrot:** Non-collinear ordered states, skyrmion crystals, complex
  antiferromagnets, multi-sublattice helimagnets.

Minimal setup example
---------------------

.. code-block:: text

   ip_mode       Q                 ! Spin-spiral minimization in initial phase
   mode          Q                 ! Evaluate best-q spiral in measurement phase
   qpoints       F                 ! Read explicit q-points
   qm_qvec       0.00  0.00  0.25
                 0.00  0.00  0.30
                 0.00  0.00  0.35
   qm_svec       1.0  0.0  0.0     ! Spiral in-plane direction
   qm_nvec       0.0  0.0  1.0     ! Spiral normal
   qm_type       H                 ! Helical (n perpendicular to q)
   qm_relax      Y
   qm_relax_mode H
   qm_relax_steps 200
   qm_relax_temp 1e-3
   qm_exclude    1                 ! Keep atom type 2 fixed
                 2
