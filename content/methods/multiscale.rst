Multiscale atomistic-continuum simulations
===========================================

.. _input-multiscale:

Overview
--------

**Multiscale simulation** combines discrete **atomistic regions** with continuous
**finite-difference (FD) regions**, enabling efficient simulation of systems with
both atomic precision and macroscopic scales. This is particularly useful for:

* **Defects and textures**: Skyrmions, domain walls, dislocations in large samples
* **Interface phenomena**: Grain boundaries, heterostructure interfaces
* **Long-range effects**: Spin torques, stray fields over large distances
* **Stress engineering**: Mechanical deformation on continuum scales combined with
  atomistic detail

The multiscale approach reduces computational cost by treating distant regions
via coarse-graining and finite-difference approximations, while maintaining
atomic resolution in regions of interest.

Physical framework
------------------

**Spatial decomposition:**

The simulation domain is partitioned into regions with different treatment:

1. **Atomistic region** (:math:`\Omega_a`): Full atomic resolution with exchange,
   DM, and anisotropy interactions.

2. **Finite-difference (continuum) region** (:math:`\Omega_c`): Continuous magnetic
   field :math:`\mathbf{m}(\mathbf{r})` evolved via finite differences on a box mesh.
   Governed by continuous Landau-Lifshitz equation with exchange coefficients.

3. **Fully coarse-grained region** (:math:`\Omega_{cg}`): Uniform magnetization (no
   dynamics). Acts as boundary condition for atomistic+FD regions.

4. **Partially coarse-grained region** (:math:`\Omega_{pcg}`): Mixed atomic-continuum
   treatment for smooth gradients at interfaces.

5. **Damping band** (:math:`\Omega_d`): Interface region with position-dependent
   damping for stable coupling (suppresses non-physical waves).

**Magnetic equation of motion:**

The atomistic region evolves as standard spin dynamics:

.. math::

   \frac{d\mathbf{m}_i}{dt} = -\gamma_0 \mathbf{m}_i \times \left(\mathbf{H}_{\text{eff}} + \alpha \frac{d\mathbf{m}_i}{dt}\right),

where :math:`\mathbf{H}_{\text{eff}}` includes exchange, DM, anisotropy, and
**interpolated moments from the continuum region**.

The continuum region satisfies the Landau-Lifshitz equation:

.. math::

   \frac{\partial \mathbf{m}}{\partial t} = -\gamma_0 \mathbf{m} \times (\nabla^2 \mathbf{m} + \mathbf{K})

with exchange length :math:`\ell_{\text{ex}} = \sqrt{A/K}` and magnetic stiffness
:math:`A` derived from continuum exchange coefficients.

**Coupling via interpolation:**

* **Atom-to-continuum**: Atomic moments interpolate into FD grid via smooth weighting.
* **Continuum-to-atom**: FD moments interpolate back to atomic sites via local
  interpolation operator.

This two-way coupling is computed on each timestep.

Multiscale configuration file
-----------------------------

Multiscale setup is controlled via a **multiscale configuration file** (`.conf`),
specified in the main `inpsd.dat` via:

::

   multiscale_setup_file    multiscale.conf

The configuration file defines:

* Spatial dimensions, universe size, periodic boundaries
* Atomistic regions, continuum regions, damping band widths
* Unit cell definition and atomic sublattice
* Exchange and DM interaction files for atomistic and continuum
* Finite-difference discretization (box mesh size)
* Interpolation and damping band parameters

**Example multiscale configuration file structure:**

.. code-block:: text

   ! Geometry
   dimension               2                    ! 1D, 2D, or 3D
   periodic_boundary       F   T   F            ! Periodic in x, y, z
   universe_size           180.0   58.888   1.0
   finitediff_boxes        90      30       1   ! FD mesh: 90x30x1 boxes

   ! Regions
   coarse_grained_width        0.0
   part_coarse_grained_width   2.0
   damping_band_width          5.0
   padding_width               1.1

   ! Atomic lattice
   atom_lattice_spacing        1.0
   unitcell_atoms              1.0   1.732   0.0
   1   0.25   0.433   0.0   1.0   1.0   0.0   0.0   ! (type, x, y, z, mx, my, mz, mag)

   ! Exchange interactions
   exchange_atoms
   1 1 /  0.5   0.866  0.0  /  0.418  ! (type1, type2 / vector / Jij)
   
   exchange_coarse
   ! (coefficients for continuum region)

   ! Dzyaloshinsky-Moriya interactions
   dm_atoms
   1 1 /  0.5   0.866  0.0  /  0.057  -0.033  0.0  ! DM vector

   dm_coarse
   ! (continuum DM tensor)

Configuration keywords
----------------------

**Geometry and simulation box**

dimension
   Spatial dimension: 1, 2, or 3. Controls how many spatial coordinates are used.

periodic_boundary
   Three logical flags (F/T) for periodic boundaries in x, y, z.

universe_size
   Three reals (Å): Size of simulation box in x, y, z directions.

**Finite-difference discretization**

finitediff_boxes
   Three integers: Number of FD boxes in x, y, z. Determines spatial resolution
   of continuum region. Each box has size :math:`\approx \text{universe\_size} / \text{boxes}`.

atom_lattice_spacing
   Real (Å): Lattice constant or inter-atomic spacing. Used for scaling atomic
   interactions to continuum.

**Regional decomposition**

coarse_grained_width
   Real (Å): Width of fully coarse-grained region (uniform magnetization).
   Set to 0 to disable. Often used as artificial boundary condition.

part_coarse_grained_width
   Real (Å): Width of partially coarse-grained transition region between atomistic
   and fully coarse-grained. Smooth gradient transition.

damping_band_width
   Real (Å): Width of damping band interface (position-dependent damping for stability).

padding_width
   Real (Å): Width of padding region around atomistic atoms. Ensures continuum
   mesh sufficient coverage.

**Atomic unit cell and lattice**

unitcell_atoms
   Three reals: Lattice vectors in Å defining unit cell.

   Then, one line per atom in unit cell:
   
   ::

      type   x   y   z   mx   my   mz   magnitude

   * ``type``: Atomic species (integer)
   * ``x, y, z``: Fractional coordinates in unit cell
   * ``mx, my, mz``: Initial moment direction
   * ``magnitude``: Moment magnitude in :math:`\mu_B`

**Atomistic-region shapes (optional)**

atomistic_shape
   Defines geometric shapes of atomistic regions. Can specify multiple shapes.

   ::

      box    center_x   center_y   center_z   /   size_x   size_y   size_z
      cylinder   center_x   center_y   center_z   /   radius   height   axis

hole_shapes
   Defines shapes to be **excluded** from atomistic region (carved out).

**Magnetic interactions: Exchange**

exchange_atoms
   Exchange interactions within atomistic region. Each line:

   ::

      type1   type2   /   dx   dy   dz   /   Jij

   where ``(dx, dy, dz)`` is the displacement vector and :math:`J_{ij}` is the
   exchange constant (meV). Multiple vectors define different neighbors.

exchange_coarse
   Exchange coefficients for continuum region. Format varies; typically
   isotropic exchange :math:`A` (meV·Å²) or tensor components.

**Magnetic interactions: Dzyaloshinsky-Moriya (DM)**

dm_atoms
   DM vectors in atomistic region:

   ::

      type1   type2   /   dx   dy   dz   /   Dx   Dy   Dz

   where :math:`\mathbf{D} = (D_x, D_y, D_z)` (meV).

dm_coarse
   DM contribution to continuum Hamiltonian. Typically specified as
   :math:`3 \times 3` tensor or vector components.

continuum_dm
   Alternative specification for continuum DM as tensor:

   ::

      continuum_dm   D_xx   D_xy   D_xz   /   D_yx   ...   ...   /   ...   ...   D_zz

continuum_exchange_coef
   Three reals or single real: Exchange stiffness :math:`A` (meV·Å²) or
   directional components :math:`(A_x, A_y, A_z)` for anisotropic exchange.

continuum_moment_magnitude
   Real (:math:`\mu_B`): Moment magnitude in continuum region.

**Advanced options**

link_error_tolerance
   Real (default 1e-5): Tolerance for interpolation weight normalization.
   Lower → more accurate interpolation, higher cost.

damping_band_window_size
   Three integers: Window size (in FD boxes) for damping band interpolation.
   Controls smoothness of damping transition.

damping_band_strength
   Real (default 0): Multiplicative factor for damping in band. Set to 0
   for no damping band suppression.

**Spin-transfer torque (optional)**

stt_vector
   Three reals: STT (spin-polarized current) vector direction and magnitude.

stt_window_size
   Three reals (Å): Spatial extent of STT application region.

**Anisotropy (optional)**

anisotropy_atoms
   (Similar format to exchange) Specifies local anisotropy axis and energy.

anisotropy_coarse
   Anisotropy in continuum region.

**Zones and moment regions (optional)**

moment_regions
   Define spatial regions with different moment magnitudes or initialization.

atom_zones
   Define atomic regions for site-specific parameters.

Input file keywords (inpsd.dat)
--------------------------------

The following keywords in ``inpsd.dat`` control multiscale operation:

ipmode
   Must be ``MS`` or ``N`` (no initial phase). Multiscale is incompatible
   with other initial phase modes (``F``, ``Q``, ``H``, ``R``).

SDEalgh
   Algorithm: Only ``1`` (Midpoint) and ``5`` (Depondt) supported in multiscale.
   Others will cause error.

multiscale_setup_file
   Filename of multiscale configuration file (see above).

do_prnmultiscale
   Print multiscale output files (Y/N): ``msregions.<simid>.out``,
   ``interface.<simid>.out``, ``dband.<simid>.out``, etc.

Output files
------------

When ``do_prnmultiscale Y``:

msregions.<simid>.out
   Index ranges of atoms in each region:

   ::

      # fully coarse grained atoms
      1    10
      # partially coarse grained atoms
      11   50
      # real atoms
      51   500
      # damped atoms
      501  520
      # non-damped atoms
      521  600

coord.<simid>.out
   Atomic coordinates from multiscale setup (3 columns: x, y, z).

interface.<simid>.out
   Interpolation weights between atomistic and continuum regions.

dband.<simid>.out
   Damping band interpolation weights.

exchange.<simid>.out
   Exchange interactions extracted/generated from multiscale configuration.

dm.<simid>.out
   DM interactions extracted from multiscale configuration.

gradlink.<simid>.out
   Spin-transfer torque gradient links (if STT enabled).

Physical examples
-----------------

**Skyrmion with stress field:**

Simulate a skyrmion on atomic scale while the far-field stress is computed
continuously via FD. Atomistic region: 60×60 Å (high resolution). Continuum:
180×180 Å. Damping band: 10 Å transition.

Configuration:
::

   dimension            2
   universe_size        180.0   180.0   1.0
   finitediff_boxes     90      90       1
   atomistic_shape
   box   90   90  0.5  /  60   60   1
   coarse_grained_width       0
   part_coarse_grained_width  10
   damping_band_width         10
   continuum_exchange_coef    0.627
   continuum_dm               0.0   0.099   0 / -0.099   0   0 / 0  0  0

**Domain wall dynamics:**

1D wall pinned by disorder, driven by magnetic field. Atomistic 100 nm window.
Continuum :math:`1\,\mu\mathrm{m}`. Record wall velocity and detailed spin texture.

**Interface and grain boundary:**

2D heterostructure with two ferromagnetic layers. Atomistic resolution at
interface (interatomic interactions matter). Continuum in bulk regions.

Computational considerations
-----------------------------

**Cost:**

* Multiscale adds ~2–5× overhead vs. pure atomistic (extra interpolation, FD solver)
* Cost scales as :math:`\sim N_{\text{atoms}} + N_{\text{FD boxes}}^{d}`
* For 2D/3D: FD cost dominates if box mesh is very fine

**Memory:**

* Atomic arrays: standard
* FD arrays: :math:`N_{\text{FD boxes}} \times 3` (one moment per box)
* Interpolation weights: :math:`N_{\text{atoms}} \times N_{\text{neighbor boxes}}`

**Stability:**

* Damping band is essential for stable coupling
* STT window size should be smaller than atomistic region
* Link error tolerance: default 1e-5 usually sufficient

**Timestep:**

Typically 0.1–1 fs (same as pure atomistic). Stability limited by FD spatial
discretization and time integration scheme.

Related keywords and cross-references
-------------------------------------

- ``ipmode``: Initial phase mode; must be ``MS`` or ``N`` with multiscale
  (see :doc:`../input/simulation`)
- ``SDEalgh``: Spin dynamics algorithm; only ``1`` and ``5`` compatible
  (see :doc:`../input/simulation`)
- ``temp``: Temperature for stochastic forces (applies to atomistic atoms only)
- ``damping``: Gilbert damping (applies to atomistic atoms; damping band modulates)
- ``do_prnmultiscale``: Print multiscale diagnostic output files

References
==========

See the centralized :doc:`../references` for full bibliographic entries:

- [Evans2014]_ - Atomistic spin model simulations of magnetic nanomaterials
- [Mendez2020]_ - Multiscale approach for magnetization dynamics: unraveling exotic magnetic states of matter
