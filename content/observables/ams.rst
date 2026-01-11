.. _label-ams:

========================================================================================
Adiabatic Magnon Spectra (AMS), Linear Spin Wave Theory (LSWT), and Topological Magnons
========================================================================================

Overview
========

UppASD implements **Linear Spin Wave Theory (LSWT)** for calculating magnon band structures in magnetic systems. This powerful capability enables:

* **Magnon dispersions**: Eigenfrequencies of collective spin excitations across the Brillouin zone
* **Density of states**: Thermodynamic and spectroscopic properties from magnon spectra
* **Topological analysis**: Chern numbers, Berry curvature, and thermal Hall conductivity

The implementation follows the **Tóth-Lake formalism** [TothLake2015]_ and extensions for non-collinear magnetism, supporting:

* Ferromagnetic (FM) and antiferromagnetic (AFM) systems
* Collinear and non-collinear magnetic configurations
* Heisenberg exchange, DMI, anisotropies, and external fields
* Spin spiral ground states with arbitrary ordering vectors

.. note::
   LSWT is valid at **low temperatures** (:math:`T \ll T_c`) where the Holstein-Primakoff approximation is accurate [Holstein1940]_. For finite-temperature dynamics and magnon-magnon interactions, use atomistic spin dynamics instead [Halilov1998]_.

Physical Background
===================

Holstein-Primakoff Transformation
----------------------------------

LSWT expands the spin operators around the classical ground state using bosonic creation (:math:`a_i^\dagger`) and annihilation (:math:`a_i`) operators:

.. math::

   S_i^+ &= \sqrt{2S} \sqrt{1 - \frac{a_i^\dagger a_i}{2S}} \, a_i \approx \sqrt{2S} \, a_i \\
   S_i^- &= \sqrt{2S} \, a_i^\dagger \sqrt{1 - \frac{a_i^\dagger a_i}{2S}} \approx \sqrt{2S} \, a_i^\dagger \\
   S_i^z &= S - a_i^\dagger a_i

The linearization (:math:`a_i^\dagger a_i \ll 2S`) leads to a quadratic Hamiltonian suitable for diagonalization.

Dynamical Matrix Construction
------------------------------

The effective Hamiltonian in momentum space takes the form:

.. math::

   \mathcal{H} = \sum_{\mathbf{q}} \begin{pmatrix} a_{\mathbf{q}}^\dagger & a_{-\mathbf{q}} \end{pmatrix}
   \begin{pmatrix}
   \mathcal{A}(\mathbf{q}) - \mathcal{C} & \mathcal{B}(\mathbf{q}) \\
   \mathcal{B}^\dagger(\mathbf{q}) & \mathcal{A}(-\mathbf{q}) - \mathcal{C}
   \end{pmatrix}
   \begin{pmatrix} a_{\mathbf{q}} \\ a_{-\mathbf{q}}^\dagger \end{pmatrix}

where:

* :math:`\mathcal{A}(\mathbf{q})`: Fourier-transformed exchange interactions
* :math:`\mathcal{B}(\mathbf{q})`: Magnon pairing terms (non-zero for AFM/non-collinear systems)
* :math:`\mathcal{C}`: On-site energy contributions (anisotropies, external fields)

The matrix elements include:

.. math::

   \mathcal{A}_{\mu\nu}(\mathbf{q}) &= \sqrt{m_\mu m_\nu} \, \mathbf{u}_\mu^\dagger \cdot \mathbf{J}_{\mu\nu}(\mathbf{q}) \cdot \mathbf{u}_\nu^* \\
   \mathcal{B}_{\mu\nu}(\mathbf{q}) &= \sqrt{m_\mu m_\nu} \, \mathbf{u}_\mu^\dagger \cdot \mathbf{J}_{\mu\nu}(\mathbf{q}) \cdot \mathbf{u}_\nu \\
   \mathcal{C}_{\mu\mu} &= \sum_\lambda m_\lambda \, \mathbf{v}_\mu^\dagger \cdot \mathbf{J}_{\mu\lambda}(0) \cdot \mathbf{v}_\lambda

Here :math:`\mathbf{u}_\mu, \mathbf{v}_\mu` are local coordinate vectors perpendicular and parallel to the equilibrium moment direction, and :math:`m_\mu` is the moment magnitude.

The Fourier-transformed interaction tensor is:

.. math::

   \mathbf{J}_{\mu\nu}(\mathbf{q}) = \sum_{\mathbf{R}} \mathbf{J}_{\mu\nu}(\mathbf{R}) \, e^{-i\mathbf{q}\cdot\mathbf{R}}

Colpa Diagonalization
----------------------

For a quadratic bosonic Hamiltonian :math:`\mathcal{H} = \Psi^\dagger h \Psi` with

.. math::

   h = \begin{pmatrix} A & B \\ B^* & A^* \end{pmatrix}, \quad
   g = \begin{pmatrix} I & 0 \\ 0 & -I \end{pmatrix}

the generalized eigenvalue problem :math:`h \mathbf{v} = \omega \, g \mathbf{v}` yields magnon eigenfrequencies :math:`\omega_n(\mathbf{q})`.

UppASD implements the **Cholesky-based Colpa method** [Colpa1978]_:

1. Ensure :math:`h` is positive-definite (add small diagonal :math:`\epsilon` if needed)
2. Cholesky decomposition: :math:`h = K^\dagger K`
3. Diagonalize :math:`K g K^\dagger` with eigenvalues :math:`\pm\omega_n`
4. Construct Bogoliubov transformation from eigenvectors

For systems where Cholesky fails (e.g., near phase transitions), a fallback **generalized eigenvalue solver** (LAPACK ``zggev``) is used.

Topological Magnons
-------------------

Berry Curvature
^^^^^^^^^^^^^^^

The Berry curvature :math:`\Omega_n(\mathbf{q})` of magnon band :math:`n` is calculated from the **link variables** between neighboring q-points:

.. math::

   U_{n}(\mathbf{q}, \mathbf{q}') = \frac{\langle u_n(\mathbf{q}) | u_n(\mathbf{q}') \rangle}{|\langle u_n(\mathbf{q}) | u_n(\mathbf{q}') \rangle|}

where :math:`|u_n(\mathbf{q})\rangle` are the normalized eigenvectors. For a plaquette with corners :math:`\mathbf{q}, \mathbf{q}+\delta_x, \mathbf{q}+\delta_x+\delta_y, \mathbf{q}+\delta_y`:

.. math::

   F_n(\mathbf{q}) = \log\left[ U_n(\mathbf{q}, \mathbf{q}+\delta_x) U_n(\mathbf{q}+\delta_x, \mathbf{q}+\delta_x+\delta_y) U_n(\mathbf{q}+\delta_x+\delta_y, \mathbf{q}+\delta_y) U_n(\mathbf{q}+\delta_y, \mathbf{q}) \right]

The Berry curvature is :math:`\Omega_n(\mathbf{q}) = \text{Im}[F_n(\mathbf{q})]`.

Chern Number
^^^^^^^^^^^^

The **Chern number** (topological invariant) for band :math:`n` is the integral of Berry curvature over the Brillouin zone:

.. math::

   C_n = \frac{1}{2\pi} \int_{\text{BZ}} \Omega_n(\mathbf{q}) \, d^2\mathbf{q} \approx \frac{1}{2\pi} \sum_{\text{plaquettes}} \text{Im}[F_n]

Integer Chern numbers (:math:`C_n = 0, \pm 1, \pm 2, \ldots`) classify topologically distinct magnon bands. Non-zero Chern numbers indicate:

* **Topological magnon insulators**: Gapped bulk with protected edge modes
* **Thermal Hall effect**: Transverse heat current under temperature gradient

Thermal Magnon Conductivity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The thermal Hall conductivity :math:`\kappa_{xy}` (W/K) is calculated from:

.. math::

   \kappa_{xy} = -\frac{k_B^2 T}{\hbar (2\pi)^2} \sum_{n} \int_{\text{BZ}} c_2(n_B[\omega_n(\mathbf{q})]) \, \Omega_n(\mathbf{q})^2 \, d^2\mathbf{q}

where :math:`n_B(\omega) = 1/(e^{\omega/k_B T} - 1)` is the Bose-Einstein distribution and :math:`c_2(x)` is defined as:

.. math::

   c_2(x) = (1+x) \left[\log\frac{1+x}{x}\right]^2 - [\log x]^2 - 2\text{Li}_2(-x)

with :math:`\text{Li}_2` the dilogarithm function.

Input Parameters
================

AMS/LSWT Calculation
--------------------

.. list-table::
   :widths: 20 15 15 50
   :header-rows: 1

   * - Keyword
     - Type
     - Default
     - Description
   * - ``do_ams``
     - Character(1)
     - 'N'
     - | Calculate AMS (Adiabatic Magnon Spectra) using collinear LSWT.
       | **Y**: Perform calculation
       | **N**: Skip
   * - ``do_diamag``
     - Character(1)
     - 'N'
     - | Calculate magnon dispersions using non-collinear LSWT formalism.
       | **Y**: Full non-collinear LSWT
       | **N**: Skip
       | **Note**: Use for non-collinear ground states or spin spirals
   * - ``nc_qvect``
     - Real(3)
     - 0 0 0
     - Ordering wave vector :math:`\mathbf{q}_0` for non-collinear calculations (reduced coordinates)
   * - ``nc_nvect``
     - Real(3)
     - 0 0 1
     - Spin rotation axis for non-collinear LSWT (Cartesian unit vector)
   * - ``nc_eps``
     - Real
     - :math:`10^{-6}`
     - | Diagonal offset for ensuring positive-definite dynamical matrix.
       | Increase (e.g., :math:`10^{-4}`) if Cholesky decomposition fails near instabilities.

Magnon Density of States
-------------------------

.. list-table::
   :widths: 20 15 15 50
   :header-rows: 1

   * - Keyword
     - Type
     - Default
     - Description
   * - ``do_magdos``
     - Character(1)
     - 'N'
     - | Calculate magnon density of states (DOS).
       | **Y**: Calculate and output DOS
       | **N**: Skip
       | **A**: Read DOS from file (useful for multiphase simulations)
   * - ``magdos_sigma``
     - Real
     - 1.0
     - Gaussian broadening :math:`\sigma` for DOS (meV). Typical: 0.5-2.0 meV
   * - ``magdos_freq``
     - Integer
     - 200
     - Number of frequency bins in DOS histogram
   * - ``magdos_lfreq``
     - Integer
     - 0
     - | Low-frequency cutoff (meV). 
       | **0**: Use minimum eigenfrequency
   * - ``magdos_hfreq``
     - Integer
     - 0
     - | High-frequency cutoff (meV).
       | **0**: Use maximum eigenfrequency
   * - ``magdosfile``
     - String
     - ''
     - Input file for reading DOS (used when ``do_magdos='A'``)
   * - ``magdos_rasamples``
     - Integer
     - 1
     - Number of random samples for random alloy DOS averaging

Topological Analysis
--------------------

.. list-table::
   :widths: 20 15 15 50
   :header-rows: 1

   * - Keyword
     - Type
     - Default
     - Description
   * - ``do_chern``
     - Character(1)
     - 'N'
     - | Calculate Chern numbers and Berry curvature.
       | **Y**: Perform topological analysis
       | **N**: Skip
   * - ``Nx``, ``Ny``, ``Nz``
     - Integer
     - 100, 100, 1
     - | Reciprocal space grid dimensions for Chern number integration.
       | Typical: 50-200 points per direction (convergence test required).
       | **Note**: Increase for accurate Chern numbers near topological phase transitions.
   * - ``Chern_qvect``
     - Real(3)
     - 0 0 0
     - Center of the Brillouin zone for Chern integration (typically :math:`\Gamma`-point)

Q-point Specification
---------------------

Magnon dispersions require a list of **q-points** (momentum space sampling). See :doc:`correlations` for details on q-point generation and dynamic structure factor output.

.. list-table::
   :widths: 20 15 15 50
   :header-rows: 1

   * - Keyword
     - Type
     - Default
     - Description
   * - ``qpoints``
     - Character(1)
     - 'C'
     - | Q-point generation mode:
       | **F**: Read from file (``qfile``)
       | **A**: Automatic generation (uniform grid)
       | **C**: Complete unit cell sampling
       | **D**: Automatic with symmetry reduction
   * - ``qfile``
     - String
     - './qfile'
     - | File containing q-points (when ``qpoints='F'``).
       | Format: First line = number of q-points, subsequent lines = :math:`q_x \, q_y \, q_z` (reduced coordinates)

Output Files
============

AMS Calculations (``do_ams='Y'``)
----------------------------------

1. **ams.<simid>.out**

   Magnon dispersion relations (eigenfrequencies in meV).

   Format::

      # q-point   qx       qy       qz       Band 1   Band 2   ...   Band NA
          1     0.0000   0.0000  -0.5000    12.345   12.345   ...   45.678
          2     0.0000   0.0000  -0.4900    12.456   12.456   ...   45.789
          ...

   * First column: q-point index
   * Columns 2-4: q-vector components (reduced coordinates)
   * Subsequent columns: Magnon eigenfrequencies (meV), sorted by magnitude

2. **jqams.<simid>.out**

   Alternative eigenvalues from :math:`J(\mathbf{q})` diagonalization (for comparison).

3. **evams.<simid>.out**

   Eigenvectors (magnon polarization vectors) at each q-point. Useful for analyzing magnon character.

4. **magdos.<simid>.out**

   Magnon density of states (when ``do_magdos='Y'``).

   Format::

      #          E(meV)         D(S(q,E))         Int(D(S))
            0.000000000      0.012345678      0.000000000
            0.100000000      0.023456789      0.001234567
            ...

   * Column 1: Energy (meV)
   * Column 2: DOS (states/meV)
   * Column 3: Integrated DOS (cumulative)

Non-collinear LSWT (``do_diamag='Y'``)
---------------------------------------

5. **ncams.<simid>.out**

   Magnon dispersions for q-points (phason modes).

6. **ncams+q.<simid>.out**, **ncams-q.<simid>.out**

   Magnon dispersions for :math:`\mathbf{q} + \mathbf{q}_0` and :math:`\mathbf{q} - \mathbf{q}_0` (spin spiral satellites).

7. **ncsqw.<simid>.out**, **ncsqw_intensity.<simid>.out**

   Dynamical structure factor :math:`S(\mathbf{q}, \omega)` tensorial components and scattering intensity.

Topological Analysis (``do_chern='Y'``)
----------------------------------------

8. **chern.<simid>.out**

   Chern numbers for all magnon bands and thermal Hall conductivity.

   Format::

      Band number ->   1   2   3   4   5   6   ...
      Ch_Number   ->   0   0   1  -1   0   0   ...
      Ch_Number+Q ->   0   0   1  -1   0   0   ...
      Ch_Number-Q ->   0   0  -1   1   0   0   ...
      Magnon Thermal Conductivity in W/K ->  1.234567890123456E-08

   * Row 1: Band indices
   * Row 2: Chern numbers for phason bands
   * Rows 3-4: Chern numbers for :math:`\pm\mathbf{q}_0` satellites (for non-collinear systems)
   * Row 5: Thermal Hall conductivity :math:`\kappa_{xy}` (W/K) at specified temperature

9. **bphase.<simid>.out**, **bphase+q.<simid>.out**, **bphase-q.<simid>.out**

   Berry phase :math:`\text{Im}[F_n(\mathbf{q})]` at each plaquette for all bands.

   Format::

      Band #          qx          qy          qz       Band1       Band2       ...
          ...

Practical Examples
==================

Example 1: Ferromagnetic Heisenberg Chain
------------------------------------------

Calculate magnon dispersion along the chain direction (:math:`\Gamma \to X`).

**inpsd.dat**::

   simid     HeisChain
   ncell     1  1  100
   BC        0  0  P

   posfile   ./posfile
   exchange  ./jfile
   momfile   ./momfile

   Mensemble 1
   Initmag   3

   ip_mode   S
   ip_nphase 1
   20000 1.0e-3 1e-16 4.0

   mode      S
   temp      1.0e-3
   damping   0.0010
   Nstep     20000
   timestep  1.0e-15

   qpoints   F
   qfile     ./qfile

   do_ams    Y
   do_magdos Y

**qfile** (101 q-points along chain)::

   101
   0.0  0.0  -0.50
   0.0  0.0  -0.49
   0.0  0.0  -0.48
   ...
   0.0  0.0   0.48
   0.0  0.0   0.49
   0.0  0.0   0.50

**Expected Output**: Cosine dispersion :math:`\omega(q) = 2JS[1 - \cos(qa)]` where :math:`J` is the exchange constant, :math:`S` the spin magnitude, and :math:`a` the lattice constant.

Physical Interpretation
^^^^^^^^^^^^^^^^^^^^^^^^

* **Acoustic mode**: :math:`\omega(0) = 0` (Goldstone mode from broken :math:`\text{SO}(3)` symmetry)
* **Maximum at zone boundary**: :math:`\omega_{\max} = 4JS` for :math:`q = \pi/a`
* **Quadratic near** :math:`\Gamma`: :math:`\omega(q) \approx JSa^2 q^2` (long-wavelength limit)

Example 2: Antiferromagnetic Square Lattice with DMI
-----------------------------------------------------

Compute topological magnons in a 2D AFM with Dzyaloshinskii-Moriya interaction.

**inpsd.dat**::

   simid     AFM_DMI_2D
   ncell     20  20  1
   BC        P   P   0

   posfile   ./posfile
   exchange  ./jfile
   dm        ./dmfile
   momfile   ./momfile

   Mensemble 1
   Initmag   3

   ip_mode   S
   ip_nphase 1
   50000 1.0 1e-16 4.0

   mode      S
   temp      0.001
   damping   0.01
   Nstep     50000
   timestep  1.0e-15

   qpoints   D          # Automatic with symmetry
   qfile     ./qfile

   do_ams    Y
   do_magdos Y
   magdos_sigma 0.5
   magdos_freq  500

   do_chern  Y
   Nx        100
   Ny        100
   Nz        1

**dmfile** (example DMI vectors)::

   1  1  1  0.5  0.0  0.0    # DMI along x-bonds
   1  2  1  0.0  0.5  0.0    # DMI along y-bonds

**Expected Output**: Non-trivial Chern numbers (:math:`C_n = \pm 1`) for acoustic magnon branches, indicating topological protection. Magnon band gaps open due to DMI.

Physical Interpretation
^^^^^^^^^^^^^^^^^^^^^^^^

* **Band inversion**: DMI induces avoided crossings at high-symmetry points
* **Topological edge states**: Chern numbers predict chiral edge magnons (not directly visible in bulk calculation but detectable via nanoribbon calculations)
* **Thermal Hall effect**: :math:`\kappa_{xy} \neq 0` even without external magnetic field

Example 3: Non-collinear Spin Spiral
-------------------------------------

Calculate magnon spectrum for a conical spin spiral with ordering vector :math:`\mathbf{q}_0 = (0, 0, 0.3)`.

**inpsd.dat**::

   simid     SpinSpiral_3D
   ncell     10  10  20
   BC        P   P   P

   posfile   ./posfile
   exchange  ./jfile
   momfile   ./momfile

   Mensemble 1
   Initmag   3

   ip_mode   S
   ip_nphase 1
   30000 1.0 1e-16 4.0

   mode      S
   temp      0.01
   damping   0.01
   Nstep     30000
   timestep  1.0e-15

   qpoints   F
   qfile     ./qfile

   do_diamag Y
   nc_qvect  0.0  0.0  0.3
   nc_nvect  0.0  0.0  1.0
   nc_eps    1.0e-5

   do_magdos Y
   magdos_sigma 1.0

**Expected Output**: Three output files (``ncams.<simid>.out``, ``ncams+q.<simid>.out``, ``ncams-q.<simid>.out``) corresponding to phason and satellite modes. The phason mode (``ncams.<simid>.out``) has :math:`\omega(\mathbf{q}_0) \approx 0` (soft mode at ordering vector).

Physical Interpretation
^^^^^^^^^^^^^^^^^^^^^^^^

* **Phason branch**: Goldstone mode associated with helical ordering
* **Amplitudon branches**: Gapped modes from transverse fluctuations
* **Folded Brillouin zone**: Effective BZ reduced due to incommensurate ordering

Example 4: Random Alloy - FeCo
-------------------------------

Magnon DOS for a random Fe₀.₅Co₀.₅ alloy with chemical disorder.

**inpsd.dat**::

   simid     FeCo_Alloy
   ncell     10  10  10
   BC        P   P   P

   posfile   ./posfile
   exchange  ./jfile
   momfile   ./momfile
   do_ralloy 1
   Nchmax    2

   Mensemble 1
   Initmag   3

   ip_mode   S
   ip_nphase 1
   20000 1.0e-3 1e-16 4.0

   mode      S
   temp      0.01
   damping   0.01
   Nstep     20000
   timestep  1.0e-15

   qpoints   C          # Full cell sampling
   qfile     ./qfile

   do_ams    Y
   do_magdos Y
   magdos_sigma 2.0
   magdos_rasamples 10   # Average over 10 random configurations

**Expected Output**: Broadened magnon DOS reflecting chemical disorder. Multiple magnon branches with varying exchange stiffness.

Physical Interpretation
^^^^^^^^^^^^^^^^^^^^^^^^

* **DOS broadening**: Disorder induces magnon lifetime effects (captured phenomenologically)
* **Configurational averaging**: Multiple random alloy samples reduce statistical noise
* **Effective medium**: Average DOS reflects mean-field properties of alloy

Troubleshooting
===============

Common Issues
-------------

**Problem 1: Negative or imaginary eigenfrequencies**

Symptoms::

   Warning in diamag: non-positive definite matrix in zpotrf

**Causes**:

* Magnetic configuration not in equilibrium (residual forces)
* Numerical instability near phase transition or soft modes
* Insufficient ``nc_eps`` for Cholesky stability

**Solutions**:

1. **Improve initial relaxation**:
   
   Increase ``ip_nphase`` iterations::

      ip_nphase 2
      50000 5.0 1e-16 4.0    # High-temp relaxation
      50000 0.1 1e-16 4.0    # Low-temp relaxation

2. **Increase diagonal offset**:
   
   ::

      nc_eps 1.0e-4    # Default: 1e-6

3. **Use Monte Carlo initialization**:
   
   ::

      ip_mode M
      ip_mcanneal 1
      ip_nphase 1
      50000 10.0 0.01 4.0

**Problem 2: Chern numbers are non-integer**

Symptoms::

   Ch_Number   ->   0   0   0.842  -0.937   0   0

**Causes**:

* Insufficient q-space sampling (grid too coarse)
* Berry curvature singularities at high-symmetry points
* Temperature effects (thermal fluctuations)

**Solutions**:

1. **Increase grid resolution**:
   
   ::

      Nx 200
      Ny 200

   Test convergence: Calculate Chern numbers for ``Nx = 50, 100, 150, 200`` and verify integer convergence.

2. **Check ground state quality**:
   
   Ensure magnetic configuration is fully relaxed (visualize spin structure).

3. **Lower temperature**:
   
   ::

      temp 0.001    # LSWT requires T << Tc

**Problem 3: Magnon DOS has negative values**

Symptoms::

   DOS contains small negative values near zero energy

**Causes**:

* Negative eigenfrequencies (see Problem 1)
* Insufficient Gaussian broadening

**Solutions**:

1. Fix negative eigenfrequencies (see Problem 1 solutions)

2. Increase broadening::

      magdos_sigma 2.0

**Problem 4: AMS output contains only acoustic modes**

Symptoms::

   All magnon branches start from zero at Gamma point

**Causes**:

* System has continuous symmetry (FM without anisotropy)
* Missing gap-opening interactions (anisotropy, external field, DMI)

**Solutions**:

1. **Add uniaxial anisotropy**:
   
   In ``kfile``::

      1  1  0.01  0.0  0.0  1.0    # 0.01 meV easy-axis anisotropy

2. **Apply external field**:
   
   In ``inpsd.dat``::

      hfield 0.0  0.0  1.0    # 1 T along z

**Problem 5: Cholesky decomposition fails frequently**

Symptoms::

   Multiple "non-positive definite matrix" warnings across many q-points

**Causes**:

* System near magnetic instability (soft modes)
* Strong frustration (competing interactions)
* Non-collinear ground state with ``do_ams='Y'`` (should use ``do_diamag='Y'``)

**Solutions**:

1. **Switch to non-collinear formalism**:
   
   ::

      do_ams N
      do_diamag Y
      nc_qvect 0.0 0.0 0.0
      nc_nvect 0.0 0.0 1.0

2. **Use robust solver**:
   
   The fallback ``zggev`` solver automatically activates when Cholesky fails.

3. **Increase stabilization**:
   
   ::

      nc_eps 1.0e-3

Convergence Tests
-----------------

**Q-point Sampling**

Magnon DOS and Chern numbers require dense sampling:

* **DOS**: Test convergence with ``qpoints='C'`` for increasing system sizes (``ncell``). Typical: 20×20×20 cells minimum.
* **Chern numbers**: Test ``Nx, Ny`` convergence explicitly. Integer Chern numbers should stabilize by 100×100 grid.

Example convergence script::

   for Nx in 50 75 100 125 150; do
       sed -i "s/Nx .*/Nx $Nx/" inpsd.dat
       sed -i "s/Ny .*/Ny $Nx/" inpsd.dat
       uppasd > log_Nx${Nx}.txt
       grep "Ch_Number" chern.HeisAFM.out >> chern_convergence.txt
   done

**Broadening Parameter**

The ``magdos_sigma`` controls DOS smoothness:

* Too small: Noisy, spiky DOS
* Too large: Over-smoothed, missing fine structure

Optimal choice: :math:`\sigma \approx 0.1 \times` (bandwidth). Typical: 0.5-2.0 meV.

Performance Considerations
--------------------------

**Memory Scaling**

LSWT memory usage: :math:`\mathcal{O}(N_A^2 \times N_q)` where :math:`N_A` = atoms per cell, :math:`N_q` = q-points.

Large systems (>100 atoms/cell) may require:

* Reduced q-point sampling (use symmetry reduction: ``qpoints='D'``)
* Sparse matrix techniques (future feature)

**CPU Scaling**

Diagonalization dominates: :math:`\mathcal{O}(N_A^3 \times N_q)` operations.

Parallelization:

* OpenMP: Q-point loop parallelized (export ``OMP_NUM_THREADS=8``)
* MPI: Not yet supported for LSWT (use OpenMP on single node)

**Typical Timings** (Intel Xeon, 16 cores):

* 1D chain (1 atom/cell, 100 q-points): <1 second
* 2D square (1 atom/cell, 10000 q-points): ~10 seconds
* 3D bulk (4 atoms/cell, 1000 q-points): ~60 seconds
* Topological analysis (100×100 grid): ~5 minutes

Physical Validation
====================

Analytical Benchmarks
---------------------

**Ferromagnetic Heisenberg Chain**

Exact dispersion:

.. math::

   \omega(q) = 2JS[1 - \cos(qa)]

where :math:`J` = exchange (meV), :math:`S` = spin, :math:`a` = lattice constant (Å).

**Square Lattice FM**

.. math::

   \omega(q_x, q_y) = 4JS[2 - \cos(q_x a) - \cos(q_y a)]

Compare with ``ams.<simid>.out`` output.

**AFM with Sublattice Magnetization**

Two magnon branches:

.. math::

   \omega_\pm(\mathbf{q}) = \sqrt{A(\mathbf{q})^2 - B(\mathbf{q})^2}

where :math:`A, B` are matrix elements from AMS calculation.

Comparison with Experiment
---------------------------

**Inelastic Neutron Scattering (INS)**

The magnetic cross-section is:

.. math::

   \frac{d^2\sigma}{d\Omega dE} \propto |F(\mathbf{Q})|^2 \sum_{\alpha\beta} \left(\delta_{\alpha\beta} - \frac{Q_\alpha Q_\beta}{Q^2}\right) S^{\alpha\beta}(\mathbf{q}, \omega)

UppASD output ``ncsqw_intensity.<simid>.out`` includes the polarization factor :math:`(\delta_{\alpha\beta} - Q_\alpha Q_\beta / Q^2)`.

**Spin-Polarized EELS**

High-resolution electron energy loss spectroscopy can resolve magnon dispersions in thin films. Compare peak positions in ``ams.<simid>.out`` with experimental dispersion curves.

Limitations and Extensions
==========================

Known Limitations
-----------------

1. **Linear approximation**: LSWT fails for large-amplitude fluctuations (:math:`T \sim T_c`, low-dimensional systems at finite :math:`T`)
2. **Magnon-magnon interactions**: Anharmonic effects (magnon decay, thermal broadening) not included
3. **Quantum corrections**: :math:`1/S` corrections (quantum fluctuations) require higher-order spin-wave theory
4. **Finite-T effects**: Zero-temperature formalism only (thermal renormalization of parameters needed)

When to Use Atomistic Spin Dynamics Instead
--------------------------------------------

* **Finite temperature** (:math:`T > 0.1 T_c`): Magnon-magnon scattering, temperature-dependent linewidths
* **Non-perturbative regimes**: Large DMI, strong frustration, non-collinear excitations
* **Metastable states**: Skyrmions, domain walls (topology-driven dynamics)
* **Long-time dynamics**: Magnon transport, magnon-phonon coupling

Future Features
---------------

Planned extensions:

* **Magnon-phonon coupling**: Spin-lattice dynamics (joint DOS, hybrid excitations)
* **Disorder effects**: Alloy-induced magnon scattering (beyond configurational averaging)
* **Surface/edge states**: Real-space Green's function methods for topological edge magnons
* **Time-resolved spectroscopy**: Pump-probe magnon dynamics after laser excitation

References
----------

See the centralized :doc:`../references` for full bibliographic entries:

- [TothLake2015]_ - Linear spin wave theory for single-Q incommensurate magnetic structures
- [Colpa1978]_ - Diagonalization of the quadratic boson Hamiltonian
- [Holstein1940]_ - Holstein-Primakoff transformation and ferromagnetic field dependence
- [Halilov1998]_ - Adiabatic spin dynamics from spin-density-functional theory


See Also
========

* :doc:`correlations` (Q-point generation and dynamic structure factor)
* :doc:`correlations` (Dynamical spin correlation :math:`S(\mathbf{q}, \omega)`)
* :doc:`stiffness` (Spin-wave stiffness and micromagnetic parameters)
* :doc:`../input/simulation` (Temperature control and equilibration)

For computational methods, see:

* **LAPACK Documentation**: Linear algebra routines (``zheev``, ``zpotrf``, ``zggev``)
* **SpinW Software**: Alternative LSWT package with GUI [Toth et al., *J. Phys.: Condens. Matter* (2015)]
