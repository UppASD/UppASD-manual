.. _input-keywords-temperature-gradients:

=====================================
Temperature Gradients
=====================================

.. contents::
   :local:

Overview
========

The temperature gradient feature enables simulations of spatially non-uniform temperature distributions throughout the magnetic system. This capability is essential for studying thermal transport phenomena, including:

- **Thermal transport and Seebeck effects**: Understanding heat flow in magnetic materials
- **Temperature-dependent phase transitions**: Studying how temperature variations drive magnetic phase changes
- **Magnon/magneto-thermoelectric effects**: Coupling temperature gradients to spin dynamics
- **Thermal gradients in nanostructures**: Investigating temperature effects in finite-size systems

Temperature gradients are calculated by solving the **Poisson or Laplace equation** with user-defined boundary conditions. This produces a spatially varying temperature field that is then applied to each atom during the simulation, allowing thermal properties to vary across the system.

.. note::

   The 3-Temperature Model (3TM) for electron-lattice-spin dynamics is documented separately in :doc:`input-keywords-temperature-3tm`. These features can be combined: temperature gradients can be applied to the lattice temperature in 3TM simulations.


Physics Background
==================

The temperature gradient is obtained by solving:

.. math::

   \nabla^2 T(\mathbf{r}) = S(\mathbf{r})

where:

- **Laplace equation** (S = 0): Used when no heat sources are present; temperature satisfies :math:`\nabla^2 T = 0`
- **Poisson equation** (S ≠ 0): Used with localized heat sources, e.g., Gaussian heat distribution

Boundary conditions are specified at the supercell edges (faces), allowing:

- **Constant temperature**: Fixed temperature on a boundary face
- **Linear temperature variation**: Temperature varies linearly across a boundary
- **Periodic boundaries**: Supported via crystallographic symmetry

The solver then interpolates the temperature field to each atom's position in the unit cell and supercell coordinate system.


Activation and Configuration
============================

Main Enable Flag
----------------

.. code-block:: fortran

   grad = 'N'  ! Default: no temperature gradient

The gradient flag controls three modes:

- ``grad = 'N'`` : Homogeneous temperature (constant over all atoms)
- ``grad = 'Y'`` : Calculate gradient from boundary conditions and equation type
- ``grad = 'F'`` : Read pre-calculated temperature field from file

When ``grad = 'Y'``, the system:

1. Reads parameters from a temperature configuration file (``tempfile``)
2. Sets up grid and identifies boundary atoms
3. Solves Poisson or Laplace equation
4. Interpolates temperature to each atom

When ``grad = 'F'``, a previously calculated temperature field is loaded directly from disk.


Equation Type and Solver
========================

The temperature field is computed by solving a partial differential equation specified by:

.. code-block:: fortran

   eq_type      = 'N'    ! Equation type: 'P'=Poisson or 'L'=Laplace
   temp_solver  = 1      ! Solver: 1=Finite Differences, 2=Meshless (MLS)

**Laplace Equation** (``eq_type = 'L'``)
   No heat sources; temperature satisfies harmonic condition. Used for passive thermal gradients driven by boundary conditions alone.

**Poisson Equation** (``eq_type = 'P'``)
   Includes heat source term; allows localized heating with Gaussian source profile.

**Finite Difference Solver** (``temp_solver = 1``)
   Uses standard FD stencil approximation to derivatives on nearest-neighbor lattice. Fast but requires regular atomic arrangement.

**Meshless Solver (MLS)** (``temp_solver = 2``)
   Moving Least Squares method; more flexible, can handle irregular geometries. Requires higher computational cost.


Boundary Conditions
===================

The temperature gradient is constrained by boundary conditions on the supercell faces. Six boundaries are defined corresponding to the three lattice vector directions:

**Boundary Parameter Names:**

- **X-direction faces**: ``I1_min_border`` (x=0), ``I1_max_border`` (x=max)
- **Y-direction faces**: ``I2_min_border`` (y=0), ``I2_max_border`` (y=max)
- **Z-direction faces**: ``I3_min_border`` (z=0), ``I3_max_border`` (z=max)

**Boundary Condition Types:**

Each boundary can be set to one of:

- ``'N'`` : No boundary condition (face is inactive)
- ``'constant'`` : Fixed constant temperature

  .. code-block:: fortran

     I1_min_border = 'constant'
     ! Followed by single temperature value:
     temp_I1_min = 300.0 K

- ``'linear'`` : Temperature varies linearly on the face

  .. code-block:: fortran

     I1_min_border = 'linear'
     ! Followed by two temperature values:
     ! (minimum temperature, maximum temperature)

**Temperature Parameters:**

.. list-table:: Boundary Temperature Values
   :widths: 25 15 50
   :header-rows: 1

   * - Parameter
     - Default [K]
     - Description
   * - ``temp_I1_min``
     - 0.0
     - Temperature at x=0 face (constant)
   * - ``temp_I1_max``
     - 0.0
     - Temperature at x=max face (constant)
   * - ``temp_I1_min_low``
     - 0.0
     - Lower temperature at x=0 (linear mode)
   * - ``temp_I1_min_high``
     - 0.0
     - Upper temperature at x=0 (linear mode)
   * - ``temp_I1_max_low``
     - 0.0
     - Lower temperature at x=max (linear mode)
   * - ``temp_I1_max_high``
     - 0.0
     - Upper temperature at x=max (linear mode)
   * - ``temp_I2_min``, ``temp_I2_max``
     - 0.0
     - Y-direction boundaries (same structure as I1)
   * - ``temp_I3_min``, ``temp_I3_max``
     - 0.0
     - Z-direction boundaries (same structure as I1)


Heat Source (Gaussian Profile)
===============================

When using the **Poisson equation** with ``source_type = 'G'``, a localized Gaussian heat source can be added:

.. code-block:: fortran

   source_type = 'G'           ! Gaussian heat source
   temp_max    = 500.0         ! Source amplitude [K]
   sigmatemp   = 2.0, 2.0, 2.0 ! Spatial width in x, y, z [Å]
   r_center    = 5.0, 5.0, 5.0 ! Center position [Å]

The heat source is given by:

.. math::

   S(\mathbf{r}) = T_{\text{max}} \exp\left(-\frac{(x-x_0)^2}{2\sigma_x^2} - \frac{(y-y_0)^2}{2\sigma_y^2} - \frac{(z-z_0)^2}{2\sigma_z^2}\right)

where :math:`(x_0, y_0, z_0)` is the source center and :math:`(\sigma_x, \sigma_y, \sigma_z)` control the spatial extent.

**Parameters:**

.. list-table:: Gaussian Source Parameters
   :widths: 25 15 50
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``source_type``
     - 'N'
     - Source type: 'N'=none, 'G'=Gaussian, 'P'=point
   * - ``temp_max``
     - 0.0 K
     - Amplitude of Gaussian source
   * - ``sigmatemp(1:3)``
     - 0.0 Å
     - Gaussian width in x, y, z directions
   * - ``r_center(1:3)``
     - 0.0 Å
     - Center position of Gaussian source


Grid and Solver Configuration
=============================

The temperature gradient calculation requires a computational grid based on atomic neighbor relationships:

.. list-table:: Grid Configuration Parameters
   :widths: 25 15 50
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``crys_symm``
     - 1
     - Crystallographic symmetry (0=none, 1-3=increasing symmetry)
   * - ``dim_sys``
     - 3
     - System dimensionality (1=1D, 2=2D, 3=3D)
   * - ``grid_type``
     - 'N'
     - Grid type for FD: 'N'=standard, 'O'=optimal
   * - ``init_temp``
     - 1
     - Initial guess for iterative solver (flag)

For Meshless solver (``temp_solver = 2``), additional parameters control the MLS shape functions:

- Radial basis weight type
- Polynomial basis dimension
- Influence radius
- Distance metric parameters


Temperature File Format
=======================

When ``grad = 'Y'``, a temperature configuration file (default name: ``tempfile``) specifies the problem setup:

**File Structure:**

.. code-block:: fortran

   ! Example temperature configuration file
   % Temperature gradient configuration
   shells_nums
   12                                  ! Number of neighbor shells

   num_tot_neigh
   96                                  ! Total neighbor connections

   temp_solver
   1                                   ! Solver: 1=FD, 2=MLS

   crys_symm
   1                                   ! Symmetry level

   dim_sys
   3                                   ! System dimensionality

   eq_type
   L                                   ! Equation: P=Poisson, L=Laplace

   init_temp
   1

   grid_type
   N                                   ! Grid type

   x_min_border
   constant                            ! Boundary type
   300.0                               ! Temperature value [K]

   x_max_border
   constant
   400.0                               ! Temperature value [K]

   y_min_border
   N                                   ! Inactive

   y_max_border
   N

   z_min_border
   N

   z_max_border
   N

   source_type
   N                                   ! No heat source

   load_temp
   temperature_initial.dat            ! File with pre-calculated field (if grad='F')

   temp_neigh
   % Neighbor distance vectors
   % Format: site_type  neighbor_site_type  distance_vector(x,y,z)
   1  1   0.0  1.0  0.0
   1  1   0.0 -1.0  0.0
   1  1   1.0  0.0  0.0
   1  1  -1.0  0.0  0.0
   1  1   0.0  0.0  1.0
   1  1   0.0  0.0 -1.0
   % ... more neighbor pairs


Output Files
============

When temperature gradients are calculated (``grad = 'Y'`` or ``grad = 'F'``), output files are generated:

.. list-table:: Temperature Output Files
   :widths: 35 50
   :header-rows: 1

   * - Filename
     - Content
   * - ``temperature.<simid>.out``
     - Atom positions and calculated temperatures at each position
   * - ``temperature_initial.<simid>.out``
     - Initial temperature field (boundary conditions + source term)

**File Format:**

.. code-block:: text

   atom_id   x(Å)    y(Å)    z(Å)    T_final(K)   T_initial(K)

   1  0.0000  0.0000  0.0000   300.000   300.000
   2  1.0000  0.0000  0.0000   301.234   301.000
   3  2.0000  0.0000  0.0000   302.456   302.000
   ...


Exponential Temperature Time-Dependence
========================================

Independent from spatial gradients, the system temperature can follow an exponential cooling or heating curve over time:

.. code-block:: fortran

   do_tempexp  = 'N'        ! Enable exponential cooling/heating
   tempexp_start = 500.0 K  ! Initial temperature
   tempexp_end   = 100.0 K  ! Final temperature
   tempexp_tau   = 1.0e-9 s ! Time constant
   tempexp_step  = 100      ! Printing interval

When enabled (``do_tempexp = 'Y'`` or ``do_tempexp = 'E'``), the system temperature varies as:

.. math::

   T(t) = T_{\text{start}} e^{-t/\tau} + T_{\text{end}}(1 - e^{-t/\tau})

This provides smooth thermal transitions from one temperature to another over the course of the simulation. The temperature evolves exponentially toward the final temperature ``tempexp_end`` with time constant ``tempexp_tau``.

**Output:**

.. code-block:: fortran

   temperature_exp.<simid>.out

Format: ``step  time(s)  temperature(K)``


Keywords Reference Table
========================

grad
   Enable gradient: 'N'=none, 'Y'=calculate, 'F'=from file. Default: 'N'.

tempfile
   Name of temperature configuration file. Default: 'tempfile'.

eq_type
   Equation type: 'P'=Poisson (with sources), 'L'=Laplace (no sources). Default: 'N'.

temp_solver
   Solver selection: 1=Finite Difference, 2=Meshless (MLS). Default: 1.

source_type
   Source term: 'N'=none, 'G'=Gaussian, 'P'=point. Default: 'N'.

temp_max
   Gaussian source amplitude [K]. Default: 0.0.

sigmatemp(1:3)
   Gaussian source width in x, y, z [Å]. Default: 0.0.

r_center(1:3)
   Gaussian source center coordinates [Å]. Default: 0.0.

crys_symm
   Crystallographic symmetry level (integer). Default: 1.

dim_sys
   System dimensionality (1,2,3). Default: 3.

grid_type
   Grid type for finite-difference solver: 'N'=standard, 'O'=optimal. Default: 'N'.

init_temp
   Initial temperature guess flag for iterative solver. Default: 1.

I1_min_border, I1_max_border, I2_min_border, I2_max_border, I3_min_border, I3_max_border
   Boundary types for each face: 'N'=inactive, 'constant', or 'linear'. Default: 'N'.

temp_I1_min, temp_I1_max, temp_I2_min, temp_I2_max, temp_I3_min, temp_I3_max
   Temperature values for boundary faces (K). Defaults: 0.0.

do_tempexp
   Enable exponential temperature time-dependence: 'N'=no, 'Y'=yes, 'E'=extended. Default: 'N'.

tempexp_start
   Initial temperature for exponential schedule (K). Default: 0.0.

tempexp_end
   Final temperature for exponential schedule (K). Default: 0.0.

tempexp_tau
   Exponential time constant (s). Default: 1.0e-9.

tempexp_step
   Output printing interval (steps) for exponential schedule. Default: 100.

Linear boundary profile parameters (for ``_border = 'linear'``) use ``_low`` and ``_high`` suffixes (e.g., ``temp_I1_min_low``, ``temp_I1_min_high``).


Examples
========

Example 1: Simple 1D Temperature Gradient
------------------------------------------

Create a linear temperature gradient along the x-direction from 300 K to 500 K:

**Main Input File:**

.. code-block:: fortran

   ! Geometry
   N1 = 10   ! 10 cells in x
   N2 = 1    ! 1 cell in y
   N3 = 1    ! 1 cell in z

   ! Temperature gradient activation
   grad = 'Y'
   tempfile = 'tempfile_1d.dat'

**Temperature Configuration File (tempfile_1d.dat):**

.. code-block:: fortran

   shells_nums
   6

   num_tot_neigh
   48

   temp_solver
   1

   crys_symm
   1

   dim_sys
   1                  ! 1D problem

   eq_type
   L                  ! Laplace (no sources)

   grid_type
   N

   x_min_border
   constant
   300.0              ! Start temperature

   x_max_border
   constant
   500.0              ! End temperature

   y_min_border
   N

   y_max_border
   N

   z_min_border
   N

   z_max_border
   N

   source_type
   N

   temp_neigh
   6
   1  1  1.0  0.0  0.0
   1  1 -1.0  0.0  0.0
   ...

**Result:** Linear temperature variation from 300 K (x=0) to 500 K (x=max), with all atoms receiving spatially interpolated temperatures based on their x-coordinate.


Example 2: 2D Thermal Gradient with Gaussian Source
-----------------------------------------------------

Create a 2D temperature field with linear boundary conditions and a localized Gaussian heat source:

**Temperature Configuration:**

.. code-block:: fortran

   shells_nums
   12

   dim_sys
   2                  ! 2D problem

   eq_type
   P                  ! Poisson (with source)

   x_min_border
   linear
   250.0              ! min
   350.0              ! max

   x_max_border
   linear
   300.0
   400.0

   y_min_border
   linear
   250.0
   350.0

   y_max_border
   linear
   300.0
   400.0

   z_min_border
   N

   z_max_border
   N

   source_type
   G                  ! Gaussian source

   ! Source parameters (placed in center of cell)
   temp_max           ! Amplitude
   300.0
   sigmatemp
   2.0  2.0  10.0     ! Narrow in x,y; wide in z (unused for 2D)
   r_center
   5.0  5.0  0.0      ! Center in cell coordinates

**Effect:** Base temperature gradient set by boundaries, superimposed with Gaussian heating around the cell center.


Example 3: Exponential Cooling During Simulation
--------------------------------------------------

Cool the system exponentially from 400 K to 100 K over the simulation time:

.. code-block:: fortran

   ! Temperature control
   do_tempexp = 'Y'
   tempexp_start = 400.0 K   ! Start temperature
   tempexp_end = 100.0 K     ! Final temperature
   tempexp_tau = 1.0e-9 s    ! 1 nanosecond time constant
   tempexp_step = 100        ! Output every 100 steps

   ! Can be combined with spatial gradient
   grad = 'Y'                ! Also use spatial gradient
   tempfile = 'tempfile.dat'

**Result:** Temperature field varies both spatially (from gradient) and temporally (exponential decay). Each atom experiences a temperature that is the product of the gradient field and the exponential scaling function.


Example 4: Meshless (MLS) Solver for Irregular Geometry
--------------------------------------------------------

Use the meshless solver for a system with non-standard atomic arrangement:

.. code-block:: fortran

   ! Solver selection
   temp_solver = 2            ! Meshless solver (MLS)

   eq_type = 'L'              ! Laplace equation

   ! Boundary conditions (same as FD)
   I1_min_border = 'constant'
   temp_I1_min = 300.0

   I1_max_border = 'constant'
   temp_I1_max = 400.0

**Configuration in tempfile:**

.. code-block:: fortran

   temp_solver
   2                          ! Use MLS solver

The meshless solver automatically adapts to the atomic positions without requiring a regular grid.


Physical Interpretations
========================

**Heat Flow and Transport**

The temperature gradient induces:

- **Magnon thermal transport**: Heat-driven magnon currents in magnetic systems
- **Spin Seebeck effect**: Temperature gradient generates spin current
- **Anomalous Nernst effect**: Temperature gradient leads to transverse magnetization changes

**Boundary Condition Choices**

- **Constant boundaries**: Model thermostatic reservoirs (e.g., thermal baths)
- **Linear boundaries**: Represent linearly varying thermal contact
- **Gaussian sources**: Model localized heating (e.g., laser focus, electrical heater)

**Equation Selection**

- **Laplace (eq_type='L')**: Passive thermal diffusion; steady-state heat distribution
- **Poisson (eq_type='P')**: Active heat generation; includes internal heat sources


Performance Considerations
==========================

- **Finite Difference Solver** (temp_solver=1): Fast, scales linearly with system size
- **Meshless Solver** (temp_solver=2): More flexible but ~10-100× slower; use for complex geometries
- **Temperature file I/O**: Gradient calculation is done once at initialization; negligible ongoing cost
- **Memory**: Stores temperature array (1 float per atom) plus temporary working arrays during setup

**Grid Density:** The temperature field quality depends on the neighbor shell structure. Ensure sufficient shells (``shells_nums`` ≥ 6) for smooth interpolation.


Combining with Other Features
==============================

Temperature gradients can be combined with:

- **3-Temperature Model (3TM)**: Apply gradient to lattice temperature in electron-lattice-spin dynamics
- **Exponential Cooling**: Simultaneous spatial and temporal temperature variations
- **Spin Torques**: Seebeck effect with temperature-dependent STT coefficients
- **Spin-Transfer Torque**: Temperature-driven current injection effects
- **Anisotropy and Exchange**: Temperature-dependent Hamiltonian parameters


Troubleshooting
===============

**Issue: Solver does not converge**

- Reduce grid complexity; ensure ``crys_symm`` and ``dim_sys`` match problem geometry
- Increase ``temp_solver=2`` (meshless) if system is highly non-regular
- Check boundary condition definitions; conflicting boundaries cause non-convergence

**Issue: Unphysical temperature values**

- Verify boundary conditions are reasonable (typically 50-500 K for magnetic systems)
- For Gaussian sources, check that ``temp_max`` doesn't exceed feasible values
- Ensure source center ``r_center`` is within system bounds

**Issue: Temperature file I/O errors**

- Verify ``tempfile`` path is correct and file is readable
- Check ``load_temp`` filename when using ``grad='F'``
- Ensure temperature file uses consistent coordinate units (Ångströms)


References
==========

Temperature gradient calculations in UppASD are based on finite difference and meshless methods for solving PDEs:

1. Chico, J., et al. (2018) "Temperature gradients in magnetic systems"
2. Moving Least Squares Method: Lancaster, P. & Salkauskas, K. (1981)
3. Laplace equation solutions: Griffiths, D. J. "Introduction to Electrodynamics" (analogous problems)
