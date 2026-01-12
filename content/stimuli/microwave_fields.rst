.. _microwave_fields:

Microwave Fields
================

The microwave field module provides comprehensive functionality for applying time-dependent and spatially-varying electromagnetic fields to simulate realistic experimental conditions in magnetic systems.

Overview
--------

The microwave field module (``MicroWaveField``) supports multiple types of electromagnetic pulses and fields:

1. **Global Monochromatic Microwave Field** - Uniform sinusoidal oscillating field applied to the entire system
2. **Site-Dependent Monochromatic Field** - Spatially varying sinusoidal field with different amplitudes at different sites
3. **Frequency-Broadened Gaussian Microwave Field** - Monochromatic field with Gaussian frequency envelope (FT-limited pulse)
4. **Spatially Gaussian-Shaped Frequency-Broadened Field** - Frequency-broadened field with additional Gaussian spatial distribution
5. **Static Gaussian Pulse** - Spatially localized static magnetic pulse with Gaussian profile
6. **Moving Gaussian Pulse** - Gaussian-shaped static field that moves along a predefined trajectory
7. **Moving Circular Pulse** - Circular-shaped static field moving along a trajectory
8. **Moving Cubic Pulse** - Cubic-shaped static field moving along a trajectory
9. **Moving Microwave Fields** - Time-dependent Gaussian/circular/cubic shaped fields moving along trajectories

Global Monochromatic Microwave Field
-------------------------------------

**Purpose**: Apply a uniform sinusoidal oscillating magnetic field to the entire system.

**Input Parameters**:

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``mwf``
     - Character (Y/P/I/S/W/N)
     - Enable monochromatic microwave field. Y=enabled, P=phase lag, I=instantaneous, S=site-dependent, W=weighted, N=disabled
   * - ``mwfampl``
     - Real
     - Amplitude of the microwave field (in mT)
   * - ``mwffreq``
     - Real
     - Frequency of the microwave oscillation (in GHz)
   * - ``mwfdir``
     - Real(3)
     - Direction vector of the microwave field (x, y, z components)
   * - ``mwf_pulse_time``
     - Integer
     - Duration of the pulse in time steps
   * - ``mwf_site_file``
     - String
     - File path for site-dependent field amplitudes (only for S/W modes)
   * - ``site_phase``
     - Character (Y/N)
     - Enable site-dependent phase modulation
   * - ``site_phase_file``
     - String
     - File path for site-dependent phase values

**Output Parameters**:

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``prn_mwf``
     - Character (Y/N)
     - Enable output of microwave field values
   * - ``mwf_step``
     - Integer
     - Sampling interval for field output
   * - ``mwf_buff``
     - Integer
     - Buffer size for field data storage

**Output File Format**:

When ``prn_mwf=Y``, the file ``mwf.simid.out`` is written containing:

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{step}`
     - :math:`B_x`
     - :math:`B_y`
     - :math:`B_z`

where :math:`B_i` are the field components in units matching the input amplitude.

**Input File Format** (for site-dependent fields):

The ``mwf_site_file`` should contain atom indices and relative amplitudes:

.. code-block:: text

    # Atom_index  Amplitude_multiplier
    1     1.0
    2     0.8
    3     0.6
    ...

Frequency-Broadened Gaussian Microwave Field
---------------------------------------------

**Purpose**: Apply a monochromatic field with Gaussian temporal envelope for FT-limited pulses.

**Input Parameters**:

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``mwf_gauss``
     - Character (Y/P/S/W/N)
     - Enable frequency-broadened microwave field. Y=enabled, P=phase lag, S=site-dependent, W=weighted, N=disabled
   * - ``mwf_gauss_ampl``
     - Real
     - Amplitude of the frequency-broadened field
   * - ``mwf_gauss_freq``
     - Real
     - Central frequency of the Gaussian envelope (in GHz)
   * - ``mwf_gauss_time_sigma``
     - Real
     - Temporal width (sigma) of the Gaussian envelope (in time units)
   * - ``mwf_gauss_dir``
     - Real(3)
     - Direction vector (x, y, z components)
   * - ``mwf_gauss_pulse_time``
     - Integer
     - Duration of the pulse in time steps
   * - ``mwf_gauss_site_file``
     - String
     - File path for site-dependent field amplitudes (S/W modes)

**Output Parameters**:

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``prn_mwf_gauss``
     - Character (Y/N)
     - Enable output of frequency-broadened field
   * - ``mwf_gauss_step``
     - Integer
     - Sampling interval
   * - ``mwf_gauss_buff``
     - Integer
     - Buffer size

**Output File Format**:

When ``prn_mwf_gauss=Y``, the file ``mwf_gauss.simid.out`` is written in the same format as monochromatic fields.

Spatially Gaussian-Shaped Frequency-Broadened Field
----------------------------------------------------

**Purpose**: Frequency-broadened field with Gaussian spatial localization around multiple center points.

**Input Parameters**:

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``mwf_gauss_spatial``
     - Character (Y/P/N)
     - Enable spatially gaussian-shaped frequency-broadened field
   * - ``mwf_gauss_spatial_ampl``
     - Real
     - Amplitude of the spatial Gaussian field
   * - ``mwf_gauss_spatial_freq``
     - Real
     - Central frequency (in GHz)
   * - ``mwf_gauss_spatial_time_sigma``
     - Real
     - Temporal Gaussian width
   * - ``mwf_gauss_spatial_space_sigma``
     - Real(3)
     - Spatial Gaussian widths (sigma_x, sigma_y, sigma_z) in Angstroms
   * - ``mwf_gauss_spatial_pulse_time``
     - Integer
     - Pulse duration in time steps
   * - ``mwf_gauss_spatial_site_file``
     - String
     - File path defining spatial centers

**Spatial Center File Format**:

The ``mwf_gauss_spatial_site_file`` should contain center coordinates:

.. code-block:: text

    # Center_x  Center_y  Center_z
    0.0    0.0    0.0
    5.0    5.0    5.0
    ...

Static Gaussian Pulse
----------------------

**Purpose**: Spatially localized static magnetic pulse with Gaussian profile.

**Input Parameters**:

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``do_gauss``
     - Character (Y/P/N)
     - Enable static Gaussian pulse
   * - ``gauss_spatial_ampl``
     - Real
     - Amplitude of the Gaussian pulse
   * - ``gauss_spatial_sigma``
     - Real(3)
     - Spatial widths (sigma_x, sigma_y, sigma_z) in Angstroms
   * - ``gauss_pulse_time``
     - Integer
     - Duration of the pulse in time steps
   * - ``gauss_site_file``
     - String
     - File path defining Gaussian center locations

**Output Parameters**:

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``prn_gauss``
     - Character (Y/N)
     - Enable output of Gaussian field
   * - ``gauss_step``
     - Integer
     - Sampling interval
   * - ``gauss_buff``
     - Integer
     - Buffer size

Moving Gaussian Pulse (Static)
------------------------------

**Purpose**: Spatially localized static Gaussian pulse that moves along a predefined trajectory.

**Input Parameters**:

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``mov_gauss``
     - Character (Y/P/N)
     - Enable moving Gaussian pulse
   * - ``mov_gauss_file``
     - String
     - File path containing trajectory coordinates
   * - ``mov_gauss_step``
     - Integer
     - Time steps between trajectory updates
   * - ``mov_gauss_pulse_time``
     - Integer
     - Total pulse duration
   * - ``mov_gauss_ampl``
     - Real
     - Amplitude of the moving pulse
   * - ``mov_gauss_space_sigma``
     - Real(3)
     - Spatial widths (sigma_x, sigma_y, sigma_z)

**Trajectory File Format**:

The ``mov_gauss_file`` should contain center coordinates at different times:

.. code-block:: text

    # Time_step  Center_x  Center_y  Center_z
    0      0.0    0.0    0.0
    100    1.0    1.0    0.0
    200    2.0    2.0    0.0
    ...

**Output Parameters**:

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``prn_mov_gauss``
     - Character (Y/N)
     - Enable output
   * - ``mov_gauss_pstep``
     - Integer
     - Sampling interval
   * - ``mov_gauss_buff``
     - Integer
     - Buffer size

Moving Circular Pulse
---------------------

**Purpose**: Circular-shaped static field moving along a trajectory.

**Input Parameters**:

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``mov_circle``
     - Character (Y/P/N)
     - Enable moving circular pulse
   * - ``mov_circle_file``
     - String
     - File path containing trajectory
   * - ``mov_circle_step``
     - Integer
     - Time steps between updates
   * - ``mov_circle_pulse_time``
     - Integer
     - Total pulse duration
   * - ``mov_circle_ampl``
     - Real
     - Amplitude
   * - ``mov_circle_radius``
     - Real
     - Radius of the circular region (in Angstroms)

**Output Parameters**:

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``prn_mov_circle``
     - Character (Y/N)
     - Enable output
   * - ``mov_circle_pstep``
     - Integer
     - Sampling interval
   * - ``mov_circle_buff``
     - Integer
     - Buffer size

Moving Cubic Pulse
------------------

**Purpose**: Cubic-shaped static field moving along a trajectory.

**Input Parameters**:

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``mov_square``
     - Character (Y/P/N)
     - Enable moving cubic pulse
   * - ``mov_square_file``
     - String
     - File path containing trajectory
   * - ``mov_square_step``
     - Integer
     - Time steps between updates
   * - ``mov_square_pulse_time``
     - Integer
     - Total pulse duration
   * - ``mov_square_ampl``
     - Real
     - Amplitude
   * - ``mov_square_dimensions``
     - Real(3)
     - Dimensions (dx, dy, dz) of the cubic region (in Angstroms)

**Output Parameters**:

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``prn_mov_square``
     - Character (Y/N)
     - Enable output
   * - ``mov_square_pstep``
     - Integer
     - Sampling interval
   * - ``mov_square_buff``
     - Integer
     - Buffer size

Moving Microwave Fields
-----------------------

The module supports three variants of moving time-dependent (microwave) fields with the same spatial shapes as static fields:

Moving Gaussian Microwave Field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Input Parameters**:

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``mwf_mov_gauss``
     - Character (Y/P/N)
     - Enable moving Gaussian microwave field
   * - ``mwf_mov_gauss_file``
     - String
     - File path containing trajectory
   * - ``mwf_mov_gauss_step``
     - Integer
     - Time steps between trajectory updates
   * - ``mwf_mov_gauss_pulse_time``
     - Integer
     - Total pulse duration
   * - ``mwf_mov_gauss_ampl``
     - Real
     - Amplitude
   * - ``mwf_mov_gauss_freq``
     - Real
     - Oscillation frequency (in GHz)
   * - ``mwf_mov_gauss_time_sigma``
     - Real
     - Temporal Gaussian envelope width
   * - ``mwf_mov_gauss_space_sigma``
     - Real(3)
     - Spatial widths (sigma_x, sigma_y, sigma_z)

**Output Parameters**:

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``prn_mwf_mov_gauss``
     - Character (Y/N)
     - Enable output
   * - ``mwf_mov_gauss_pstep``
     - Integer
     - Sampling interval
   * - ``mwf_mov_gauss_buff``
     - Integer
     - Buffer size

Moving Circular Microwave Field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Input Parameters**:

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``mwf_mov_circle``
     - Character (Y/P/N)
     - Enable moving circular microwave field
   * - ``mwf_mov_circle_file``
     - String
     - File path containing trajectory
   * - ``mwf_mov_circle_step``
     - Integer
     - Time steps between updates
   * - ``mwf_mov_circle_pulse_time``
     - Integer
     - Total pulse duration
   * - ``mwf_mov_circle_ampl``
     - Real
     - Amplitude
   * - ``mwf_mov_circle_freq``
     - Real
     - Oscillation frequency
   * - ``mwf_mov_circle_time_sigma``
     - Real
     - Temporal Gaussian width
   * - ``mwf_mov_circle_radius``
     - Real
     - Radius of the circular region

**Output Parameters**:

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``prn_mwf_mov_circle``
     - Character (Y/N)
     - Enable output
   * - ``mwf_mov_circle_pstep``
     - Integer
     - Sampling interval
   * - ``mwf_mov_circle_buff``
     - Integer
     - Buffer size

Moving Cubic Microwave Field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Input Parameters**:

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``mwf_mov_square``
     - Character (Y/P/N)
     - Enable moving cubic microwave field
   * - ``mwf_mov_square_file``
     - String
     - File path containing trajectory
   * - ``mwf_mov_square_step``
     - Integer
     - Time steps between updates
   * - ``mwf_mov_square_pulse_time``
     - Integer
     - Total pulse duration
   * - ``mwf_mov_square_ampl``
     - Real
     - Amplitude
   * - ``mwf_mov_square_freq``
     - Real
     - Oscillation frequency
   * - ``mwf_mov_square_time_sigma``
     - Real
     - Temporal Gaussian width
   * - ``mwf_mov_square_dimensions``
     - Real(3)
     - Dimensions (dx, dy, dz) of the cubic region

**Output Parameters**:

.. list-table::
   :widths: 25 15 50
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``prn_mwf_mov_square``
     - Character (Y/N)
     - Enable output
   * - ``mwf_mov_square_pstep``
     - Integer
     - Sampling interval
   * - ``mwf_mov_square_buff``
     - Integer
     - Buffer size

Flag Values
-----------

**Field Enable Flags** (mwf, mwf_gauss, etc.):

- ``Y`` - Enabled with standard behavior
- ``P`` - Enabled with phase lag correction
- ``I`` - Instantaneous field (monochromatic only)
- ``S`` - Site-dependent field using site file
- ``W`` - Weighted site-dependent field
- ``N`` - Disabled (default)

**Print Flags** (prn_mwf, prn_gauss, etc.):

- ``Y`` - Enable output to file
- ``N`` - Disable output (default)

Output Files
------------

The module generates output files when the respective ``prn_*`` flags are set to ``Y``. All output files follow the naming convention:

- Monochromatic field: ``mwf.simid.out``
- Frequency-broadened field: ``mwf_gauss.simid.out``
- Static Gaussian: ``gauss.simid.out``
- Moving Gaussian: ``mov_gauss.simid.out``
- Moving circular: ``mov_circle.simid.out``
- Moving cubic: ``mov_square.simid.out``
- Moving Gaussian microwave: ``mwf_mov_gauss.simid.out``
- Moving circular microwave: ``mwf_mov_circle.simid.out``
- Moving cubic microwave: ``mwf_mov_square.simid.out``
- Gaussian spatial frequency-broadened: ``mwf_gauss_spatial.simid.out``

Units and Conventions
---------------------

- **Amplitude**: mT (millitesla)
- **Frequency**: GHz
- **Spatial dimensions**: Angstroms (Å)
- **Time**: Simulation time steps
- **Sigma parameters**: Real space coordinates (Angstroms) for spatial, time steps for temporal

Example Configuration
---------------------

A typical configuration file section for microwave fields might look like:

.. code-block:: text

    # Monochromatic microwave field
    mwf              Y
    mwfampl          100.0
    mwffreq          10.0
    mwfdir           1.0  0.0  0.0
    mwf_pulse_time   10000
    prn_mwf          Y
    mwf_step         100
    mwf_buff         100

    # Frequency-broadened Gaussian microwave field
    mwf_gauss        Y
    mwf_gauss_ampl   150.0
    mwf_gauss_freq   10.5
    mwf_gauss_time_sigma  50.0
    mwf_gauss_dir    0.0  1.0  0.0
    mwf_gauss_pulse_time  5000
    prn_mwf_gauss    Y
    mwf_gauss_step   100
    mwf_gauss_buff   50

    # Moving Gaussian pulse (static)
    mov_gauss        Y
    mov_gauss_file   mov_gauss_traj.dat
    mov_gauss_step   100
    mov_gauss_pulse_time  8000
    mov_gauss_ampl   200.0
    mov_gauss_space_sigma  2.0  2.0  2.0
    prn_mov_gauss    Y
    mov_gauss_pstep  100
    mov_gauss_buff   80

Notes
-----

- Multiple field types can be enabled simultaneously and will be superposed
- The direction vectors should typically be normalized or scaled appropriately
- Site-dependent fields require corresponding input files with proper formatting
- Trajectory files for moving fields must contain sufficient data points for the entire simulation
- The temporal Gaussian envelope (``*_time_sigma``) creates FT-limited pulses with frequency bandwidth inversely proportional to the temporal width
- Spatial Gaussian widths (``*_space_sigma``) define the extent of field localization in 3D space
