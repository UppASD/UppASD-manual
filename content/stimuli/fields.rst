Time- and space-dependent fields
================================

Parameters for applying time and space-dependent magnetic fields
-----------------------------------------------------------------

The UppASD code supports a comprehensive suite of **time and space-dependent magnetic field
profiles** that can be applied to the system during simulations. These include pulsed
fields with various temporal shapes, monochromatic and broadened microwave fields with
spatial and temporal modulation, moving field profiles, and demagnetization effects.

--------------------------------------------------
Field overview and use cases
--------------------------------------------------

Time and space-dependent fields enable:

- **Magnetic field pulses**: Ultra-short field transients with various shapes (square, exponential, Gaussian, polynomial decay)
- **Microwave excitation**: Monochromatic, broadband, and amplitude-modulated microwave fields at frequencies up to THz
- **Spatial field inhomogeneity**: Gaussian, circular, and square-shaped field profiles localized in space
- **Moving field patterns**: Dynamically positioned fields that traverse the sample (e.g., field-driven domain wall motion)
- **Demagnetization effects**: Back-action of the sample's magnetization on applied fields via shape anisotropy
- **Frequency broadening**: Spectral richness for simulating laboratory pulse shapes with finite bandwidth

Applications include:

- **Spin dynamics**: Response of magnetic systems to femtosecond laser-induced magnetization changes
- **Microwave ferromagnetic resonance (FMR)**: Frequency-dependent resonance behavior
- **Domain wall motion**: Driving and tracking magnetic domain wall propagation
- **Magnon dynamics**: Creating and studying coherent spin wave excitations
- **Thermally-assisted switching**: Field-assisted magnetization reversal at elevated temperatures

--------------------------------------------------
Static and pulsed field overview
--------------------------------------------------

The external field is calculated as the superposition of:

1. **Global static field** ``hfield``: Uniform field applied everywhere
2. **Site-dependent fields**: Field variations by atomic type or position
3. **Pulsed field**: Time-dependent field envelope ``bpulsefield(t)``
4. **Microwave fields**: Oscillating fields at specified frequencies
5. **Spatial fields**: Position-dependent field envelopes (Gaussian, circular, square)
6. **Moving fields**: Spatially and temporally varying field profiles
7. **Demagnetization field**: Self-consistent back-action of magnetization

The total field at atom ``i`` at time ``t`` is:

.. math::

   \mathbf{B}(i,t) = \mathbf{B}_0 + \mathbf{B}_{site}(i) + \mathbf{B}_{\text{pulse}}(t) + \mathbf{B}_{\text{MW}}(i,t) + \mathbf{B}_{\text{spatial}}(i,t) + \mathbf{B}_{\text{demag}}(t)

--------------------------------------------------
Magnetic field pulses
--------------------------------------------------

Field pulses are time-dependent field envelopes defined by the parameter ``do_bpulse``.

do_bpulse
   Enable magnetic field pulse. Values: ``0``=disabled, ``1``=exponential, ``2``=Gaussian,
   ``3``=polynomial decay, ``4``=square pulse, ``5``=site-dependent static, ``6``=atom-dependent static.
   Default: ``0``.

bpulsefile
   File containing pulse parameters (amplitude, time window, shape details). Required if ``do_bpulse > 0``.
   Default: ``bpulsefile``.

**Field pulse file format** (``bpulsefile``):

.. code-block:: text

   ! Header
   B0_x B0_y B0_z         ! Initial amplitude (Tesla)
   dt_pulse                ! Time step for pulse (seconds)
   nstep_pulse             ! Number of pulse time steps
   npar                    ! Number of parameters for pulse shape
   par1 par2 ... par_npar  ! Shape-dependent parameters

**Pulse shapes:**

1. **Exponential pulse** (``do_bpulse = 1``): Smooth onset and offset with constant plateau

   .. math::

      B(t) = B_0 \begin{cases}
      A \cdot \exp\left(\frac{t - t_{\text{start}}}{\tau_1}\right) & t \le t_{\text{start}} \\
      A & t_{\text{start}} < t < t_{\text{end}} \\
      A \cdot \exp\left(\frac{t_{\text{end}} - t}{\tau_2}\right) & t \ge t_{\text{end}}
      \end{cases}

   Parameters: ``[t_rise, t_start, t_end, t_fall, B_rise, B_max]``

2. **Gaussian pulse** (``do_bpulse = 2``): Symmetric smooth pulse

   .. math::

      B(t) = B_0 \cdot A \cdot \exp\left(-\frac{(t - t_{\text{center}})^2}{2\sigma^2}\right)

   Parameters: ``[t_center, \sigma, B_peak, B_sustain]``

3. **Polynomial decay** (``do_bpulse = 3``): Power-law decay envelope

   .. math::

      B(t) = B_0 \cdot B \cdot (t - t_0)^n \exp(-\alpha t)

   Parameters: ``[t_0, t_max, n, \alpha, B_rise, B_max]``

4. **Square pulse** (``do_bpulse = 4``): Instantaneous rise and fall

   .. math::

      B(t) = B_0 \begin{cases}
      0 & t < t_{\text{start}} \\
      A & t_{\text{start}} \le t \le t_{\text{end}} \\
      0 & t > t_{\text{end}}
      \end{cases}

   Parameters: ``[dummy, t_start, t_end, dummy, dummy, A]``

--------------------------------------------------
Monochromatic microwave fields
--------------------------------------------------

Oscillating fields at a specific frequency can represent FMR excitation or
resonant spin-wave coupling:

mwf
   Enable monochromatic microwave field. Modes: ``Y``=global, ``S``=site-dependent, ``W``=weighted site,
   ``P``=printing, ``I``=printing with intensity, ``N``=disabled. Default: ``N``.

mwfampl
   Amplitude of the monochromatic microwave field (Tesla). Default: ``0.0`` T.

mwffreq
   Frequency of the monochromatic microwave field (GHz). Default: ``0.0`` GHz.

mwfdir
   Direction vector of the microwave field: ``mwfdir <x> <y> <z>``. Default: ``1.0 0.0 0.0``.

mwf_site_file
   File with atom indices (one per line) specifying where the microwave field applies for site-dependent modes.

mwf_pulse_time
   Duration in simulation steps for the monochromatic microwave field; 0 means applied for the entire simulation.

The monochromatic field is evaluated as:

.. math::

   \mathbf{B}_{\text{MW}}(t) = B_0 \cos(2\pi f t) \hat{\mathbf{n}}

where :math:`B_0` is amplitude, :math:`f` is frequency, and :math:`\hat{\mathbf{n}}` is direction.

--------------------------------------------------
Frequency-broadened microwave fields
--------------------------------------------------

For more realistic microwave pulses with finite spectral width:

mwf_gauss
   Enable Gaussian-broadened microwave field. Modes analogous to ``mwf``. Default: ``N``.

mwf_gauss_ampl
   Amplitude (Tesla). Default: ``0.0``.

mwf_gauss_freq
   Center frequency (GHz). Default: ``0.0``.

mwf_gauss_time_sigma
   Frequency width (GHz). Default: ``1.0``.

mwf_gauss_dir
   Direction vector. Default: ``1.0 0.0 0.0``.

mwf_gauss_site_file
   Atom index file for site-dependent broadened microwave fields.

mwf_gauss_pulse_time
   Duration in steps; 0 for entire simulation.

The frequency spectrum is generated via inverse Fourier transform to obtain time-domain
waveform with spectral content defined by Gaussian envelope.

--------------------------------------------------
Static Gaussian-shaped fields
--------------------------------------------------

Spatially localized field profiles for studying local excitations:

do_gauss
   Enable static Gaussian-shaped spatial field profile (Y=yes, P=printing, N=disabled). Default: ``N``.

gauss_spatial_ampl
   Amplitude (Tesla). Default: ``0.0``.

gauss_spatial_sigma
   Spatial width ``<x> <y> <z>`` in Å. Default: ``1.0 1.0 1.0``.

gauss_site_file
   File specifying Gaussian center position.

gauss_pulse_time
   Duration in steps; 0 means entire simulation.

The spatial profile is:

.. math::

   B_{\text{gauss}}(\mathbf{r}) = B_0 \exp\left(-\frac{(x-x_0)^2}{2\sigma_x^2} - \frac{(y-y_0)^2}{2\sigma_y^2} - \frac{(z-z_0)^2}{2\sigma_z^2}\right)

--------------------------------------------------
Moving field profiles
--------------------------------------------------

Fields that translate through the sample to drive domain walls or other localized excitations:

**Moving Gaussian field:**

mov_gauss
   Enable moving Gaussian field (Y=yes, P=printing, N=disabled). Default: ``N``.

mov_gauss_ampl
   Amplitude (Tesla). Default: ``0.0``.

mov_gauss_space_sigma
   Spatial extent ``<x> <y> <z>`` in Å. Default: ``1.0 1.0 1.0``.

mov_gauss_file
   Trajectory waypoint file (Cartesian coords, one per line).

mov_gauss_step
   Steps between position updates. Default: ``100``.

mov_gauss_pulse_time
   Total duration in steps; 0 means entire simulation.

**Moving circular field:**

mov_circle
   Enable moving circular field (Y=yes, P=printing, N=disabled). Default: ``N``.

mov_circle_ampl
   Amplitude (Tesla). Default: ``0.0``.

mov_circle_radius
   Radius in Å. Default: ``1.0``.

mov_circle_file
   Trajectory file for circular motion.

mov_circle_step
   Steps between position updates. Default: ``100``.

**Moving square field:**

mov_square
   Enable moving square/cubic field (Y=yes, P=printing, N=disabled). Default: ``N``.

mov_square_ampl
   Amplitude (Tesla). Default: ``0.0``.

mov_square_dimensions
   Side lengths ``<x> <y> <z>`` in Å. Default: ``1.0 1.0 1.0``.

mov_square_file
   Trajectory file for cubic motion.

--------------------------------------------------
Moving microwave fields
--------------------------------------------------

Combining oscillating microwave and static spatial/moving profiles:

mwf_mov_gauss
   Moving microwave Gaussian field (Y=yes, P=printing, N=disabled). Oscillating field with moving Gaussian envelope. Default: ``N``.

mwf_mov_gauss_ampl
   Amplitude (Tesla). Default: ``0.0``.

mwf_mov_gauss_freq
   Oscillation frequency (GHz). Default: ``0.0``.

mwf_mov_gauss_time_sigma
   Frequency broadening (GHz). Default: ``1.0``.

mwf_mov_gauss_space_sigma
   Spatial extent ``<x> <y> <z>`` in Å. Default: ``1.0 1.0 1.0``.

mwf_mov_gauss_file
   Trajectory file. Default: ``mwf_mov_gauss_file``.

mwf_mov_gauss_step
   Position update interval (steps). Default: ``100``.

Similarly, ``mwf_mov_circle`` and ``mwf_mov_square`` provide moving microwave versions of
circular and square field profiles with all analogous parameters.

**Spatial Gaussian microwave field:**

mwf_gauss_spatial
   Enable spatially-modulated broadened microwave field (Y=yes, P=printing, N=disabled). Combines Gaussian frequency broadening with spatial Gaussian envelope. Default: ``N``.

mwf_gauss_spatial_ampl
   Amplitude (Tesla). Default: ``0.0``.

mwf_gauss_spatial_freq
   Center frequency (GHz). Default: ``0.0``.

mwf_gauss_spatial_time_sigma
   Frequency width (GHz). Default: ``1.0``.

mwf_gauss_spatial_space_sigma
   Spatial extent ``<x> <y> <z>`` in Å. Default: ``1.0 1.0 1.0``.

--------------------------------------------------
Demagnetization field
--------------------------------------------------

Self-consistent back-action of the system's magnetization on applied fields via shape anisotropy:

demag
   Enable demagnetization field calculation (Y=yes, N=no). Default: ``N``.

demag1, demag2, demag3
   Include demagnetization in x, y, z directions (Y=yes, N=no). Defaults: ``N``.

demagvol
   Effective volume for demagnetization calculations. Default: ``1.0``.

The demagnetization field is calculated as:

.. math::

   \mathbf{B}_{\text{demag}} = -\frac{\mu_0}{\text{Vol}} \left(\mathbf{M}_{\text{total}} \odot \hat{\mathbf{n}}\right)

where :math:`\mathbf{M}_{\text{total}}` is the total magnetization, :math:`\odot` is
element-wise multiplication with the direction mask, and the negative sign indicates
the field opposes magnetization (destabilizing in the magnetized direction).

For a thin film with normal ``z``:
- Only ``demag3='Y'`` is needed (out-of-plane demagnetization)
- :math:`B_{\text{demag}} = -\frac{\mu_0}{V} M_z`

For a bulk sample, all three components may be relevant depending on magnetization direction.

--------------------------------------------------
Output and printing options
--------------------------------------------------

Several flags control whether field data is printed to files:

prn_mwf
   Print monochromatic microwave field to output (Y/N). Default: ``N``.

prn_gauss
   Print static Gaussian field to output (Y/N). Default: ``N``.

prn_mwf_gauss
   Print broadened microwave field to output (Y/N). Default: ``N``.

prn_mov_gauss
   Print moving Gaussian field to output (Y/N). Default: ``N``.

prn_mwf_mov_gauss
   Print moving microwave Gaussian to output (Y/N). Default: ``N``.

When enabled, field data is written per time step for analysis and visualization.

--------------------------------------------------
Site-dependent field files
--------------------------------------------------

Several field types support site-dependent specification via external files:

**Format: ``mwf_site_file``, ``mwf_gauss_site_file``, etc.**

.. code-block:: text

   <atom_index_1>
   <atom_index_2>
   ...

Each line specifies one atom (1-indexed) that receives the field. Atoms not listed
receive zero field (for site-dependent modes S/W).

**Format: Trajectory files** (``mov_gauss_file``, ``mov_circle_file``, etc.)

.. code-block:: text

   <x_1> <y_1> <z_1>
   <x_2> <y_2> <z_2>
   ...

Each line specifies a waypoint in Cartesian coordinates (Ångströms). The field is
interpolated between waypoints over successive position-update intervals.

--------------------------------------------------
Example: Microwave ferromagnetic resonance
--------------------------------------------------

Simulating FMR by applying a monochromatic microwave field:

**inpsd.dat:**

.. code-block:: text

   !! Global static field
   hfield           0.1  0.0  0.0
   
   !! Monochromatic microwave field
   mwf              Y
   mwfampl          0.001
   mwffreq          10.0
   mwfdir           0.0  1.0  0.0
   mwf_pulse_time   0

With ``mwfampl = 0.001`` T and ``mwffreq = 10`` GHz, this simulates small-angle FMR
excitation perpendicular to the applied field. The frequency-dependent absorption
can be studied by sweeping ``mwffreq`` or ``hfield`` across resonance.

--------------------------------------------------
Example: Femtosecond laser-induced switching
--------------------------------------------------

Simulating ultrafast magnetization dynamics from laser pulse:

**inpsd.dat:**

.. code-block:: text

   !! Pulse shape (exponential with rapid rise/fall)
   do_bpulse        1
   bpulsefile       laser.pulse
   
   !! Demagnetization due to sample shape
   demag            Y
   demag3           Y
   demagvol         100.0

**File: laser.pulse**

.. code-block:: text

   ! Laser-induced field parameters
   0.0  0.0  0.5    ! B0: 0.5 Tesla in z-direction
   0.1              ! dt_pulse: 0.1 fs time step
   100              ! nstep: 100 steps = 10 fs total
   6                ! 6 parameters for exponential shape
   0.05 0.5 9.95 0.5 0.1 1.0

This applies a 0.5 T field pulse in the z-direction with exponential envelope,
simulating ultrafast laser-induced (de)magnetization.

--------------------------------------------------
Example: Domain wall motion with moving field
--------------------------------------------------

Driving a domain wall using a moving field profile:

**inpsd.dat:**

.. code-block:: text

   !! Moving Gaussian field
   mov_gauss        Y
   mov_gauss_ampl   0.05
   mov_gauss_space_sigma  5.0  5.0  5.0
   mov_gauss_file   trajectory.txt
   mov_gauss_step   10
   mov_gauss_pulse_time 0

**File: trajectory.txt**

.. code-block:: text

   0.0   0.0   10.0
   10.0  0.0   10.0
   20.0  0.0   10.0
   30.0  0.0   10.0
   40.0  0.0   10.0

This moves a Gaussian field along the x-direction at z = 10 Å with spatial width 5 Å.
The field pushes the domain wall forward. Adjusting ``mov_gauss_ampl`` and ``mov_gauss_step``
controls wall velocity and possible depinning phenomena.

--------------------------------------------------
Technical considerations
--------------------------------------------------

**Frequency sampling:**
- Microwave fields are evaluated at each time step, so maximum frequency accessible is limited by Nyquist: :math:`f_{\max} = 1 / (2 \Delta t)`
- For ``timestep = 1 fs`` (typical), :math:`f_{\max} \approx 0.5` PHz (infrared regime)
- Typical GHz microwave: :math:`f = 10` GHz is well-resolved

**Spatial resolution:**
- Field profiles are evaluated at atom positions
- Gaussian sigma should be >~ 2× atomic spacing for smooth modulation
- Smaller sigma values give sharper field gradients

**Computational cost:**
- Pulsed and microwave fields add minimal cost (simple evaluations)
- Moving fields require position tracking and interpolation
- Demagnetization field requires global magnetization reduction (~O(Natom))

**Units:**
- Magnetic field: Tesla (T)
- Frequency: Gigahertz (GHz) (internally converted to rad/s)
- Time: Femtoseconds (fs)
- Length: Ångströms (Å)
- Volume: Ų for demagnetization

--------------------------------------------------

References
----------

See the centralized :doc:`../references` for full bibliographic entries:

- [Pereiro2014]_ - Topological properties in magnetic topological insulators

See also
--------

- :doc:`../input/system` (system geometry for demagnetization)
- :doc:`../input/hamiltonian` (static exchange interactions)
- :doc:`../input/core-files` (core file formats: ``posfile``, ``momfile``, etc.)
