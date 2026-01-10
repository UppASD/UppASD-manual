Correlation functions
=====================

Parameters for measuring spin correlations and dynamical structure factors
--------------------------------------------------------------------------

The UppASD code can calculate **space- and time-displaced correlation functions**
and their Fourier transforms to obtain the **dynamical structure factor**
:math:`S(\mathbf{q},\omega)` [Bergman2010]_, which describes magnon dispersions
and magnetic excitations in the simulated system. Additional correlation types
(displacement, velocity, angular momentum) can also be sampled for multiphysics
simulations.

--------------------------------------------------
Spin correlation functions: definition and uses
--------------------------------------------------

The **connected time- and space-displaced spin correlation function** is:

.. math::

   C^k(\mathbf{r}-\mathbf{r'},t) = \langle m^k_{\mathbf{r}}(t) m^k_{\mathbf{r'}}(0) \rangle - \langle m^k_{\mathbf{r}}(t) \rangle \langle m^k_{\mathbf{r'}}(0) \rangle

where :math:`k \in \{x,y,z\}` is a Cartesian component, and angular brackets denote
ensemble averages. The **dynamical structure factor** :math:`S^k(\mathbf{q},\omega)`
is obtained by Fourier transforming both space and time:

.. math::

   S^k(\mathbf{q},\omega) = \frac{1}{\sqrt{2\pi}N} \sum_{\mathbf{r},\mathbf{r'}} e^{i\mathbf{q}\cdot(\mathbf{r}-\mathbf{r'})} \int_{-\infty}^{\infty} e^{i\omega t} C^k(\mathbf{r}-\mathbf{r'},t) dt

This function describes:

- **Magnon dispersion relations** :math:`\omega(\mathbf{q})` (peaks in :math:`S(\mathbf{q},\omega)`)
- **Magnon linewidths** (widths of peaks, related to magnon lifetimes)
- **Density of magnon states** (integrated over :math:`\mathbf{q}`)
- **Scattering intensities** in neutron or x-ray experiments [Squires1978]_

The **static correlation function** (equal-time correlation :math:`C^k(\mathbf{r},0)`)
and its Fourier transform :math:`S^k(\mathbf{q})` reveal:

- **Magnetic ordering vectors** (wavevector :math:`\mathbf{q}_0` at which :math:`S(\mathbf{q})` peaks)
- **Correlation lengths** (width of :math:`S(\mathbf{q})` near the peak)
- **Phase transitions** (changes in ordering as a function of temperature)

.. tip::

   **Input keywords:** For a complete list of all correlation measurement parameters
   (``do_sc``, ``do_sr``, ``sc_step``, ``qpoints``, etc.), see the comprehensive
   reference in :doc:`input-keywords-observables`.

--------------------------------------------------
Sampling strategy and resolution
--------------------------------------------------

Proper sampling of :math:`S(\mathbf{q},\omega)` requires careful choice of:

1. **Temporal sampling period** ``sc_step``: Controls the maximum frequency accessible

   .. math::

      \omega_{\max} = \frac{\pi}{\text{timestep} \times \text{sc_step}}

2. **Number of temporal samples** ``sc_nstep``: Controls frequency resolution

   .. math::

      \Delta \omega = \frac{2\omega_{\max}}{N_{\text{step}}} = \frac{2\pi}{N_{\text{step}} \times \text{timestep} \times \text{sc_step}}

3. **Total sampling time**: Determines lowest frequency accessed

   .. math::

      t_{\max} = \text{sc_nstep} \times \text{timestep} \times \text{sc_step}

   .. math::

      \omega_{\min} = \frac{2\pi}{t_{\max}} = \frac{2\pi}{N_{\text{step}} \times \text{timestep} \times \text{sc_step}}

4. **q-point mesh**: Determines spatial resolution and Brillouin zone coverage

Spin correlation keywords
-------------------------

do_sc
   Spin correlation sampling mode. Choose **Y** for full :math:`S(\mathbf{q},\omega)` and :math:`S(\mathbf{q})` sampling, **Q** for frequency-dependent :math:`S(\mathbf{q},\omega)` only, **T** for time-dependent :math:`S(\mathbf{q},t)` only, **C** for static :math:`S(\mathbf{q})` only, and **N** to disable (default).

do_sr
   Sample static correlation :math:`G(\mathbf{r})` in real space directly (Y=yes, N=no). When enabled, :math:`C(\mathbf{r},0)` is computed without requiring a q-point mesh.

sc_step
   Number of MD/MC steps between temporal correlation samples; controls maximum frequency :math:`\omega_{\max} = \pi / (\text{timestep} \times \text{sc_step})`. Default: **10**.

sc_nstep
   Number of temporal samples collected for each correlation measurement. Controls frequency resolution :math:`\Delta \omega = 2\pi / (N_{\text{step}} \times \text{timestep} \times \text{sc_step})`. Default: **100**.

sc_sep
   Number of steps between independent correlation measurements (used when ``do_sc = C``). Default: **1**.

qpoints
   Q-point mesh generation method: **F**=read Cartesian coordinates from file, **G**=read from ``qpoints.dat``, **A**=automatic generation (not implemented), **C**=full reciprocal cell, **D**=read fractional coordinates, **B**=read with weights. Default: **F**.

qfile
   Name of file containing q-point vectors. Default: **qpoints**.

do_sc_proj
   Sublattice-projected spin correlation (Y=yes, C=static only, N=no). Computes :math:`S(\mathbf{q},\omega)` per atomic type. Default: **N**.

do_sc_projch
   Chemical-species-projected spin correlation (Y=yes, C=static only, N=no). Similar to ``do_sc_proj`` but groups by chemical species. Default: **N**.

do_qt_traj
   Print time-dependent equal-time correlation :math:`S(\mathbf{q},t)` to file (Y=yes, N=no). Only for ``do_sc = C``. Default: **N**.

do_sc_local_axis
   Project correlations onto local quantization axis (Y=yes, N=no). Measures longitudinal/transverse components (:math:`S_\parallel`, :math:`S_\perp`) instead of Cartesian components. Default: **N**.

sc_local_axis_mix
   Mixing parameter for updating the local quantization axis (relevant when ``do_sc_local_axis = Y``). Range 0.0–1.0; default **0.05**.

sc_window_fun
   Windowing function for time-domain signal before Fourier transform: **1**=rectangular, **2**=Hann, **3**=Hamming, **4**=Blackman-Harris. Default: **1**.

do_sc_dosonly
   Compute magnon density of states (DOS) only, skip full :math:`S(\mathbf{q},\omega)` output (Y=yes, N=no). Default: **N**.

do_connected
   Include only the connected part of the correlation (Y=yes, N=no). Default: **Y**.

Alternative correlation types
------------------------------

In addition to spin correlations, UppASD can measure correlations of other observables:

do_uc
   Lattice displacement correlation (Y=yes, N=no). For spin-lattice coupled systems, measures ionic displacement–displacement correlations. Default: **N**.

do_ur
   Sample displacement correlation in real space directly (Y=yes, N=no). Default: **N**.

do_vc
   Velocity correlation of ions (Y=yes, N=no). Useful for phonon spectroscopy studies. Default: **N**.

do_vr
   Sample velocity correlation in real space directly (Y=yes, N=no). Default: **N**.

do_lc
   Angular momentum correlation (Y=yes, N=no). Measures orbital angular momentum fluctuations in materials with strong spin-orbit coupling. Default: **N**.

do_lr
   Sample angular momentum correlation in real space directly (Y=yes, N=no). Default: **N**.

Energy-resolved sampling
-------------------------

As an alternative to explicit temporal sampling, the energy resolution can be specified directly:

sc_emax
   Maximum energy for :math:`S(\mathbf{q},\omega)` sampling (in mRy). If specified, this value replaces the frequency range computed from ``sc_step`` and ``timestep``. Default: **0** (use ``sc_step`` instead).

sc_eres
   Energy resolution for :math:`S(\mathbf{q},\omega)` (in mRy). If specified, ``sc_nstep`` is computed as :math:`N_{\text{step}} = \lfloor E_{\max} / E_{\text{res}} \rfloor`. Default: **0** (use ``sc_nstep`` instead).

When both ``sc_emax`` and ``sc_eres`` are provided, ``sc_step`` and ``sc_nstep`` are automatically calculated from the Nyquist-Shannon theorem.

Adiabatic Magnon Spectra (AMS)
--------------------------------

The **Adiabatic Magnon Spectra** provide a fast estimate of magnon dispersions from linear spin wave theory, computed directly from the Hamiltonian exchange interactions without requiring long-time simulations.

do_ams
   Enable Adiabatic Magnon Spectra calculation (Y=yes, N=no). Computes magnon dispersions from the Hamiltonian for collinear magnetic structures. Default: **N**.

do_magdos
   Compute magnon density of states from AMS (Y=yes, N=no, A=read from file). If **Y**, the DOS is integrated over the Brillouin zone after computing :math:`\omega(\mathbf{q})`. Default: **N**.

magdos_freq
   Number of frequency points for MDOS calculation. Typical: **200**. Default: **200**.

magdos_sigma
   Gaussian broadening (in meV) applied to MDOS peaks. Typical: 10–50 meV. Default: **30.0** meV.

Autocorrelation sampling
-------------------------

For certain applications (e.g., studying specific defect sites or impurities), **local autocorrelations** can be measured on individual atoms:

do_autocorr
   Enable autocorrelation sampling at selected sites (Y=yes, N=no). Requires a file specifying atom indices. Default: **N**.

acfile
   File specifying atom indices for autocorrelation measurements (one index per line). Used only if ``do_autocorr = Y``. Default: **autocorr.dat**.

--------------------------------------------------
Q-point file formats
--------------------------------------------------

**Format for Cartesian coordinates** (``qpoints = F`` or ``qpoints = G``):

.. code-block:: text

   <number of q-points>
   <q_x_1>  <q_y_1>  <q_z_1>
   <q_x_2>  <q_y_2>  <q_z_2>
   ...

**Format for fractional coordinates** (``qpoints = D``):

.. code-block:: text

   <number of q-points>
   <q1>  <q2>  <q3>
   ...

where :math:`\mathbf{q} = q_1 \mathbf{b}_1 + q_2 \mathbf{b}_2 + q_3 \mathbf{b}_3` with
:math:`\mathbf{b}_i` the reciprocal lattice vectors.

**Format with weights** (``qpoints = B`` or ``qpoints = I``):

.. code-block:: text

   <number of q-points>
   <q_x_1>  <q_y_1>  <q_z_1>  <weight_1>
   <q_x_2>  <q_y_2>  <q_z_2>  <weight_2>
   ...

Weights are used in DOS calculations to account for symmetry-reduced Brillouin zone sampling.

--------------------------------------------------
Output files
--------------------------------------------------

When correlation sampling is enabled, the code produces:

- ``sqw.<simid>.out``: Dynamic structure factor :math:`S(\mathbf{q},\omega)` (full 3D array)
- ``sqt.<simid>.out``: Time-dependent structure factor :math:`S(\mathbf{q},t)` (for ``do_sc = T`` or ``Y``)
- ``sq.<simid>.out``: Static structure factor :math:`S(\mathbf{q})` (for ``do_sc = C`` or ``Y``)
- ``dos.<simid>.out``: Magnon density of states (if enabled)
- ``ams.<simid>.out``: Adiabatic Magnon Spectra from AMS calculation
- Real-space files (if ``do_sr = Y``): Spatial correlations :math:`G(\mathbf{r})`

--------------------------------------------------
Example: Computing magnon dispersion
--------------------------------------------------

Typical ``inpsd.dat`` settings for computing the magnon dispersion relation:

.. code-block:: text

   !! Enable full S(q,omega) sampling
   do_sc        Y
   
   !! Frequency sampling: 100 points, spacing 10 MD steps
   sc_nstep     100
   sc_step      10
   
   !! Q-point mesh from file
   qpoints      F
   qfile        qpoints.txt
   
   !! Optional: subtypes resolved
   do_sc_proj   Y
   
   !! Windowing to reduce spectral leakage
   sc_window_fun 2

With ``timestep = 1e-15`` s, this gives:

- :math:`\omega_{\max} = \pi / (10^{-15} \times 10) = 3.14 \times 10^{14}` rad/s ≈ 50 meV
- :math:`\Delta \omega = 2 \times 50 / 100 = 1.0` meV (frequency resolution)

--------------------------------------------------
Example: Static correlation and phase diagram
--------------------------------------------------

To track magnetic ordering across a phase transition:

.. code-block:: text

   !! Sample only static S(q)
   do_sc        C
   do_sr        Y
   
   !! Frequent measurements to track transition
   sc_sep       1
   do_qt_traj   Y
   
   !! q-point mesh covering full BZ
   qpoints      C
   qfile        qpoints.txt
   
   !! Optional: monitor local quantization axis
   do_sc_local_axis  Y
   sc_local_axis_mix 0.1

--------------------------------------------------
Example: Fast magnon spectra from AMS
--------------------------------------------------

For quick estimates without long-time simulations:

.. code-block:: text

   !! Enable fast AMS calculation
   do_ams       Y
   do_magdos    Y
   
   !! High frequency resolution for DOS
   magdos_freq  400
   magdos_sigma 20.0
   
   !! Q-point mesh for dispersion
   qpoints      C
   qfile        qpoints.txt

References and further reading
------------------------------

See :doc:`references` for bibliography entries used in this page.

See also
--------

- :doc:`input-keywords-autocorrelation` (local autocorrelation sampling)
- :doc:`input-keywords-montecarlo` (Monte Carlo sampling methods)
- :doc:`input-keywords-system` (system setup and q-points)
