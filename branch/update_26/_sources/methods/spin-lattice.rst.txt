Spin-lattice dynamics (SLD) simulations
=======================================

.. _input-spinlattice:

Overview
--------

**Spin-lattice dynamics (SLD)** extends pure spin dynamics (SD) by treating the
atomic lattice (ionic coordinates) as dynamical degrees of freedom. This allows
simulation of coupled spin and lattice evolution, including:

* Thermal expansion and contraction from magnetic excitations
* Magnetoelastic coupling: how spin textures modulate lattice deformations
* Magnetostriction: inverse magnetoelastic effect (magnetic field → strain)
* Ultrafast demagnetization with lattice heating (combined with 3TM)
* Phonon-magnon interactions mediating spin dynamics

The SLD module couples spin-magnetic Hamiltonians with lattice-vibrational
Hamiltonians through **spin-lattice interaction terms** (ML and MML).

.. tip::

   **Input keywords:** For a complete list of lattice dynamics measurement parameters
   (``do_lavrg``, ``lavrg_step``, ``do_ltottraj``, ``do_ld``, ``do_velrsc``, etc.),
   see the comprehensive reference in :doc:`../input/observables`.

Physical Hamiltonian
--------------------

The complete SLD Hamiltonian comprises:

.. math::

   \mathcal{H}_{\text{SLD}} = \mathcal{H}_{\text{magnetic}} + \mathcal{H}_{\text{lattice}} + \mathcal{H}_{\text{SL}},

where:

**Magnetic Hamiltonian** (:math:`\mathcal{H}_{\text{mag}}`)
   Standard spin-only terms (exchange, anisotropy, applied field):

   .. math::

      \mathcal{H}_{\text{mag}} = -\sum_{ij} J_{ij}\mathbf{m}_i\cdot\mathbf{m}_j
      - \sum_i K_i (m_i^z)^2 - \sum_i \mathbf{H}_{\text{ext}}\cdot\mathbf{m}_i.

**Lattice (ionic) Hamiltonian** (:math:`\mathcal{H}_{\text{lat}}`)
   Purely ionic forces and elastic energy:

   .. math::

      \mathcal{H}_{\text{lat}} = \sum_i \frac{1}{2} m_{\text{ion}} v_i^2 + \mathcal{V}_{\text{LL}},

   where :math:`\mathcal{V}_{\text{LL}}` is the **lattice-lattice (LL) potential**:

   .. math::

      \mathcal{V}_{\text{LL}} = \frac{1}{2}\sum_{ij} \mathbf{u}_i \cdot \Phi_{ij} \cdot \mathbf{u}_j.

   Here :math:`\mathbf{u}_i` is the atomic displacement from equilibrium, and
   :math:`\Phi_{ij}` is the **force constant tensor** (second derivative of lattice potential).

**Spin-lattice coupling** (:math:`\mathcal{H}_{\text{SL}}`)
   Interaction between magnetic moments and lattice distortions:

   .. math::

      \mathcal{H}_{\text{SL}} = \sum_i \mathbf{u}_i \cdot \mathbf{F}_{i}^{\text{ML}} + \sum_{ij}
      (\mathbf{m}_i \cdot \mathbf{m}_j) \mathbf{u}_{ij} \cdot \Phi_{ij}^{\text{MML}} \cdot \mathbf{u}_{ij}.

   Two types of coupling:

   * **ML (Magneto-Lattice)**: Linear in magnetic moment, couples magnetic field to lattice strain
   * **MML (Magneto-Magnetic-Lattice)**: Quadratic in moments, modulates exchange with displacement

Equations of motion
-------------------

SLD integrates coupled second-order ODEs for atomic positions:

.. math::

   m_{\text{ion}} \ddot{\mathbf{u}}_i = -\frac{\partial \mathcal{H}}{\partial \mathbf{u}_i} - \gamma_{\text{lat}} \dot{\mathbf{u}}_i + \mathbf{F}_i^{\text{th}},

where:

* :math:`\mathbf{F}_i = -\frac{\partial \mathcal{H}}{\partial \mathbf{u}_i}` is the effective lattice force
* :math:`\gamma_{\text{lat}}` is the lattice damping (friction coefficient)
* :math:`\mathbf{F}_i^{\text{th}} = \sqrt{2 m_{\text{ion}} \gamma_{\text{lat}} k_B T} \, \boldsymbol{\xi}_i(t)`
  is the thermal Langevin force (white noise)

The effective force has three components:

.. math::

   F_i^{\alpha} = \underbrace{-\frac{\partial \mathcal{V}_{\text{LL}}}{\partial u_i^\alpha}}_{\text{LL force}} - 
   \underbrace{\frac{\partial \mathcal{V}_{\text{ML}}}{\partial u_i^\alpha}}_{\text{ML force}} - 
   \underbrace{\frac{\partial \mathcal{V}_{\text{MML}}}{\partial u_i^\alpha}}_{\text{MML force}}.

**LL force** (lattice-lattice):

   .. math::

      F_i^{\text{LL}} = -\sum_j \Phi_{ij} \cdot \mathbf{u}_j,

   standard harmonic force from neighboring ionic displacements.

**ML force** (magneto-lattice):

   .. math::

      F_i^{\text{ML}} = -\sum_j M_{ij}^{\alpha\beta} m_j^\alpha,

   where :math:`M_{ij}` is the magneto-lattice tensor coupling moment component
   :math:`\alpha` to lattice direction :math:`\beta`.

**MML force** (magneto-magnetic-lattice):

   .. math::

      F_i^{\text{MML}} = -\frac{1}{2}\sum_{ij} (\mathbf{m}_i \cdot \mathbf{m}_j) \cdot \Phi_{ij}^{\text{MML}} \cdot \mathbf{u}_{ij},

   couples exchange interaction strength modulation to lattice distortion.

Integration scheme: Generalized Langevin Verlet (GJF)
-----------------------------------------------------

UppASD uses the **generalized Langevin Verlet (GJF) integrator**, a robust
second-order scheme that preserves the canonical ensemble at constant temperature:

.. math::

   \mathbf{u}(t+\Delta t) &= \mathbf{u}(t) + \Delta t \mathbf{v}(t) + \frac{(\Delta t)^2}{2m}
   \left[\mathbf{F}(t) - \gamma \mathbf{v}(t) + \mathbf{F}_{\text{th}}(t)\right],\\
   \mathbf{v}(t+\Delta t) &= \frac{\mathbf{u}(t+\Delta t) - \mathbf{u}(t-\Delta t)}{2\Delta t}.

The scheme is implemented via:

1. ``u_gjfverlet``: Update positions using force at time :math:`t`
2. ``v_gjfverlet``: Update velocities using forces at :math:`t` and :math:`t+\Delta t`

Thermal noise is generated at each step by ``lattrannum`` with proper scaling by
:math:`\sqrt{T(t)}` for temperature-dependent simulations.

Interaction file formats
------------------------

Three types of lattice interactions are specified via external files:

**LL file (lattice-lattice force constants)**

   Format: :math:`i, j, \Phi_{xx}, \Phi_{xy}, \Phi_{xz}, \Phi_{yy}, \Phi_{yz}, \Phi_{zz}`

   Each line specifies a pair :math:`(i,j)` and the 6 independent components of
   the symmetric force constant tensor (J/Ų). Typically read from Phonopy or
   ab initio lattice dynamics calculations.

   Example:
   ::

      1    1    100.0    0.0      0.0      100.0   0.0      100.0
      1    2    -50.0    0.5      0.0      -50.0   0.0      -50.0

**ML file (magneto-lattice coupling)**

   Format: :math:`i, j, M_{xx}, M_{xy}, M_{xz}, M_{yx}, M_{yy}, M_{yz}, M_{zx}, M_{zy}, M_{zz}`

   A 3×3 tensor :math:`M_{ij}` coupling moment of atom :math:`j` to lattice displacement at atom :math:`i`.
   Units: meV/Å. Can be fitted to magnetostriction data or computed from first principles.

   Example:
   ::

      1    1    2.5      0.0      0.0      0.0      2.5      0.0      0.0      0.0      2.5
      1    2    1.2      0.1      0.0      0.1      1.2      0.05     0.0      0.05     1.2

**MML file (magneto-magnetic-lattice coupling)**

   Format: :math:`i, j, \Phi^{\text{MML}}_{xx}, \Phi^{\text{MML}}_{xy}, \ldots`

   6×6 symmetric tensor modulating exchange (second-order magnon-phonon coupling).
   Rare in practice; typically only diagonal (isotropic exchange modulation).

Keyword reference
-----------------

Activation and basic settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

do_ld
   Enable lattice dynamics (Y/*N*, default N). Must be set to ``Y`` to activate LD/SLD.

phonfile
   File containing initial lattice configuration (atomic positions). Only needed if
   ``initlatt`` requires external data.

Lattice interactions
^^^^^^^^^^^^^^^^^^^^

ll
   Enable lattice-lattice (LL) interactions. Followed by filename containing force constant tensors.
   Mandatory for LD/SLD. File format: atom indices and 6-component force constant tensor.

ml
   Enable magneto-lattice (ML) coupling. Followed by filename containing ML tensors.
   Recommended for SLD. File format: atom indices and 3×3 ML coupling tensor.

mml
   Enable magneto-magnetic-lattice (MML) coupling. Followed by filename. Advanced feature.
   Modulates exchange via lattice displacement.

lll
   Enable three-body lattice interactions (anharmonic). Advanced; rarely used.

ll_scale
   Rescaling factor for all LL force constants (float, default 1.0). Useful for
   parametric studies or matching to experiments.

ll_phonopy
   Read LL force constants from Phonopy output (phonopy.yaml). Followed by filename.
   Automatically parses Hessian matrix from Phonopy calculations.

ll_phonopycoordfile
   Coordinate file for Phonopy read. Format: atomic positions and types for Phonopy parsing.

Damping and temperature control
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

lattdamp
   Lattice damping (friction) coefficient :math:`\gamma_{\text{lat}}` (fs⁻¹, default 0.0).
   Controls thermal relaxation of lattice. Typical range: 0.001–0.1 fs⁻¹.

do_velrsc
   Enable velocity rescaling thermostat (Y/*N*). Rescales velocities to maintain constant T.

velrsc_step
   Interval for velocity rescaling (integer steps, default 100).

velrsc_taut
   Time constant for rescaling :math:`\tau` (fs, default 100). Higher → gentler rescaling.

Initialization
^^^^^^^^^^^^^^^

initlatt
   Lattice initialization mode:

   * ``1``: From external phonfile
   * ``2``: From lattice restart file (lattrestartfile)
   * ``3``: Equilibrium (zero displacement and velocity)
   * ``4``: Perturbed (random displacements for phonon generation)

initexc
   Excitation amplitude for perturbed initialization (Å, default 0.01).

lattrestartfile
   File with previous lattice configuration (displacements and velocities).
   For restarting from checkpoint.

lattroteul
   Euler angle for lattice rotation (radians). Rotates entire lattice before simulation.

lattrotang
   Rotation angles as 3-vector. Alternative to lattroteul for Cartesian specification.

Center-of-mass control
^^^^^^^^^^^^^^^^^^^^^^

do_set_avrgp0
   Enforce zero average momentum (Y/*N*). Prevents center-of-mass drift.

do_set_avrgu0
   Enforce zero average displacement (Y/*N*). Fixes center-of-mass position.

Measurements and output
^^^^^^^^^^^^^^^^^^^^^^^

do_lavrg
   Lattice averages (Y/*N*). Outputs ``lattavrg.<simid>.out`` with time-averaged
   displacements, velocities, kinetic energy.

do_proj_lavrg
   Type-projected (sublattice) lattice averages (Y/*N*).

do_projch_lavrg
   Chemical-species projected averages (Y/*N*).

lavrg_step
   Interval for sampling lattice averages (integer steps, default 100).

lavrg_buff
   Buffer size before writing (integer, default 10).

do_ltottraj
   Full lattice trajectories (Y/*N*). Outputs ``lattices.<simid>.out`` with atomic
   displacements :math:`\mathbf{u}_i(t)` at every step. Large file!

ltottraj_step
   Sampling interval for trajectories (integer steps, default 1000).

Harmonic corrections (advanced)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

do_n3
   Include third-order (cubic anharmonic) corrections (Y/*N*). Advanced;
   requires cubic force constant tensors.

Numerical options (advanced)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

mml_diag
   MML coupling model: 0 (off), 1 (diagonal), 2 (full tensor).

mml_ene_opt
   Energy optimization for MML (Y/*N*). Iteratively refines MML parameters.

do_ethermfield
   Print thermal force contributions to file (Y/*N*). Debug feature.

do_prn_eeff
   Print total effective lattice field to file (Y/*N*). For verification.

Output files
------------

lattavrg.<simid>.out
   Lattice averages (if ``do_lavrg Y``). Columns:

   ::

      step   <u_x>    <u_y>    <u_z>    <v_x>    <v_y>    <v_z>    E_kin    E_pot

   Useful for tracking thermal expansion, acoustic phonons, lattice heating.

lattices.<simid>.out
   Full lattice trajectories (if ``do_ltottraj Y``). Format per step:

   ::

      # mstep   Natom
      iatom    u_x       u_y       u_z       v_x       v_y       v_z

   Large file; use for detailed analysis or visualization.

Example setup: Magnetic-thermal expansion
------------------------------------------

.. code-block:: text

   ! Spin-lattice dynamics: magnetic system with magnetostriction
   do_ld               Y

   ! Lattice interactions from ab initio (Phonopy)
   ll_phonopy          phonopy.yaml
   ll_phonopycoordfile phonopy_coords.dat

   ! Magneto-lattice coupling (magnetostriction)
   ml                  ml_couplings.dat

   ! Damping and temperature
   lattdamp            0.01          ! fs⁻¹ (moderate friction)
   do_velrsc           Y
   velrsc_step         100
   velrsc_taut         100.0         ! fs

   ! Initialization: equilibrium lattice, then heat
   initlatt            3
   do_set_avrgu0       Y
   do_set_avrgp0       Y

   ! Measurements
   do_lavrg            Y
   lavrg_step          50
   do_ltottraj         N             ! Turn on only for short test runs

Example setup: Phonon-magnon coupling (SLD + 3TM)
-------------------------------------------------

.. code-block:: text

   ! Combined ultrafast demagnetization + lattice heating
   do_ld               Y
   do_3tm              Y             ! Couple with 3TM

   ! LL and ML interactions
   ll                  ll_ni.dat
   ml                  ml_ni.dat

   ! Temperature and damping
   lattdamp            0.05
   temp_spin_init      300.0
   temp_latt_init      300.0
   temp_elec_init      300.0

   ! 3TM pulse coupling lattice heating
   P_pulse             100.0
   t0_pulse            1.0e-12
   sigma_pulse         2.12e-14

   ! Output: lattice heating timescale
   do_lavrg            Y
   lavrg_step          10

Physical examples and applications
----------------------------------

**Thermal expansion:**

Magnetic moments create local magnetostriction via ML coupling. When spins
heat (disorder increases) as T rises, the coupling :math:`F^{\text{ML}} = -M \mathbf{m}`
weakens, allowing lattice to expand. Effective negative thermal expansion coefficient
near Curie point.

**Phonon softening near magnetic transitions:**

Exchange modulation through MML changes effective LL force constants dynamically.
Near the critical temperature, certain phonon branches soften, indicating proximity
to spin-disorder transition.

**Ultrafast demagnetization + lattice dynamics:**

When combined with 3TM, the hot electron bath (:math:`T_e \gg T_l`) drives fast
spin disorder. ML coupling transfers this to the lattice, causing rapid thermal
expansion and phonon heating on 100 fs timescale—observable in time-resolved
X-ray or electron diffraction.

**Magnetoelastic domains:**

SLD can simulate equilibrium magnetic and elastic domain coexistence. Competing
LL elastic energy and ML magnetoelastic energy create complex multidomain patterns
at low T.

Computational performance
--------------------------

**Cost scaling:**

- SLD adds ~2–3× overhead vs. pure SD (force calculations :math:`O(N\times N_{\text{neigh}})`)
- Lattice trajectory output (``do_ltottraj Y``) adds significant I/O; use sparingly
- Thermal noise generation (``lattdamp > 0``) adds ~10–20% CPU overhead

**Memory:**

- Additional arrays: displacements, velocities, accelerations :math:`\sim 3 \times 3 \times N_{\text{atom}}`
- Force constant tensors from LL file: :math:`N_{\text{pairs}} \times 6` floats

**Recommended timesteps:**

- LD only (``do_ld Y, do_sd N``): :math:`\Delta t = 1`–5 fs (lattice ~0.1 THz)
- SLD (``do_ld Y, do_sd Y``): :math:`\Delta t = 0.1`–1 fs (magnetic ~10 THz ≫ lattice)
- Ultra-SLD with 3TM: :math:`\Delta t \lesssim 1\ \text{fs}` (capture electron heating ~100 fs)

Related keywords and cross-references
-------------------------------------

- ``do_3tm``: Couple lattice heating to ultrafast electron dynamics (see :doc:`../stimuli/temperature-3tm`)
- ``lattdamp``: Lattice friction; compare with ``damping`` (magnetic Gilbert damping)
- ``temp``: Global temperature for stochastic forces
- ``Mensemble``: Number of replicas; important for averaging out thermal noise in lattice

References
==========

See the centralized :doc:`../references` for full bibliographic entries:

- [Hellsvik2019]_ - General method for atomistic spin-lattice dynamics with first-principles accuracy
- [Eriksson2017]_ - Atomistic spin dynamics foundations and applications
- [Ma2012]_ - Spin-lattice-electron dynamics simulations of magnetic materials
