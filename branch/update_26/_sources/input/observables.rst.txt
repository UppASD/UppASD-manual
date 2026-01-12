Observables and measurements
============================

This page is a comprehensive reference for all observable and measurement-related input keywords in UppASD. Keywords are organized by category. For detailed physics background and formulas, see the dedicated chapters: :doc:`../observables/averages`, :doc:`../observables/correlations`, :doc:`../observables/autocorrelation`, and :doc:`../observables/polarization`.

Overview
--------

Most observables are controlled by a pair of keywords:

- ``do_<observable>`` – Enables measurement (Y/N, or specialized mode)
- ``<observable>_step`` – Sampling interval in simulation steps
- ``<observable>_buff`` – Optional buffer size before file write

This design minimizes I/O overhead while maintaining fine temporal resolution.

Energy and thermodynamics
--------------------------

``plotenergy``
   Enable energy breakdown (Y=yes, *N=no*). When enabled, writes total energy and
   contributions from each term in the Hamiltonian (exchange, anisotropy, DMI, etc.)
   to ``energy.simid.out``. Default: **N**.

Average magnetization
---------------------

``do_avrg``
   Enable sampling of system-wide average magnetization and higher-order moments
   (Y=yes, *N=no*). Output: ``averages.simid.out``. See :doc:`../observables/averages` for details.
   Default: **N**.

``avrg_step``
   Sampling interval in simulation steps (integer). Default: **100**.

``avrg_buff``
   Buffer size before writing (integer). Default: **10**.

``do_proj_avrg``
   Type-projected (sublattice) or site-projected average magnetization
   (Y=by type, A=by site, *N=no*). Output: ``projavgs.simid.out``.
   Automatically enables ``do_avrg``. Default: **N**.

``do_projch_avrg``
   Chemical-species-projected averages for random alloys (Y=yes, *N=no*).
   Requires ``do_ralloy 1``. Output: ``projchavgs.simid.out``. Default: **N**.

Cumulants and thermodynamic fluctuations
----------------------------------------

``do_cumu``
   Enable cumulant sampling for phase transition analysis and thermodynamic moments
   (Y=standard, A=antiferromagnetic, *N=no*). Calculates Binder cumulant, magnetic
   susceptibility, specific heat. **Automatically enabled for Monte Carlo simulations.**
   Output: ``cumulants.simid.out``, ``cumulants.simid.json``. Default: **N**.

``cumu_step``
   Sampling interval for cumulants (integer steps). Default: **50**.

``cumu_buff``
   Buffer size (integer). Default: **10**.

``do_cumu_proj``
   Type-projected cumulants (Y=yes, *N=no*). Computes phase transition parameters
   per atom type. Output: ``cumulants_proj.simid.out``. Default: **N**.

Trajectories: full moment configurations
-----------------------------------------

``do_tottraj``
   Sample and print all magnetic moment configurations (Y=yes, *N=no*).
   Generates large ``moments.simid.out`` file. Use with caution for big systems.
   Default: **N**.

``tottraj_step``
   Sampling interval (integer steps). Default: **1000**.

``tottraj_buff``
   Buffer size (integer). Default: **10**.

``ntraj``
   Number of individual atom trajectories to sample separately (integer ≥ 0).
   Followed by ``ntraj`` input lines, each with format:
   ``atom_index traj_step traj_buff``.
   Default: **0** (no individual trajectories).

Spin correlations and structure factors
---------------------------------------

``do_sc``
   Spin correlation sampling mode (Y=full, Q=frequency only, T=time-dependent only,
   C=static only, *N=disabled*). Computes dynamical structure factor :math:`S(\mathbf{q},\omega)`
   and/or static :math:`S(\mathbf{q})`. See :doc:`../observables/correlations`.
   Default: **N**.

``do_sr``
   Real-space correlation :math:`G(\mathbf{r})` (Y=yes, *N=no*). Default: **N**.

``sc_step``
   Temporal sampling interval (integer steps). Controls max frequency.
   Default: **10**.

``sc_nstep``
   Number of temporal samples per correlation measurement (integer).
   Default: **100**.

``sc_sep``
   Steps between independent measurements for ``do_sc = C`` (integer).
   Default: **1**.

``qpoints``
   Q-point mesh specification mode: F=Cartesian from file, G=from qpoints.dat,
   C=full reciprocal cell, D=fractional coords, B=with weights.
   Default: **F**.

``qfile``
   Filename containing q-point vectors. Default: **qpoints**.

``do_sc_proj``
   Type-projected spin correlations (Y=yes, C=static only, *N=no*).
   Default: **N**.

``do_sc_projch``
   Chemical-species-projected correlations (Y=yes, C=static only, *N=no*).
   Default: **N**.

``do_qt_traj``
   Write time-dependent equal-time correlation :math:`S(\mathbf{q},t)` (Y=yes, *N=no*).
   Default: **N**.

``do_sc_local_axis``
   Use local quantization axis (Y=yes, *N=no*). Computes parallel/perpendicular
   components. Default: **N**.

``sc_local_axis_mix``
   Mixing parameter for local axis update (0.0–1.0). Default: **0.05**.

``sc_window_fun``
   Windowing function (1=rectangular, 2=Hann, 3=Hamming, 4=Blackman-Harris).
   Default: **1**.

``do_sc_dosonly``
   Magnon density-of-states only (Y=yes, *N=no*). Skips full :math:`S(\mathbf{q},\omega)`.
   Default: **N**.

``do_connected``
   Include only connected part (Y=yes, *N=no*). Default: **Y**.

``do_uc``
   Lattice displacement correlations (Y=yes, *N=no*). For SLD systems. Default: **N**.

Autocorrelation functions
--------------------------

``do_autocorr``
   Enable spin autocorrelation sampling (Y=yes, *N=no*). Measures relaxation dynamics.
   Output: ``autocorr.simid.out``. Requires external file ``acfile`` listing waiting times.
   See :doc:`../observables/autocorrelation`. Default: **N**.

``acfile``
   File containing waiting times (one per line, in simulation steps). Mandatory if ``do_autocorr Y``.

``ac_step``
   Sampling interval (integer steps). Default: **100**.

``ac_buff``
   Buffer size (integer). Default: **10**.

``do_macro_cells``
   Spatial binning for autocorrelation (Y=yes, *N=no*). Outputs per-macrocell files.
   Default: **N**.

Polarization and chirality
---------------------------

``do_pol``
   Enable ferroelectric polarization measurement (Y=yes, *N=no*).
   Computes :math:`\mathbf{P} = \gamma \sum_{i,j} \hat{\mathbf{e}}_{ij} \times (\mathbf{m}_i \times \mathbf{m}_j)`.
   Output: ``polarization.simid.out``. Requires ``do_sortcoup N``.
   See :doc:`../observables/polarization`. Default: **N**.

``max_pol_nn``
   Number of neighbors in polarization sum (integer). Default: **6**.

``pol_step``
   Sampling interval (integer steps). Default: **100**.

``pol_buff``
   Buffer size (integer). Default: **10**.

``do_chir``
   Enable scalar chirality measurement (Y=yes, *N=no*). Measures local spin
   chirality :math:`\boldsymbol{\chi}_{ij} = \mathbf{m}_i \times \mathbf{m}_j`.
   Default: **N**.

Stiffness and micromagnetic parameters
--------------------------------------

``do_stiffness``
   Calculate spin-wave stiffness and exchange constant (Y=yes, *N=no*).
   Computes effective continuum parameters from atomistic simulations.
   Default: **N**.

``eta_min``
   Convergence parameter lower bound for stiffness (float). Default: **6**.

``eta_max``
   Convergence parameter upper bound for stiffness (float). Default: **12**.

``alat``
   Lattice constant (in meters) for stiffness normalization (float).
   Required if ``do_stiffness Y``.

Lattice dynamics observables (SLD)
----------------------------------

For spin-lattice dynamics (``do_ld Y``), additional observables track lattice dynamics:

``do_lavrg``
   Lattice averages (Y=yes, *N=no*). Outputs ``lattavrg.<simid>.out`` with
   displacements, velocities, kinetic energy. Default: **N**.

``do_proj_lavrg``
   Type-projected lattice averages (Y=yes, *N=no*). Default: **N**.

``do_projch_lavrg``
   Chemical-species-projected lattice averages (Y=yes, *N=no*). Default: **N**.

``lavrg_step``
   Sampling interval (integer steps). Default: **100**.

``lavrg_buff``
   Buffer size (integer). Default: **10**.

``do_ltottraj``
   Full lattice trajectories (Y=yes, *N=no*). Outputs ``lattices.<simid>.out``.
   **Generates large files.** Default: **N**.

``ltottraj_step``
   Sampling interval (integer steps). Default: **1000**.
