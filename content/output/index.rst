Output files
============

Depending on the settings chosen in the input file, UppASD prints out a varying number of output files. These all share the suffix *.simid.out* where *simid* is the simulation handle defined in the input file.

Simulation and Hamiltonian output
---------------------------------

**coord.simid.out** is written if ``do_prnstruct`` is switched on. Prints out the coordinates of each moment in the system in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{iatom}`
     - :math:`x`
     - :math:`y`
     - :math:`z`

**struct.simid.out** is written if ``do_prnstruct`` is switched on. Prints out the neighbour coupling list for the exchange couplings of the system, including Cartesian coordinates of exchange couplings grouped into coordination shells. Handy for checking if the system is set up correctly. *Warning:* this file might be very large for a realistic system, be mindful of that.

For scalar exchange interactions (default), the format is:

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{iatom}`
     - :math:`\text{jatom}`
     - :math:`\text{itype}`
     - :math:`\text{jtype}`
     - :math:`r_{ij}^x`
     - :math:`r_{ij}^y`
     - :math:`r_{ij}^z`
     - :math:`J_{ij}`
     - :math:`|r_{ij}|`

For tensor exchange interactions (when exchange coupling matrices are used), the format includes all nine tensor components:

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{iatom}`
     - :math:`\text{jatom}`
     - :math:`\text{itype}`
     - :math:`\text{jtype}`
     - :math:`r_{ij}^x`
     - :math:`r_{ij}^y`
     - :math:`r_{ij}^z`
     - :math:`J_{ij}^{xx}`
     - :math:`J_{ij}^{xy}`
     - :math:`J_{ij}^{xz}`
     - :math:`J_{ij}^{yx}`
     - :math:`J_{ij}^{yy}`
     - :math:`J_{ij}^{yz}`
     - :math:`J_{ij}^{zx}`
     - :math:`J_{ij}^{zy}`
     - :math:`J_{ij}^{zz}`
     - :math:`|r_{ij}|`

**dmdata.simid.out** is written if the DM interaction is defined. Prints out the DM coupling for each atom in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{iatom}`
     - :math:`\text{jatom}`
     - :math:`\text{itype}`
     - :math:`\text{jtype}`
     - :math:`r_{ij}^x`
     - :math:`r_{ij}^y`
     - :math:`r_{ij}^z`
     - :math:`D_{ij}^x`
     - :math:`D_{ij}^y`
     - :math:`D_{ij}^z`
     - :math:`|r_{ij}|`

**aniso1.simid.out** is written if the anisotropy is defined. Prints the anisotropy parameters for each atom in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`e_x`
     - :math:`e_y`
     - :math:`e_z`
     - :math:`k_1`
     - :math:`k_2`

where the first three entries are the direction of the anisotropy axis.

**biqdmdata.simid.out** is written if the effective quadratic DM interaction is defined. Prints out the effective quadratic DM coupling for each atom.

**bqdata.simid.out** is written if the bq interaction is defined. Prints out the bq coupling for each atom.

**dmstruct.simid.out** is written if the DM interaction is defined and ``do_prnstruct`` is switched on. Prints out the coupling list for the DM vector of the system. Similar to the data presented in **struct.simid.out**.

**inp.simid.out** extensive output of the values assigned to global variables after reading ``inpsd.dat`` and accompanying files.

**pddata.simid.out** is written if anisotropic exchange interaction pd interaction is defined. Prints out the effective pd couplings for each atom.


Measured observables
--------------------

**averages.simid.out** is written if measurement phase is run in SD mode. Prints out the average magnetization as a function of simulation time, in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{step}`
     - :math:`m_x`
     - :math:`m_y`
     - :math:`m_z`
     - :math:`m`
     - :math:`\sigma(m)`

where :math:`step` is the simulation time expressed in terms of the number of time steps, :math:`m_x`, :math:`m_y` and :math:`m_z` are the components of the intensive average magnetization *i.e.*, :math:`m_x=\frac{1}{N}\sum_i m_{x,i}`), :math:`m=\sqrt{m_x^2+m_y^2+m_z^2}`, and so on. :math:`\sigma(m)` is the standard deviation of :math:`m` when the number of ensembles is larger than one.

.. \vindex{cumulants.simid.out} \index{Binder cumulant} \index{Susceptibility} \index{Specific heat}

**cumulants.simid.out** prints out the running time averages of the intensive magnetization and its higher order moments, in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{step}`
     - :math:`\langle M \rangle`
     - :math:`\langle M \rangle^2`
     - :math:`\langle M \rangle^4`
     - :math:`U_4`
     - :math:`\chi`
     - :math:`C_V`

where, brackets denote time averaged quantities and :math:`U_4=1-\frac{1 \langle M \rangle^4}{3 \langle M \rangle^2}` is the fourth order Binder cumulant, useful for estimating transition temperatures [Binder2009]_, :math:`\chi` is the magnetic susceptibility, and :math:`C_V` is the heat capacity.

**afmcumulants.simid.out** is written if the ``do_cumu_proj`` flag is switched on. Prints out the running time averages of the antiferromagnetic order parameter and its higher order moments, in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{step}`
     - :math:`\langle L \rangle`
     - :math:`\langle L^2 \rangle`
     - :math:`\langle L^4 \rangle`
     - :math:`U_4^L`
     - :math:`\chi`
     - :math:`C_V`

where :math:`U_4^L` is the Binder cumulant for the AFM order parameter :math:`L`, :math:`\chi` is the susceptibility, and :math:`C_V` is the heat capacity.

**projcumulants.simid.out** is written if the ``do_cumu_proj`` flag is switched on. Prints out the running time averages of the intensive magnetization and its higher order moments, projected to each atom type, in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{step}`
     - :math:`\text{type}`
     - :math:`\langle M \rangle`
     - :math:`\langle M^2 \rangle`
     - :math:`\langle M^4 \rangle`
     - :math:`U_4`
     - :math:`\chi`

where the `type` column indicates the atom type index.

.. **mcinitial.simid.out** is written if initial phase is set to MC mode. Prints out the final iterations of the MC initial phase.
.. , in the format
.. %\begin{equation}\nonumber
.. %  mcstep, m, U_4, \chi
.. %\end{equation}
.. %\noindent where $\chi$ is the magnetic susceptibility. This is useful for checking whether or not the initial run has thermalized before entering the measurement stage.

.. %\subsubsection*{mcmeasure.simid.out}
.. %Is written if measurement phase is set to MC mode. Prints out the quantities measured in MC mode, using the same format used for \rfilename{mcinitial.simid.out}.

**moments.simid.out** is written if the ``do_tottraj`` flag is switched on. Prints the configuration of all magnetic moments at regular interval in time in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{step}`
     - :math:`\text{iatom}`
     - :math:`m_{x,i}`
     - :math:`m_{y,i}`
     - :math:`m_{z,i}`
     - :math:`m_i`

note that this file is very large. It is useful for creating animations of the evolution in time of the magnetic configuration of the system. Printed for ensemble nr 1.

**polarization.simid.out** prints out the average ferroelectric polarization as a function of simulation time, in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{step}`
     - :math:`p_x`
     - :math:`p_y`
     - :math:`p_z`
     - :math:`p`
     - :math:`\sigma(p)`

where :math:`step` is the simulation time expressed in terms of the number of time steps, :math:`p_x`, :math:`p_y` and :math:`p_z` are the components of the intensive average polarization (*i.e*, :math:`p_x=\frac{1}{N}\sum_i p_{x,i}`,\ldots, ) and :math:`p=\sqrt{p_x^2+p_y^2+p_z^2}`. :math:`\sigma(p)` is the standard deviation of :math:`p` when the number of ensembles is larger than one.

**projavgs.simid.out** is written if the ``do_proj_avrg`` flag is switched on. Prints out the same thermodynamic averages printed in ``averages.simid.out``, but projected to each atom type sublattice, in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{step}`
     - :math:`\text{type}`
     - :math:`\langle M \rangle`
     - :math:`\sigma(M)`
     - :math:`\langle M \rangle_x`
     - :math:`\langle M \rangle_y`
     - :math:`\langle M \rangle_z`

where the `type` column indicates the atom type index.

**projchavgs.simid.out** is written if the ``do_proj_avrg`` flag is switched on. Prints out the same thermodynamic averages printed in ``averages.simid.out``, but projected to each chemical type sublattice, in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{step}`
     - :math:`\text{type}`
     - :math:`\langle M \rangle`
     - :math:`\sigma(M)`
     - :math:`\langle M \rangle_x`
     - :math:`\langle M \rangle_y`
     - :math:`\langle M \rangle_z`
     - :math:`\sum\langle M \rangle`

where the `type` column indicates the chemical type index.

**restart.simid.out** the magnetic configuration of the system at a specific point in time. Can be used as input when the ``initmag`` flag is set to 4. The format is:

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{ensemble}`
     - :math:`\text{iatom}`
     - :math:`|m_i|`
     - :math:`m_{x,i}`
     - :math:`m_{y,i}`
     - :math:`m_{z,i}`

The first line contains the time step, followed by all magnetic moments for each ensemble.

**trajectory.simid.out** the trajectory as a function of time step for an individual magnetic moment in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{step}`
     - :math:`\text{iatom}`
     - :math:`m_{x,i}`
     - :math:`m_{y,i}`
     - :math:`m_{z,i}`
     - :math:`m_m`
     - :math:`m_i`

If the number ``ntraj`` is also defined to be greater than 1, the code prints out ntraj files named ``trajectory.simid.XXX.Y.out`` where XXX is the trajectory number (padded to 3 digits) and Y is the ensemble number.

**sq.simid.out** is written if the ``do_sc`` flag is Y or C. Prints out the static correlation function in reciprocal space :math:`S(q)` in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{nq}`
     - :math:`q_x`
     - :math:`q_y`
     - :math:`q_z`
     - :math:`S^x(\mathbf{q})`
     - :math:`S^y(\mathbf{q})`
     - :math:`S^z(\mathbf{q})`
     - :math:`S(\mathbf{q})`

**sqt0.simid.out** is written if the ``do_sc`` flag is Y or C and the ``do_qt_traj`` flag is Y. Prints out the trajectory in time of the equal time correlation function :math:`S(q)` in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{step}`
     - :math:`\text{nq}`
     - :math:`q_x`
     - :math:`q_y`
     - :math:`q_z`
     - :math:`S^x(\mathbf{q})`
     - :math:`S^y(\mathbf{q})`
     - :math:`S^z(\mathbf{q})`
     - :math:`S(\mathbf{q})`

**sra.simid.out** is written if the ``do_sc`` flag is Y or C. Prints out the static correlation function in real space :math:`S(r)` in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`|r|`
     - :math:`S^x(r)`
     - :math:`S^y(r)`
     - :math:`S^z(r)`
     - :math:`S(r)`

**sqt.simid.out** is written if the ``do_sc`` flag is switched on. Prints out the time-resolved structure factor :math:`S(q,t)` in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{nstep}`
     - :math:`\text{nq}`
     - :math:`\operatorname{Re}[S^x(\mathbf{q},t)]`
     - :math:`\operatorname{Im}[S^x(\mathbf{q},t)]`
     - :math:`\operatorname{Re}[S^y(\mathbf{q},t)]`
     - :math:`\operatorname{Im}[S^y(\mathbf{q},t)]`
     - :math:`\operatorname{Re}[S^z(\mathbf{q},t)]`
     - :math:`\operatorname{Im}[S^z(\mathbf{q},t)]`

.. %This file can be very large.

**projsqt.simid.out** is written if the ``do_sc_proj`` flag is switched on. Prints out the same information printed in ``sqt.simid.out``, but projected to each atom type present in the system.

.. %This file can be very large.

**sqw.simid.out** is written if the ``do_sc`` flag is switched on. Prints out the frequency-resolved dynamic structure factor :math:`S(q,\omega)` in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{nq}`
     - :math:`q_x`
     - :math:`q_y`
     - :math:`q_z`
     - :math:`\text{nstep}`
     - :math:`S^x(\mathbf{q},\omega)`
     - :math:`S^y(\mathbf{q},\omega)`
     - :math:`S^z(\mathbf{q},\omega)`
     - :math:`S(\mathbf{q},\omega)`

.. %This file can be very large.

**projsqw.simid.out** is written if the ``do_sc_proj`` flag is switched on. Prints out the same information printed in ``sqw.simid.out``, but projected to each atom type present in the system.

.. %This file can be very large.

**swdos.simid.out** is written if the ``do_sc`` flag is switched on. Prints out the :math:`S(q,\omega)` 'density of states' as a function of energy, in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`E` (meV)
     - :math:`D(S^x(q,E))`
     - :math:`D(S^y(q,E))`
     - :math:`D(S^z(q,E))`
     - :math:`D(|S(q,E)|)`

where :math:`D` denotes the density of states obtained by integrating over all q-points. For scalar correlation functions, only two columns are printed: :math:`E` and :math:`D(S(q,E))`.

**sknumber.simid.out** is written if the ``skyno`` flag is Y. Prints out the skyrmion number as a function of time, in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{step}`
     - :math:`N_{\text{sk}}`
     - :math:`\langle N_{\text{sk}} \rangle`
     - :math:`\sigma(N_{\text{sk}})`

where :math:`N_{\text{sk}}` is the instantaneous skyrmion number, :math:`\langle N_{\text{sk}} \rangle` is the cumulative average, and :math:`\sigma(N_{\text{sk}})` is the standard deviation.

**cmass_skynum.simid.out** is written if the ``do_skyno_cmass`` flag is Y. Prints out the center of mass position of the skyrmion charge density, in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{step}`
     - :math:`r_{\text{sk},x}`
     - :math:`r_{\text{sk},y}`
     - :math:`r_{\text{sk},z}`

where :math:`\mathbf{r}_{\text{sk}}` is the center of mass coordinate of the topological charge.

**proj_sknumber.XX.simid.out** is written if the ``do_proj_skyno`` flag is Y. Prints out the skyrmion number projected onto each atom type, where XX is the atom type index. Format is the same as ``sknumber.simid.out``.

**dens_skynum.simid.out** is written if the ``do_skyno_den`` flag is Y. Prints out the site-dependent skyrmion number density, in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{step}`
     - :math:`\text{site}`
     - :math:`n_{\text{sk},i}`

where :math:`n_{\text{sk},i}` is the local skyrmion number density at site :math:`i`.

**totenergy.simid.out** is written if the ``plotenergy`` flag is switched on. Prints out the total energy of the system as a function of time step, in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{step}`
     - :math:`E_{\text{tot}}`
     - :math:`E_{\text{exc}}`
     - :math:`E_{\text{ani}}`
     - :math:`E_{\text{DM}}`
     - :math:`E_{\text{PD}}`
     - :math:`E_{\text{BiqDM}}`
     - :math:`E_{\text{BQ}}`
     - :math:`E_{\text{dip}}`
     - :math:`E_{\text{Zeeman}}`
     - :math:`E_{\text{LSF}}`
     - :math:`E_{\text{chir}}`
     - :math:`E_{\text{ring}}`
     - :math:`E_{\text{SA}}`

where the columns represent total energy, exchange, anisotropy, DM interaction, pseudo-dipolar, biquadratic DM, biquadratic, dipolar, Zeeman (external field), lattice spin fluctuation, chiral, ring exchange, and shape anisotropy contributions respectively. The unit for all entries is mRy per atom.

**stdenergy.simid.out** is written if the ``plotenergy`` flag is switched on and ``Mensemble > 1``. Prints out the standard deviation of the energy contributions across ensembles, using the same format as ``totenergy.simid.out``.

.. only: false

  **magnon_curr.simid.out** is written if the ``do_currents`` flag is switched on. Prints out the magnon current density at each site, in the format
  
  .. list-table::
     :widths: auto
     :header-rows: 0
     :class: borderless centered
  
     * - :math:`\text{step}`
       - :math:`\text{site}`
       - :math:`j_{m,x}`
       - :math:`j_{m,y}`
       - :math:`j_{m,z}`
       - :math:`|\mathbf{j}_m|^2`
  
  where :math:`\mathbf{j}_m` is the magnon current density vector.
  
  **heat_curr.simid.out** is written if the ``do_currents`` flag is switched on. Prints out the heat current density at each site, using the same format as ``magnon_curr.simid.out``.
  
  **heat_curr2.simid.out** is written if the ``do_currents`` flag is switched on. Prints out an alternative formulation of the heat current density at each site, using the same format as ``magnon_curr.simid.out``.
  
  **psi_data.simid.out** is written if the ``do_currents`` flag is switched on. Prints out the complex order parameter :math:`\psi` at each site for each ensemble, in the format
  
  .. list-table::
     :widths: auto
     :header-rows: 0
     :class: borderless centered
  
     * - :math:`\text{step}`
       - :math:`\text{site}`
       - :math:`\text{ensemble}`
       - :math:`|\psi|`
       - :math:`\text{arg}(\psi)`
       - :math:`\text{Re}(\psi)`
       - :math:`\text{Im}(\psi)`


Lattice dynamics observables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**lattaverages.simid.out** is written if the lattice dynamics module is active. Prints out the average ionic displacements and velocities as a function of simulation time, in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{step}`
     - :math:`u_x`
     - :math:`u_y`
     - :math:`u_z`
     - :math:`u`
     - :math:`\sigma(u)`
     - :math:`v_x`
     - :math:`v_y`
     - :math:`v_z`
     - :math:`v`
     - :math:`\sigma(v)`

where :math:`\mathbf{u}` is the average displacement vector and :math:`\mathbf{v}` is the average velocity vector.

**lattmomenta.simid.out** is written if the lattice dynamics module is active. Prints out the average ionic momenta and angular momenta, in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{step}`
     - :math:`p_x`
     - :math:`p_y`
     - :math:`p_z`
     - :math:`p`
     - :math:`\sigma(p)`
     - :math:`L_x`
     - :math:`L_y`
     - :math:`L_z`
     - :math:`L`
     - :math:`\sigma(L)`

where :math:`\mathbf{p}` is the momentum vector and :math:`\mathbf{L}` is the angular momentum vector.

**lattenergy.simid.out** is written if the lattice dynamics module is active. Prints out the lattice energy components, in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{step}`
     - :math:`E_{\text{ld-pot}}`
     - :math:`E_{\text{sd-pot}}`
     - :math:`E_{\text{sld-pot}}`
     - :math:`E_{\text{tot-pot}}`
     - :math:`E_{\text{kin}}`
     - :math:`E_{\text{tot}}`
     - :math:`T_{\text{ion}}`

where the columns represent lattice-derived potential, spin-derived potential, spin-lattice coupling potential, total potential, kinetic energy, total energy, and ionic temperature respectively.

**lattcumenergy.simid.out** is written if the lattice dynamics module is active. Prints out cumulative lattice energy statistics including heat capacity, in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{step}`
     - :math:`\langle E_{\text{pot}} \rangle`
     - :math:`\langle E_{\text{kin}} \rangle`
     - :math:`\langle E_{\text{tot}} \rangle`
     - :math:`\sigma(E_{\text{tot}})`
     - :math:`C_V`
     - :math:`T`

where brackets denote time-averaged quantities, :math:`C_V` is the heat capacity, and :math:`T` is the temperature.

**lattprojavgs.simid.out** is written if the lattice dynamics module is active and ``do_proj_avrg`` is switched on. Prints out the same lattice dynamics averages printed in ``lattaverages.simid.out``, but projected to each atom type.

**lattprojchavrgs.simid.out** is written if the lattice dynamics module is active and ``do_proj_avrg`` is switched on. Prints out the same lattice dynamics averages printed in ``lattaverages.simid.out``, but projected to each chemical type.

**disp.simid.out** is written if the ``do_tottraj`` flag is switched on for lattice dynamics. Prints out the displacement and velocity trajectories for all atoms, in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{step}`
     - :math:`\text{iatom}`
     - :math:`u_x`
     - :math:`u_y`
     - :math:`u_z`
     - :math:`v_x`
     - :math:`v_y`
     - :math:`v_z`
     - :math:`L_x`
     - :math:`L_y`
     - :math:`L_z`

**gdisp.simid.out** is written if the ``do_tottraj`` flag is switched on for lattice dynamics. Prints out the global (lab frame) coordinates and velocities, in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{step}`
     - :math:`\text{iatom}`
     - :math:`x`
     - :math:`y`
     - :math:`z`
     - :math:`v_x`
     - :math:`v_y`
     - :math:`v_z`
     - :math:`L_x`
     - :math:`L_y`
     - :math:`L_z`

where the coordinates include both initial positions and displacements (global coordinates), :math:`\mathbf{v}` is the velocity vector, and :math:`\mathbf{L}` is the angular momentum vector.

**rcoord.simid.out** is written if the ``do_tottraj`` flag is switched on for lattice dynamics. Prints out the real-space coordinates including displacements, in the format

.. list-table::
   :widths: auto
   :header-rows: 0
   :class: borderless centered

   * - :math:`\text{iatom}`
     - :math:`x`
     - :math:`y`
     - :math:`z`

**disptraj.simid.T.E.out** is written if lattice dynamics trajectories are sampled (similar to ``trajectory.simid.out`` for magnetic moments). Prints selected atomic displacement trajectories where T is the trajectory number and E is the ensemble number.


**References**

See the centralized :doc:`../references` for full bibliographic entries:

- [Binder2009]_ - Guide to Monte Carlo simulation in statistical physics

