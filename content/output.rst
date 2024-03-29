Output files
============

Depending on the settings chosen in the input file, UppASD prints out a varying number of output files. These all share the suffix *.simid.out* where *simid* is the simulation handle defined in the input file.

Simulation and Hamiltonian output
---------------------------------

**aniso1.simid.out** is written if the anisotropy is defined. Prints the anisotropy parameters for each atom in the format

.. math::

  aniso1.simid.out:   e_x, e_y, e_z, k_1, k_2

where the three first entries are the direction of the anisotropy axis.

**biqdmdata.simid.out** is written if the effective quadratic DM interaction is defined. Prints out the effective quadratic DM coupling for each atom.

**bqdata.simid.out** is written if the bq interaction is defined. Prints out the bq coupling for each atom.

**coord.simid.out** is written if ``do_prnstruct`` is switched on. Prints out the coordinates of each moment in the system.

**dmdata.simid.out** is written if the DM interaction is defined. Prints out the DM coupling for each atom.

**dmstruct.simid.out** is written if the DM interaction is defined and ``do_prnstruct`` is switched on. Prints out the coupling list for the DM vector of the system. Similar to the data presented in **struct.simid.out**.

**inp.simid.out** extensive output of the values assigned to global variables after reading ``inpsd.dat`` and accompanying files.

**pddata.simid.out** is written if anisotropic exchange interaction pd interaction is defined. Prints out the effective pd couplings for each atom.

**struct1.simid.out** is written if ``do_prnstruct`` is switched on. Prints out the neighbour coupling list for the exchange couplings of the system. Handy for checking if the system is set up correctly. *Warning* this file might be very large for a realistic system, be mindful of that.

**struct.simid.out** is written if ``do_prnstruct`` is switched on. Prints out the Cartesian coordinates of the exchange couplings. The entries are grouped into coordination shells. Handy for checking if the system is set up correctly. *Warning* this file might be very large for a realistic system, be mindful of that.


Measured observables
--------------------

**averages.simid.out** is written if measurement phase is run in SD mode. Prints out the average magnetization as a function of simulation time, in the format

.. math::

  averages.simid.out:   step, m_x, m_y, m_z, m, \sigma(m)

where :math:`step` is the simulation time expressed in terms of the number of time steps, :math:`m_x`, :math:`m_y` and :math:`m_z` are the components of the intensive average magnetization *i.e.*, :math:`m_x=\frac{1}{N}\sum_i m_{x,i}`), :math:`m=\sqrt{m_x^2+m_y^2+m_z^2}`, and so on. :math:`\sigma(m)` is the standard deviation of :math:`m` when the number of ensembles is larger than one.

.. \vindex{cumulants.simid.out} \index{Binder cumulant} \index{Susceptibility} \index{Specific heat}

**cumulants.simid.out}** prints out the running time averages of the intensive magnetization and its higher order moments, in the format

.. math::

  cumulants.simid.out:   step, \langle M \rangle, \langle M \rangle^2, \langle M \rangle^4, U_4, \chi, C_v

where, brackets denote time averaged quantities and :math:`U_4=1-\frac{1 \langle M \rangle^4}{3 \langle M \rangle^2}` is the fourth order Binder cumulant, useful for estimating transition temperatures [Binder2009]_, :math:`\chi` is the magnetic susceptibility, and :math:`C_V` is the heat capacity.

.. **mcinitial.simid.out** is written if initial phase is set to MC mode. Prints out the final iterations of the MC initial phase.
.. , in the format
.. %\begin{equation}\nonumber
.. %  mcstep, m, U_4, \chi
.. %\end{equation}
.. %\noindent where $\chi$ is the magnetic susceptibility. This is useful for checking whether or not the initial run has thermalized before entering the measurement stage.

.. %\subsubsection*{mcmeasure.simid.out}
.. %Is written if measurement phase is set to MC mode. Prints out the quantities measured in MC mode, using the same format used for \rfilename{mcinitial.simid.out}.

**moments.simid.out** is written if the ``do_tottraj`` flag is switched on. Prints the configuration of all magnetic moments at regular interval in time in the format

.. math::

  moments.simid.out:   step, atom, m_{x,i}, m_{y,i}, m_{z,i}, m_{i}

note that this file is very large. It is useful for creating animations of the evolution in time of the magnetic configuration of the system. Printed for ensemble nr 1.

**polarization.simid.out** prints out the average ferroelectric polarization as a function of simulation time, in the format

.. math::

  polarization.simid.out:   step, p_x, p_y, p_z, p, \sigma(p)

where :math:`step` is the simulation time expressed in terms of the number of time steps, :math:`p_x`, :math:`p_y` and :math:`p_z` are the components of the intensive average polarization (*i.e*, :math:`p_x=\frac{1}{N}\sum_i p_{x,i}`,\ldots, ) and :math:`p=\sqrt{p_x^2+p_y^2+p_z^2}`. :math:`\sigma(p)` is the standard deviation of :math:`p` when the number of ensembles is larger than one.

**projavs.simid.out** is written of the ``do_proj_avrg`` flag is switched on. Prints out the same thermodynamic averages printed in ``averages.simid.out``, but projected to each atom type sublattice. The format is also identical to ``averages.simid.out``, except for the addition of a column indicating the sublattice.

**restart.simid.out** the magnetic configuration of the system at a specific point in time. Can be used as input when the ``initmag`` flag is set to 4.

**trajectory.simid.out** the trajectory as a function of time step for an individual magnetic moment on format

.. math::

  trajectory.simid.out:   step, atom nr, m_{x,i},m_{y,i},m_{z,i},m_m,m_i

if the number ``ntraj`` is also defined to be greater than 1, the code prints out ntraj files named ``trajectory.simid.ntraj.ensemblnr.out``.

**sq.simid.out}** is written if the ``do_sc`` flag is C. Prints out the static correlation function in reciprocal space :math:`S(q)` in the format

.. math::
   
  sq.simid.out:   nq,q_x,q_y,q_z,S^x(\mathbf{q}),S^y(\mathbf{q}),S^z(\mathbf{q}),S(\mathbf{q})

**sqt0.simid.out** is written if the ``do_sc`` flag is C and the ``do_qt_traj`` flag is Y. Prints out the trajectory in time of the equal time correlation function :math:`S(q)` in the format

.. math::

  sqt0.simid.out:step,nq,q_x,q_y,q_z,S^x(\mathbf{q}),S^y(\mathbf{q}),S^z(\mathbf{q}),S(\mathbf{q})

**sra.simid.out** is written if the ``do_sc`` flag is C. Prints out the static correlation function in real space :math:`S(r)` in the format

.. math::

  sra.simid.out:   |r|,S^x(r),S^y(r),S^z(r),S(r)

**sqt.simid.out** is written of the ``do_sc`` flag is switched on. Prints out the time-resolved structure factor :math:`S(q,t)` in the format

.. math::

  sqt.simid.out:   nstep, nq, \mathrm{Re}[S^x(\mathbf{q},t)], \mathrm{Im}[S^x(\mathbf{q},t)], \mathrm{Re}[S^y(\mathbf{q},t)], \mathrm{Im}[S^y(\mathbf{q},t)], \mathrm{Re}[S^z(\mathbf{q},t)], \mathrm{Im}[S^z(\mathbf{q},t)]

.. %This file can be very large.

**projsqt.simid.out** is written of the ``do_sc_proj`` flag is switched on. Prints out the same information printed in ``sqt.simid.out``, but projected to each atom type present in the system.

.. %This file can be very large.

**sqw.simid.out** is written of the ``do_sc`` flag is switched on. Prints out the frequency-resolved dynamic structure factor :math:`S(q,\omega)` in the format

.. math::

  sqw.simid.out:   nq, q_x, q_y, q_z, nstep, S^x(\mathbf{q},\omega), S^y(\mathbf{q},\omega), S^z(\mathbf{q},\omega), S(\mathbf{q},\omega)

.. %This file can be very large.

**projsqw.simid.out** is written of the ``do_sc_proj`` flag is switched on. Prints out the same information printed in ``sqw.simid.out``, but projected to each atom type present in the system.

.. %This file can be very large.

**swdos.simid.out** is written of the ``do_sc`` flag is switched on. Prints out the :math:`S(q,\omega)` 'density of states' as a function of energy.

**totenergy.simid.out** is written if the ``plotenergy`` flag is switched on. Prints out the total energy of the system as a function of time step.
