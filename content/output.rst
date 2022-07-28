Output files
============

Depending on the settings chosen in the input file, UppASD prints out a varying number of output files. These all share the suffix *.simid.out* where *simid* is the simulation handle defined in the input file.

Simulation and Hamiltonian output
---------------------------------

**aniso1.simid.out** is written if the anisotropy is defined. Prints the anisotropy parameters for each atom in the format::

  aniso1.simid.out:  e_x,e_y,e_z,k_1,k_2

where the three first entries are the direction of the anisotropy axis.

**biqdmdata.simid.out** is written if the effective quadratic DM interaction is defined. Prints out the effective quadratic DM coupling for each atom.

**bqdata.simid.out** is written if the bq interaction is defined. Prints out the bq coupling for each atom.

**coord.simid.out** is written if \rkeyword{do_prnstruct} is switched on. Prints out the coordinates of each moment in the system.

**dmdata.simid.out** is written if the DM interaction is defined. Prints out the DM coupling for each atom.
.. % and their strengths for each atom.

**dmstruct.simid.out** is written if the DM interaction is defined and **do_prnstruct** is switched on. Prints out the coupling list for the DM vector of the system. Similar to the data presented in **struct.simid.out**.

**inp.simid.out** extensive output of the values assigned to global variables after reading ``inpsd.dat`` and accompanying files.

**pddata.simid.out** is written if anisotropic exchange interaction pd interaction is defined. Prints out the effective pd couplings for each atom.

**struct1.simid.out** is written if ``do_prnstruct`` is switched on. Prints out the neighbour coupling list for the exchange couplings of the system. Handy for checking if the system is set up correctly. *Warning* this file might be very large for a realistic system, be mindful of that.

**struct.simid.out** is written if ``do_prnstruct`` is switched on. Prints out the Cartesian coordinates of the exchange couplings. The entries are grouped into coordination shells. Handy for checking if the system is set up correctly. *Warning* this file might be very large for a realistic system, be mindful of that.


Measured observables
--------------------

\subsubsection*{averages.simid.out} \vindex{averages.simid.out} \index{Magnetization}
Is written if measurement phase is run in SD mode. Prints out the average magnetization as a function of simulation time, in the format
\begin{fileinfo}[averages.simid.out:]
\begin{lstlisting}
  step, $m_x$, $m_y$, $m_z$, $m$, $\sigma$(m)
\end{lstlisting}
\end{fileinfo}
\noindent where $step$ is the simulation time expressed in terms of the number of time steps, $m_x$, $m_y$ and $m_z$ are the components of the intensive average magnetization (\textit{i.e.}, $m_x=\frac{1}{N}\sum_i m_{x,i}$), $m=\sqrt{m_x^2+m_y^2+m_z^2}$, and so on. $\sigma(m)$ is the standard deviation of $m$ when the number of ensembles is larger than one.
%%%%%%
\subsubsection*{cumulants.simid.out} \vindex{cumulants.simid.out} \index{Binder cumulant} \index{Susceptibility} \index{Specific heat}
Prints out the running time averages of the intensive magnetization and its higher order moments, in the format  
\begin{fileinfo}[cumulants.simid.out:]
\begin{lstlisting}
  step, $<M>$, $<M>^2$, $<M>^4$, $U_4$, $\chi$, $C_v$
\end{lstlisting}
\end{fileinfo}
 \noindent where, brackets denote time averaged quantities and $U_4=1-\frac{1 <M>^4}{3 <M>^2}$ is the fourth order `Binder' cumulant, useful for estimating transition temperatures ~\cite{Binder2009}, $\chi$ is the magnetic susceptibility, and $C_V$ is the heat capacity.
%%%%%%
%\subsubsection*{mcinitial.simid.out} 
%Is written if initial phase is set to MC mode. Prints out the final iterations of the MC initial phase, in the format
%\begin{equation}\nonumber
%  mcstep, m, U_4, \chi
%\end{equation}
%\noindent where $\chi$ is the magnetic susceptibility. This is useful for checking whether or not the initial run has thermalized before entering the measurement stage.
%%%%%%
%\subsubsection*{mcmeasure.simid.out} 
%Is written if measurement phase is set to MC mode. Prints out the quantities measured in MC mode, using the same format used for \rfilename{mcinitial.simid.out}.
%%%%%%
\subsubsection*{moments.simid.out} \vindex{moments.simid.out} \index{Trajectories}
Is written if the \rkeyword{do_tottraj} flag is switched on. Prints the configuration of all magnetic moments at regular interval in time in the format 
\begin{fileinfo}[moments.simid.out:]
\begin{lstlisting}
  step, atom, $m_{x,i}$, $m_{y,i}$, $m_{z,i}$, $m_{i}$
\end{lstlisting}
\end{fileinfo}
\noindent  Note that this file is very large. It is useful for creating animations of the evolution in time of the magnetic configuration of the system. Printed for ensemble nr 1.
%%%%%%
\subsubsection*{polarization.simid.out} \vindex{polarization.simid.out} \index{Polarization}
Prints out the average ferroelectric polarization as a function of simulation time, in the format
\begin{fileinfo}[polarization.simid.out:]
\begin{lstlisting}
  step, $p_x$, $p_y$, $p_z$, $p$, $\sigma(p)$
\end{lstlisting}
\end{fileinfo}
\noindent where $step$ is the simulation time expressed in terms of the number of time steps, $p_x$, $p_y$ and $p_z$ are the components of the intensive average polarization (\textit{i.e.}, $p_x=\frac{1}{N}\sum_i p_{x,i}$,\ldots, ) and $p=\sqrt{p_x^2+p_y^2+p_z^2}$. $\sigma(p)$ is the standard deviation of $p$ when the number of ensembles is larger than one.
%%%%%%
\subsubsection*{projavs.simid.out} \vindex{projavs.simid.out}
Is written of the \rkeyword{do_proj_avrg} flag is switched on. Prints out the same thermodynamic averages printed in \rfilename{averages.simid.out}, but projected to each atom type sublattice. The format is also identical to \rfilename{averages.simid.out}, except for the addition of a column indicating the sublattice.
%%%%%%
\subsubsection*{restart.simid.out}\vindex{restart.simid.out}
The magnetic configuration of the system at a specific point in time. Can be used as input when the \rkeyword{initmag} flag is set to 4.
%%%%%%
\subsubsection*{trajectory.simid.out} \vindex{trajectory.simid.out} \index{Trajectories}
The trajectory as a function of time step for an individual magnetic moment on format
\begin{fileinfo}[trajectory.simid.out:]
\begin{lstlisting}
  step, atom nr, $m_{x,i}$,$m_{y,i}$,$m_{z,i}$,$m_m$,$m_i$
\end{lstlisting}
\end{fileinfo}
\noindent If the number \rkeyword{ntraj} is also defined to be greater than 1, the code prints out ntraj files named \texttt{trajectory.simid.ntraj.ensemblnr.out}.
%%%%%%
\subsubsection*{sq.simid.out}\vindex{sq.simid.out} \index{Correlation function!Spatial}
Is written if the \rkeyword{do_sc} flag is C. Prints out the static correlation function in reciprocal space $S(q)$ in the format
\begin{fileinfo}[sq.simid.out:]
\begin{lstlisting}
nq,$q_x$,$q_y$,$q_z$,$S^x(\mathbf{q})$,$S^y(\mathbf{q})$,$S^z(\mathbf{q})$,$S(\mathbf{q})$
\end{lstlisting}
\end{fileinfo}
%%%%%%
\subsubsection*{sqt0.simid.out} \vindex{sqt0.simid.out} \index{Correlation function!Equal time}
Is written if the \rkeyword{do_sc} flag is C and the \rkeyword{do_qt_traj} flag is Y. Prints out the trajectory in time of the equal time correlation function $S(q)$ in the format
\begin{fileinfo}[sqt0.simid.out:]
\begin{lstlisting}
step,nq,$q_x$,$q_y$,$q_z$,$S^x(\mathbf{q})$,$S^y(\mathbf{q})$,$S^z(\mathbf{q})$,$S(\mathbf{q})$
\end{lstlisting}
\end{fileinfo}
%%%%%%
\subsubsection*{sra.simid.out} \vindex{sra.simid.out} \index{Correlation function!Spatial}
Is written if the \rkeyword{do_sc} flag is C. Prints out the static correlation function in real space $S(r)$ in the format
\begin{fileinfo}[sra.simid.out:]
\begin{lstlisting}
$|r|$,$S^x(r)$,$S^y(r)$,$S^z(r)$,$S(r)$
\end{lstlisting}
\end{fileinfo}
%%%%%%
\subsubsection*{sqt.simid.out} \vindex{sqt.simid.out} \index{Correlation function!Time resolved}
Is written of the \rkeyword{do_sc} flag is switched on. Prints out the time-resolved structure factor $S(q,t)$ in the format
\begin{fileinfo}[sqt.simid.out:]
\begin{lstlisting}
  nstep, nq, $\mathrm{Re}[S^x(\mathbf{q},t)]$,   
  $\mathrm{Im}[S^x(\mathbf{q},t)]$, $\mathrm{Re}[S^y(\mathbf{q},t)]$, $\mathrm{Im}[S^y(\mathbf{q},t)]$, $\mathrm{Re}[S^z(\mathbf{q},t)]$, $\mathrm{Im}[S^z(\mathbf{q},t)]$
\end{lstlisting}
\end{fileinfo}
%\noindent This file can be very large.
%%%%%%
\subsubsection*{projsqt.simid.out} \vindex{projsqt.simid.out}\index{Correlation function!Projected}
Is written of the \rkeyword{do_sc_proj} flag is switched on. Prints out the same information printed in \rfilename{sqt.simid.out}, but projected to each atom type present in the system. %This file can be very large.
%%%%%%
\subsubsection*{sqw.simid.out}  \vindex{sqw.simid.out} \index{Correlation function!Dynamical structure factor}
Is written of the \rkeyword{do_sc} flag is switched on. Prints out the frequency-resolved structure factor $S(q,\omega)$ in the format
\begin{fileinfo}[sqw.simid.out: \\]
\begin{lstlisting}
  nq, $q_x$, $q_y$, $q_z$, nstep, $S^x(\mathbf{q},\omega)$, $S^y(\mathbf{q},\omega)$, $S^z(\mathbf{q},\omega)$, $S(\mathbf{q},\omega)$
\end{lstlisting}
\end{fileinfo}
%\noindent This file can be very large.
%%%%%%
\subsubsection*{projsqw.simid.out} \vindex{projsqw.simid.out} \index{Correlation function!Projected}
Is written of the \rkeyword{do_sc_proj} flag is switched on. Prints out the same information printed in \rfilename{sqw.simid.out}, but projected to each atom type present in the system. %This file can be very large.
%%%%%%
\subsubsection*{swdos.simid.out} \vindex{swdos.simid.out} \index{Density of states}
Is written of the \rkeyword{do_sc} flag is switched on. Prints out the $S(q,\omega)$ `density of states' as a function of energy.
%%%%%%
\subsubsection*{totenergy.simid.out} \vindex{totenergy.simid.out} \index{Energy}
Is written if the \rkeyword{plotenergy} flag is switched on. Prints out the total energy of the system as a function of time step.
