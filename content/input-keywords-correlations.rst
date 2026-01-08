inpsd.dat keywords: correlations
=================================

Parameters for measuring of correlation functions
-------------------------------------------------
.. %spin wave sampling and correlations

.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------------------------------------------------+
|  do_sc        |   Flag to determine if spin correlations should be analysed (Q= :math:`S(\mathbf{q},\omega)`, *N=no*), |
|               |   C= :math:`G(r)`). Setting this flag to Q or C measures space- and time-displaced correlation         |
|               |   functions. The spatial time dependent correlation function $C(\mathbf{r},t)$ is defined as           |
+---------------+--------------------------------------------------------------------------------------------------------+

.. math::

  C^k (\mathbf{r}-\mathbf{r'},t) = \langle m^k_{\mathbf{r}}(t) m^k_{\mathbf{r'}}(0) \rangle - \langle m^k_{\mathbf{r}}(t) \rangle \langle m^k_{\mathbf{r'}}(0) \rangle,
  \label{eqn:cf}

where the angular brackets signify an ensemble average and :math:`k` the Cartesian component. The dynamical structure factor is then obtained by Fourier transforming :math:`C(\mathbf{r},t)` as
  
.. math::

  S^k(\mathbf{q},\omega) = \frac{1}{\sqrt{2\pi}N} \sum_{\mathbf{r},\mathbf{r'}} e^{i\mathbf{q}\cdot(\mathbf{r}-\mathbf{r'})} \int_{-\infty}^{\infty} e^{i\omega t} C^k (\mathbf{r}-\mathbf{r'},t) dt,

and this function describes the energy dispersion for excited spin waves present in the simulated system [Bergman2010]_. If the flag is set to C, the static correlation function :math:`G(\mathbf{r})` and its Fourier transform :math:`S(\mathbf{q})` are measured. By locating the maximum of :math:`S(\mathbf{q})`, the ordering vector of the simulated system can be determined. In this case it is important to have a ``qfile`` that includes :math:`\mathbf{q}` -vectors in the whole Brillouin zone.

.. %By default, both $S(\mathbf{q},\omega)$ and $S(\mathbf{q},t)$ are written to files but if only one of these correlation functions is %wanted, a selective printing can be obtained by giving the values \rkeyword{do_sc} = W or T, instead of the normal choice of %\rkeyword{do_sc} = Y. 

In order to obtain a useful :math:`S(\mathbf{q},\omega)` measurement, it is important to understand the sampling of the function that is determined by ``sc_nstep``, ``sc_step``, and ``timestep``.

.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------------------------------------------------+
|  do_sc_proj   |    Flag to determine if type projected spin correlation should be analyzed (Q=yes, C=G(r), *N=no*)     |
+---------------+--------------------------------------------------------------------------------------------------------+
|  do_sc_projc  |    Flag to determine if chemical type projected spin correlation should be analyzed of random alloys   |
|               |    (Q=yes, C=G(r), *N=no*).                                                                            |
+---------------+--------------------------------------------------------------------------------------------------------+
|  do_qt_traj   |    Flag to determine if the time evolution of the equal time spin correlation :math:`S(\mathbf{q})`    |
|               |    should be written to file (Y=yes, *N=no*).                                                          |
+---------------+--------------------------------------------------------------------------------------------------------+

This works only if ``do_sc C``. The function :math:`S(\mathbf{q})` is sampled every ``sc_sep`` time step and can give insight in the phase transitions in systems with more than one magnetic order parameter. Suggested use is to first determine the magnetic phase diagram and the associated ordering vectors by sampling :math:`S(\mathbf{q})` (as described above).
The order parameters can then be specified in a ``qpoints`` file and followed in simulations where the systems is driven out of equilibrium by an external perturbation in form of an applied magnetic field, a heat pulse or a two-magnon Raman scattering excitation.

.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------------------------------------------------+
|  sc_mode      |    Flag to determine when to transform the spin correlations (0=in memory, 1=in scratch file,          |
|               |    *2=on the fly*). Options 0 and 1 generate enormous files.                                           |
+---------------+--------------------------------------------------------------------------------------------------------+
|  sc_nstep     |    Number of steps to sample. This number sets the resolution of time/frequency based correlation      |
|               |    functions by deciding the number of measured times/frequencies to include in the calculation.       |
+---------------+--------------------------------------------------------------------------------------------------------+
|  sc_step      |    Number of time steps between each sampling. This number determines the time/frequency range over    |
|               |    which correlation functions are measured. The mininum sample time is given by                       |
|               |    ``timestep`` * ``sc_step`` and the maximal sampling time is then determinded by                     |
|               |    ``sc_nstep`` * ``timestep`` * ``sc_step``. The minimal/maximal frequencies are then determined by   |
|               |    the inverse of the maximal/minimal sampling time.                                                   |
+---------------+--------------------------------------------------------------------------------------------------------+
|  sc_sep       |    Number of time steps between the start of subsequent spin correlation measurements.                 |
+---------------+--------------------------------------------------------------------------------------------------------+
|  qpoints      |    Flag for for generation of q-point mesh necessary for :math:`S(\mathbf{q},\omega)` calculations.    |
|               |    (F=external file with Cartesian coordinates}, A=automatic, C=full cell, P=extended plane spanned by |
|               |    the first and third reciprocal lattice vector, D=external file with direct coordinates).            |
+---------------+--------------------------------------------------------------------------------------------------------+
|  sc_window_fun|    Choice of windowing function for the Fourier transforms used in :math:`S(\mathbf{q},\omega)`        |
|               |    calculations (*1=box*, 2=Hann, 3=Hamming, 4=Blackman-Harris).                                       |
+---------------+--------------------------------------------------------------------------------------------------------+
|  do_ams       |    Spin wave dispersion from the Fourier transform av the exchange interactions, so-called Adiabatic   |
|               |    Magnon Spectra (AMS) (Y=yes, *N=no*). This version only handles AMS in collinear magnetic           |
|               |    structures but it is very fast and can therefore be a good option for comparison with the full      |
|               |    dynamical spectra. If ``do_ams Y`` then one must provide a qfile just as in the case of             |
|               |    :math:`S(\mathbf{q},\omega)`.                                                                       |
+---------------+--------------------------------------------------------------------------------------------------------+
|  do_magdos    |    Magnon density of states (MDOS) from AMS (Y=yes, *N=no*, A=read from file).                         |
+---------------+--------------------------------------------------------------------------------------------------------+
|  magdos_freq  |    Number of frequencies in MDOS calculation from AMS. Around 200 is recommended.                      |
+---------------+--------------------------------------------------------------------------------------------------------+
|  magdos_sigma |    Gaussian broadening (in meV) in MDOS calculation from AMS (around 30 is recommended).               |
+---------------+--------------------------------------------------------------------------------------------------------+
|  do_autocorr  |    Flag to enable autocorrelation sampling (Y=yes, *N=no*).                                            |
+---------------+--------------------------------------------------------------------------------------------------------+
|  acfile       |    External file containing waiting times for the autocorrelation measurements.                        |
+---------------+--------------------------------------------------------------------------------------------------------+

.. %\litem{sc_navrg} Number of spin correlation measurements to average over.

.. \litem{do_sc_local_axis} Modify the sampling for $S(q,\omega)$ so that $S^\bot$ and $S^\parallel$ are sampled instead of $S^x$, $S^y$, $S^z$. This normally improves the simulated spectra for ferromagnets but should be used with care since it can, if misused, suppress low-level excitations. (\emph{Y},N)

.. \litem{sc_local_axis_mix} Determines the rate of updating the local quantization axis used when \rkeyword{do_sc_local_axis}=Y. Values larger than zero can be useful if there are unwanted fluctuations such as global rotations of the whole systems, which can happen for in particular for finite systems such as clusters.
