Core input files
================

inpsd.dat
---------

A file with the hardcoded name ``inpsd.dat`` is the main input file necessary to run UppASD. Contained in this file are also the names of the files containing the exchange interactions, the atomic positions, and the atomic moments. Although the names of these files are arbitrary, in this manual they are referred to as the ``jfile``, ``posfile`` and ``momfile``, respectively. Other optional files containing information such as the uniaxial anisotropy and the Dzyaloshinskii-Moriya vectors may also be included, as described below.

The input format is keyword based. The code is programmed to search for given keywords, and then read in the values that follow. If no keyword is given, a default value is set. As an example of a standard ``inpsd.dat`` file layout, input for Fe in body centered cube (bcc) structure is shown below (as found in the examples directory). More advanced examples like supercells and random alloys follows later, but let's keep things simple for now::

  simid bccFe100                                    
  ncell  12              12              12  System size            
  BC     P               P               P   Boundary conditions (0=vacuum, P=periodic)
  cell  -0.5000000000    0.5000000000    0.5000000000
         0.5000000000   -0.5000000000    0.5000000000
         0.5000000000    0.5000000000   -0.5000000000
  Sym    0  Symmetry of exchange bonding vectors (0 for no, 1 for cubic, 2 for 2d cubic, 3 for hexagonal)
  
  posfile    ./posfile
  momfile    ./momfile
  exchange   ./jASD1
  anisotropy ./kfile
  
  do_ralloy  0
  Mensemble  1
  tseed      4499
  maptype    1
  
  SDEalgh    1                         SDE solver: 1=midpoint, 2=heun, 3=heun3, 4=Heun_proper, 5=Depondt
  Initmag    3                         Initial config of moments (1=random, 2=cone, 3=spec., 4=file)
  #restartfile ./restart.bccFe100.out
  
  ip_mode     M                        Initial phase parameters
  ip_mcanneal 1                        --
  10000 300 1.00e-16 0.3               --
  
  mode     M                           S=SD, M=MC
  temp     300                         Measurement phase parameters
  mcNstep  50000
  Nstep    50000
  damping  0.1
  timestep 1.0e-16
  
  do_avrg        Y                     Measure averages
  
  do_cumu        Y
  cumu_step      50
  cumu_buff      10
  
  do_tottraj     N                     Measure moments
  tottraj_step   1000
  
  plotenergy     1
  
  do_sc          C
  do_ams         Y
  do_magdos      Y
  magdos_freq    200
  magdos_sigma   30
  qpoints        C
  
  do_stiffness   Y
  eta_max        12
  eta_min        6
  alat           2.83e-10
  
  
While the meaning of most of the entries in this particular example may be obvious, each input field will be described later in this manual. In short, the input will perform a Monte Carlo simulation at T=300 K allowing 10 000 steps to reach equilbrium and then 50 000 steps to measure observables. However, it is also clear that more information than in ``inpsd.dat`` are required and must be read in from external files in order for the system to be fully defined. These are:

- ``posfile``
- ``momfile``
- ``exchange``

posfile
-------

The positions of the atoms in the unit cell are given in basis vector coordinates.

.. or in Cartesian coordinates.

While these can be listed directly in the ``inpsd.dat`` file, it is typically more convenient to give them in an external file, ``posfile``. For the example above the positions are given as::

  1 1   0.000000  0.000000  0.000000

The first entry indicates the *site number*, whereas the second one indicates the *atom type*. The concept of atom type is central when setting up UppASD simulations because every ``type`` is associated with a particular set of exchange couplings. In this specific case there is only one atom type (and one site), namely Fe with atomic position at the origin.

.. %In the case of random alloy, two extra columns are required for atom component (third column) and its concentration (fourth %column). In the case of a binary 30-70 AB alloy in the B2 structure, the corresponding \rfilename{posfile} looks like:

.. %\begin{fBox} \index{Random alloy}
.. %\begin{Verbatim} 
.. %1 1 1  0.30     0.000000  0.000000  0.000000
.. %1 1 2  0.70     0.000000  0.000000  0.000000
.. %2 1 1  0.30     0.500000  0.500000  0.500000
.. %2 1 2  0.70     0.500000  0.500000  0.500000
.. %\end{Verbatim} 
.. %\end{fBox}


momfile
-------

This file lists the magnetic moments of the atoms in the unit cell. Also, if the ``initmag`` entry is set to 3, the initial direction of the moments is read from this file. For bcc Fe::

  1 1 2.2459 0.0 0.0 1.0 

The first entry indicates the site number (same as first column in posfile), the second entry the chemical type (always 1 for non-random systems), and the third entry indicates the magnitude of the magnetic moment in :math:`\mu_{\mathrm{B}}`, as calculated or estimated from an electronic structure calculation or similar. The last three entries indicate the initial :math:`x`, :math:`y`, and :math:`z` components of the moment (assuming ``initmag`` is set to 3).

.. %For random alloy, magnetic moment of each type is needed. For binary a AB alloy (like Fe-Ni) in the B2 structure, the corresponding %\rfilename{momfile} :
..
.. %\begin{fBox} \index{Random alloy}
.. %\begin{Verbatim}
.. %1 1 2.23 1.0 0.0 0.0
.. %1 2 0.60 1.0 0.0 0.0 
.. %2 1 2.23 1.0 0.0 0.0
.. %2 2 0.60 1.0 0.0 0.0 
.. %\end{Verbatim}  
.. %\end{fBox}


.. _exchange:

exchange
--------

This file lists the exchange couplings within the system. The content and length of this file depends on the symmetry of the system, and the number of atom types present. If no symmetry is used, *i.e.* ``sym 0`` (as in example), all exchange interactions within each interaction shell must be specified. The first shell of the bcc lattice contains 8 interactions, so for ``sym 0`` 8 eight couplings need to be specified in the exchange file. If symmetry is used, then only one interaction in each shell is specified and the program will automatically find the other couplings within the shell depending on the crystal symmetry. For the present Fe example using maptype 1, the first line reads::

  1 1 -0.500 -0.500 -0.500 1.359407144 0.866

The first two entries indicate the sites, which corresponds to the types that one whishes to map, :math:`i` and :math:`j`. In this case as both atoms have the same type, one can indicate the interactions between atoms in site 1 and 2, as 1-1. An example on how to deal with more atom types in the unit cell will be presented later in the manual.

The third, fourth and fifth entries specify the interaction vector between the atoms and depending on choice of the maptype, it has different meaning. Using maptype 1, the vector is specified in Cartesian coordinates. If the SPR-KKR software is used, this corresponds to columns eight, nine and ten in the exchange parameter outfile.
If instead maptype 2 is used, the coordination vector is specified in lattice coordinates and the first line in jfile modifies to::

  1 1 -1 -1 -1 1.359407144 0.866

Once again taking SPR-KKR as an example, that corresponds to columns five, six and seven in the exchange parameters outfile.
The sixth entry in jfile is the exchange energy in mRy and the last entry (not read and optional) is the distance between atoms.

These files together with the inpsd.dat forms the minimal set that is required to run a full ASD or MC simulation. Optionally, there are plenty other external files that may be used for more specific applications and features.  

.. %In systems with more than one atom type, the exchange parameters need to be defined between all atoms (\textit{i.e.} for 2 atoms %the interactions are between the 1-1, 1-2, 2-1 and 2-2 atoms). Note also that if no symmetry is assumed  (\rkeyword{sym} is set to 0), %then the $J_{ij}$  parameters  need to be specified for all neighbours, even those belonging to the same coordination shell.
.. 
.. %In case of random alloy, the \rfilename{jfile} has two additional columns specifying interactions between atom types, like A-A, A-B, B-%A and B-B interactions for a binary alloy. For a B2 binary alloy (model system with only NN interactions), the corresponding %\rfilename{jfile} looks like:
.. 
.. %\begin{fBox} \index{Random alloy}
.. %\begin{Verbatim}
.. %1 1 1 1   0.5 0.5 0.5   1.00000
.. %1 1 1 2   0.5 0.5 0.5   0.50000
.. %1 1 2 1   0.5 0.5 0.5   0.50000
.. %1 1 2 2   0.5 0.5 0.5   0.20000
.. %\end{Verbatim} 
.. %\end{fBox}
