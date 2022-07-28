Input files
===========


inpsd.dat
---------

A file with the hardcoded name ``inpsd.dat`` is the main input file necessary to run UppASD. Contained in this file are also the names of the files containing the exchange interactions, the atomic positions, and the atomic moments. Although the names of these files is arbitrary, in this manual they are referred to as the ``jfile``, ``posfile`` and ``momfile``, respectively. Other optional files containing information such as the uniaxial anisotropy and the Dzyaloshinskii-Moriya vectors may also be included, as described below.

The input format is keyword based. The code is programmed to search for given keywords, and then read in the values that follow. If no keyword is given, a default value is set. As an example of a standard ``inpsd.dat`` file layout, input for a Fe in bcc lattice is shown below (as found in the examples directory). More advanced examples like supercells and random alloys follows later, but let's keep things simple for now::

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
  
  
While the meaning of most of the entries in this particular example may be obvious, each input field will be described later in this manual. In short, the input will perform a Monte Carlo simulation at T=300 K allowing 10 000 steps to reach equilbrium and then 50 000 steps to measure observables. However,  it is also clear that more information than in ``inpsd.dat`` are required and must be read in from external files in order for the system to be fully defined. These are:
  

posfile
-------

The positions of the atoms in the unit cell are given in basis vector coordinates.
.. or in Cartesian coordinates.
While these can be listed directly in the ``inpsd.dat`` file, it is typically more convenient to give them in an external file, ``posfile``. For the example above the positions are given as::

  1 1   0.000000  0.000000  0.000000

The first entry indicates the *site number*, whereas the second one indicates the *atom type*. The concept of atom type is central when setting up UppASD simulations because every ``type`` is associated with a particular set defined exchange couplings. In this specific case there is only one atom type (and one site), namely Fe with atomic position at the origin. 

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

The first entry indicates the site number (same as first column in posfile), the second entry the chemical type (always 1 for non-random systems), and the third entry indicates the magnitude of the magnetic moment in :math:`\mu_{\mathrm{B}}`, as calculated or estimated from an electronic structure calculation or similar. The last three entries indicate the initial $x$, $y$, and $z$ components of the moment (assuming  \rkeyword{initmag} is set to 3).

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

exchange
--------

This file lists the exchange couplings within the system. The content and length of this file depends on the symmetry of the system, and the number of atom types present. If no symmetry is used, *i.e.* sym 0 (as in example), all exchange interactions within each interaction shell must be specified. For the bcc lattice, that means that first shell contain 8 interactions and so forth. If symmetry is used, then only one interaction in each shell is specified and the program will automatically found the others within the shell depending on the crystal symmetry. For the present Fe example using maptype 1, the first line reads::

  1 1 -0.500 -0.500 -0.500 1.359407144 0.866

The first two entries indicate the sites, which corresponds to the types that one whishes to map, $i$ and $j$. In this case as both atoms have the same type, one can indicate the interactions between atoms in site 1 and 2, as 1-1, an example on how to deal with more atom types in the unit cell will be presented afterwards. 

The third, fourth and fifth entries specify the interaction  vector between the atoms and depending on choice of the maptype, it has different meaning.  Using maptype 1, the vector is specified in carteisan coodinates.  If the SPR-KKR software is used, that is directly columns eight, nine and ten in the exchange parameter outfile. 
If instead maptype 2 is used, the coordination vector is put in basis coordinates and the first line in jfile modifies to::

  1 1 -1 -1 -1 1.359407144 0.866

Once again taking SPR-KKR as an example, that corresponds to columns five, six and seven in the exchange parameters outfile.
The sixth entry in jfile is the exchange energy in mRy and last entry (not read and optional) is the distance between atoms.

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
.. 

Random alloys and more than one atom in the cell
------------------------------------------------

When the system of interest contains more than one atom in the cell and/or have some chemical disorder, then the setup naturally becomes slightly more complicated. The necessary modifications in the input files are here demonstrated using the system FeCo (found in the examples directory) with the composition 50-50 an example. We are using two different setups for the system, either using an ordered "supercell" or as a random binary alloy. 
 
**FeCo supercell (B2)**
An ordered structure of FeCo with 50-50 composition can be represented in the B2 (CsCl) crystal structure which is a simple cubic lattice with two basis atoms. The Fe atoms occupy the corners and Co atoms the center. In the inpsd.dat input file, the Bravais lattice vectors needs to be specified for the simple cubic lattice::

  cell   1.00000   0.00000   0.00000         
         0.00000   1.00000   0.00000
         0.00000   0.00000   1.00000

The next step is to specify the basis, i.e. the positions of the Fe and Co atoms. We recommend as in the previous example of Fe, always use a separate file (posfile). Fe and Co occupy two different sites in the cell and are of different atom types, the posfile then takes the form::
 
  1 1   0.000000  0.000000  0.000000
  2 2   0.500000  0.500000  0.500000

First line, denotes the Fe atom that is site number 1 and atom type 1 (first and second column) with position 0 0 0 (corner). Second line is the corresponding information for the Co atom that has position 0.5 0.5 0.5 (center in cell).  Once the atom type numbers are set in this file, it will carry over the information in the other files which then needs to be consistent. Now when the simulation cell is set up, we need to specify the magnetic moments and then the (exchange) interactions between them. Starting with magnetic moments, the corresponding momfile::

  1 1 2.7207 0.0 0.0 1.0 
  2 1 1.7202 0.0 0.0 1.0

Once again, the first line specifies the Fe (with site number 1 and chemical type 1) with moment 2.7207 Bohr (from a DFT calculation) and initial moment direction along the z-direction (initmag 3). The second line specifies the same information but for site number 2, i.e. Co that has moment 1.7202 Bohr from calculation. Now both the cell and magnetic moments on each site are specified, what is left to do is the specification of exchange interactions between the moments. From experience, this is the most crucial part in the setup and most easily to get it wrong. The full jfile in the example is longer than specified here (due to the lack of symmetry), here we only show one of the nearest neighbour interactions. We have Fe and Co moments in the cell, a Fe moment could interact with other Fe (Fe-Fe) or with Co (Fe-Co). Vice versa, a Co moment could interact with Fe (Co-Fe) or with other Co (Co-Co). To be complete, we need to specify all the interactions, i.e. Fe-Fe, Fe-Co, Co-Fe and Co-Co interactions. The jfile (using maptype 2) then contains the following blocks::

  1 1  0  0 -1   0.031818272 1.000
  1 2  0  0  0   1.839624404 0.866
  2 1  0  0  0   1.839624404 0.866
  2 2  0  0 -1   0.059966387 1.000

Remember that the ``types`` of atoms that the exchange interactions is valid for, are given in the first two columns of the jfile which specify the ``sites`` :math:`i` and :math:`j`.  The sites correspond to the information given in the posfile. First line then specifies a Fe-Fe interaction, second line Fe-Co, third line Co-Fe and fourth line Co-Co. 

**FeCo random alloy**
UppASD has the capability to deal with chemical disorder in one or several sublattices of a system. Taking Fe-Co as example, it is natually occuring in the bcc lattice (for Co concentrations less  than :math:`\approx 70\%` with random arrangement of the Fe and Co atoms. Internally within the program, a supercell is created with the target composition set by the user. The required input files needs some modifications that are discussed here. First of all, the flag do_ralloy in the inpsd.dat file needs to set to 1. Then, as ususal, the Bravais lattice needs to be specified and in this case we are using the primitive bcc lattice with its lattice vectors::

  cell         -0.5000000    0.5000000    0.5000000
                0.5000000   -0.5000000    0.5000000
                0.5000000    0.5000000   -0.5000000

So far, the setup is not any different from a non-random system. However, the position file looks a bit different. Now we have two chemical types (Fe and Co), each with a certain concentration, that are both situated on the same sublattice::

  1 1 1  0.500   0.000000  0.000000  0.000000
  1 1 2  0.500   0.000000  0.000000  0.000000

Compare to non-random systems, the posfile now has two additional columns. The third column specify the chemical type (Fe or Co), each with its concentration (fourth column). The concentrations do not need to add up to 1, if smaller then the system becomes diluted with random voids (vacancies) in it. In the present example, Fe (chemical type 1) and Co (chemical type 2) both have 50\% concentration. Next, we need to specify the magnetic moments on each sublattice and for each chemical type. The corresponding momfile::

  1 1 2.4850 0.0 0.0 1.0
  1 2 1.7041 0.0 0.0 1.0

The first column always specifies the site number (same as column 1 in the posfile) and column 2 specifies the chemical type (same as column 3 in the posfile). In the example, the first line corresponds to Fe moment and second line the Co moment. The only remaining part is the specification of exchange interactions. Somewhat similar to the FeCo B2 example, we have four distinct set of exchange interactions (Fe-Fe,Fe-Co,Co-Fe and Co-Co), however in this case all interactions are taking place within the same sublattice. A subset of the jfile (first shell) has the following shape (maptype 2)::

  1 1 1 1 -1 -1 -1 1.970049732 0.866
  1 1 1 2 -1 -1 -1 1.947329604 0.866
  1 1 2 1 -1 -1 -1 1.947329604 0.866
  1 1 2 2 -1 -1 -1 1.238957583 0.866

The first and second columns are the same as the jfile for non random systems and specifies the \textbf{sites} $i$ and $j$ and thus their corresponding atomic (sublattice) \textbf{types}. In this case, we only have one sublattice so it is 1 for all interactions. The third and fourth columns specifies the chemical types of the atoms on that particular sublattice and from top to bottom in this example that means Fe-Fe, Fe-Co,Co-Fe and Co-Co interactions.  


Input Entries
-------------

The following entries are currently implemented in UppASD. Where applicable, the default entry setting is underlined.
.. this is subset of the more relevant flags available for inpsd.dat


System parameters
^^^^^^^^^^^^^^^^^

.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------------------------------------------------+
|  simid        |    The 8 character long simulation id. All output files will include the ``simid`` as a label.         |
+---------------+--------------------------------------------------------------------------------------------------------+
|  cell         |    The three lattice vectors describing the cell.                                                      |
+---------------+--------------------------------------------------------------------------------------------------------+
|  ncell        |    Number of repetitions of the cell in each of the lattice vector directions.                         |
+---------------+--------------------------------------------------------------------------------------------------------+
|  bc           |    Boundary conditions (P=periodic, 0=free).                                                           |
+---------------+--------------------------------------------------------------------------------------------------------+
|  natoms       |    Number of atoms in one cell. (Not needed if a ``posfile`` is provided)                              |
+---------------+--------------------------------------------------------------------------------------------------------+
|  ntypes       |    Number of types atoms in one cell. (Not needed if a ``posfile`` is provided)                        |
+---------------+--------------------------------------------------------------------------------------------------------+
|  posfile      |    External file for the positions of the atoms in one cell, with the site number and type of the atom.|
+---------------+--------------------------------------------------------------------------------------------------------+
|  momfile      |    External file describing the magnitudes and directions of magnetic moments.                         |
+---------------+--------------------------------------------------------------------------------------------------------+
|  posfiletype  |    Flag to change between \emph{C=Cartesian} or D=direct coordinates in posfile.                       |
+---------------+--------------------------------------------------------------------------------------------------------+
|  set_landeg   |    Flag for assigning different values of the gyromagnetic factor for the moments. Set to 0 by default.|
+---------------+--------------------------------------------------------------------------------------------------------+


Hamiltonian parameters
^^^^^^^^^^^^^^^^^^^^^^

.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------------------------------------------------+
|  exchange        |    External file for Heisenberg exchange couplings on the form                                      |
+---------------+--------------------------------------------------------------------------------------------------------+

.. math::

   \mathcal{H}_{\mathrm{XC}} = - \sum_{i\neq j}J_{ij}  \mathbf{e}_i \cdot \mathbf{e}_j ,\label{XC_ham}

where :math:`J_{ij}` is the Heisenberg exchange interaction between atoms :math:`i` and :math:`j`. For an example of the file, see the description in Sec.\ref{fxc}.

.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------------------------------------------------+
|  dm        |    External file for Dzyaloshinskii-Moriya (DM) exchange couplings on the form                            |
+---------------+--------------------------------------------------------------------------------------------------------+

.. math::

  \mathcal{H}_{\mathrm{DM}} = - \sum_{i\neq j}\mathbf{D}_{ij}  \cdot \left(\mathbf{e}_i \times \mathbf{e}_j\right),

where :math:`\mathbf{D}_{ij}` is the DM vector. The format is similar to that of the exchange file, *i.e.* in a 2d square lattice it may look something like::

  1 1  1.0000  0.0000 0.0000  0.30000  0.00000  0.00000
  1 1 -1.0000  0.0000 0.0000 -0.30000 -0.00000 -0.00000
  1 1  0.0000  1.0000 0.0000  0.00000  0.30000  0.00000
  1 1  0.0000 -1.0000 0.0000 -0.00000 -0.30000 -0.00000

The first two entries specify site numbers in the chemical unit cell. The third to fifth entries specify the vector :math:`\mathbf{r}_{ij}` in terms of the lattice vectors, and the final three entries specify the DM vector :math:`\mathbf{D}_{ij}`.
.. %Note that in this case the \rkeyword{maptype} flag has been set to 2 in the input file.

.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------------------------------------------------+
|  pd        |    External file for anisotropic symmetric exchange coupling on the form                                  |
+---------------+--------------------------------------------------------------------------------------------------------+

.. math::

  \mathcal{H}_{\mathrm{ani}} = -\sum_{i\neq j} \sum_{\alpha,\beta}\mathbf{J}_{ij}^{\alpha \beta} \ m_i^{\alpha} m_j^{\beta},

where :math:`\mathbf{J}_{ij}^{\alpha \beta}` are the pd couplings and indices :math:`\alpha` and :math:`\beta` refer to the :math:`x`, :math:`y`, and :math:`z`-components of the spins. The format is similar to that of the exchange file. An example file for anisotropic symmetric exchange (here for ``maptype=1`` and ``posfiletype=D``) is::

  1  1   0.25  0.00 -0.25  0.00 -0.01  0.00  0.00  0.00  0.00  

The first two entries indicate the site number and the type of atom, respectively. The third, fourth and fifth entries specify the coordination shell in direct coordinates. The sixth to eleventh entry specify the coupling strength in order :math:`J^{xx}`, :math:`J^{yy}`, :math:`J^{zz}`, :math:`J^{xy}` :math:`(=J^{yx})`, :math:`J^{xz}` :math:`(=J^{zx})`, :math:`J^{yz}` :math:`(=J^{zy})`.

+---------------+--------------------------------------------------------------------------------------------------------+
|  bq       |    External file for biquadratic exchange coupling on the form                                             |
+---------------+--------------------------------------------------------------------------------------------------------+

.. math::

  \mathcal{H}_{\mathrm{bq}} = -\sum_{i\neq j}B_{ij} \left( \mathbf{e}_i\cdot\mathbf{e}_j \right)^2.

The format is identical to that of the ``exchange`` file discussed above, with the values for the exchange couplings :math:`J_{ij}` replaced by the biquadratic exchange couplings :math:`B_{ij}`.

+---------------+--------------------------------------------------------------------------------------------------------+
|  biqdm    |    External file for effective quadratic Dzyaloshinskii-Moriya coupling                                    |
+---------------+--------------------------------------------------------------------------------------------------------+

.. \footnote{For a motivation of this coupling, see Giovannetti \textit{et al.}, Phys. Rev. Lett. \textbf{106}, 026401 (2011)}, which takes the form
   
.. math::

  \mathcal{H}_{\mathrm{biqdm}} = -\sum_{i\neq j}F_{ij} \left( \mathbf{e}_i\times\mathbf{e}_j \right)^2.

The format is identical to that of the ``exchange`` file discussed above, with the values for the exchange couplings :math:`J_{ij}` replaced by the quadratic effective Dzyaloshinskii-Morya exchange coupling :math:`F_{ij}`.

+---------------+--------------------------------------------------------------------------------------------------------+
|  do_tensor  |    Tensorial exchange coupling                                                                           |
+---------------+--------------------------------------------------------------------------------------------------------+

This switch allows the exchange data to be read in according to the tensorial representation of the Heisenberg Hamiltonian, as implemented in the Vienna-Budapest SKKR code.~\cite{Udvardi2003} In this case, the exchange Hamiltonian is defined as

.. math::

  \mathcal{H}_{\mathrm{Tens}} = \sum_{i,j} \mathbf{e}_i \mathcal{J}_{ij} \mathbf{e}_j.

Here, :math:`\mathcal{J}_{ij}=-J_{ij}\mathcal{I} + \mathcal{J}^S_{ij} +  \mathcal{J}^A_{ij}` is a :math:`3 \times 3` tensor (in which :math:`\mathcal{I}` is the unit matrix), the trace of which is equal to the exchange constant as defined in Eq.~(\ref{exchange}) by~\cite{Udvardi2003}. 
.. %
.. %\begin{equation}
.. % J_{ij} = \frac{1}{3} \mathrm{Tr}(\mathcal{J}_{ij}).
.. %\end{equation}
.. %
In this formalism, the anti-symmetric part of the tensor are proportional to the components of the DM vector :math:`\mathbf{D}_{ij}` in Eq.~(\ref{DM_ham}), as :math:`D_{ij}^x=\frac{1}{2}(J_{ij}^{yz}-J_{ij}^{zy})`, :math:`D_{ij}^y=\frac{1}{2}(J_{ij}^{xz}-J_{ij}^{zx})` and :math:`D_{ij}^z=\frac{1}{2}(J_{ij}^{xy}-J_{ij}^{yx})`. In order to define the first shell of exchange parameters in bcc Fe using this formalism, the exchange file would be changed to look as follows::

  0 0 1 2 0.00134 0.0 0.0 0.0 0.00134 0.0 0.0 0.0 0.00134
  0 0 2 1 0.00134 0.0 0.0 0.0 0.00134 0.0 0.0 0.0 0.00134

*NB*: ``maptype`` must be set to 2 in order to use the tensorial format. In addition, since SKKR prints the exchange in Ry, UppASD reads this input in Ry and not in mRy as usual. 

+---------------+--------------------------------------------------------------------------------------------------------+
|  anisotropy    |     External file for anisotropy strengths and directions.                                            |
+---------------+--------------------------------------------------------------------------------------------------------+

The single-ion, or uniaxial, anisotropy is defined as
   
.. math::

  \mathcal{H}^{\mathrm{U}}_{\mathrm{ani}} = \sum_i K_1^{\mathrm{U}} (\mathbf{e}_i\cdot\mathbf{e}_i)^2 + K_2^{\mathrm{U}} (\mathbf{e}_i\cdot\mathbf{e}_i)^4,

where :math:`K_1` and :math:`K_2` are the strength of the linear and four-fold term along an axis with direction :math:`\mathbf{e}_i`. In a cubic system, one must also define the so-called cubic anisotropy, given by

.. math::

  \mathcal{H}^{\mathrm{C}}_{\mathrm{ani}} = \sum_i K_1^{\mathrm{C}} (m_{i,x}^2m_{i,y}^2 + m_{i,y}^2m_{i,z}^2 + m_{i,z}^2m_{i,x}^2 ) + K_2^{\mathrm{C}} m_{i,x}^2 m_{i,y}^2 m_{i,z}^2,

where :math:`(m_x, m_y, m_z)=\mathbf{m}`. UppASD is able to read in either Eq.~(\ref{uniaxial}) or Eq.~(\ref{cubic}), or even both. For bcc Fe, a plausible ``kfile`` might be::

  1   2   -0.020    0.000    0.0    1.0    0.0    0.1
  2   2   -0.020    0.000    0.0    1.0    0.0    0.1   

The first entry lists the atom number, whereas the second entry indicates if the uniaxial (\texttt{1}), cubic (\texttt{2}) or both (\texttt{7}) anisotropies are to be mounted. The second and third entries list the strength of :math:`K_1` and :math:`K_2`, respectively. The fifth to seventh entries indicate the components of the vector :math:`\mathbf{e}_i`. Finally, in the instance of the second entry being set to 7, the final entry indicates the ratio between  :math:`K^{\mathrm{U}}_{\mathrm{ani}}` and  :math:`K^{\mathrm{C}}_{\mathrm{ani}}`.

+---------------+------------------------------------------------------------------------------------------------------------------------+
|  sym      |    Flag to determine the assumed symmetry of the system (\emph{0=none}, 1=cubic, 2=2d cubic (in :math:`xy` plane),         |
|           |    3=hexagonal).                                                                                                           |
+---------------+------------------------------------------------------------------------------------------------------------------------+

It is also possible to provide symmetry operations manually. This is done by setting ``sym`` to 4 and then create an additional input file
``sym.mat`` containing the number of symmetry operations followed by the operations in matrix form. 
An example of {\it sym.mat} for only inversion symmetry can look like::

  2 
    1.0000   0.0000  0.0000
    0.0000   1.0000  0.0000
    1.0000   0.0000  1.0000
   -1.0000   0.0000  0.0000
    0.0000  -1.0000  0.0000
    0.0000   0.0000 -1.0000

Do not forget the identity operation when using custom symmetry operations. The symmetry operations act on ``exchange``, ``bq``, ``pd`` couplings, but not on ``dm`` or ``biqdm`` couplings. Note that the ``sym`` flag only concerns how the program will treat the exchange couplings, it does thus not have to reflect the proper symmetry of the simulated system. *I.e*, if the exchange interactions given in ``posfile`` are not symmetry reduced, then ``sym`` should be set to :math:`0` even if the system has more symmetry than the identity symmetry.

+---------------+------------------------------------------------------------------------------------------------------------------------+
|  maptype  |    Flag that determines how the coordinates for the different exchange couplings are given.                                |
+---------------+------------------------------------------------------------------------------------------------------------------------+

For \emph{1=coordinates} the coordinates are given in Cartesian or direct coordinates (see ``posfiletype``). For 2 the coordinates of a coupling vector are implicitly given by specifying that the coupling links atom :math:`i` with atom :math:`j` (for an example, see ``dm``).

+---------------+------------------------------------------------------------------------------------------------------------------------+
|  do_prnstruct  |    Flag to print lattice structure (\emph{0=off}/1=on/2=print only coordinates).                                      |
+---------------+------------------------------------------------------------------------------------------------------------------------+

Useful for checking if the system geometry and couplings are correctly set up.

.. %\litem{do_dip} Flag for enabling dipole-dipole interactions (\emph{0=off}/1=on).



General simulation parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------------------------------------------------+
|  do_ralloy    |    Flag to set if a random alloy is being simulated (*0=off*/1=on).                                    |
+---------------+--------------------------------------------------------------------------------------------------------+
|  aunits       |    Implement atomic units, *i.e.* set :math:`k_B`, :math:`\hbar`, ... :math:`=1` (Y/\emph{N}). If this |
|               |    is switched on, the \rkeyword{timestep} in SD mode should be roughly 0.1:math:`J_{ij}`.             |
+---------------+--------------------------------------------------------------------------------------------------------+
|  sdealgh      |    Switch for choosing SDE solver (\emph{1=Midpoint}, 4=Heun , 5=Depondt-Mertens). The default option  |
|               |    runs the semi-implicit midpoint solver developed by Mentink \textit{et al}.~\cite{Mentink2010}.     |
|               |    In this case, as when using the Depondt-Mertens solver~\cite{Depondt2009}, the \rkeyword{timestep}  |
|               |    can be as large as  10$^{-16}$ seconds, but this should \textit{always} be checked carefully        |
+---------------+--------------------------------------------------------------------------------------------------------+
|  mensemble    |    Number of ensembles to simulate. The default value is 1, but this may be increased to improve       |
|               |    statistics, especially if investigating laterally confined systems, such as finite                  |
|               |    clusters or other low-dimensional systems.                                                          |
+---------------+--------------------------------------------------------------------------------------------------------+
|  tseed        |    Random number seed for the stochastic field simulating the fluctuations due to temperature.         |
|               |    Default value is 1.                                                                                 |
+---------------+--------------------------------------------------------------------------------------------------------+
|  do_sortcoup  |    Flag to specify if the arrays of couplings should be sorted or not (\emph{Y=yes}, N=no). Of         |
|               |    importance for sampling of polarization. Could be very slow if long range interactions.             |
+---------------+--------------------------------------------------------------------------------------------------------+


Initialization parameters
^^^^^^^^^^^^^^^^^^^^^^^^^

.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------------------------------------------------+
|  initmag      |    Switch for setting up the initial configuration of the magnetic moments (1=Random distribution,     |
|               |    2=Cone, 3=aligned along direction defined in momfile, *4=Read from restartfile*).                   |
+---------------+--------------------------------------------------------------------------------------------------------+
|  restartfile  |    External file containing stored snapshot from previous simulation (used when initmag=4).            |
|               |    The format coincides with the format of the output file ``restart.simid.out``.                      |
+---------------+--------------------------------------------------------------------------------------------------------+
|  mseed        |    Random number seed for magnetic moments if initmag=1. Set to 1 by default.                          |
+---------------+--------------------------------------------------------------------------------------------------------+
|  theta0       |    If \rkeyword{initmag}=2, the magnetic moments are randomly distributed in a cone                    |
|               |    prescribed by this angle, and \rkeyword{phi0}. Set to 0 by default.                                 |
+---------------+--------------------------------------------------------------------------------------------------------+
|  phi0         |    Cone angle for initmag=2. Set to 0 by default.                                                      |
+---------------+--------------------------------------------------------------------------------------------------------+
|  roteul       |    Perform global rotation of magnetization. Set to 0 by default.                                      |
+---------------+--------------------------------------------------------------------------------------------------------+
|  rotang       |    Euler angles describing the rotation if roteul=1.                                                   |
+---------------+--------------------------------------------------------------------------------------------------------+
|  initexc      |    Perform initial excitation of the spin system (\emph{N=none}, I=Vacancies,                          |
|               |    R=Two magnon Raman scattering).                                                                     |
+---------------+--------------------------------------------------------------------------------------------------------+
|  initconc     |    Concentration of vacancies or two magnon spin scattering.                                           |
+---------------+--------------------------------------------------------------------------------------------------------+
|  initneigh    |    eighbour index referring to the list of neighbours for Heisenberg exchange. Determines which spins  |
|               |    to swap in two magnon spin scattering.                                                              |
+---------------+--------------------------------------------------------------------------------------------------------+


Initial phase parameters
^^^^^^^^^^^^^^^^^^^^^^^^

.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------------------------------------------------+
|  ip_mode      |    Mode for initial phase run (S=SD, M=Monte Carlo, H=Heat bath Monte Carlo, \emph{N=none}).           |
+---------------+--------------------------------------------------------------------------------------------------------+
|  ip_temp      |    Temperature for initial phase run if Monte Carlo (ip_mode=M or H).                                  |
+---------------+--------------------------------------------------------------------------------------------------------+
|  ip_hfield    |    External applied field (in units of Tesla) for initial phase run, given in Cartesian coordinates,   |
|               |    *e.g.* ``hfield   1.0   0.0   0.0``.                                                                |
+---------------+--------------------------------------------------------------------------------------------------------+
|  ip_mcnstep   |    Number of Monte Carlo sweeps (MCS) over the system if ip_mode=M or H.                               |
+---------------+--------------------------------------------------------------------------------------------------------+
|  ip_damping   |    Damping parameter $\alpha$ for SD initial phase. Default value is 0.05.                             |
+---------------+--------------------------------------------------------------------------------------------------------+
|  ip_nphase    |    Number of initial phases to be done with SD.                                                        |
+---------------+--------------------------------------------------------------------------------------------------------+

This must be followed by ``ip_nphase`` lines containing number of steps, temperature, timestep and damping for each phase. An example (for an initialization with the temperature decreasing from 300 K to 10 K) can look like::

  ip_nphase 3
  20000 300.0 1.0d-16  0.1
  20000 100.0 1.0d-16  0.1
  30000 010.0 1.0d-16  0.1

.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------------------------------------------------+
|  ip_mcanneal  |     Number of initial phases to be done with MC.                                                       |
+---------------+--------------------------------------------------------------------------------------------------------+

This must be followed by ``ip_mcanneal`` lines containing number of steps and temperature for each phase. An example (for an initialization with the temperature decreasing from 300 K to 10 K) can look like::

  ip_mcanneal 3
  20000 300.0
  20000 100.0 
  30000 010.0 


Measurement phase parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------------------------------------------------+
|  mode         |    Mode for measurement phase run (\emph{S=SD}, M=Monte Carlo, H=Heat bath Monte Carlo).               |
+---------------+--------------------------------------------------------------------------------------------------------+
|  temp         |    Temperature for measurement phase.                                                                  |
+---------------+--------------------------------------------------------------------------------------------------------+
|  hfield       |    External applied field (in units of Tesla) for measurement phase.                                   |
+---------------+--------------------------------------------------------------------------------------------------------+
|  mcnstep      |    Number of Monte Carlo sweeps (MCS) over the system if mode=M or H.                                  |
+---------------+--------------------------------------------------------------------------------------------------------+
|  damping      |    Damping parameter $\alpha$ for SD measurement phase. Default value is 0.05.                         |
+---------------+--------------------------------------------------------------------------------------------------------+
|  timestep     |    Time step between SD iterations. Unless ``aunits Y``, this should typically be set to a value       |
|               |    between :math:`10^{-17}` and :math:`10^{-15}` seconds, depending on the system and SDE solver.      |
+---------------+--------------------------------------------------------------------------------------------------------+
|  relaxtime    |    Relaxation time in LLG+I equation (if sdealgh=11).                                                  |
+---------------+--------------------------------------------------------------------------------------------------------+
|  set_bpulse   |    Add magnetic field pulse ``0=no``, :math:`1-4` for different shapes)                                |
+---------------+--------------------------------------------------------------------------------------------------------+


Parameters for measuring of observables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Typically the measurement of each observable is controlled by two parameters in a combination as follows; ``do_observable`` that enables the measurement and ``observable_step`` that determines the frequency of the measurements. Here the ``observable`` should be replaced by the internal name of the wanted quantity i.e. ``do_avrg`` and ``avrg_step`` for the average magnetization.

.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------------------------------------------------+
|  plotenergy   |    Flag to enable the calculation of the energy of the system projected to the different components of |
|               |    the Hamiltonian. ``{0=off}/1=on``)                                                            .     |
+---------------+--------------------------------------------------------------------------------------------------------+
|  do_avrg      |    Sample and print average magnetization, and its higher order moments. ``Y/N``                       |
+---------------+--------------------------------------------------------------------------------------------------------+
|  do_proj_avrg |    Sample and print type (*i.e*. sublattice) projected average moments. (``Y/N/A``).                   |
+---------------+--------------------------------------------------------------------------------------------------------+
| do_projch_avrg|    Sample and print chemical (*i.e.*} sublattice) projected average moments (``Y/N/A``).               |
+---------------+--------------------------------------------------------------------------------------------------------+
|  avrg_step    |    Add magnetic field pulse (\emph{0=no}, $1-4$ for different shapes)                                  |
+---------------+--------------------------------------------------------------------------------------------------------+
|  avrg_buff    |    Number of samplings of averages to buffer between printing to file. Set to 10 by default.           |
+---------------+--------------------------------------------------------------------------------------------------------+
|  do_cumu      |    Sample cumulants (Y/N). Automatically enabled for Monte Carlo simulations.                          |
+---------------+--------------------------------------------------------------------------------------------------------+
|  cumu_step    |    Number of time steps between sampling of cumulants. Set to 25 by default.                           |
+---------------+--------------------------------------------------------------------------------------------------------+
|  do_tottraj   |    Sample and print all trajectories (moments) in the system. (Y/N). Generates the (rather large)      |
|               |    ``moments.simid.out`` file.                                                                         |
+---------------+--------------------------------------------------------------------------------------------------------+
|  tottraj_step |    Number of time steps between samplings of moments. Set to 1000 by default.                          |
+---------------+--------------------------------------------------------------------------------------------------------+
|  ntraj        |    Number of individual trajectories to sample and print. Followed by ``ntraj`` lines describing atoms |
|               |    to sample, time step between samples and steps to buffer. Set to 0 by default.                      |
+---------------+--------------------------------------------------------------------------------------------------------+
|  do_pol       |    Sample and print average ferroelectric polarization (Y/N) according to the expression               |
|               |    :math:`P\propto \gamma\sum_{i,j}\hat{\mathbf{e}}_{ij}\times(\mathbf{m}_i\times\mathbf{m}_j)`.       |
|               |    Uses the neighbour lists set up for exchange but here the sum is performed up to a threshold        |
|               |    ``max_pol_nn``. For this construction to work, it is important to set the flag ``do_sortcoup N`.    |
+---------------+--------------------------------------------------------------------------------------------------------+
|  max_pol_nn   |    Number of neighbours to use when evaluating the polarization.                                       |
+---------------+--------------------------------------------------------------------------------------------------------+
|  pol_step     |    Number of time steps between sampling of polarization averages. Set to 100 by default.              |
+---------------+--------------------------------------------------------------------------------------------------------+
|  pol_buff     |  Number of samplings of polarization averages to buffer between printing to file. Set to 10 by default.|
+---------------+--------------------------------------------------------------------------------------------------------+
|  do_stiffness |    Calculation of spin-wave stiffness (and tensor) and micromagnetic exchange constant (Y/N).          |
+---------------+--------------------------------------------------------------------------------------------------------+
|  eta_min      |    Lowest value of auxiliary convergence parameter in stiffness calculation (recommended around 6-8)   |
+---------------+--------------------------------------------------------------------------------------------------------+
|  eta_max      |    Largest value of auxiliary convergence parameter in stiffness calculation (recommended around 10-12)|
+---------------+--------------------------------------------------------------------------------------------------------+
|  alat         |    Lattice constant (in m) for calculation of exchange stiffness                                       |
+---------------+--------------------------------------------------------------------------------------------------------+


Parameters for measuring of correlation functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. %spin wave sampling and correlations

.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------------------------------------------------+
|  do_sc        |   Flag to determine if spin correlations should be analysed (Q=S($\mathbf{q}$,$\omega$), \emph{N=no},  |
|               |   C=G(r)). Setting this flag to Q or C measures space- and time-displaced correlation functions.       |
|               |   The spatial time dependent correlation function $C(\mathbf{r},t)$ is defined as                      |
+---------------+--------------------------------------------------------------------------------------------------------+

.. math::

  C^k (\mathbf{r}-\mathbf{r'},t) = \langle m^k_{\mathbf{r}}(t) m^k_{\mathbf{r'}}(0) \rangle - \langle m^k_{\mathbf{r}}(t) \rangle \langle m^k_{\mathbf{r'}}(0) \rangle,
  \label{eqn:cf}

where the angular brackets signify an ensemble average and $k$ the Cartesian component. The dynamical structure factor is then obtained by Fourier transforming $C(\mathbf{r},t)$ as
  
.. math::

  S^k(\mathbf{q},\omega) = \frac{1}{\sqrt{2\pi}N} \sum_{\mathbf{r},\mathbf{r'}} e^{i\mathbf{q}\cdot(\mathbf{r}-\mathbf{r'})} \int_{-\infty}^{\infty} e^{i\omega t} C^k (\mathbf{r}-\mathbf{r'},t) dt,

and this function describes the energy dispersion for excited spin waves present in the simulated system.~\cite{Bergman2010}. If the flag is set to C, the static correlation function :math:`G(\mathbf{r})` and its Fourier transform :math:`S(\mathbf{q})` are measured. By locating the maximum of :math:`S(\mathbf{q})`, the ordering vector of the simulated system can be determined. In this case it is important to have a ``qfile`` that includes :math:`\mathbf{q}`-vectors in the whole Brillouin zone.
.. %By default, both $S(\mathbf{q},\omega)$ and $S(\mathbf{q},t)$ are written to files but if only one of these correlation functions is %wanted, a selective printing can be obtained by giving the values \rkeyword{do_sc} = W or T, instead of the normal choice of %\rkeyword{do_sc} = Y. 
In order to obtain a useful $S(\mathbf{q},\omega)$ measurement, it is important to understand the sampling of the function that is determined by ``sc_nstep``, ``sc_step``, and ``timestep``.

.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------------------------------------------------+
|  do_sc_proj   |    Flag to determine if type projected spin correlation should be analyzed (Q=yes, C=G(r),\emph{N=no}).|
+---------------+--------------------------------------------------------------------------------------------------------+
|  do_sc_projc  |    Flag to determine if chemical type projected spin correlation should be analyzed of random alloys   |
|               |    (Q=yes, C=G(r),\emph{N=no}).                                                                        |
+---------------+--------------------------------------------------------------------------------------------------------+
|  do_qt_traj   |    Flag to determine if the time evolution of the equal time spin correlation $S(\mathbf{q})$ should be|
|               |    written to file (Y=yes,\emph{N=no}).                                                                |
+---------------+--------------------------------------------------------------------------------------------------------+

This works only if ``do_sc C``. The function :math:`S(\mathbf{q})` is sampled every ``sc_sep`` time step and can give insight in the phase transitions in systems with more than one magnetic order parameter. Suggested use is to first determine the magnetic phase diagram and the associated ordering vectors by sampling :math:`S(\mathbf{q})` (as described above). The 
 order parameters can then be specified in a ``qpoints`` file and followed in simulations where the systems is driven out of equilibrium by an external perturbation in form of an applied magnetic field, a heat pulse or a two-magnon Raman scattering excitation.
 
.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------------------------------------------------+
|  sc_mode      |    Flag to determine when to transform the spin correlations (0=in memory, 1=in scratch file,          |
|               |    \emph{2=on the fly}). Options 0 and 1 generate enormous files.                                      |
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
|               |    calculations (\emph{1=box}, 2=Hann, 3=Hamming, 4=Blackman-Harris).                                  |
+---------------+--------------------------------------------------------------------------------------------------------+
|  do_ams       |    Spin wave dispersion from the Fourier transform av the exchange interactions, so-called Adiabatic   |
|               |    Magnon Spectra (AMS) (Y=yes, \emph{N=no}). This version only handles AMS in collinear magnetic      |
|               |    structures but it is very fast and can therefore be a good option for comparison with the full      |
|               |    dynamical spectra. If ``do_ams Y`` then one must provide a qfile just as in the case of             |
|               |    :math:`S(\mathbf{q},\omega)`.                                                                       |
+---------------+--------------------------------------------------------------------------------------------------------+
|  do_magdos    |    Magnon density of states (MDOS) from AMS (Y=yes, \emph{N=no}, A=read from file).                    |
+---------------+--------------------------------------------------------------------------------------------------------+
|  magdos_freq  |    Number of frequencies in MDOS calculation from AMS. Around 200 is recommended.                      |
+---------------+--------------------------------------------------------------------------------------------------------+
|  magdos_sigma |    Gaussian broadening (in meV) in MDOS calculation from AMS (around 30 is recommended).               |
+---------------+--------------------------------------------------------------------------------------------------------+
|  do_autocorr  |    Flag to enable autocorrelation sampling (Y=yes, \emph{N=no}).                                       |
+---------------+--------------------------------------------------------------------------------------------------------+
|  acfile       |    External file containing waiting times for the autocorrelation measurements.                        |
+---------------+--------------------------------------------------------------------------------------------------------+

.. %\litem{sc_navrg} Number of spin correlation measurements to average over.

.. \litem{do_sc_local_axis} Modify the sampling for $S(q,\omega)$ so that $S^\bot$ and $S^\parallel$ are sampled instead of $S^x$, $S^y$, $S^z$. This normally improves the simulated spectra for ferromagnets but should be used with care since it can, if misused, suppress low-level excitations. (\emph{Y},N)

.. \litem{sc_local_axis_mix} Determines the rate of updating the local quantization axis used when \rkeyword{do_sc_local_axis}=Y. Values larger than zero can be useful if there are unwanted fluctuations such as global rotations of the whole systems, which can happen for in particular for finite systems such as clusters.





