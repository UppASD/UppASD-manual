System and geometry
===================

UppASD features more than 300 keywords for the ``inpsd.dat`` file. In the following some of the keywords are described. Where applicable, the default value for the keyword is underlined.

.. this is subset of the more relevant flags available for inpsd.dat

System parameters
-----------------

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
|  posfiletype  |    Flag to change between *C=Cartesian* or *D=direct* coordinates in posfile.                          |
+---------------+--------------------------------------------------------------------------------------------------------+
|  set_landeg   |    Flag for assigning different values of the gyromagnetic factor for the moments. Set to 0 by default.|
+---------------+--------------------------------------------------------------------------------------------------------+
|  alat         |    Lattice parameter (m). Only needed for spin-torques, dipole treatment and 3TM.                      |
+---------------+--------------------------------------------------------------------------------------------------------+
