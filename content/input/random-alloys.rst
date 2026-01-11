Random alloys and multi-atom cells
==================================

When the system of interest contains more than one atom in the cell and/or have some chemical disorder, then the setup naturally becomes slightly more complicated. The necessary modifications in the input files are here demonstrated using the system FeCo (found in the examples directory) with the composition 50-50 an example. We are using two different setups for the system, either using an ordered "supercell" or as a random binary alloy. 
 
**FeCo supercell (B2)**
An ordered structure of FeCo with 50-50 composition can be represented in the B2 (CsCl) crystal structure which is a simple cubic lattice with two basis atoms. The Fe atoms occupy the corners and Co atoms the center. In the inpsd.dat input file, the Bravais lattice vectors needs to be specified for the simple cubic lattice::

  cell   1.00000   0.00000   0.00000         
         0.00000   1.00000   0.00000
         0.00000   0.00000   1.00000

The next step is to specify the basis, i.e. the positions of the Fe and Co atoms. We recommend as in the previous example of Fe, always use a separate file (posfile). Fe and Co occupy two different sites in the cell and are of different atom types, the posfile then takes the form::
 
  1 1   0.000000  0.000000  0.000000
  2 2   0.500000  0.500000  0.500000

First line, denotes the Fe atom that is site number 1 and atom type 1 (first and second column) with position 0 0 0 (corner). Second line is the corresponding information for the Co atom that has position 0.5 0.5 0.5 (center in cell). Once the atom type numbers are set in this file, it will carry over the information in the other files which then needs to be consistent. Now when the simulation cell is set up, we need to specify the magnetic moments and then the (exchange) interactions between them. Starting with magnetic moments, the corresponding momfile::

  1 1 2.7207 0.0 0.0 1.0 
  2 1 1.7202 0.0 0.0 1.0

Once again, the first line specifies the Fe (with site number 1 and chemical type 1) with moment 2.7207 :math:`\mu_{\mathrm{B}}` (from a DFT calculation) and initial moment direction along the z-direction (``initmag`` 3). The second line specifies the same information but for site number 2, i.e. Co that has moment 1.7202 :math:`\mu_{\mathrm{B}}` from calculation.

Now that both the cell and magnetic moments on each site are specified, what is left to do is the specification of exchange interactions between the moments. From experience, this is the most crucial part in the setup, and most easily to get it wrong. The full jfile in the example is longer than specified here (due to the lack of symmetry), here we only show one of the nearest neighbour interactions. We have Fe and Co moments in the cell, a Fe moment could interact with other Fe (Fe-Fe) or with Co (Fe-Co). Vice versa, a Co moment could interact with Fe (Co-Fe) or with other Co (Co-Co). To be complete, we need to specify all the interactions, i.e. Fe-Fe, Fe-Co, Co-Fe and Co-Co interactions. The jfile (using ``maptype`` 2) then contains the following blocks::

  1 1  0  0 -1   0.031818272 1.000
  1 2  0  0  0   1.839624404 0.866
  2 1  0  0  0   1.839624404 0.866
  2 2  0  0 -1   0.059966387 1.000

Remember that the ``types`` of atoms that the exchange interactions is valid for, are given in the first two columns of the jfile which specify the ``sites`` :math:`i` and :math:`j`.  The sites correspond to the information given in the posfile. First line then specifies a Fe-Fe interaction, second line Fe-Co, third line Co-Fe and fourth line Co-Co. 

**FeCo random alloy**
UppASD has the capability to deal with chemical disorder in one or several sublattices of a system. Taking Fe-Co as example, it is natually occuring in the bcc lattice (for Co concentrations less  than :math:`\approx 70\%` with random arrangement of the Fe and Co atoms. Internally within the program, a supercell is created with the target composition set by the user. The required input files needs some modifications that are discussed here. First of all, the flag do_ralloy in the inpsd.dat file needs to be set to 1. Then, as ususal, the Bravais lattice needs to be specified and in this case we are using the primitive bcc lattice with its lattice vectors::

  cell         -0.5000000    0.5000000    0.5000000
                0.5000000   -0.5000000    0.5000000
                0.5000000    0.5000000   -0.5000000

So far, the setup is not any different from a non-random system. However, the position file looks a bit different. Now we have two chemical types (Fe and Co), each with a certain concentration, that are both situated on the same sublattice::

  1 1 1  0.500   0.000000  0.000000  0.000000
  1 1 2  0.500   0.000000  0.000000  0.000000

Compare to non-random systems, the posfile now has two additional columns. The third column specify the chemical type (Fe or Co), each with its concentration (fourth column). The concentrations do not need to add up to 100%, if smaller then the system becomes diluted with random voids (vacancies) in it. In the present example, Fe (chemical type 1) and Co (chemical type 2) both have 50\% concentration. Next, we need to specify the magnetic moments on each sublattice and for each chemical type. The corresponding momfile::

  1 1 2.4850 0.0 0.0 1.0
  1 2 1.7041 0.0 0.0 1.0

The first column always specifies the site number (same as column 1 in the posfile) and column 2 specifies the chemical type (same as column 3 in the posfile). In the example, the first line corresponds to Fe moment and second line the Co moment. The only remaining part is the specification of exchange interactions. Somewhat similar to the FeCo B2 example, we have four distinct set of exchange interactions (Fe-Fe,Fe-Co,Co-Fe and Co-Co), however in this case all interactions are taking place within the same sublattice. A subset of the jfile (first shell) has the following shape (``maptype`` 2)::

  1 1 1 1 -1 -1 -1 1.970049732 0.866
  1 1 1 2 -1 -1 -1 1.947329604 0.866
  1 1 2 1 -1 -1 -1 1.947329604 0.866
  1 1 2 2 -1 -1 -1 1.238957583 0.866

The first and second columns are the same as the jfile for non random systems and specifies the *sites* :math:`i` and :math:`j` and thus their corresponding atomic (sublattice) *types*. In this case, we only have one sublattice so it is 1 for all interactions. The third and fourth columns specifies the chemical types of the atoms on that particular sublattice and from top to bottom in this example that means Fe-Fe, Fe-Co,Co-Fe and Co-Co interactions.
