inpsd.dat keywords: Hamiltonian parameters
==========================================

Hamiltonian parameters
----------------------

.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------------------------------------------------+
|  exchange        |    External file for Heisenberg exchange couplings on the form                                      |
+---------------+--------------------------------------------------------------------------------------------------------+

.. math::

   \mathcal{H}_{\mathrm{XC}} = - \sum_{i\neq j}J_{ij}  \mathbf{e}_i \cdot \mathbf{e}_j ,\label{XC_ham}

where :math:`J_{ij}` is the Heisenberg exchange interaction between atoms :math:`i` and :math:`j`. For an example of the file, see :ref:`exchange`.

.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------------------------------------------------+
|  dm        |    External file for Dzyaloshinskii-Moriya (DM) exchange couplings on the form                            |
+---------------+--------------------------------------------------------------------------------------------------------+

.. math::

  \mathcal{H}_{\mathrm{DM}} = - \sum_{i\neq j}\mathbf{D}_{ij} \cdot \left(\mathbf{e}_i \times \mathbf{e}_j\right),

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

This switch allows the exchange data to be read in according to the tensorial representation of the Heisenberg Hamiltonian, as implemented in the Vienna-Budapest SKKR code [Udvardi2003]_ . In this case, the exchange Hamiltonian is defined as

.. math::

  \mathcal{H}_{\mathrm{Tens}} = \sum_{i,j} \mathbf{e}_i \mathcal{J}_{ij} \mathbf{e}_j.

Here, :math:`\mathcal{J}_{ij}=-J_{ij}\mathcal{I} + \mathcal{J}^S_{ij} +  \mathcal{J}^A_{ij}` is a :math:`3 \times 3` tensor (in which :math:`\mathcal{I}` is the unit matrix), the trace of which is equal to the exchange constant as defined in [Udvardi2003]_. 

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

.. Eq.~(\ref{uniaxial}) or Eq.~(\ref{cubic}), or even both.

where :math:`(m_x, m_y, m_z)=\mathbf{m}`. UppASD is able to read in either uniaxial or cubic anisotropy, or both of them. For bcc Fe, a plausible ``kfile`` might be::

  1   2   -0.020    0.000    0.0    1.0    0.0    0.1
  2   2   -0.020    0.000    0.0    1.0    0.0    0.1   

The first entry lists the atom number, whereas the second entry indicates if the uniaxial (\texttt{1}), cubic (\texttt{2}) or both (\texttt{7}) anisotropies are to be mounted. The second and third entries list the strength of :math:`K_1` and :math:`K_2`, respectively. The fifth to seventh entries indicate the components of the vector :math:`\mathbf{e}_i`. Finally, in the instance of the second entry being set to 7, the final entry indicates the ratio between  :math:`K^{\mathrm{U}}_{\mathrm{ani}}` and  :math:`K^{\mathrm{C}}_{\mathrm{ani}}`.

+---------------+------------------------------------------------------------------------------------------------------------------------+
|  sym      |    Flag to determine the assumed symmetry of the system (*0=none*, 1=cubic, 2=2d cubic (in :math:`xy` plane),              |
|           |    3=hexagonal).                                                                                                           |
+---------------+------------------------------------------------------------------------------------------------------------------------+

It is also possible to provide symmetry operations manually. This is done by setting ``sym`` to 4 and then create an additional input file
``sym.mat`` containing the number of symmetry operations followed by the operations in matrix form. 
An example of ``sym.mat`` for only inversion symmetry can look like::

  2 
    1.0000   0.0000  0.0000
    0.0000   1.0000  0.0000
    1.0000   0.0000  1.0000
   -1.0000   0.0000  0.0000
    0.0000  -1.0000  0.0000
    0.0000   0.0000 -1.0000

Do not forget the identity operation when using custom symmetry operations. The symmetry operations act on ``exchange``, ``bq``, ``pd`` couplings, but not on ``dm`` or ``biqdm`` couplings. Note that the ``sym`` flag only concerns how the program will treat the exchange couplings, it does thus not have to reflect the proper symmetry of the simulated system. Thus, if the exchange interactions given in ``posfile`` are not symmetry reduced, then ``sym`` should be set to :math:`0` even if the system has more symmetry than the identity symmetry.

+---------------+------------------------------------------------------------------------------------------------------------------------+
|  maptype  |    Flag that determines how the coordinates for the different exchange couplings are given.                                |
+---------------+------------------------------------------------------------------------------------------------------------------------+

For *1=coordinates* the coordinates are given in Cartesian or direct coordinates (see ``posfiletype``). For 2 the coordinates of a coupling vector are implicitly given by specifying that the coupling links atom :math:`i` with atom :math:`j` (for an example, see ``dm``).

+---------------+------------------------------------------------------------------------------------------------------------------------+
|  do_prnstruct  |    Flag to print lattice structure (*0=off*/1=on/2=print only coordinates).                                           |
+---------------+------------------------------------------------------------------------------------------------------------------------+

Useful for checking if the system geometry and couplings are correctly set up.

.. %\litem{do_dip} Flag for enabling dipole-dipole interactions (\emph{0=off}/1=on).
