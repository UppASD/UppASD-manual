.. _label-stiffness:

Spin Stiffness and Micromagnetic Parameters
===========================================

The stiffness calculation feature extracts micromagnetic parameters from atomistic exchange interactions, bridging the gap between **atomistic spin dynamics** and **micromagnetic continuum models**. This enables:

- **Spin-wave stiffness** :math:`D` calculation for ferromagnets
- **Exchange stiffness** :math:`A` (micromagnetic exchange constant) extraction
- **DMI spiralization** tensor :math:`D_0` for systems with Dzyaloshinskii-Moriya interaction
- **Mean-field critical temperature** :math:`T_c^{\text{MFA}}` estimation
- **Random alloy analysis** with site-resolved stiffness and :math:`T_c` values

These parameters are essential for:

- **Connecting atomistic and continuum descriptions**: Directly comparing atomistic Heisenberg exchange with micromagnetic models
- **Spin-wave spectroscopy interpretation**: Relating measured spin-wave dispersions to atomistic parameters
- **Domain wall studies**: Predicting domain wall widths from exchange and anisotropy
- **Skyrmion stability**: Understanding skyrmion sizes from DMI and exchange balance


Physics Background
------------------

Spin-Wave Stiffness (Pajda Formalism)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The formalism follows **Pajda et al., PRB 64, 174402 (2001)**, which defines the spin-wave stiffness for ferromagnetic systems as:

.. math::

   D = \frac{2}{3} \sum_{ij} J_{ij} r_{ij}^2

where:

- :math:`J_{ij}` is the exchange interaction between sites :math:`i` and :math:`j`
- :math:`r_{ij}` is the distance vector between the sites
- The sum runs over all exchange pairs weighted by distance squared

This formulation accounts for multi-sublattice systems and random alloys by constructing a **stiffness matrix** :math:`D_{ab}^{\alpha\beta}` where :math:`a,b` label sublattices and :math:`\alpha,\beta` denote spatial directions (x, y, z).

**Extension to Anisotropic Systems:**

For non-cubic systems, the full tensorial form is calculated:

.. math::

   D_{\alpha\beta} = 2 \sum_{ij} J_{ij} r_{ij}^\alpha r_{ij}^\beta

This captures directional variations in spin-wave propagation.

Exchange Stiffness (Micromagnetic Constant)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The micromagnetic **exchange stiffness** :math:`A` relates to the spin-wave stiffness via:

.. math::

   A = \frac{D \cdot M_s}{4}

where :math:`M_s` is the saturation magnetization. Units:

- :math:`D`: meV·Å²
- :math:`M_s`: μ_B/m³
- :math:`A`: pJ/m (typical values: 5-20 pJ/m for transition metal ferromagnets)

This parameter appears in the **micromagnetic energy density**:

.. math::

   E_{\text{ex}} = A \int |\nabla \mathbf{m}|^2 \, dV

where :math:`\mathbf{m}` is the normalized magnetization direction.

DMI Spiralization Tensor
~~~~~~~~~~~~~~~~~~~~~~~~~

For systems with **Dzyaloshinskii-Moriya interaction** (DMI), the spiralization tendency is quantified by:

.. math::

   D_{0}^{\alpha\beta} = \sum_{ij} \mathbf{D}_{ij}^\alpha \cdot \mathbf{r}_{ij}^\beta

where :math:`\mathbf{D}_{ij}` is the DMI vector between sites :math:`i,j`. This tensor determines:

- **Spiral wavelength**: :math:`\lambda = 4\pi D / D_0`
- **Skyrmion size**: Proportional to :math:`D/D_0`
- **Helical pitch**: Intrinsic pitch of chiral spin textures

Units: meV·Å (energy × length).

Convergence and Eta Parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The stiffness calculation employs an exponential damping factor to ensure convergence:

.. math::

   D(\eta) = \sum_{ij} J_{ij} r_{ij}^2 \exp(-\eta r_{ij})

The parameter :math:`\eta` (eta) weights near-neighbor vs. far-neighbor contributions. The final stiffness is extrapolated to :math:`\eta \to 0` using:

1. **Rational polynomial fitting** (Padé approximants): More robust for oscillatory behavior
2. **Least-squares (LSQ) fitting**: Standard polynomial fit to :math:`D(\eta)` vs. :math:`\eta`

Users specify an eta range [``eta_min``, ``eta_max``] for extrapolation. Typical values: ``eta_min=10``, ``eta_max=20``.

Activation and Configuration
-----------------------------

Main Enable Flags
~~~~~~~~~~~~~~~~~

.. code-block:: fortran

   do_stiffness = 'N'       ! Enable stiffness calculation
   do_dm_stiffness = 'N'    ! Enable DMI spiralization calculation
   prn_J0_matrix = 'N'      ! Print full J0 exchange matrix

**Calculation Timing:**

Stiffness calculations are performed **once during initialization** (not dynamically during time evolution). They extract ground-state properties from the Hamiltonian.

**When to Use:**

- ``do_stiffness = 'Y'``: Always enable for ferromagnets when extracting micromagnetic parameters
- ``do_dm_stiffness = 'Y'``: Enable when DMI is present (``dm_file`` specified) to study chiral systems
- ``prn_J0_matrix = 'Y'``: For debugging or detailed analysis of exchange matrix structure

Convergence Parameters
~~~~~~~~~~~~~~~~~~~~~~

The eta convergence range controls the extrapolation to :math:`\eta = 0`:

.. list-table:: Convergence Parameters
   :widths: 25 15 50
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ``eta_min``
     - int
     - Minimum eta index for extrapolation (typically 6-12)
   * - ``eta_max``
     - int
     - Maximum eta index for extrapolation (typically 12-24)

**Physical Meaning:**

- Smaller ``eta_min``, ``eta_max``: Include more distant neighbors (slower decay)
- Larger values: Focus on nearest neighbors (faster decay)

**Recommended Settings:**

- **Simple cubic/BCC/FCC**: ``eta_min=10``, ``eta_max=20``
- **Complex structures**: ``eta_min=6``, ``eta_max=12``
- **Random alloys**: ``eta_min=6``, ``eta_max=12`` (more sampling needed)

The code internally uses :math:`\eta_{\text{eff}} = 0.1 \times \text{eta}` as the damping parameter in Ångströms.

Output Files and Results
------------------------

Main Output File
~~~~~~~~~~~~~~~~

When ``do_stiffness = 'Y'``, results are written to:

.. code-block:: fortran

   micro_ASD.<simid>.out

This file contains comprehensive micromagnetic information:

**File Structure:**

.. code-block:: text

   **************** MICROMAGNETIC INFORMATION *********************
     Unit cell volume         :  1.234E-29  m^3
     Saturation magnetization :  1234.56  MA/m
     Anisotropy density       :  0.245  MJ/m^3
     Anisotropy per atom      :  0.456  meV
     Domain wall width        :  12.34  nm
     ETA MIN :     10
     ETA MAX :     20
   ****************************************************************
   
   Exchange stiffness at T=0 K:    15.234  pJ/m
   **************** EXCHANGE STIFFNESS MATRIX [pJ/m] **************
         14.123        0.000        0.000
          0.000       14.123        0.000
          0.000        0.000       15.678
   ****************************************************************
   
   Spin wave stiffness from Jijs:   123.45 ±   0.12  meVÅ^2
   **************** SPIN WAVE STIFFNESS MATRIX [meVA^2] ***********
        122.345        0.000        0.000
          0.000      122.345        0.000
          0.000        0.000      125.678
   ****************************************************************
   
   Exchange stiffness at T=0 K (LSQ):    15.123  pJ/m
   [LSQ matrices follow similar format]
   
   **************** J0 VECTOR [meV] *******************************
        1      123.456
        2      234.567
   ****************************************************************
   
   Tc-MFA from stiffness :     1234.5 K

DMI Spiralization Output
~~~~~~~~~~~~~~~~~~~~~~~~

When ``do_dm_stiffness = 'Y'``, additional DMI information is appended:

.. code-block:: text

   DM matrix (q,n)   D_{ij} x r_{ij}
   **************** DMI SPIRALIZATION MATRIX [meVA] ***************
               r_x         r_y         r_z
    D_x      12.345      0.000      0.000
    D_y       0.000     12.345      0.000
    D_z       0.000      0.000     15.678
   ****************************************************************
   
   ******* DMI SPIRALIZATION MATRIX (LSQ fit) [meVA] **************
   [LSQ fit values]
   ****************************************************************
   
   **************** DMI SPIRALIZATION MATRIX [mJ/m^2] *************
   [Converted to micromagnetic units]
   ****************************************************************
   
   **************** SPIRAL WAVELENGTH MATRIX (STIFF/DMI) [A] ******
        123.456        999.999        999.999
        999.999        123.456        999.999
        999.999        999.999        156.789
   ****************************************************************

The spiral wavelength :math:`\lambda = 4\pi D_{\alpha\beta} / D_0^{\alpha\beta}` predicts intrinsic helical pitch.

Random Alloy Output
~~~~~~~~~~~~~~~~~~~~

For random alloy simulations (``do_ralloy = 1``), site-resolved properties are calculated:

.. code-block:: text

   *** Supercell averaged properties *******************************
   *****************************************************************
   Exchange stiffness:      12.345  pJ/m
   Spin wave  stiffness:   123.456  meVÅ^2
   Tc-MFA:                 1234.56  K
   J0 VECTOR [meV]:
        1      123.456
        2      234.567
   
   *** Site resolved properties *************************************
    atom       A_xc[pJ/m] D_xc[meVÅ^2] Tc-MFA[K]
        1        12.123     121.234     1200.12
        2        12.567     125.678     1268.45
      ...

This shows how local chemical environment affects stiffness and :math:`T_c`.


J0 Matrix Printing
------------------

When ``prn_J0_matrix = 'Y'``, the full exchange matrix is printed:

.. code-block:: text

   **************** J0 MATRIX [meV] *******************************
        1      1     123.456
        1      2      12.345
        2      1      12.345
        2      2     234.567
   ****************************************************************

Format: ``sublattice_i  sublattice_j  J0_value[meV]``

Terminal Output
~~~~~~~~~~~~~~~

During calculation, progress is reported:

.. code-block:: text

   Calculate the ferromagnetic stiffness done
   Calculate the DMI spiralization done
   Calculate the ferromagnetic stiffness of random alloy done
   Tc-MFA from stiffness :     1234.5 K


**Keywords Reference Table**

.. list-table:: Stiffness Calculation Keywords
   :widths: 25 12 12 50
   :header-rows: 1

   * - Keyword
     - Type
     - Default
     - Description
   * - ``do_stiffness``
     - char
     - 'N'
     - Enable stiffness calculation: 'N'=off, 'Y'=on
   * - ``do_dm_stiffness``
     - char
     - 'N'
     - Enable DMI spiralization: 'N'=off, 'Y'=on
   * - ``prn_J0_matrix``
     - char
     - 'N'
     - Print full J0 exchange matrix: 'N'=off, 'Y'=on
   * - ``eta_min``
     - int
     - 0
     - Minimum eta index for convergence extrapolation
   * - ``eta_max``
     - int
     - 0
     - Maximum eta index for convergence extrapolation

Examples
--------

Example 1: BCC Iron - Basic Stiffness Calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculate exchange stiffness and spin-wave stiffness for BCC Fe:

**Input File (inpsd.dat):**

.. code-block:: fortran

   ! System specification
   simid  bccFe
   cell   2.86 2.86 2.86
   BC     P P P
   ncell  10 10 10

   ! Moments
   Natom  1
   NA     1
   momfile momfile.dat

   ! Exchange
   exchange jfile.dat

   ! Stiffness calculation
   do_stiffness Y
   eta_min 12
   eta_max 18

**Moment File (momfile.dat):**

.. code-block:: fortran

   1  0.0  0.0  1.0  2.2

**Exchange File (jfile.dat):**

.. code-block:: fortran

   1  1   1   1.0  0.0  0.0   21.5
   1  1   2   2.0  0.0  0.0   -2.1
   ...

**Expected Output (micro_ASD.bccFe.out):**

.. code-block:: text

   Exchange stiffness at T=0 K:    16.234  pJ/m
   Spin wave stiffness from Jijs:   281.45 ±   0.34  meVÅ^2
   Tc-MFA from stiffness :      1043.2 K

   EXCHANGE STIFFNESS MATRIX [pJ/m]
        16.234        0.000        0.000
         0.000       16.234        0.000
         0.000        0.000       16.234

**Interpretation:** Isotropic stiffness (cubic symmetry), typical Fe values.

Example 2: FePt with DMI Spiralization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculate both exchange and DMI stiffness for a chiral ferromagnet:

**Input File:**

.. code-block:: fortran

   simid  FePt_DMI
   cell   3.85 3.85 3.85
   BC     P P P
   ncell  8 8 8

   ! System
   Natom  2
   NA     2
   momfile momfile.dat

   ! Exchange and DMI
   exchange jfile.dat
   dm       dmfile.dat

   ! Stiffness with DMI
   do_stiffness Y
   do_dm_stiffness Y
   eta_min 10
   eta_max 16

**DM File (dmfile.dat):**

.. code-block:: fortran

   1  2   1   0.0  0.0  1.0   0.5  0.0  0.0
   1  2   2   1.0  0.0  0.0   0.0  0.5  0.0
   ...

**Expected Output:**

.. code-block:: text

   Exchange stiffness at T=0 K:    12.567  pJ/m
   
   DMI SPIRALIZATION MATRIX [meVA]
               r_x         r_y         r_z
    D_x       3.456      0.000      0.000
    D_y       0.000      3.456      0.000
    D_z       0.000      0.000      0.000
   
   SPIRAL WAVELENGTH MATRIX (STIFF/DMI) [A]
        145.678        999.999        999.999
        999.999        145.678        999.999
        999.999        999.999        999.999

**Interpretation:** 

- Skyrmion-hosting system with spiral wavelength ~14.6 nm
- Uniaxial DMI (D_z components active)
- Skyrmion diameter proportional to wavelength


Example 3: Random Alloy (FeCo) - Site-Resolved Analysis
--------------------------------------------------------

Analyze local variations in random FeCo alloy:

**Input File:**

.. code-block:: fortran

   simid  FeCo_random
   cell   2.85 2.85 2.85
   BC     P P P
   ncell  4 4 4

   ! Random alloy setup
   Natom  64
   NA     1
   do_ralloy 1
   Nchmax 2
   
   ! Chemical types
   momfile momfile.dat
   exchange jfile.dat
   
   ! Stiffness
   do_stiffness Y
   eta_min 6
   eta_max 12

**Expected Output:**

.. code-block:: text

   *** Supercell averaged properties *******************************
   Exchange stiffness:      14.567  pJ/m
   Spin wave  stiffness:   234.567  meVÅ^2
   Tc-MFA:                 1156.78  K
   
   *** Site resolved properties *************************************
    atom       A_xc[pJ/m] D_xc[meVÅ^2] Tc-MFA[K]
        1        14.123     230.234     1140.12
        2        14.890     238.678     1172.45
        3        13.987     228.123     1135.89
      ...

**Interpretation:**

- ~3% variation in local stiffness due to chemical disorder
- Fe-rich sites: slightly lower stiffness
- Co-rich sites: slightly higher stiffness and :math:`T_c`


Example 4: Anisotropic System (HCP Cobalt)
-------------------------------------------

Calculate stiffness tensor for hexagonal system:

**Input File:**

.. code-block:: fortran

   simid  hcpCo
   cell   2.51 2.51 4.07  ! a, a, c lattice parameters
   BC     P P P
   ncell  8 8 6

   Natom  2
   NA     2
   momfile momfile.dat
   exchange jfile.dat
   anisotropy kfile.dat

   ! Stiffness
   do_stiffness Y
   eta_min 10
   eta_max 18

**Expected Output:**

.. code-block:: text

   Anisotropy density       :  1.234  MJ/m^3
   Anisotropy per atom      :  0.145  meV
   Domain wall width        :  8.67  nm
   
   EXCHANGE STIFFNESS MATRIX [pJ/m]
        18.234        0.000        0.000
         0.000       18.234        0.000
         0.000        0.000       16.123
   
   SPIN WAVE STIFFNESS MATRIX [meVA^2]
        305.678        0.000        0.000
          0.000      305.678        0.000
          0.000        0.000      270.123

**Interpretation:**

- In-plane (x,y) stiffness ~13% higher than c-axis (z)
- Anisotropic domain wall widths depending on propagation direction
- Domain wall width calculated from :math:`\Delta = \sqrt{A/K}` where K is anisotropy density

Physical Interpretations and Applications
------------------------------------------

Domain Wall Width
~~~~~~~~~~~~~~~~~

For systems with uniaxial anisotropy, the Bloch/Néel domain wall width is:

.. math::

   \Delta = \sqrt{\frac{A}{K}}

where :math:`K` is the anisotropy density (MJ/m³). Typical values:

- **Soft magnets** (low K): :math:`\Delta \sim 10\text{-}100` nm
- **Hard magnets** (high K): :math:`\Delta \sim 1\text{-}10` nm

Spin-Wave Dispersion
~~~~~~~~~~~~~~~~~~~~

The spin-wave frequency in the long-wavelength limit follows:

.. math::

   \omega(q) \approx \frac{2D}{\hbar} q^2

Relating measured spin-wave stiffness to atomistic exchange:

- **Inelastic neutron scattering**: Measures :math:`D` experimentally
- **Atomistic validation**: Compare calculated :math:`D` with experiment
- **Anisotropic propagation**: Tensor components give direction-dependent dispersion

Skyrmion Size
~~~~~~~~~~~~~~

For skyrmions stabilized by DMI:

.. math::

   R_{\text{sky}} \propto \frac{D}{D_0}

Typical values:

- **Fe/Ir interfaces**: :math:`D_0 \sim 3` meV·Å, :math:`D \sim 100` meV·Å² → :math:`R \sim 10\text{-}20` nm
- **Bulk MnSi**: :math:`D_0 \sim 1` meV·Å, :math:`D \sim 50` meV·Å² → :math:`R \sim 50\text{-}100` nm

Mean-Field Tc Estimation
~~~~~~~~~~~~~~~~~~~~~~~~~

The printed :math:`T_c^{\text{MFA}}` provides an upper bound for the critical temperature:

.. math::

   T_c^{\text{MFA}} = \frac{2}{3k_B} \lambda_{\text{max}}(J_0)

where :math:`\lambda_{\text{max}}` is the largest eigenvalue of the J0 exchange matrix. Typical accuracy:

- **3D systems**: Overestimates by 20-40%
- **2D systems**: Overestimates by 50-100%
- **Random alloys**: Site-resolved :math:`T_c` shows local ordering tendencies

Comparison with Experiment
---------------------------

**Spin-Wave Stiffness D:**

- **Neutron scattering**: Directly measures D from low-q dispersion
- **Brillouin light scattering**: Probes surface spin waves
- **Typical agreement**: Within 10-20% for first-principles exchange parameters

**Exchange Stiffness A:**

- **Domain wall imaging**: Infer A from domain wall width measurements
- **Micromagnetic simulations**: Use calculated A as input
- **Typical agreement**: Within 20-30% (more indirect measurement)

**DMI Spiralization D0:**

- **Spin spiral wavelength**: Measured via neutron diffraction or LEEM
- **Skyrmion size**: Observed with magnetic imaging (MFM, STXM)
- **Typical agreement**: Within factor of 2 (strong correlation effects)

Performance Considerations
--------------------------

**Computational Cost:**

- **Calculation time**: Negligible (~seconds) for systems up to NA=10 sublattices
- **Memory**: Scales as :math:`O(N_A^2 \times \eta_{\text{max}})` for stiffness matrices
- **Bottleneck**: Eigenvalue decomposition of NA×NA matrices (repeated for each eta)

**Convergence:**

- **Eta range**: Wider range (larger ``eta_max - eta_min``) improves extrapolation accuracy
- **Fitting method**: Rational polynomial typically more stable than LSQ for oscillatory systems
- **Random alloys**: Require smaller eta (faster decay) due to disorder

Troubleshooting
---------------

**Issue: Large discrepancy between rational and LSQ fits**

- **Cause**: Oscillatory behavior in :math:`D(\eta)` due to crystallographic frustration
- **Solution**: Increase ``eta_max`` to capture more oscillations; use rational fit (more robust)

**Issue: Negative stiffness values**

- **Cause**: Dominant antiferromagnetic interactions (system not ferromagnetic)
- **Solution**: Verify ground state is ferromagnetic; stiffness calculation assumes FM ordering

**Issue: Tc-MFA far from experimental values**

- **Expected**: MFA overestimates :math:`T_c` by 20-50% in 3D systems
- **Not a bug**: Use as qualitative trend, not quantitative prediction
- **Monte Carlo**: Run ``do_mc='Y'`` for more accurate :math:`T_c` estimation

**Issue: DMI spiralization gives unrealistic wavelengths**

- **Cause**: Small DMI values lead to large :math:`\lambda = 4\pi D/D_0`
- **Check**: Verify DMI parameters in dmfile; ensure DMI is significant compared to exchange
- **Physical**: Some systems have weak DMI → no stable spirals/skyrmions

**Issue: Site-resolved Tc shows unphysical variations (>50%)**

- **Cause**: Insufficient statistics in random alloy sampling
- **Solution**: Increase supercell size (ncell); average over multiple configurations

Combining with Other Features
------------------------------

Stiffness calculations are typically used alongside:

- **Spin-wave calculations** (``do_ams='Y'``): Compare calculated D with spin-wave analysis
- **GNEB calculations** (``do_gneb='Y'``): Use A to estimate energy barriers
- **Monte Carlo** (``do_mc='Y'``): Validate Tc-MFA against finite-temperature simulations
- **Topology analysis** (``do_tottraj='Y'``): Relate skyrmion size to D/D0 ratio

References
----------

The stiffness calculation implementation in UppASD is based on:

1. **Pajda, M., Kudrnovský, J., Turek, I., Drchal, V., & Bruno, P.**
   "Ab initio calculations of exchange interactions, spin-wave stiffness constants, and Curie temperatures of Fe, Co, and Ni"
   *Physical Review B* **64**, 174402 (2001)
   - Original formalism for spin-wave stiffness from atomistic exchange

2. **Udvardi, L., Szunyogh, L., Palotás, K., & Weinberger, P.**
   "First-principles relativistic study of spin waves in thin magnetic films"
   *Physical Review B* **68**, 104436 (2003)
   - Extension to thin films and surfaces

3. **Eriksson, O., Bergman, A., Bergqvist, L., & Hellsvik, J.**
   "Atomistic Spin Dynamics: Foundations and Applications"
   *Oxford University Press* (2017)
   - Comprehensive treatment of atomistic-to-micromagnetic mapping

4. **Hellsvik, J., Skubic, B., Nordström, L., Sanyal, B., Eriksson, O., Nordblad, P., & Warnicke, P.**
   "Dynamics of diluted magnetic semiconductors from atomistic spin-dynamics simulations: Mn-doped GaAs as a case study"
   *Physical Review B* **78**, 144419 (2008)
   - Random alloy stiffness calculations

5. **Beg, M., Pepper, R. A., & Fangohr, H.**
   "User interfaces for computational science: A domain specific language for OOMMF embedded in Python"
   *AIP Advances* **7**, 056025 (2017)
   - Micromagnetic exchange constant usage and conventions
