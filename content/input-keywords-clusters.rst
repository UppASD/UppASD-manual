inpsd.dat keywords: cluster analysis
==============================================================

Cluster analysis keywords for grouping spins or atoms into clusters.
--------------------------------------------------------------------------------

The UppASD code can perform simple cluster analyses for post-processing and
coarse-grained measurements. This page documents the input keywords used to
control cluster detection and reporting.

cluster_method
   Method used to define clusters. Typical values: ``distance`` (geometric cutoff), ``bond`` (bond connectivity), ``label`` (read clusters from file). Default: ``distance``.

cluster_cutoff
   Distance cutoff (in Å) used when ``cluster_method = distance``. Two sites belong to the same cluster if their separation is less than this value. Default: **3.0**.

cluster_file
   Path to a file containing pre-defined cluster assignments (used when ``cluster_method = label``). Format: one integer cluster id per atom, matching the system indexing used in the input files.

cluster_min_size
   Minimum number of atoms/spins for a cluster to be reported. Clusters smaller than this value are ignored. Default: **2**.

prn_clusters
   Print cluster statistics to output (Y=yes, N=no). If enabled, the code writes cluster sizes, centroids, and a list of member indices. Default: **N**.

save_clusters
   Save cluster assignments to disk (Y=yes, N=no). When enabled, writes ``clusters_<simid>.dat`` with one cluster id per line. Default: **N**.

See also
--------

- :doc:`input-keywords-system` (system setup and atom indices)
- :doc:`input-keywords-observables` (writing per-site outputs)

inpsd.dat keywords: cluster embedding
======================================

Parameters for embedding impurity clusters within a host system
---------------------------------------------------------------

The UppASD code provides functionality to **embed magnetic impurity clusters**
within a host system, allowing users to study local properties of defects, interfaces,
and heterostructures without breaking full system symmetries. The cluster method enables
different magnetic moments, exchange interactions, Dzyaloshinskii-Moriya (DM) interactions,
anisotropies, and gyromagnetic ratios inside the cluster compared to the host
[Chico2014]_.

--------------------------------------------------
Cluster embedding: concept and use cases
--------------------------------------------------

The cluster embedding approach is particularly useful for:

- **Impurity atoms or clusters**: Studying magnetic defects embedded in a magnetic host
- **Interface magnetism**: Analyzing magnetic properties at heterostructure interfaces
- **Alloyed systems**: Modeling chemical disorder with distinct magnetic properties
- **Magnetic inhomogeneities**: Systems where different regions have fundamentally different exchange coupling
- **Multi-phase structures**: Interfaces between ferromagnetic and antiferromagnetic regions

The method works by:

1. Defining a **cluster unit cell** with its own lattice vectors and atomic basis
2. Specifying the number of **repetitions** of this unit cell in each direction
3. Providing separate Hamiltonian parameters (exchange, DM, anisotropy) for cluster atoms
4. Mapping cluster atoms to specific positions in the host supercell
5. Overwriting host Hamiltonian with cluster values for mapped atoms

This allows arbitrary spatial distribution of impurities while treating them with
a potentially very different magnetic structure than the host.

--------------------------------------------------
Cluster geometry and structure
--------------------------------------------------

The cluster is defined by:

- **Unit cell**: A cluster basis of :math:`N_A^{(c)}` atoms and lattice vectors :math:`\mathbf{C}_1^{(c)}, \mathbf{C}_2^{(c)}, \mathbf{C}_3^{(c)}`
- **Supercell**: Repetitions :math:`N_1^{(c)} \times N_2^{(c)} \times N_3^{(c)}` of the cluster unit cell
- **Total atoms**: :math:`N_{\text{atom}}^{(c)} = N_A^{(c)} \times N_1^{(c)} \times N_2^{(c)} \times N_3^{(c)}`
- **Mapping**: Each cluster atom is mapped to a specific position in the host via ``index_clus``

The cluster lattice is **independent** of the host lattice, allowing flexible positioning
and orientation of the cluster within the host supercell.

--------------------------------------------------
Input file parameters for cluster embedding
--------------------------------------------------

.. warning:: The original tabular grid has been replaced with a definition-style keyword list to
   avoid fragile grid-table parsing errors. See the keywords below for parameter descriptions.

do_cluster
   Enable cluster embedding functionality (Y=yes, N=no). When enabled, the code reads separate cluster data files and maps cluster atoms to the host system. Default: **N**.

posfile_clus
   File containing cluster unit cell atomic positions and lattice vectors. Format depends on ``posfiletype``: Cartesian (C) or direct/fractional (D). Default: **posfile_clus**.

momfile_clus
   File containing cluster magnetic moments (magnitude, direction, Landé if present). One line per atom type in the cluster unit cell. Default: **momfile_clus**.

jfile_clus
   Exchange coupling files for the cluster (same format as host ``jfile``). May be multiple files for different shells. Default: **jfile_clus**.

dmfile_clus
   Dzyaloshinskii-Moriya (DM) vectors for the cluster (site1 site2 [chem_i chem_j] DM_x DM_y DM_z). Only read if ``do_dm=1``. Default: **dmfile_clus**.

safile_clus
   Scalar anisotropic exchange (DM-like) for cluster (same format as ``dmfile_clus``). Default: **safile_clus**.

kfile_clus
   Magnetic anisotropy data for cluster atoms. Format: site type K_1 K_2 K_3 K_4 K_5 K_6. Read if ``do_anisotropy_clus=1``. Default: **kfile_clus**.

N1_clus, N2_clus, N3_clus
   Number of cluster unit cell repetitions in x, y, z directions (defines the cluster supercell). Defaults: **1**.

NA_clus
   Number of atoms in cluster unit cell (must match positions in ``posfile_clus``). Default: **1**.

C1_clus, C2_clus, C3_clus
   Cluster unit cell lattice vectors (Cartesian). Example: ``C1_clus 3.0 0.0 0.0``.

bc_clus
   Boundary conditions for cluster: **P**=periodic, **O**=open. Default: **P**.

do_anisotropy_clus
   Read anisotropy data for cluster atoms from ``kfile_clus`` (1=yes, 0=no). Default: **0**.

mult_axis_clus
   Multiple anisotropy axes for cluster atoms (Y=yes, N=no). Default: **N**.

random_anis_clus
   Randomize anisotropy directions in cluster (Y=yes, N=no). Default: **N**.

random_anis_den_clus
   Density of random anisotropy (0.0-1.0). Used if ``random_anis_clus=Y``. Default: **0.0**.

do_fixed_mom
   Fixed moment calculation (Y=yes, N=no). When enabled, atoms marked with ``fixed_flag=1`` in restart files remain frozen. Default: **N**.

--------------------------------------------------
Cluster position file format
--------------------------------------------------

**File: ``posfile_clus``** (Cartesian coordinates, ``posfiletype=C``):

.. code-block:: text

   <site> <type> <chem_type> <concentration> <x> <y> <z>

where:
  - ``site``: Atom number in the cluster unit cell (1 to NA_clus)
  - ``type``: Atom type (integer; types are renumbered: cluster types = host_types + NT_host)
  - ``chem_type``: Chemical species (for random alloys)
  - ``concentration``: Occupancy probability (0.0-1.0)
  - ``x, y, z``: Cartesian coordinates (Ångströms)

**File: ``posfile_clus``** (Direct/fractional coordinates, ``posfiletype=D``):

.. code-block:: text

   <site> <type> <chem_type> <concentration> <q1> <q2> <q3>

where :math:`\mathbf{r} = q_1 \mathbf{C}_1^{(c)} + q_2 \mathbf{C}_2^{(c)} + q_3 \mathbf{C}_3^{(c)}`.

**File: ``momfile_clus``** (magnetic moments):

.. code-block:: text

   <site> <chem_type> <magnitude> <e_x> <e_y> <e_z> [<Landeg>] [<induced_flag>]

where:
  - ``site``: Atom number in cluster unit cell
  - ``chem_type``: Chemical species
  - ``magnitude``: Magnetic moment magnitude :math:`|\mathbf{m}|` (in :math:`\mu_B`)
  - ``e_x, e_y, e_z``: Unit direction vector (will be normalized automatically)
  - ``Landeg``: Gyromagnetic ratio (optional, only if ``set_landeg=1``)
  - ``induced_flag``: Whether moment is induced (1) or permanent (0) (optional, for LSF)

--------------------------------------------------
Hamiltonian files for clusters
--------------------------------------------------

**File: ``jfile_clus``** (Exchange couplings):

For non-random alloys (``do_ralloy=0``):

.. code-block:: text

   <site_i> <site_j> [<distance_x> <distance_y> <distance_z>] <J>

For random alloys (``do_ralloy=1``):

.. code-block:: text

   <site_i> <site_j> <chem_i> <chem_j> [<distance_x> <distance_y> <distance_z>] <J>

When ``do_jtensor=1``, the coupling is a :math:`3 \times 3` tensor:

.. code-block:: text

   <site_i> <site_j> <J_xx> <J_xy> <J_xz> <J_yx> <J_yy> <J_yz> <J_zx> <J_zy> <J_zz>

**File: ``dmfile_clus``** (DM vectors):

.. code-block:: text

   <site_i> <site_j> [<chem_i> <chem_j>] <DM_x> <DM_y> <DM_z>

The DM interaction contributes:

.. math::

   \mathcal{H}_{\text{DM}} = \mathbf{D}_{ij} \cdot (\mathbf{S}_i \times \mathbf{S}_j)

**File: ``kfile_clus``** (Anisotropy):

.. code-block:: text

   <site> <aniso_type> <K_1> <K_2> <K_3> <K_4> <K_5> <K_6>

where ``aniso_type`` specifies the anisotropy form:
  - **0**: No anisotropy
  - **1**: Uniaxial (easy axis): :math:`\mathcal{H} = -K (S_z)^2`
  - **2**: Cubic: :math:`\mathcal{H} = -K (S_x^4 + S_y^4 + S_z^4)`

The six parameters are axis components and coupling constants.

--------------------------------------------------
Cluster atom index mapping
--------------------------------------------------

The code maintains a **mapping array** ``index_clus`` that connects each cluster atom
to its position in the host system. This is determined by:

1. **Cluster generation**: A regular lattice of cluster unit cells is created
2. **Host overlap**: The generated cluster positions are compared to host atom positions
3. **Nearest-neighbor matching**: Each cluster atom is mapped to the nearest host atom within tolerance
4. **Index array**: ``index_clus(i_clus) = i_host`` stores the mapping

Once mapping is complete:

- Host atom types at mapped positions are renumbered: :math:`\text{type}_{\text{host}} + N_T^{\text{host}}`
- Host Hamiltonian values (exchange, DM, anisotropy) are **overwritten** with cluster values
- Unmapped cluster atoms expand the host system size

--------------------------------------------------
Interaction handling: overwrite strategy
--------------------------------------------------

When cluster embedding is active, the code replaces host Hamiltonian parameters with cluster
values for all atoms identified as part of the cluster:

1. **Exchange couplings**: For each cluster atom and its neighbors, the exchange tensor is
   replaced with the cluster-specified value. New interactions not in the host are added.

2. **DM interactions**: Cluster DM vectors completely replace host DM vectors for mapped pairs.

3. **Anisotropy**: Cluster anisotropy constants and axes replace host values.

4. **Gyromagnetic ratio**: Cluster-specific ``Landeg`` overrides host value for each atom.

This allows **sharp interfaces** between host and cluster with potentially very different
magnetic properties.

--------------------------------------------------
Fixed moment calculations
--------------------------------------------------

The **fixed moment approach** allows selective freezing of magnetic moments:

do_fixed_mom
   Enable fixed moment mode (Y=yes, N=no). When active, uses restart files that include a ``fixed_flag`` for each atom. Atoms with ``fixed_flag=1`` are frozen while others evolve normally. Default: **N**.

In fixed moment mode:

- Atoms with ``fixed_flag=1`` have their magnetic moments frozen (no dynamics)
- Atoms with ``fixed_flag=0`` evolve normally
- Useful for **studying excitations in a static background**: e.g., magnons above a frozen skyrmion
- Compatible with cluster embedding: entire cluster can be frozen or individual atoms

The restart/configuration file format includes this flag:

.. code-block:: text

   <step_number> (header)
   <iatom> <itype> <mmom> <e_x> <e_y> <e_z> <fixed_flag> [<induced_flag>]
   ...

When ``initmag=4`` (restart-based initialization), the fixed_flag values are read from file.
Otherwise, they are determined from ``inp_fixed_mom_flag`` in the input structure.

--------------------------------------------------
Example: Impurity cluster in a host
--------------------------------------------------

Setup for embedding a 2×2×2 cluster of 4 Fe atoms at positions (0,0,0) and fractional
translations within a larger host system:

**inpsd.dat settings:**

.. code-block:: text

   ! Cluster embedding
   do_cluster       Y
   posfile_clus     cluster_pos.txt
   momfile_clus     cluster_mom.txt
   jfile_clus       cluster_J.txt
   N1_clus          2
   N2_clus          2
   N3_clus          2
   NA_clus          1
   C1_clus          3.0  0.0  0.0
   C2_clus          0.0  3.0  0.0
   C3_clus          0.0  0.0  3.0

**File: cluster_pos.txt** (1 atom in unit cell):

.. code-block:: text

   1 1 1 1.0 0.0 0.0 0.0

**File: cluster_mom.txt**:

.. code-block:: text

   1 1 2.5 0.0 0.0 1.0

(2.5 μ_B pointing along z)

**File: cluster_J.txt** (nearest-neighbor in-cluster coupling):

.. code-block:: text

   1 1 1 0 1.0 0.0 0.0 -5.0
   1 1 1 0 0.0 1.0 0.0 -5.0
   1 1 1 0 0.0 0.0 1.0 -5.0

This creates a 2×2×2 cluster of 8 Fe atoms with strong FM coupling within the cluster,
which is embedded at the location where host atoms are found.

--------------------------------------------------
Example: Interface magnetism with fixed moment
--------------------------------------------------

Studying interface excitations between ferromagnetic (FM) and antiferromagnetic (AFM) clusters:

**inpsd.dat:**

.. code-block:: text

   do_cluster       Y
   do_fixed_mom     Y
   posfile_clus     interface_pos.txt
   momfile_clus     interface_mom.txt
   jfile_clus       interface_J.txt
   N1_clus          1
   N2_clus          1
   N3_clus          4
   NA_clus          2
   C1_clus          3.0  0.0  0.0
   C2_clus          0.0  3.0  0.0
   C3_clus          0.0  0.0  3.0

**File: interface_pos.txt** (2 atoms per unit cell):

.. code-block:: text

   1 1 1 1.0 0.0 0.0 0.0
   2 2 1 1.0 0.5 0.5 0.0

(Two sublattices at (0,0,0) and (0.5,0.5,0) in fractional coordinates)

**File: interface_mom.txt**:

.. code-block:: text

   1 1 2.0 0.0 0.0  1.0 1
   2 1 2.0 0.0 0.0 -1.0 1

(Both sublattices frozen: fixed_flag=1)

**File: interface_J.txt**:

.. code-block:: text

   1 2 1 0 1 -8.0
   1 2 1 0 0 -8.0
   1 2 1 0 0 -8.0
   2 1 1 0 -1 -8.0
   2 1 1 0 0 -8.0
   2 1 1 0 0 -8.0

This AFM alignment (opposite moments, AFM coupling) remains static while allowing other
atoms in the system to evolve.

--------------------------------------------------
Technical notes and limitations
--------------------------------------------------

**Symmetry breaking:**
  The cluster method explicitly breaks crystal and magnetic symmetries. Consequently:
  - The ``do_reduced`` option is not compatible with cluster embedding
  - Full supercell treatment is required
  - Simulation times may be longer due to reduced symmetry

**Boundary effects:**
  - Cluster atoms at boundaries may have fewer neighbors than bulk atoms
  - The mapping tolerance (0.005 Ångströms) may need adjustment for mismatched lattices
  - Atoms not successfully mapped expand the host system

**Performance:**
  - Cluster embedding adds minimal computational overhead during simulation
  - Most cost is in initial setup and neighbor list generation
  - For very small clusters or large host cells, overhead is negligible

**LSF (Local Spin Frame) compatibility:**
  - When ``do_lsf=Y``, cluster momfiles must include configuration numbers as first column
  - Multiple configurations can be specified for LSF-based sampling

**Random alloys:**
  - Cluster atoms can have chemical disorder (``do_ralloy=1``)
  - Chemical types and concentrations are specified in ``posfile_clus``
  - Different exchange couplings can apply to each chemical combination

--------------------------------------------------
References and further reading
--------------------------------------------------

For the cluster embedding methodology and applications:

.. [Chico2014] J. Chico, L. E. F. Foa Torres, R. H. Barco, O. Eriksson, and A. Bergman, "Unraveling novel quantum critical behavior in a two-dimensional Heisenberg antiferromagnet," *Phys. Rev. B* **90**, 064425 (2014).

See also
--------

- :doc:`input-keywords-system` (host system definition)
- :doc:`input-core-files` (core file formats: ``posfile``, ``momfile``, etc.)
- :doc:`input-keywords-hamiltonian` (exchange coupling definitions)
- :doc:`input-keywords-observables` (site-resolved outputs and per-site writing)