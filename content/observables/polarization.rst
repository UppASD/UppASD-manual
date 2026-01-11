Polarization and chirality measurements
========================================

.. _input-polarization:

Overview
--------

UppASD computes spin-driven **ferroelectric polarization** and **local chirality**
from the inverse Dzyaloshinskii-Moriya (DM) mechanism and spin-current models.
The polarization :math:`\mathbf{P}` is proportional to the vector spin chirality
:math:`\sum_{i,j}\hat{\mathbf{e}}_{ij}\times(\mathbf{m}_i\times\mathbf{m}_j)`,
where :math:`\hat{\mathbf{e}}_{ij}` is the normalized bond vector connecting
neighboring spins. This quantity is crucial for multiferroic systems where
magnetic order induces electric polarization.

Physical mechanism and formula
------------------------------

The spin-induced polarization originates from relativistic spin-orbit coupling
and inverse DM interactions. For a pair of spins :math:`\mathbf{m}_i` and
:math:`\mathbf{m}_j` separated by :math:`\mathbf{r}_{ij}=\mathbf{r}_j-\mathbf{r}_i`,
the local contribution is:

.. math::

   \mathbf{p}_{ij} \propto \hat{\mathbf{e}}_{ij} \times (\mathbf{m}_i \times \mathbf{m}_j),

where :math:`\hat{\mathbf{e}}_{ij} = \mathbf{r}_{ij}/|\mathbf{r}_{ij}|`. The
total polarization is the sum over all bonds (typically nearest neighbors):

.. math::

   \mathbf{P} = \gamma \sum_{i,j} \hat{\mathbf{e}}_{ij} \times (\mathbf{m}_i \times \mathbf{m}_j),

with :math:`\gamma` a material-dependent coupling constant (absorbed into the
definition here). This form captures cycloidal and helical spin spirals that
break inversion symmetry.

**Key points:**

- :math:`\mathbf{P}` is zero for collinear ferromagnets and antiferromagnets
- Non-zero for non-collinear structures (spirals, skyrmions, conical phases)
- Sign of :math:`\mathbf{P}` flips under time reversal or spatial inversion
- Magnitude scales with spin canting angle and DM coupling strength

.. tip::

   **Input keywords:** For a complete list of all polarization and chirality
   measurement parameters (``do_pol``, ``max_pol_nn``, ``pol_step``, ``do_chir``, etc.),
   see the comprehensive reference in :doc:`../input/observables`.

Mathematical derivation
-----------------------

Starting from the cross product identity, we expand the scalar triple product:

.. math::

   \mathbf{p}_{ij} = \hat{\mathbf{e}}_{ij} \times (\mathbf{m}_i \times \mathbf{m}_j).

Using the vector identity :math:`\mathbf{a}\times(\mathbf{b}\times\mathbf{c})=
\mathbf{b}(\mathbf{a}\cdot\mathbf{c}) - \mathbf{c}(\mathbf{a}\cdot\mathbf{b})`,
we can write:

.. math::

   \mathbf{m}_i \times \mathbf{m}_j = \begin{vmatrix}
   \hat{\mathbf{x}} & \hat{\mathbf{y}} & \hat{\mathbf{z}} \\
   m_i^x & m_i^y & m_i^z \\
   m_j^x & m_j^y & m_j^z
   \end{vmatrix},

and then form the triple product with :math:`\hat{\mathbf{e}}_{ij}`. In
component form (as implemented in ``buffer_pol``):

.. math::

   t_x &= m_i^y m_j^z - m_i^z m_j^y, \\
   t_y &= m_i^z m_j^x - m_i^x m_j^z, \\
   t_z &= m_i^x m_j^y - m_i^y m_j^x,

   p_x &= e_{ij}^y t_z - e_{ij}^z t_y, \\
   p_y &= e_{ij}^z t_x - e_{ij}^x t_z, \\
   p_z &= e_{ij}^x t_y - e_{ij}^y t_x.

The code sums these contributions over all neighbor pairs (up to ``max_pol_nn``)
and all atoms.

Local chirality: scalar vs vector
---------------------------------

**Vector chirality** (or spin chirality) is the antisymmetric part of the spin
correlation between neighbors:

.. math::

   \boldsymbol{\chi}_{ij} = \mathbf{m}_i \times \mathbf{m}_j.

This measures the local rotation sense of the magnetization. A finite
:math:`\boldsymbol{\chi}` indicates broken mirror symmetry.

**Scalar chirality** is the projection of :math:`\boldsymbol{\chi}` onto a third
spin direction (triangular plaquette). In UppASD, the local chirality
(``s_cross_s``) is the bond-wise vector :math:`\mathbf{m}_i\times\mathbf{m}_j`,
summed over neighbors of each atom.

Implemented in ``measure_chirality(Natom, Mensemble, emomM, max_no_neigh, nlist, nlistsize)``.

**Physical interpretation:**

- Non-zero :math:`\boldsymbol{\chi}` signals non-coplanar spin textures
- Related to topological Hall effect and Berry curvature in electronic structure
- For triangular lattices, the scalar chirality over a plaquette
  :math:`\mathbf{m}_1\cdot(\mathbf{m}_2\times\mathbf{m}_3)` is directly tied to
  the topological charge density

Local polarization: site-resolved
----------------------------------

Local polarization :math:`\mathbf{p}_i` at site :math:`i` is the sum of
polarization contributions from all bonds :math:`ij` connected to :math:`i`:

.. math::

   \mathbf{p}_i = \sum_{j\in\text{neighbors}(i)} \hat{\mathbf{e}}_{ij} \times (\mathbf{m}_i \times \mathbf{m}_j).

This provides a real-space map of polarization density, useful for visualizing
domain walls, skyrmion cores, or spiral modulations.

Implemented in ``measure_local_pol(Natom, Mensemble, emomM, max_no_neigh, nlist, nlistsize)``.

**Applications:**

- Identify regions of strong ferroelectric response
- Correlate with topological structures (skyrmions induce radial :math:`\mathbf{p}`)
- Analyze magnetoelectric coupling at interfaces or defects

Neighbor list and cutoff considerations
---------------------------------------

Polarization calculations sum over neighbors in the exchange neighbor list up to
index ``max_pol_nn``. This cutoff controls which bonds contribute:

- **Short cutoff** (nearest neighbors only): Captures local spiral twists but
  may miss long-range correlations
- **Long cutoff** (many shells): Includes distant pairs, increasing computational
  cost and potentially introducing spurious contributions

**Recommendation:** Set ``max_pol_nn`` to the number of nearest-neighbor shells
relevant for the magnetic structure (typically 6-12 for 2D lattices, 1-2 for
simple spirals).

**Important flag:** Use ``do_sortcoup N`` to prevent reordering of the neighbor
list, which would break the correspondence between geometric bonds
(:math:`\hat{\mathbf{e}}_{ij}`) and exchange neighbors. Without this, the
polarization formula is undefined.

Keyword reference
-----------------

do_pol
   Enable polarization measurement (Y/*N*). Activates calculation of total
   :math:`\mathbf{P}` and writes ``polarization.<simid>.out``.

pol_step
   Number of time steps between polarization samples (integer, default 100).

pol_buff
   Number of samples to buffer before writing to file (integer, default 10).

do_loc_pol
   Local (site-resolved) polarization (Y/*N*). Outputs :math:`\mathbf{p}_i` to
   ``prestart.<simid>.out`` for visualization.

do_chiral
   Chirality measurement (Y/*N*). Computes and outputs local vector chirality
   :math:`\mathbf{m}_i\times\mathbf{m}_j` summed over neighbors. Writes
   ``crestart.<simid>.out``.

max_pol_nn
   Maximum number of neighbors for polarization sum (integer, default 6). Must
   not exceed the exchange neighbor list size. Larger values include more shells
   but increase computational cost.

Output files and analysis
-------------------------

polarization.<simid>.out
   Time-series of ensemble-averaged polarization. Format:

   ::

      step   <P>_x   <P>_y   <P>_z   <P>   sigma(P)

   where :math:`\langle P\rangle = \sqrt{\langle P_x\rangle^2 + \langle P_y\rangle^2
   + \langle P_z\rangle^2}` and :math:`\sigma(P)` is the standard deviation over
   ensembles.

prestart.<simid>.out
   Site-resolved local polarization :math:`\mathbf{p}_i` at the final time step.
   Format per atom: ``step``, ``iatom``, :math:`p_x`, :math:`p_y`, :math:`p_z`,
   :math:`|\mathbf{p}_i|`.

crestart.<simid>.out
   Site-resolved local chirality :math:`\boldsymbol{\chi}_i` (sum of
   :math:`\mathbf{m}_i\times\mathbf{m}_j` over neighbors). Format identical to
   ``prestart.<simid>.out``.

Example setup
-------------

**Basic polarization measurement:**

.. code-block:: text

   ! Ensure neighbor list is not reordered
   do_sortcoup  N

   ! Polarization settings
   do_pol       Y
   pol_step     50
   pol_buff     20
   max_pol_nn   12         ! Include up to 2nd neighbors

**Local polarization and chirality:**

.. code-block:: text

   do_sortcoup  N
   do_pol       Y
   do_loc_pol   Y          ! Output site-resolved polarization
   do_chiral    Y          ! Output site-resolved chirality
   pol_step     100
   max_pol_nn   6

Physical examples and use cases
-------------------------------

**Cycloidal spin spirals:**

In a cycloidal spiral with :math:`\mathbf{q}\parallel\hat{\mathbf{z}}` and
rotation in the :math:`xy`-plane, the polarization is perpendicular to both
:math:`\mathbf{q}` and the spiral plane:

.. math::

   \mathbf{m}(z) = \cos(qz)\,\hat{\mathbf{x}} + \sin(qz)\,\hat{\mathbf{y}}
   \quad\Rightarrow\quad \mathbf{P}\parallel\hat{\mathbf{z}}.

**Néel skyrmions:**

Skyrmions with radial magnetization profiles induce radial or tangential
:math:`\mathbf{p}_i` patterns, depending on the Bloch vs Néel character.
Summing over the skyrmion gives a net :math:`\mathbf{P}=0` for isolated
skyrmions (topological cancellation), but non-zero at domain walls or
lattice edges.

**Multiferroic domain walls:**

At a 180° domain wall separating opposite spiral helicities, the local
polarization :math:`\mathbf{p}_i` flips sign, creating a ferroelectric domain
structure co-located with the magnetic wall.

Computational performance and accuracy
--------------------------------------

**Performance:**

- Polarization calculation is :math:`O(N\times N_{\text{neigh}})` per time step
- Negligible overhead compared to energy/force calculations for typical
  ``max_pol_nn`` values (≤12)
- Local polarization (``do_loc_pol``) adds :math:`O(N)` storage but minimal CPU cost

**Numerical accuracy:**

- Results are exact for the discrete lattice (no gradient approximations)
- Sensitive to the quality of the magnetic structure (equilibration, convergence)
- Ensemble averaging (``Mensemble`` > 1) reduces thermal fluctuations in
  :math:`\mathbf{P}`

**Improving results:**

- Use fine lattices to resolve short-wavelength spirals
- Ensure sufficient equilibration before sampling (especially in Monte Carlo)
- Cross-check with analytical models for simple spirals (e.g., single-:math:`q`
  helices should give :math:`|\mathbf{P}|\propto|\sin\theta|` for canting
  angle :math:`\theta`)

Related keywords and cross-references
-------------------------------------

- ``do_sortcoup``: Must be ``N`` for polarization (see :doc:`../input/hamiltonian`)
- ``max_no_neigh``: Defines exchange neighbor list size
- ``do_skyno``, ``do_skyno_den``: Topological measurements often correlate with
  local polarization (see :doc:`topology`)
- Spin dynamics: polarization dynamics can track magnon-driven ferroelectric
  switching

References
==========

See the centralized :doc:`../references` for full bibliographic entries:

- [Katsura2005]_ - Spin current and magnetoelectric effect in noncollinear magnets
- [Mostovoy2006]_ - Ferroelectricity in spiral magnets
- [Tokura2010]_ - Multiferroics with spiral spin orders
