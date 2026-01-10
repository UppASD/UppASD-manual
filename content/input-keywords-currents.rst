Spin torques
============

Parameters for current-induced spin torques
--------------------------------------------

The UppASD code supports several forms of current-induced torques acting on
atomic magnetic moments, including [RalphStiles2008]_, [Manchon2019]_:

- :ref:`Volume spin-transfer torques (STT) of Zhang–Li type <torque-zhang-li>`
- :ref:`Slonczewski-type STT from a fixed spin polarization <torque-slonczewski>`
- :ref:`Spin–orbit torques (SOT) in a general phenomenological form <torque-sot>`
- :ref:`Spin Hall effect (SHE) induced torques in thin-film geometries <torque-she>`

All torque contributions are implemented in **effective-field form** and enter
the Landau–Lifshitz–Gilbert (LLG) equation through the total effective field
acting on each magnetic moment [Meo2023]_.

.. important::

   Since the spin and orbital torques depend on current densities, it is important to provide the correct lattice parameter for the system.
   See the keyword ``alat`` in :doc:`system keywords <input-keywords-system>`.

-------------------------------------------------
Canonical current density representation
-------------------------------------------------

Internally, all current-driven torque mechanisms are expressed in terms of a
**canonical charge-current density vector**

.. math::

   \mathbf{j} \equiv \text{jdens}

with units of **A/m²**.

This vector defines both:

- the **magnitude** of the applied current, :math:`|\mathbf{j}|`, and
- its **direction**, :math:`\hat{\mathbf{j}} = \mathbf{j}/|\mathbf{j}|`.

Optionally, a site-resolved current density
:math:`\mathbf{j}_i` may be used for spatially inhomogeneous currents.

-------------------------------------------------
Legacy current input (jvec)
-------------------------------------------------

For backward compatibility, the legacy input keyword ``jvec`` may be used.
If ``jdens`` is not provided but ``jvec`` (or ``jvecfile``) is present, the code
automatically converts ``jvec`` into a physical current density using the
internal conversion factor ``stt_dens_conv``.

In this case:

.. math::

   \mathbf{j} = \text{stt_dens_conv} \times \mathbf{j}_{\text{vec}}

A note is printed at startup indicating that legacy input has been used and
showing the resulting current density in physical units.

.. _torque-zhang-li:

Zhang–Li spin-transfer torque (volume STT)
-------------------------------------------------

The Zhang–Li torque is enabled by setting [ZhangLi2004]_

.. tabularcolumns:: |l|l|
 
+---------------+--------------------------------------------------------------+
|  stt          |  Type of spin-transfer torque (A=Zhang–Li, S=Slonczewski,    |
|               |  *N=none*).                                                  |
+---------------+--------------------------------------------------------------+

For ``stt = A``, the magnetization dynamics includes the terms

.. math::

   \frac{\partial \mathbf{m}}{\partial t}
   \supset
   -(\mathbf{u}\cdot\nabla)\mathbf{m}
   + \beta\,\mathbf{m}\times(\mathbf{u}\cdot\nabla)\mathbf{m},

where

.. math::

   \mathbf{u} = \frac{P\mu_B}{e M_s}\,\mathbf{j}.

In the atomistic implementation:

- the **strength** of the torque scales linearly with the current magnitude
- the gradient operator is evaluated along the current direction, and the overall
   prefactor depends on the current magnitude and material parameters.

.. _torque-slonczewski:

Slonczewski spin-transfer torque
--------------------------------------------------

A fixed-polarization spin-transfer torque is enabled by setting [Slonczewski1996]_

.. tabularcolumns:: |l|l|

+---------------+--------------------------------------------------------------+
|  stt          |  Type of spin-transfer torque (S=Slonczewski,                |
|               |  A=Zhang–Li, *N=none*).                                      |
+---------------+--------------------------------------------------------------+

This torque represents spin injection from a reference magnetic layer with a
prescribed polarization direction :math:`\mathbf{p}`.

The effective-field contribution is written in the form

.. math::

   \mathbf{B}_{\mathrm{STT}}
   \propto
   \mathbf{p}
   + \mathbf{m}\times\mathbf{p},

which, when inserted into the LLG equation, generates the conventional
Slonczewski field-like and damping-like spin-transfer torques

.. math::

   \mathbf{T}_{\mathrm{FL}} \propto \mathbf{m}\times\mathbf{p}, \qquad
   \mathbf{T}_{\mathrm{DL}} \propto \mathbf{m}\times(\mathbf{m}\times\mathbf{p}).

The strength of the torque scales linearly with the magnitude of the applied
charge current density :math:`|\mathbf{j}|`.

-------------------------------------------------
Polarization specification
-------------------------------------------------

The fixed polarization direction may be specified as

- ``stt_pol_vec`` – Uniform polarization direction :math:`\mathbf{p}` (Cartesian components).
- ``stt_site_pol`` – Flag to enable site-dependent polarization (Y=yes, *N=no*).
- ``stt_site_file`` – External file containing site-resolved polarization vectors.

If a site-dependent polarization is not specified, the uniform vector
``stt_pol_vec`` is applied to all magnetic atoms.

-------------------------------------------------
Notes
-------------------------------------------------

- Fixed STT does not depend on spatial gradients of the magnetization.
- It does not require spin–orbit coupling.
- It is distinct from both Zhang–Li STT and SHE/SOT.
- The same canonical current density ``jdens`` is used to control its magnitude.

.. _torque-sot:

General spin–orbit torque (SOT)
--------------------------------------------------

The general SOT is enabled by ``do_sot`` (Y=yes, *N=no*) [GambardellaMiron2011]_, [Manchon2019]_, [Haney2013]_.

The SOT is written in effective-field form as

.. math::

   \mathbf{B}_{\mathrm{SOT}}
   =
   -(\tau_{\mathrm{FL}} - \lambda\,\tau_{\mathrm{DL}})\,\mathbf{P}
   -(\tau_{\mathrm{DL}} + \lambda\,\tau_{\mathrm{FL}})
   \left(\mathbf{m}\times\mathbf{P}\right),

where:

- :math:`\mathbf{P}` is the spin polarization direction,
- :math:`\tau_{\mathrm{FL}}` and :math:`\tau_{\mathrm{DL}}` are **effective-field
  amplitudes**,
- :math:`\lambda` is the local Gilbert damping parameter.

The quantities ``sot_field`` and ``sot_damping`` correspond to
:math:`\tau_{\mathrm{FL}}` and :math:`\tau_{\mathrm{DL}}`, respectively.

**Important:**  
``sot_field`` and ``sot_damping`` are **effective fields**, not torques.
They generate field-like and damping-like torques through the LLG equation.

.. tabularcolumns:: |l|l|

+----------------+-------------------------------------------------------------+
|  sot_field     |  Field-like SOT effective field amplitude (internal field   |
|                |  units).                                                    |
+----------------+-------------------------------------------------------------+
|  sot_damping   |  Damping-like SOT effective field amplitude (internal field |
|                |  units).                                                    |
+----------------+-------------------------------------------------------------+

.. _torque-she:

Spin Hall effect (SHE) torque
--------------------------------------------------

The SHE torque is enabled by ``do_she`` (Y=yes, *N=no*) [Liu2012]_, [Manchon2019]_, [Meo2023]_.

The SHE torque uses the canonical current density ``jdens`` and generates a
spin polarization direction :math:`\boldsymbol{\sigma}` according to one of
two modes:

1. **Geometry-based mode (recommended for thin films)**  
   If a surface/interface normal vector ``she_n_vec`` is provided,

   .. math::

      \boldsymbol{\sigma}
      =
      \frac{\hat{\mathbf{n}}\times\hat{\mathbf{j}}}
           {|\hat{\mathbf{n}}\times\hat{\mathbf{j}}|}.

2. **Direct polarization mode**  
   If ``she_n_vec`` is not provided, the direction of ``jdens`` itself is
   interpreted as the spin polarization direction:

   .. math::

      \boldsymbol{\sigma} = \hat{\mathbf{j}}.

A message is printed at startup indicating which mode is used.

The resulting effective field generates both field-like and damping-like
spin–orbit torques through the LLG equation.

-------------------------------------------------
Units and conventions
-------------------------------------------------

- ``jdens``: charge current density, **A/m²**
- ``jvec``: legacy input, converted internally to A/m²
- ``sot_field``, ``sot_damping``: effective fields (same units as other fields
  in UppASD)
- Polarization vectors (``stt_pol_vec``, ``sot_pol_vec``,
  :math:`\boldsymbol{\sigma}`): unit vectors

-------------------------------------------------
Notes on usage
-------------------------------------------------

- The general SOT and SHE torques are **independent mechanisms**; either or both
  may be enabled.
- The SOT amplitudes may be specified phenomenologically without reference to
  current density.
- The SHE torque provides a geometry-based current-to-torque mapping and is
  intended for thin-film bilayer systems.
- For thick samples, care should be taken when using interfacial torque models.

