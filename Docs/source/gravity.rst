.. _ch:gravity:

*******
Gravity
*******


Introduction
============

Gravity Types
--------------------

Castro can incorporate gravity as a constant, monopole approximation,
or a full Poisson solve. To enable gravity in the code, set::

    USE_GRAV = TRUE

in the ``GNUmakefile``, and then turn it on in the inputs file
via ``castro.do_grav`` = 1. If you want to incorporate a point mass
(through ``castro.point_mass``), you must have::

    USE_POINTMASS = TRUE

in the ``GNUmakefile``.

There are currently three options for how gravity is calculated,
controlled by setting ``gravity.gravity_type``. The options are
``ConstantGrav``, ``PoissonGrav``, or ``MonopoleGrav``.
Again, these are only relevant if ``USE_GRAV =
TRUE`` in the ``GNUmakefile`` and ``castro.do_grav`` = 1 in the inputs
file. If both of these are set then the user is required to specify
the gravity type in the inputs file or the program will abort.

.. note:: make sure you have set the ``problem::center[]`` variable
   appropriately for you problem.  This can be done by directly
   setting it in the ``problem_initialize()`` function.


Integration Strategy
--------------------

Castro uses subcycling to integrate levels at different timesteps.
The gravity algorithm needs to respect this to obtain full accuracy.
When self-gravity is computed via a multigrid solve
(``gravity.gravity_type = PoissonGrav``), we correct for this (though
not for other gravity types). At coarse-fine interfaces, the stencil
used in the Laplacian understands the coarse-fine interface and is
different than the stencil used in the interior.

There are two types of solves that we discuss with AMR:

-  *composite solve* : This is a multilevel solve, starting at
   a coarse level (usually level 0) and solving for the potential on
   all levels up to the finest level.

-  *level solve* : This solves for the potential only on
   a particular level. Finer levels are ignored. At coarse-fine
   interfaces, the data from the coarse levels acts as Dirichlet
   boundary conditions for the current-level-solve.

The overall integration strategy is as follows, and is similar to
the discussion in :cite:`castro_I`. Briefly:

-  At the beginning of a simulation, we do a multilevel composite
   solve (if ``gravity.no_composite`` = 0).

   We also do a multilevel composite solve after each regrid.

-  The old-time gravity on the coarse level is defined based on
   this composite solve, but we also do a level solve on the coarse
   level, and use it to determine the difference between the composite
   solve and the level solve, and store that in a MultiFab.

-  After the hydro advance on the coarse level, we do another level
   solve, and use the (level solve - compositive solve) as a lagged
   predictor of how much we need to add back to that level solve to get
   an effective new-time value for phi on the coarse level, and that’s
   what defines the phi used for the new-time gravity

-  Then we do the fine grid timestep(s), each using the same
   strategy

-  At an AMR synchronization step across levels (see Section
   :ref:`sec:amr_synchronization` for a
   description of when these synchronizations occur), if we’re
   choosing to synchronize the gravitational field across levels
   (``gravity.no_sync`` = 0) we then do a solve starting from the coarse
   grid that adjusts for the mismatch between the fine-grid phi and
   the coarse-grid phi, as well as the mismatch between the fine-grid
   density fluxes and the coarse-grid density fluxes, and add the
   resulting sync solve phi to both the coarse and the fine level

   Thus, to within the gravity error tolerance, you get the same final
   result as if you had done a full composite solve at the end of the
   timestep (assuming ``gravity.no_sync`` = 0).

If you do ``gravity.no_composite`` = 1, then you never do a full
multilevel solve, and the gravity on any level is defined only by the
solve on that level. The only time this would be appropriate is if
the fine level(s) cover essentially all of the mass on the grid for
all time.

Controls
--------

-  For the full Poisson solver
   (``gravity.gravity_type`` = ``PoissonGrav``), the behavior
   of the full Poisson solve / multigrid solver is controlled by
   ``gravity.no_sync`` and ``gravity.no_composite``.

-  For isolated boundary conditions, and when
   ``gravity.gravity_type`` = ``PoissonGrav``, the parameters
   ``gravity.max_multipole_order`` and
   ``gravity.direct_sum_bcs`` control the accuracy of
   the Dirichlet boundary conditions. These are described in
   Section `2.3.2 <#sec-poisson-3d-bcs>`__.

-  For ``MonopoleGrav``, in 1D we must have ``coord_sys`` = 2, and in
   2D we must have ``coord_sys`` = 1.

The following parameters apply to gravity
solves:

-  ``gravity.gravity_type`` : how should we calculate gravity?
   Can be ``ConstantGrav``, ``PoissonGrav``, or ``MonopoleGrav``

-  ``gravity.const_grav`` : if ``gravity.gravity_type`` =
   ``ConstantGrav``, set the value of constant gravity (default: 0.0)

-  ``gravity.no_sync`` : ``gravity.gravity_type`` =
   ``PoissonGrav``, do we perform the “sync solve"? (0 or 1; default: 0)

-  ``gravity.no_composite`` : if gravity.gravity_type
   = ``PoissonGrav``, whether to perform a composite solve (0 or 1;
   default: 0)

-  ``gravity.max_solve_level`` : maximum level to solve
   for :math:`\phi` and :math:`\mathbf{g}`; above this level, interpolate from
   below (default: ``MAX_LEV``-1)

-  ``gravity.abs_tol`` : if ``gravity.gravity_type`` = ``PoissonGrav``,
   this is the absolute tolerance for the Poisson solve. You can
   specify a single value for this tolerance (or do nothing, and get a
   reasonable default value), and then the absolute tolerance used by
   the multigrid solve is ``abs_tol`` :math:`\times 4\pi G\,
   \rho_{\text{max}}` where :math:`\rho_{\text{max}}` is the maximum
   value of the density on the domain. On fine levels, this absolute
   tolerance is multiplied by :math:`(\mathtt{ref\_ratio})^2` to account
   for the change in the absolute scale of the Laplacian operator. You
   can also specify an array of values for ``abs_tol``, one for each
   possible level in the simulation, and then the scaling by
   :math:`(\mathtt{ref\_ratio})^2` is not applied.

-  ``gravity.rel_tol`` : if ``gravity.gravity_type`` = ``PoissonGrav``,
   this is the relative tolerance for the Poisson solve. By default it
   is zero. You can specify a single value for this tolerance and it
   will apply on every level, or you can specify an array of values
   for ``rel_tol``, one for each possible level in the
   simulation. This replaces the old parameter ``gravity.ml_tol``.

-  ``gravity.max_multipole_order`` : if ``gravity.gravity_type`` =
   ``PoissonGrav``, this is the max :math:`\ell` value to use for
   multipole BCs (must be :math:`\geq 0`; default: 0)

-  ``gravity.direct_sum_bcs`` : if ``gravity.gravity_type`` =
   ``PoissonGrav``, evaluate BCs using exact sum (0 or 1; default: 0)

-  ``gravity.drdxfac`` : ratio of dr for monopole gravity
   binning to grid resolution

The follow parameters affect the coupling of hydro and gravity:

-  ``castro.do_grav`` : turn on/off gravity

-  ``castro.moving_center`` : do we recompute the center
   used for the multipole gravity solver each step?

-  ``castro.point_mass`` : point mass at the center of the star
   (must be :math:`\geq 0`; default: 0.0)

Note that in the following, ``MAX_LEV`` is a hard-coded parameter
in ``Source/Gravity.cpp`` which is currently set to 15. It
determines how many levels can be tracked by the ``Gravity`` object.

Types of Approximations
=======================

``ConstantGrav``
----------------

Gravity can be defined as constant in direction and magnitude,
defined in the inputs file by::

   gravity.const_grav = -9.8

for example, to set the gravity to have magnitude :math:`9.8` in the
negative :math:`y`-direction if in 2D, negative :math:`z`-direction if in 3-D.
The actual setting is done in Gravity.cpp as::

     grav.setVal(const_grav, AMREX_SPACEDIM-1, 1, ng);

Note that at present we do not fill the gravitational potential
:math:`\phi` in this mode; it will be set to zero.

Note: ``ConstantGrav`` can only be used along a Cartesian direction
(vertical for 2D axisymmetric).

.. _sec-monopole-grav:

``MonopoleGrav``
----------------

``MonopoleGrav`` integrates the mass distribution on the grid
in spherical shells, defining an enclosed mass and uses this
to compute the gravitational potential and acceleration in a
spherically-symmetric fashion.

-  In 1D spherical coordinates we compute

   .. math:: g(r) = -\frac{G M_{\rm enclosed}}{ r^2}

   where :math:`M_{\rm enclosed}` is calculated from the density at
   the time of the call.

   For levels above the coarsest level we define the extent of that
   level’s radial arrays as ranging from the center of the star (:math:`r=0`)
   to the cell at that level farthest away from the origin. If there
   are gaps between fine grids in that range then we interpolate the
   density from a coarser level in order to construct a continuous
   density profile. We note that the location of values in the density
   profile and in the gravitational field exactly match the location of
   data at that level so there is no need to interpolate between points
   when mapping the 1D radial profile of :math:`g` back onto the original
   grid.

-  In 2D or 3D we compute a 1D radial average of density and use
   this to compute gravity as a one-dimensional integral, then
   interpolate the gravity vector back onto the Cartesian grid
   cells. At the coarsest level we define the extent of the 1D arrays
   as ranging from the center of the star to the farthest possible
   point in the grid (plus a few extra cells so that we can fill ghost
   cell values of gravity). At finer levels we first define a single
   box that contains all boxes on that fine level, then we interpolate
   density from coarser levels as needed to fill the value of density
   at every fine cell in that box. The extent of the radial array is
   from the center of the star to the *nearest* cell on one of the
   faces of the single box. This ensures that all cells at that
   maximum radius of the array are contained in this box.

   We then average the density onto a 1D radial array. We note that
   there is a mapping from the Cartesian cells to the radial array and
   back; unlike the 1D case this requires interpolation. We use
   quadratic interpolation with limiting so that the interpolation
   does not create new maxima or minima.

   The default resolution of the radial arrays at a level is the grid
   cell spacing at that level, i.e., :math:`\Delta r = \Delta x`.
   For increased accuracy, one can define ``gravity.drdxfac`` as a number
   greater than :math:`1` (:math:`2` or :math:`4` are recommended) and
   the spacing of the radial array will then satisfy :math:`\Delta x /
   \Delta r =` drdxfac.  Individual Cartesian grid cells are
   subdivided by drdxfac in each coordinate direction for the
   purposing of averaging the density, and the integration that
   creates :math:`g` is done at the finer resolution of the new
   :math:`\Delta r`.

   Note that the center of the star is defined in the subroutine
   ``probinit`` and the radius is computed as the distance from that
   center.

   .. note:: there is an additional correction at the corners in
             ``make_radial_grav`` that accounts for the volume in a shell
             that is not part of the grid.

What about the potential in this case? when does
``make_radial_phi`` come into play?

``PoissonGrav``
---------------

The most general case is a self-induced gravitational field,

.. math:: \mathbf{g}(\mathbf{x},t) = \nabla \phi

where :math:`\phi` is defined by solving

.. math::
   \mathbf{\Delta} \phi = 4 \pi G \rho
   :label: eq:Self Gravity

We only allow ``PoissonGrav`` in 2D or 3D because in 1D, computing
the monopole approximation in spherical coordinates is faster and more
accurate than solving the Poisson equation.

Poisson Boundary Conditions: 2D
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In 2D, if boundary conditions are not periodic in both directions, we
use a monopole approximation at the coarsest level. This involves
computing an effective 1D radial density profile (on level = 0 only),
integrating it outwards from the center to get the gravitational
acceleration :math:`\mathbf{g}`, and then integrating :math:`g`
outwards from the center to get :math:`\phi` (using :math:`\phi(0) =
0` as a boundary condition, since no mass is enclosed at :math:`r =
0`). For more details, see Section `2.2 <#sec-monopole-grav>`__.

.. _sec-poisson-3d-bcs:

Poisson Boundary Conditions: 3D
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following describes methods for doing isolated boundary
conditions. The best reference for Castro’s implementation of this
is :cite:`katz:2016`.

-  **Multipole Expansion**

   In 3D, by default, we use a multipole expansion to estimate the value
   of the boundary conditions. According to, for example, Jackson’s
   *Classical Electrodynamics* (with the corresponding change to
   Poisson’s equation for electric charges and gravitational
   ”charges”), an expansion in spherical harmonics for :math:`\phi` is

   .. math::
      \phi(\mathbf{x}) = -G\sum_{l=0}^{\infty}\sum_{m=-l}^{l} \frac{4\pi}{2l + 1} q_{lm} \frac{Y_{lm}(\theta,\phi)}{r^{l+1}},
      :label: spherical_harmonic_expansion

   The origin of the coordinate system is taken to be the ``center``
   variable, that must be declared and stored in the ``probdata``
   module in your project directory. The validity of the expansion used
   here is based on the assumption that a sphere centered on
   ``center``, of radius approximately equal to the size of half the
   domain, would enclose all of the mass. Furthermore, the lowest order
   terms in the expansion capture further and further departures from
   spherical symmetry. Therefore, it is crucial that ``center`` be
   near the center of mass of the system, for this approach to achieve
   good results.

   The multipole moments :math:`q_{lm}` can be calculated by expanding the
   Green’s function for the Poisson equation as a series of spherical
   harmonics, which yields

   .. math:: q_{lm} = \int Y^*_{lm}(\theta^\prime, \phi^\prime)\, {r^\prime}^l \rho(\mathbf{x}^\prime)\, d^3x^\prime. \label{multipole_moments_original}

   Some simplification of :eq:`spherical_harmonic_expansion` can
   be achieved by using the addition theorem for spherical harmonics:

   .. math::

      \begin{aligned}
        &\frac{4\pi}{2l+1} \sum_{m=-l}^{l} Y^*_{lm}(\theta^\prime,\phi^\prime)\, Y_{lm}(\theta, \phi) = P_l(\text{cos}\, \theta) P_l(\text{cos}\, \theta^\prime) \notag \\
        &\ \ + 2 \sum_{m=1}^{l} \frac{(l-m)!}{(l+m)!} P_{l}^{m}(\text{cos}\, \theta)\, P_{l}^{m}(\text{cos}\, \theta^\prime)\, \left[\text{cos}(m\phi)\, \text{cos}(m\phi^\prime) + \text{sin}(m\phi)\, \text{sin}(m\phi^\prime)\right].\end{aligned}

   Here the :math:`P_{l}^{m}` are the associated Legendre polynomials and the
   :math:`P_l` are the Legendre polynomials. After some algebraic
   simplification, the potential outside of the mass distribution can be
   written in the following way:

   .. math:: \phi(\mathbf{x}) \approx -G\sum_{l=0}^{l_{\text{max}}} \left[Q_l^{(0)} \frac{P_l(\text{cos}\, \theta)}{r^{l+1}} + \sum_{m = 1}^{l}\left[ Q_{lm}^{(C)}\, \text{cos}(m\phi) + Q_{lm}^{(S)}\, \text{sin}(m\phi)\right] \frac{P_{l}^{m}(\text{cos}\, \theta)}{r^{l+1}} \right].

   The modified multipole moments are:

   .. math::

      \begin{aligned}
        Q_l^{(0)}   &= \int P_l(\text{cos}\, \theta^\prime)\, {r^{\prime}}^l \rho(\mathbf{x}^\prime)\, d^3 x^\prime \\
        Q_{lm}^{(C)} &= 2\frac{(l-m)!}{(l+m)!} \int P_{l}^{m}(\text{cos}\, \theta^\prime)\, \text{cos}(m\phi^\prime)\, {r^\prime}^l \rho(\mathbf{x}^\prime)\, d^3 x^\prime \\
        Q_{lm}^{(S)} &= 2\frac{(l-m)!}{(l+m)!} \int P_{l}^{m}(\text{cos}\, \theta^\prime)\, \text{sin}(m\phi^\prime)\, {r^\prime}^l \rho(\mathbf{x}^\prime)\, d^3 x^\prime.\end{aligned}

   Our strategy for the multipole boundary conditions, then, is to pick
   some value :math:`l_{\text{max}}` that is of sufficiently high order to
   capture the distribution of mass on the grid, evaluate the discretized
   analog of the modified multipole moments for :math:`0 \leq l \leq
   l_{\text{max}}` and :math:`1 \leq m \leq l`, and then directly compute the
   value of the potential on all of the boundary zones. This is
   ultimately an :math:`\mathcal{O}(N^3)` operation, the same order as the
   monopole approximation, and the wall time required to calculate the
   boundary conditions will depend on the chosen value of
   :math:`l_{\text{max}}`.

   The number of :math:`l` values calculated is controlled by
   ``gravity.max_multipole_order`` in your inputs file. By default, it
   is set to ``0``, which means that a monopole approximation is
   used. There is currently a hard-coded limit of
   :math:`l_{\text{max}} = 50`. This is because the method used to
   generate the Legendre polynomials is not numerically stable for
   arbitrary :math:`l` (because the polynomials get very large, for
   large enough :math:`l`).

-  **Direct Sum**

   Up to truncation error caused by the discretization itself, the
   boundary values for the potential can be computed exactly by a direct
   sum over all cells in the grid. Suppose I consider some ghost cell
   outside of the grid, at location :math:`\mathbf{r}^\prime \equiv (x^\prime,
   y^\prime, z^\prime)`. By the principle of linear superposition as
   applied to the gravitational potential,

   .. math:: \phi(\mathbf{r}^\prime) = \sum_{\text{ijk}} \frac{-G \rho_{\text{ijk}}\, \Delta V_{\text{ijk}}}{\left[(x - x^\prime)^2 + (y - y^\prime)^2 + (z - z^\prime)^2\right]^{1/2}},

   where :math:`x = x(i)`, :math:`y = y(j)` and :math:`z = z(k)` are
   constructed in the usual sense from the zone indices. The sum here
   runs over every cell in the physical domain (that is, the
   calculation is :math:`\mathcal{O}(N^3)` for each boundary
   cell). There are :math:`6N^2` ghost cells needed for the Poisson
   solve (since there are six physical faces of the domain), so the
   total cost of this operation is :math:`\mathcal{O}(N^5)` (this only
   operates on the coarse grid, at present). In practice, we use the
   domain decomposition inherent in the code to implement this solve:
   for the grids living on any MPI task, we create six :math:`N^2`
   arrays representing each of those faces, and then iterate over
   every cell on each of those grids, and compute their respective
   contributions to all of the faces. Then, we do a global reduce to
   add up the contributions from all cells together. Finally, we place
   the boundary condition terms appropriate for each grid onto its
   respective cells.

   This is quite expensive even for reasonable sized domains, so this
   option is recommended only for analysis purposes, to check if the
   other methods are producing accurate results. It can be enabled by
   setting ``gravity.direct_sum_bcs`` = 1 in your inputs file.

Point Mass
----------

Pointmass gravity works with all other forms of gravity, it is not a
separate option. Since the Poisson equation is linear in potential
(and its derivative, the acceleration, is also linear), the point mass
option works by adding the gravitational acceleration of the point
mass onto the acceleration from whatever other gravity type is under
in the simulation.

.. note:: The point mass may have a mass < 0

A useful option is ``point_mass_fix_solution``. If set to 1, then it
takes all zones that are adjacent to the location of the center
variable and keeps their density constant. Any changes in density that
occur after a hydro update in those zones are reset, and the mass
deleted is added to the pointmass. (If there is expansion, and the
density lowers, then the point mass is reduced and the mass is added
back to the grid). This calculation is done in
``pointmass_update()`` in ``Castro_pointmass.cpp``.

GR correction
=============

In the cases of compact objects or very massive stars, the general
relativity (GR) effect starts to play a role [2]_. First, we consider
the hydrostatic equilibrium due to effects of GR then derive
GR-correction term for Newtonian gravity.  The correction term is
applied to the monopole approximation only when ``USE_GR = TRUE`` is
set in the ``GNUmakefile``.

The formulae of GR-correction here are based on
:cite:`grbk1`. For detailed physics, please refer to
:cite:`grbk2`. For describing very strong gravitational
field, we need to use Einstein field equations

.. math::
   R_{ik}-\frac{1}{2}g_{ik}R=\frac{\kappa}{c^{2}}T_{ik} \quad , \quad
   \kappa=\frac{8\pi G}{c^{2}}\quad ,
   :label: field

where :math:`R_{ik}` is the Ricci tensor, :math:`g_{ik}` is the metric
tensor, :math:`R` is the Riemann curvature, :math:`c` is the speed of
light and :math:`G` is gravitational constant. :math:`T_{ik}` is the
energy momentum tensor, which for ideal gas has only the non-vanishing
components :math:`T_{00}` = :math:`\varrho c^2` , :math:`T_{11}` =
:math:`T_{22}` = :math:`T_{33}` = :math:`P` ( contains rest mass and
energy density, :math:`P` is pressure). We are interested in
spherically symmetric mass distribution. Then the line element
:math:`ds` for given spherical coordinate :math:`(r, \vartheta,
\varphi)` has the general form

.. math::

   \label{metric}
     ds^{2} = e^{\nu}c^{2}dt^{2}-e^{\lambda}dr^{2}-r^{2}(d\vartheta^{2}+\sin^{2}
     \vartheta d\varphi) \quad ,

with :math:`\nu = \nu(r)`, :math:`\lambda = \lambda(r)`. Now we can
put the expression of :math:`T_{ik}` and :math:`ds` into :eq:`field`,
then field equations can be reduced to 3 ordinary differential
equations:

.. math::
   \frac{\kappa P}{c^{2}} =
      e^{-\lambda}\left (\frac{\nu^{\prime}}{r}+\frac{1}{r^{2}} \right )-\frac{1}{r^{2}}
      \quad ,
   :label: diff1

.. math::
   \frac{\kappa P}{c^{2}} =
     \frac{1}{2}e^{-\lambda}\left (\nu^{\prime\prime}+\frac{1}{2}{\nu^{\prime}}^{2}+\frac{\nu^
       {\prime}-\lambda^{\prime}}{r}
      -\frac{\nu^{\prime}\lambda^{\prime}}{2} \right ) \quad ,
   :label: diff2

.. math::
   \kappa \varrho =
     e^{-\lambda}\left (\frac{\lambda^{\prime}}{r}-\frac{1}{r^{2}}\right )+\frac{1}{r^{2}} \quad ,
   :label: diff3

where primes means the derivatives with respect to :math:`r`. After
multiplying with :math:`4\pi r^2`, :eq:`diff3` can be
integrated and yields

.. math::

   \label{gmass1}
     \kappa m = 4\pi r (1-e^{-\lambda}) \quad ,

the :math:`m` is called “gravitational mass” inside r defined as

.. math::

   \label{gmass2}
     m = \int_{0}^{r}4\pi r^{2}  \varrho dr\quad .

For the :math:`r = R`, :math:`m` becomes the mass :math:`M` of the
star. :math:`M` contains not only the rest mass but the whole energy
(divided by :math:`c^2`), that includes the internal and gravitational
energy. So the :math:`\varrho = \varrho_0 +U/c^2` contains the whole
energy density :math:`U` and rest-mass density
:math:`\varrho_0`. Differentiation of :eq:`diff1` with respect to
:math:`r` gives :math:`P = P^{\prime}(\lambda,\lambda^{\prime},
\nu,\nu^{\prime},r)`, where
:math:`\lambda,\lambda^{\prime},\nu,\nu^{\prime}` can be eliminated by
:eq:`diff1`, :eq:`diff2`, :eq:`diff3`. Finally we reach
*Tolman-Oppenheinmer-Volkoff (TOV)* equation for hydrostatic
equilibrium in general relativity:

.. math::
   \frac{dP}{dr} = -\frac{Gm}{r^{2}}\varrho \left (1+\frac{P}{\varrho
       c^{2}}\right )\left (1+\frac{4\pi r^3 P}{m c^{2}}\right ) \left (1-\frac{2Gm}{r c^{2}} \right)^{-1} \quad .
   :label: tov

For Newtonian case :math:`c^2 \rightarrow  \infty`, it reverts to usual form

.. math::

   \label{newton}
     \frac{dP}{dr} = -\frac{Gm}{r^{2}}\varrho \quad .

Now we take effective monopole gravity as

.. math::
   \tilde{g} = -\frac{Gm}{r^{2}} (1+\frac{P}{\varrho
     c^{2}})(1+\frac{4\pi r^3 P}{m c^{2}}) (1-\frac{2Gm}{r c^{2}})^{-1}  \quad .
   :label: tov2

For general situations, we neglect the :math:`U/c^2` and potential
energy in m because they are usually much smaller than
:math:`\varrho_0`. Only when :math:`T` reaches :math:`10^{13} K`
(:math:`KT \approx m_{p} c^2`, :math:`m_p` is proton mass) before it
really makes a difference. So :eq:`tov2` can be expressed as

.. math::

   \label{tov3}
     \tilde{g} = -\frac{GM_{\rm enclosed}}{r^{2}} \left (1+\frac{P}{\varrho
       c^{2}} \right )\left (1+\frac{4\pi r^3 P}{M_{\rm enclosed} c^{2}} \right ) \left (1-\frac{2GM_{\rm enclosed}}{r c^{2}} \right )^{-1} \quad ,

where :math:`M_{enclosed}` has the same meaning as with the
``MonopoleGrav`` approximation.

Hydrodynamics Source Terms
==========================

We use a standard predictor-corrector formalism for updating the momentum.
Our first update is equal to :math:`\Delta t \times \mathbf{S}^n`, where
:math:`\mathbf{S}^n` is the value of the source terms at the old-time
(which is usually called time-level :math:`n`). At the end of the timestep,
we do a corrector step where we subtract off :math:`\Delta t / 2 \times \mathbf{S}^n`
and add on :math:`\Delta t / 2 \times \mathbf{S}^{n+1}`, so that at the end of
the timestep the source term is properly time centered.

The energy update is done in an explicitly “conservative” fashion.
We evaluate the change in energy at cell edges, using the hydrodynamical mass
fluxes, permitting total energy to be conserved (excluding possible losses at
open domain boundaries). See :cite:`katzthesis` for some more details.

.. [2]
   Note: The GR
   code and text here were contributed by Ken Chen of Univ. of
   Minnesota.
