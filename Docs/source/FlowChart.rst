*********
Flowchart
*********

Introduction
============

There are several different time-evolution methods currently
implemented in Castro. As best as possible, they share the same
driver routines and use preprocessor or runtime variables to separate
the different code paths.  These fall into two categories:

.. index:: castro.time_integration_method, USE_SIMPLIFIED_SDC, USE_TRUE_SDC

-  Strang+CTU: the Strang evolution does the burning on the
   state for :math:`\Delta t/2`, then updates the hydrodynamics using the
   burned state, and then does the final :math:`\Delta t/2` burning. No
   explicit coupling of the burning and hydro is done.  This code
   path uses the corner-transport upwind (CTU) method (the unsplit,
   characteristic tracing method of :cite:`colella:1990`).  This is the default method.

   The MHD solver uses this same driver.

-  SDC: a class of iterative methods that couples the advection and reactions
   such that each process explicitly sees the effect of the other.  We have
   two SDC implementations in Castro.

   - The "simplified SDC" method is based on the CTU hydro update.  We
     iterate over the construction of this term, using a lagged
     reaction source as inputs and do the final conservative update by
     integrating the reaction system using an ODE solver with the
     explicit advective source included in a
     piecewise-constant-in-time fastion.

   - The "true SDC" method.  This fully couples the hydro and reactions
     to either 2nd or 4th order.  This approximates the integral in
     time using a simple quadrature rule, and integrates the hydro
     explicitly and reactions implicitly to the next time node.
     Iterations allow each process to see one another and achieve
     high-order in time convergence.  This is described in :cite:`castro-sdc`.


The time-integration method used is controlled by
``castro.time_integration_method``.

  * ``time_integration_method = 0``: this is the original Castro method,
    described in :cite:`castro_I`.  This uses Strang splitting and the CTU
    hydrodynamics scheme.

  * ``time_integration_method = 1``: unused (in Castro 19.08 and
    earlier, this was a method-of-lines integration method with Strang
    splitting for reactions.)

  * ``time_integration_method = 2``: this is a full implementation of
    the spectral deferred corrections formalism, with both 2nd and 4th
    order integration implemented.  At the moment, this does not support
    multilevel domains.  Note: because of differences in the interfaces with the 
    default Strang method, you must compile with ``USE_TRUE_SDC = TRUE`` for this
    method to work (in particular, this defines ``EXTRA_THERMO`` which enables some
    additional EOS derivatives).

  * ``time_integration_method = 3``: this is the simplifed SDC method
    described above that uses the CTU hydro advection and an ODE
    reaction solve.  Note: because this requires a different set of
    state variables, you must compile with ``USE_SIMPLIFIED_SDC = TRUE`` for this
    method to work (in particular, this defines ``PRIM_SPECIES_HAVE_SOURCES``).

.. index:: USE_SIMPLIFIED_SDC, USE_TRUE_SDC

.. note::

   By default, the code is compiled for Strang-split CTU evolution
   (``time_integration_method = 0``).  Because the size of the
   different state arrays differs with the other integration schemes,
   support for them needs to be compiled in, using
   ``USE_SIMPLIFIED_SDC=TRUE`` for the simplified-SDC method
   (``time_integration_method=3``) and ``USE_TRUE_SDC=TRUE`` for the
   true SDC method (``time_integration_method = 2``).

.. note::

   MHD and radiation are currently only supported by the Strang+CTU
   evolution time integration method.

Several helper functions are used throughout:

.. index:: clean_state

-  ``clean_state``:
   There are many ways that the hydrodynamics state may become
   unphysical in the evolution. The ``clean_state()`` routine
   enforces some checks on the state. In particular, it

   #. enforces that the density is above ``castro.small_dens``

   #. enforces that the speeds in the state don't exceed ``castro.speed_limit``

   #. normalizes the species so that the mass fractions sum to 1

   #. syncs up the linear and hybrid momenta (for ``USE_HYBRID_MOMENTUM=TRUE``)

   #. resets the internal energy if necessary (too small or negative)
      and computes the temperature for all zones to be thermodynamically
      consistent with the state.

.. _flow:sec:nosdc:

Main Driver—All Time Integration Methods
========================================

This driver supports the Strang CTU integration.
(``castro.time_integration_method`` = 0)

The main evolution for a single step is contained in
``Castro_advance.cpp``, as ``Castro::advance()``. This does
the following advancement. Note, some parts of this are only done
depending on which preprocessor directives are defined at
compile-time—the relevant directive is noted in the [ ] at the start
of each step.

#. *Initialization* (``initialize_advance()``)

   This sets up the current level for advancement. The following
   actions are performend (note, we omit the actions taken for a retry,
   which we will describe later):

   -  Sync up the level information to the Fortran-side of Castro.

   -  Do any radiation initialization.

   -  Set the maximum density used for Poisson gravity tolerances.

   -  Initialize all of the intermediate storage arrays (like those
      that hold source terms, etc.).

   -  Swap the StateData from the new to old (e.g., ensures that
      the next evolution starts with the result from the previous step).

   -  Call ``clean_state``.

   -  Create the MultiFabs that hold the primitive variable information
      for the hydro solve.

   -  Zero out all of the fluxes.

   -  For true SDC, initialize the data at all time nodes (see :ref:`sec:flow_true_sdc`).

#. *Advancement*

   Call ``do_advance`` to take a single step, incorporating
   hydrodynamics, reactions, and source terms.

   For radiation-hydrodynamics, this step does the
   advective (hyperbolic) portion of the radiation update only.
   Source terms, including gravity, rotation, and diffusion are
   included in this step, and are time-centered to achieve second-order
   accuracy.

   .. index:: retry

   If ``castro.use_retry`` is set, then we subcycle the current
   step if we violated any stability criteria to reach the desired
   :math:`\Delta t`. The idea is the following: if the timestep that you
   took had a timestep that was not sufficient to enforce the stability
   criteria that you would like to achieve, such as the CFL criterion
   for hydrodynamics or the burning stability criterion for reactions,
   you can retry the timestep by setting ``castro.use_retry`` = 1 in
   your inputs file. This will save the current state data at the
   beginning of the level advance, and then if the criteria are not
   satisfied, will reject that advance and start over from the old
   data, with a series of subcycled timesteps that should be small
   enough to satisfy the criteria. Note that this will effectively
   double the memory footprint on each level if you choose to use it.
   See :ref:`ch:retry` for more details on the retry mechanism.

   .. note::

      Only Strang+CTU and simplified-SDC support retries.

#. [AUX_UPDATE] *Auxiliary quantitiy evolution*

   Auxiliary variables in Castro are those that obey a continuity
   equation (with optional sources) that are passed into the EOS, but
   not subjected to the constraint on mass fractions (summing to one).

   The advection and source terms are already dealt with in the
   main hydrodynamics advance (above step). A user-supplied routine
   ca_auxupdate can be provided here to further update these
   quantities.

#. *Radial data and [POINTMASS] point mass*

   If ``castro.spherical_star`` is set, then we average the state data
   over angles here to create a radial profile. This is then used in the
   boundary filling routines to properly set Dirichlet BCs when our domain
   is smaller than the star, so the profile on the boundaries will not
   be uniform.

   If ``castro.point_mass_fix_solution`` is set, then we
   change the mass of the point mass that optionally contributes to the
   gravitational potential by taking mass from the surrounding zones
   (keeping the density in those zones constant).

#. [RADIATION] *Radiation implicit update*

   The ``do_advance()`` routine only handled the hyperbolic
   portion of the radiation update. This step does the implicit solve
   (either gray or multigroup) to advance the radiation energies to the
   new time level. Note that at the moment, this is backward-difference
   implicit (first-order in time) for stability.

   This is handled by ``final_radiation_call()``.

#. [PARTICLES] *Particles*

   If we are including passively-advected particles, they are
   advanced in this step.

#. *Finalize*

   This cleans up at the end of a step:

   -  Update the flux registers to account for mismatches at
      coarse-fine interfaces. This cleans up the memory used during
      the step.

   -  Free any memory allocated for the level advance.


.. _sec:strangctu:

Strang+CTU Evolution
====================

``do_advance_ctu()`` in ``Castro_advance_ctu.cpp`` 

This described the flow using Strang splitting and the CTU
hydrodynamics (or MHD) method, including gravity, rotation, and
diffusion.  This integration is selected via
``castro.time_integration_method = 0``.

The system advancement: reactions, hydrodynamics, diffusion, rotation,
and gravity are all considered here.

Consider our system of equations as:

.. math:: \frac{\partial\Ub}{\partial t} = {\bf A}(\Ub) + \Rb(\Ub) + \Sb,

where :math:`{\bf A}(\Ub) = -\nabla \cdot \Fb(\Ub)`, with :math:`\Fb` the flux vector, :math:`\Rb` are the reaction
source terms, and :math:`\Sb` are the non-reaction source terms, which
includes any user-defined external sources, :math:`\Sb_{\rm ext}`. We use
Strang splitting to discretize the advection-reaction equations. In
summary, for each time step, we update the conservative variables,
:math:`\Ub`, by reacting for half a time step, advecting for a full time
step (ignoring the reaction terms), and reacting for half a time step.
The treatment of source terms complicates this a little. The actual
update, in sequence, looks like:

.. math::
   \begin{aligned}
   \Ub^\star &= \Ub^n + \frac{\dt}{2}\Rb(\Ub^n) \\
   \Ub^{n+1,(a)} &= \Ub^\star + \dt\, \Sb(\Ub^\star) \\
   \Ub^{n+1,(b)} &= \Ub^{n+1,(a)} + \dt\, {\bf A}(\Ub^\star) \\
   \Ub^{n+1,(c)} &= \Ub^{n+1,(b)} + \frac{\dt}{2}\, [\Sb(\Ub^{n+1,(b)}) - \Sb(\Ub^\star)] \\
   \Ub^{n+1}     &= \Ub^{n+1,(c)} + \frac{\dt}{2} \Rb(\Ub^{n+1,(c)})
   \end{aligned}
   :label: eq:source_correct

Note that in the first step, we add a full :math:`\Delta t` of the old-time
source to the state. This prediction ensures consistency when it
comes time to predicting the new-time source at the end of the update.
The construction of the advective terms, :math:`{\bf A(\Ub)}` is purely
explicit, and based on an unsplit second-order Godunov method. We
predict the standard primitive variables, as well as :math:`\rho e`, at
time-centered edges and use an approximate Riemann solver construct
fluxes.

At the beginning of the time step, we assume that :math:`\Ub` and the gravitational potential, :math:`\phi`, are
defined consistently, i.e., :math:`\rho^n` and :math:`\phi^n` satisfy the Poisson equation:

.. math::

   \Delta \phi^n = 4\pi G\rho^n

(see :ref:`ch:gravity` for more details about how the Poisson equation is solved.)  
Note that in
:eq:`eq:source_correct`, we can actually do some
sources implicitly by updating density first, and then momentum,
and then energy. This is done for rotating and gravity, and can
make the update more akin to:

.. math:: \Ub^{n+1,(c)} = \Ub^{n+1,(b)} + \frac{\dt}{2} [\Sb(\Ub^{n+1,(c)}) - \Sb(\Ub^n)]

If we are including radiation, then this part of the update algorithm
only deals with the advective / hyperbolic terms in the radiation update.

Here is the single-level algorithm. The goal here is to update the
``State_Type``  ``StateData`` from the old to new time (see
§ :ref:`soft:sec:statedata`). We will use the following notation
here, consistent with the names used in the code:

-  ``S_old`` is a MultiFab reference to the old-time-level
   ``State_Type`` data.

-  ``Sborder`` is a MultiFab that has ghost cells and is
   initialized from ``S_old``. This is what the hydrodynamic
   reconstruction will work from.

-  ``S_new`` is a MultiFab reference to the new-time-level
   ``State_Type`` data.

- ``old_source`` is a MultiFab reference to the old-time-level ``Source_Type`` data.

- ``new_source`` is a MultiFab reference to the new-time-level ``Source_Type`` data.


Single Step Flowchat
--------------------

In the code, the objective is to evolve the state from the old time,
``S_old``, to the new time, ``S_new``.

#. *Initialize*

   A. In ``initialize_do_advance()``, create ``Sborder``, initialized from ``S_old``

   B. If ``source_term_predictor == 1``, then aply
      ``source_corrector`` to ``sources_for_hydro``.  For Strang+CTU,
      this was the effect of initializing the source terms as:

      .. math::

         \Sb = \frac{dt}{2} (d\Sb/dt)^{n-1/2}

   C. Check for NaNs in the initial state, ``S_old``.


#. *React* :math:`\Delta t/2` [``strang_react_first_half()`` ]

   Update the solution due to the effect of reactions over half a time
   step. The integration method and system of equations used here is
   determined by a host of runtime parameters that are part of the
   Microphysics package. But the basic idea is to evolve the energy
   release from the reactions, the species mass fractions, and
   temperature through :math:`\Delta t/2`.

   Using the notation above, we begin with the time-level :math:`n` state,
   :math:`\Ub^n`, and produce a state that has evolved only due to reactions,
   :math:`\Ub^\star`.

   .. math::

      \begin{aligned}
          (\rho e)^\star &= (\rho e)^n - \frac{\dt}{2} \rho H_\mathrm{nuc} \\
          (\rho E)^\star &= (\rho E)^n - \frac{\dt}{2} \rho H_\mathrm{nuc} \\
          (\rho X_k)^\star &= (\rho X_k)^n + \frac{\dt}{2}(\rho\omegadot_k).
        \end{aligned}

   Here, :math:`H_\mathrm{nuc}` is the energy release (erg/g/s) over the
   burn, and :math:`\omegadot_k` is the creation rate for species :math:`k`.

   After exiting the burner, we call the EOS with :math:`\rho^\star`,
   :math:`e^\star`, and :math:`X_k^\star` to get the new temperature, :math:`T^\star`.

   Note that the density, :math:`\rho`, does not change via reactions in the
   Strang-split formulation.

   The reaction data needs to be valid in the ghost cells, so the reactions
   are applied to the entire patch, including ghost cells.

   After reactions, ``clean_state`` is called.

   At the end of this step, ``Sborder`` sees the effects of the
   reactions.

#. *Construct time-level n sources and apply*
   [``construct_old_gravity()``, ``do_old_sources()`` ]

   The time level :math:`n` sources are computed, and added to the
   StateData ``Source_Type``. 

   The sources that we deal with here are:

   A. sponge : the sponge is a damping term added to
      the momentum equation that is designed to drive the velocities to
      zero over some timescale. Our implementation of the sponge
      follows that of Maestro :cite:`maestro:III`

   B. external sources : users can define problem-specific sources
      in the ``ext_src_?d.f90`` file. Sources for the different
      equations in the conservative state vector, :math:`\Ub`, are indexed
      using the integer keys defined in ``meth_params_module``
      (e.g., URHO).

      This is most commonly used for external heat sources (see the
      ``toy_convect`` problem setup) for an example. But most
      problems will not use this.

   C. [``MHD``] thermal source: for the MHD system, we are including
      the "pdV" work for the internal energy equation as a source term
      rather than computing it from the Riemann problem.  This source is
      computed here for the internal energy equation.

   D. [``DIFFUSION``] diffusion : thermal diffusion can be
      added in an explicit formulation. Second-order accuracy is
      achieved by averaging the time-level :math:`n` and :math:`n+1` terms, using
      the same predictor-corrector strategy described here.

      Note: thermal diffusion is distinct from radiation hydrodynamics.

      Also note that incorporating diffusion brings in an additional
      timestep constraint, since the treatment is explicit. See
      Chapter :ref:`ch:diffusion` for more details.

   E. [``HYBRID_MOMENTUM``] angular momentum


   F. [``GRAVITY``] gravity:

      For full Poisson gravity, we solve for for gravity using:

      .. math::

         \gb^n = -\nabla\phi^n, \qquad
               \Delta\phi^n = 4\pi G\rho^n,

      The construction of the form of the gravity source for the
      momentum and energy equation is dependent on the parameter
      ``castro.grav_source_type``. Full details of the gravity
      solver are given in Chapter :ref:`ch:gravity`.


   G. [``ROTATION``] rotation

      We compute the rotational potential (for use in the energy update)
      and the rotational acceleration (for use in the momentum
      equation). This includes the Coriolis and centrifugal terms in a
      constant-angular-velocity co-rotating frame. The form of the
      rotational source that is constructed then depends on the
      parameter ``castro.rot_source_type``. More details are
      given in Chapter :ref:`ch:rotation`.

   The source terms here are evaluated using the post-burn state,
   :math:`\Ub^\star` (``Sborder``), and later corrected by using the
   new state just before the burn, :math:`\Ub^{n+1,(b)}`. This is compatible
   with Strang-splitting, since the hydro and sources takes place
   completely inside of the surrounding burn operations.

   The old-time source terms are stored in ``old_source``.

   The sources are then applied to the state after the burn,
   :math:`\Ub^\star` with a full :math:`\Delta t` weighting (this will
   be corrected later). This produces the intermediate state,
   :math:`\Ub^{n+1,(a)}` (stored in ``S_new``).

#. *Construct the hydro / MHD update* [``construct_ctu_hydro_source()``, ``construct_ctu_mhd_source()``]

   The goal is to advance our system considering only the advective
   terms (which in Cartesian coordinates can be written as the
   divergence of a flux).

   We do the hydro update in two parts—first we construct the
   advective update and store it in the hydro_source
   MultiFab, then we do the conservative update in a separate step. This
   separation allows us to use the advective update separately in more
   complex time-integration schemes.

   In the Strang-split formulation, we start the reconstruction using
   the state after burning, :math:`\Ub^\star` (``Sborder``).  For the
   CTU method, we predict to the half-time (:math:`n+1/2`) to get a
   second-order accurate method. Note: ``Sborder`` does not know of
   any sources except for reactions. 

   The method done here differs depending on whether we are doing hydro or MHD.

   A. hydrodynamics

      The advection step is complicated, and more detail is given in
      Section :ref:`Sec:Advection Step`. Here is the summarized version:

      i. Compute primitive variables.

      ii. Convert the source terms to those acting on primitive variables

      iii. Predict primitive variables to time-centered edges.

      iv. Solve the Riemann problem.

      v. Compute fluxes and advective term.

   B. MHD

      The MHD update is described in :ref:`ch:mhd`.

   To start the hydrodynamics/MHD source construction, we need to know
   the hydrodynamics source terms at time-level :math:`n`, since this
   enters into the prediction to the interface states. This is
   essentially the same vector that was computed in the previous step,
   with a few modifications. The most important is that if we set
   ``castro.source_term_predictor``, then we extrapolate the source
   terms from :math:`n` to :math:`n+1/2`, using the change from the
   previous step.

   Note: we neglect the reaction source terms, since those are already
   accounted for in the state directly, due to the Strang-splitting
   nature of this method.

   The update computed here is then immediately applied to
   ``S_new``.

#. *Clean State* [``clean_state()``]

   This is done on ``S_new``.

   After these checks, we check the state for NaNs.

#. *Update radial data and center of mass for monopole gravity*

   These quantities are computed using ``S_new``.

#. *Correct the source terms with the n+1
   contribution* [``construct_new_gravity()``, ``do_new_sources`` ]

   If we are doing self-gravity, then we first compute the updated gravitational
   potential using the updated density from ``S_new``.

   Now we correct the source terms applied to ``S_new`` so they are time-centered.
   Previously we added :math:`\Delta t\, \Sb(\Ub^\star)` to the state, when
   we really want 
   :math:`(\Delta t/2)[\Sb(\Ub^\star + \Sb(\Ub^{n+1,(b)})]` .

   We start by computing the source term vector :math:`\Sb(\Ub^{n+1,(b)})`
   using the updated state, :math:`\Ub^{n+1,(b)}`. We then compute the
   correction, :math:`(\Delta t/2)[\Sb(\Ub^{n+1,(b)}) - \Sb(\Ub^\star)]` to
   add to :math:`\Ub^{n+1,(b)}` to give us the properly time-centered source,
   and the fully updated state, :math:`\Ub^{n+1,(c)}`. 

   This correction is stored
   in the ``new_sources`` MultiFab [1]_.

   In the process of updating the sources, we update the temperature to
   make it consistent with the new state.

#. *React* :math:`\Delta t/2` [``strang_react_second_half()``]

   We do the final :math:`\dt/2` reacting on the state, begining with :math:`\Ub^{n+1,(c)}` to
   give us the final state on this level, :math:`\Ub^{n+1}`.

   This is largely the same as ``strang_react_first_half()``, but
   it does not currently fill the reactions in the ghost cells.

#. *Finalize* [``finalize_do_advance()``]

   Finalize does the following:

   A. for the momentum sources, we compute :math:`d\Sb/dt`, to use in the
      source term prediction/extrapolation for the hydrodynamic
      interface states during the next step.

   B. If we are doing the hybrid momentum algorithm, then we sync up
      the hybrid and linear momenta

A summary of which state is the input and which is updated for each of
these processes is presented below:

.. table:: update sequence of state arrays for Strang-CTU
   :align: center

   +--------------------+-----------+---------------------+---------------------+
   | *step*             | ``S_old`` | ``Sborder``         | ``S_new``           |
   +====================+===========+=====================+=====================+
   | 1. init            | input     | updated             |                     |
   +--------------------+-----------+---------------------+---------------------+
   | 2. react           |           | input / updated     |                     |
   +--------------------+-----------+---------------------+---------------------+
   | 3. old sources     |           | input               | updated             |
   +--------------------+-----------+---------------------+---------------------+
   | 4. hydro           |           | input               | updated             |
   +--------------------+-----------+---------------------+---------------------+
   | 5. clean           |           |                     | input / updated     |
   +--------------------+-----------+---------------------+---------------------+
   | 6. radial / center |           |                     | input               |
   +--------------------+-----------+---------------------+---------------------+
   | 7. correct sources |           |                     | input / updated     |
   +--------------------+-----------+---------------------+---------------------+
   | 8. react           |           |                     | input / updated     |
   +--------------------+-----------+---------------------+---------------------+


.. _sec:flow_true_sdc:

SDC Evolution
=============

The SDC evolution is selected by ``castro.time_integration_method = 2``.  It
does away with Strang splitting and instead couples the reactions and hydro
together directly.

.. note::

   At the moment, the SDC solvers do not support multilevel or AMR
   simulation.

.. note::

   The code must be compiled with ``USE_TRUE_SDC = TRUE`` to use this
   evolution type.

The SDC solver follows the algorithm detailed in :cite:`castro_sdc`.
We write our evolution equation as:

.. math::
   \frac{\partial \Ub}{\partial t} = {\bf A}(\Ub) + {\bf R}(\Ub)

where :math:`{\bf A}(\Ub) = -\nabla \cdot {\bf F}(\Ub) + {\bf S}(\Ub)`, with the 
hydrodynamic source terms, :math:`{\bf S}` grouped together with the flux divergence.

The SDC update looks at the solution a several time nodes (the number
depending on the desired temporal order of accuracy), and iteratively
updates the solution from node :math:`m` to :math:`m+1` as:

.. math::
   \begin{align}
   \avg{\Ub}^{m+1,(k+1)} = \avg{\Ub}^{m,(k+1)} &+ \Delta t \left [ \avg{{\bf A}(\Ub)}^{m,(k+1)} - \avg{{\bf A}(\Ub)}^{m,(k)} \right ] \\
                                   &+ \Delta t \left [ \avg{{\bf R}(\Ub)}^{m+1,(k+1)} - \avg{{\bf R}(\Ub)}^{m+1,(k)} \right ] \\
                                   &+ \int_{t^m}^{t^{m+1}} \left [ \avg{{\bf A}(\Ub)}^{(k)} + \avg{{\bf R}(\Ub)}^{(k)} \right ] dt
   \end{align}


.. index:: castro.sdc_order, castro.sdc_quadrature

Where :math:`k` is the iteration index.  In the SDC formalism, each
iteration gains us an order of accuracy in time, up to the order with
which we discretize the integral at the end of the above expression.
We also write the conservative state as :math:`\avg{\Ub}` to remind us
that it is the cell average and not the cell-center.  This distinction
is important when we consider the 4th order method.

In Castro, there are two parameters that together determine the number
and location of the temporal nodes, the accuracy of the integral, and
hence the overall accuracy in time: ``castro.sdc_order`` and
``castro.sdc_quadrature``. 

``castro.sdc_quadrature = 0`` uses
Gauss-Lobatto integration, which includes both the starting and ending
time in the time nodes.  This gives us the trapezoid rule for 2nd
order methods and Simpson's rule for 4th order methods.  Choosing
``castro.sdc_quadrature = 1`` uses Radau IIA integration, which includes
the ending time but not the starting time in the quadrature.


.. table:: SDC quadrature summary
   :align: center

   +--------------+---------------+---------------+-------------------+------------------+
   |``sdc_order`` |``quadrature`` |  # of         |  temporal         |  description     |
   |              |               |  time nodes   |  accuracy         |                  |
   +==============+===============+===============+===================+==================+
   |       2      |         0     |          2    |                2  | trapezoid rule   |
   +--------------+---------------+---------------+-------------------+------------------+
   |       2      |         1     |          3    |                2  | Simpson's rule   |
   +--------------+---------------+---------------+-------------------+------------------+
   |       4      |         0     |          3    |                4  | Radau 2nd order  |
   +--------------+---------------+---------------+-------------------+------------------+
   |       4      |         1     |          4    |                4  | Radau 4th order  |
   +--------------+---------------+---------------+-------------------+------------------+

The overall evolution appears as:

.. index:: k_new, A_old, A_new, R_old

#. *Initialization* (``initialize_advance``)

   We first do a ``clean_state`` on the old data (``S_old``).

   Next we create the ``MultiFab`` s that store the needed information
   at the different time nodes.  Each of the quantities below is a
   vector of size ``SDC_NODES``, whose components are the ``MultiFab``
   for that time node:


    * ``k_new`` : the current solution at this time node.

      Note that
      ``k_new[0]`` is aliased to ``S_old``, the solution at the start
      of the step, since this never changes (so long as the 0th time
      node is the start of the timestep).

    * ``A_old`` : the advective term at each time node at the old
      iteration.

    * ``A_new`` : the advective term at each time node at the current
      iteration.
    
    * ``R_old`` : the reactive source term at each time node at the old
      iteration.

#. *Advancement*

   Our iteration loop calls ``do_advance_sdc`` to update the solution through
   all the time nodes for a single iteration.

   The total number of iterations is ``castro.sdc_order`` + ``castro.sdc_extra``.

#. *Finalize*

   This clears the ``MultiFab`` s we allocated.

SDC Single Iteration Flowchart
------------------------------

.. index:: do_advance_sdc

The update through all time nodes for a single iteration is done by
``do_advance_sdc``.  The basic update appears as:

Throughout this driver we use the ``State_Type`` ``StateData`` as
storage for the current node.  In particular, we use the new time slot
in the ``StateData`` (which we refer to as ``S_new``) to allow us to
do ``FillPatch`` operations.

#. *Initialize*

   We allocate ``Sborder``.  Just like with the Strang CTU driver, we
   will use this as input into the hydrodynamics routines.

#. Loop over time nodes

   We'll use ``m`` to denote the current time node and ``sdc_iter`` to
   denote the current (0-based) iteration.  In our loop over time
   nodes, we do the following for each node:

   * Load in the starting data

     * ``S_new`` :math:`\leftarrow` ``k_new[m]``

     * ``clean_state`` on ``S_new``

     * Fill ``Sborder`` using ``S_new``

   * Construct the hydro sources and advective term

     Note: we only do this on the first time node for ``sdc_iter`` = 0, and
     we don't need to do this for the last time node on the last
     iteration.

     * Call ``do_old_sources`` filling the ``Source_Type``
       ``StateData``, ``old_source``.

     * Convert the sources to 4th order averages if needed.

     * ``sources_for_hydro`` :math:`\leftarrow` ``old_source``

     * Convert the conserved variables to primitive variables

     * Call ``construct_mol_hydro_source`` to get the advective update
       at the current time node, stored in ``A_new[m]``.
 
   * Bootstrap the first iteration.

     For the first iteration, we don't have the old iteration's
     advective and reaction terms needed in the SDC update.  So for
     the first time node (``m = 0``) on the first iteration, we do:

     * ``A_old[n]`` = ``A_old[0]``, where ``n`` loops over all time nodes.

     * Compute the reactive source using the ``m = 0`` node's state and
       store this in ``R_old[0]``.

       Then fill all other time nodes as: ``R_old[n]`` = ``R_old[0]``

    * Do the SDC update from node ``m`` to ``m+1``.

      We call ``do_sdc_update()`` to do the update in time to the next
      node.  This solves the nonlinear system (when we have reactions)
      and stores the solution in ``k_new[m+1]``.

#. Store the advective terms for the next iteration.

   Since we are done with this iteration, we do: ``A_old[n]``
   :math:`\leftarrow` ``A_new[n]``.

   We also store ``R_old`` for the next iteration.  We do this by
   calling the reaction source one last time using the data for each
   time node.

#. Store the new-time solution.

   On the last iteration, we save the solution to the ``State_Type`` ``StateData``:

   ``S_new`` :math:`\leftarrow` ``k_new[SDC_NODES-1]``

#. Call ``finalize_do_advance`` to clean up the memory.
   

Simplified-SDC Evolution
========================

The simplified SDC method uses the CTU advection solver together with
an ODE solution to update the compute advective-reacting system.  This
is selected by ``castro.time_integration_method = 3``.

We use one additional StateData type here, ``Simplified_SDC_React_Type``,
which will hold the reactive source needed by hydrodynamics.

.. note::

   The code must be compiled with ``USE_SIMPLIFIED_SDC = TRUE`` to use this
   evolution type.


We express our system as:

.. math:: \Ub_t = \mathcal{A}(\Ub) + \Rb(\Ub)

here :math:`\mathcal{A}` is the advective source, which includes both the
flux divergence and the hydrodynamic source terms (e.g. gravity):

.. math:: \mathcal{A}(\Ub) = -\nabla \cdot \Fb(\Ub) + \Sb

The simplified-SDC version of the main advance loop looks similar to the Strang CTU
version, but includes an iteration loop over the hydro, gravity, and
reaction update. So the only difference happens in step 2 of the
flowchart outlined in § \ `2 <#flow:sec:nosdc>`__. In particular this
step now proceeds as a loop over ``do_advance_ctu``.  The differences
with the Strang CTU version are highlighted below.


Note that the
radiation implicit update is not done as part of the Simplified-SDC iterations.

Simplified_SDC Hydro Advance
----------------------------

The evolution in ``do_advance`` is substantially different than the
Strang case. In particular, reactions are not evolved. Here we
summarize those differences.

#. *Initialize* [``initialize_do_advance()``]

   This is unchanged from the initialization in the CTU Strang algorithm.

#. *Construct time-level n sources and apply*
   [``construct_old_gravity()``, ``do_old_sources()``]

   Unlike the Strang case, there is no need to extrapolate source
   terms to the half-time for the prediction (the
   ``castro.source_term_predictor`` parameter), since the
   Simplified-SDC provides a natural way to approximate the
   time-centered source—we simply use the iteratively-lagged new-time
   source.  We add the corrector from the previous iteration to the
   source Multifabs before adding the current source.  The corrector
   (stored in ``source_corrector``) has the form:

   .. math::

      \Sb^\mathrm{corr} = \frac{1}{2} \left ( \Sb^{n+1,(k-1)} - S^n \right )

   where :math:`\Sb^n` does not have an iteration subscript, since we always have the
   same old time state.  

   Applying this corrector to the the source at time :math:`n`, will give
   us a source that is time-centered,

   .. math::

      {\bf S}(\Ub)^{n+1/2} = \frac{1}{2} \left ( {\bf S}(\Ub)^n + {\bf S}(\Ub)^{n+1,(k-1)} \right )

   For constructing the time-level :math:`n` source, there are no
   differences compared to the Strang algorithm.

#. *Construct the hydro update* [``construct_hydro_source()``]

   In predicting the interface states, we use an iteratively-lagged
   approximation to the reaction source on the primitive variables,
   :math:`\mathcal{I}_q^{k-1}`.  This addition is done in
   ``construct_ctu_hydro_source()`` after the source terms are
   converted to primitive variables.

   The result of this is an approximation to :math:`- [\nabla \cdot {\bf F}]^{n+1/2}` (not yet the full :math:`\mathcal{A}(\Ub)`)
   stored in ``hydro_sources``.

#. *Clean State* [``clean_state()``]

#. *Update radial data and center of mass for monopole gravity*

#. *Correct the source terms with the n+1 contribution*
   [``construct_new_gravity()``, ``do_new_sources()`` ]

#. *React* :math:`\Delta t` [``react_state()``]

   We first compute :math:`\mathcal{A}(\Ub)` using ``hydro_sources``,
   ``old_source``, and ``new_source`` via the ``sum_of_source()``
   function.  This produces an advective source of the form:
   
   .. math::

      \left [ \mathcal{A}(\Ub) \right ]^{n+1/2} = - [\nabla \cdot {\bf F}]^{n+1/2} + \frac{1}{2} (S^n + S^{n+1})

   We burn for the full :math:`\Delta t` including the advective
   update as a source, integrating

      .. math:: \frac{d\Ub}{dt} = \left [ \mathcal{A}(\Ub) \right ]^{n+1/2} + \Rb(\Ub)

   The result of evolving this equation is stored in ``S_new``.

   Note, if we do not actually burn in a zone (because we don't meet
   the thermodynamic threshold) then this step does nothing, and the
   state updated just via hydrodynamics in ``S_new`` is kept.

#. *Clean state*: This ensures that the thermodynamic state is
   valid and consistent.

#. *Construct reaction source terms*: Construct the change
   in the primitive variables due only to reactions over the
   timestep, :math:`\mathcal{I}_q^{k}`. This will be used in the next
   iteration.

#. *Finalize* [``finalize_do_advance()``]

   This differs from Strang finalization in that we do not construct
   :math:`d\Sb/dt`, but instead store the total hydrodynamical source
   term at the new time. As discussed above, this will be used in the
   next iteration to approximate the time-centered source term.

.. [1]
   The correction for gravity is slightly different since we directly compute the time-centered gravitational source term using the hydrodynamic fluxes.
