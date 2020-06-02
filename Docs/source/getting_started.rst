***************
Getting Started
***************

.. note::

   Castro has two source dependencies: `AMReX <https://github.com/AMReX-Codes/amrex>`_, the adaptive mesh
   library, and `StarKiller Microphysics <https://github.com/starkiller-astro/Microphysics>`_, the collection of equations
   of state, reaction networks, and other microphysics.  The
   instructions below describe how to get these dependencies automatically
   with Castro.

Downloading the Code
====================

Castro is maintained as a repository on GitHub, and can be obtained
via standard git clone commands. First, make sure that git
is installed on your machine—we recommend version 1.7.x or higher.


#. Clone/fork the Castro repository from the AMReX-Astro GitHub
   organization, using either HTTP access::

       git clone --recursive https://github.com/AMReX-Astro/Castro.git

   or SSH access if you have an SSH key enabled with GitHub::

       git clone --recursive git@github.com:AMReX-Astro/Castro.git

   The ``--recursive`` option to ``git clone`` is used to ensure
   that all of Castro's dependencies are downloaded. Currently this
   requirement is for the AMReX mesh refinement framework, which is
   maintained in the AMReX-Codes organization on GitHub, and the
   Microphysics repository from the starkiller-astro organization.
   AMReX adds the necessary code for the driver code for the simulation,
   while Microphysics adds the equations of state, reaction
   networks, and other microphysics needed to run Castro. If you forget
   to do a recursive clone, you can rectify the situation by running
   the following from the top-level of the Castro directory::

       git submodule update --init --recursive

   .. note::

      By default, you will be on the ``master`` branch of the source.
      Development on Castro (and its primary dependencies, AMReX and
      Microphysics) is done in the ``development`` branch, so you
      should work there if you want the latest source::

        git checkout development

      The Castro team runs nightly regression testing on the
      ``development`` branch, so bugs are usually found and fixed
      relatively quickly, but it is generally less stable than staying
      on the ``master`` branch.

#. We recommend setting the ``CASTRO_HOME`` environment
   variable to point to the path name where you have put Castro.
   Add the following to your ``.bashrc``::

       export CASTRO_HOME="/path/to/Castro/"

#. You can keep the code up to date with::

       git pull --recurse-submodules

   The recommended frequency for doing this is monthly, if you are on the
   stable ``master`` branch of the code; we issue a new release of the code
   at the beginning of each month.

#. (optional, for developers) If you prefer, you can maintain AMReX and
   Microphysics as standalone repositories rather than as git submodules.
   To do so, you can clone them from GitHub using::

       git clone https://github.com/AMReX-Codes/amrex.git
       git clone https://github.com/starkiller-astro/Microphysics.git

   or via SSH as::

       git clone git@github.com:/AMReX-Codes/amrex.git
       git clone git@github.com:/starkiller-astro/Microphysics.git

   Then, set the ``AMREX_HOME`` environment variable to point to the
   ``amrex/`` directory, and the ``MICROPHYSICS_HOME`` environment
   variable to point to the ``Microphysics/`` directory. Castro will
   look there instead of in its local ``external/`` subdirectory.

#. (optional, for developers) We use ``clang-format`` to format the code 
   in the Castro repository. To automate this, we recommend setting up a 
   pre-commit hook by copying the script in the ``Utils`` directory. You 
   will also need to change the permissions of the script to be 
   executable::

        cp Utils/scripts/pre-commit .git/hooks
        chmod +x .git/hooks/pre-commit

Building the Code
=================

In Castro each different problem setup is stored in its own
sub-directory under ``Castro/Exec/``. You build the
Castro executable in the problem sub-directory. Here we’ll
build the Sedov problem:

#. From the directory in which you checked out the Castro git repo,
   type::

       cd Castro/Exec/hydro_tests/Sedov

   This will put you into a directory in which you can run the Sedov
   problem in 1-d, 2-d or 3-d.

#. In ``Sedov/``, edit the ``GNUmakefile``, and set

   * ``DIM = 2``

     This is the dimensionality—here we pick 2-d.

   * ``COMP = gnu``

     This is the set of compilers. GNUu are a good default choice
     (this will use g++ and gfortran). You can also choose ``pgi`` and
     ``intel`` for example.

     If you want to try other compilers than the GNU suite and they
     don’t work, please let us know.

   * ``DEBUG = FALSE``

     This disables debugging checks and results in a more optimized
     executable.

   * ``USE_MPI = FALSE``

     This turns off parallelization via MPI. Set it to ``TRUE`` to build
     with MPI—this requires that you have the MPI library installed on
     your machine. In this case, the build system will need to know
     about your MPI installation. This can be done by editing the
     makefiles in the AMReX tree, but the default fallback is to look
     for the standard MPI wrappers (e.g. ``mpic++`` and ``mpif90``) to do
     the build.

#. Now type ``make``.

   The resulting executable will look something like
   ``Castro2d.gnu.ex``, which means this is a 2-d version
   of the code compiled with ``COMP = gnu``.

Running the Code
================

#. Castro takes an input file that overrides the runtime parameter defaults.
   The code is run as::

       ./Castro2d.gnu.ex inputs.2d.cyl_in_cartcoords

   This will run the 2-d cylindrical Sedov problem in Cartesian
   (:math:`x`-:math:`y` coordinates). You can see other possible
   options, which should be clear by the names of the inputs files.

#. You will notice that running the code generates directories that
   look like ``plt00000/``, ``plt00020/``, etc, and ``chk00000/``,
   ``chk00020/``, etc. These are “plotfiles” and “checkpoint”
   files. The plotfiles are used for visualization, the checkpoint
   files are used for restarting the code.

Visualization of the Results
============================

There are several options for visualizing the data. The popular VisIt
package supports the AMReX file format natively, as does the yt python
package [2]_. The standard tool used within the AMReX-community is
Amrvis, which we demonstrate here. Amrvis is available on github.

Please know that we do have a number of conversion routines to other
formats (such as matlab), but it is hard to describe them all. If you
would like to display the data in another format, please let us know
(again, asalmgren@lbl.gov) and we will point you to whatever we have


Amrvis
^^^^^^

#. Get Amrvis::

       git clone https://github.com/AMReX-Codes/Amrvis

   Then cd into ``Amrvis/``, edit the ``GNUmakefile`` there
   to set ``DIM = 2``, and again set ``COMP`` to compilers that
   you have. Leave ``DEBUG = FALSE``.

   Type ``make`` to build, resulting in an executable that
   looks like ``amrvis2d...ex``.

   If you want to build amrvis with ``DIM = 3``, you must first
   download and build volpack::

       git clone https://ccse.lbl.gov/pub/Downloads/volpack.git

   Then cd into ``volpack/`` and type ``make``.

   Note: Amrvis requires the OSF/Motif libraries and headers. If you
   don’t have these you will need to install the development version
   of motif through your package manager.  On most Linux
   distributions, the motif library is provided by the openmotif
   package, and its header files (like ``Xm.h``) are provided by
   openmotif-devel. If those packages are not installed, then use the
   package management tool to install them, which varies from
   distribution to distribution, but is straightforward.  ``lesstif``
   gives some functionality and will allow you to build the Amrvis
   executable, but Amrvis may not run properly.

   You may then want to create an alias to amrvis2d, for example::

       alias amrvis2d=/tmp/Amrvis/amrvis2d...ex

   where ``/tmp/Amrvis/amrvis2d...ex`` is the full path and name of
   the Amrvis executable.

#. Configure Amrvis:

   Copy the ``amrvis.defaults`` file to your home directory (you can
   rename it to ``.amrvis.defaults`` if you wish). Then edit the
   file, and change the palette line to point to the full
   path/filename of the ``Palette`` file that comes with Amrvis.

#. Visualize:

   Return to the ``Castro/Exec/hydro_tests/Sedov`` directory. You should
   have a number of output files, including some in the form ``pltXXXXX``,
   where XXXXX is a number corresponding to the timestep the file
   was output.

   ``amrvis2d filename`` to see a single plotfile, or ``amrvis2d -a
   plt*``, which will animate the sequence of plotfiles.

   Try playing around with this—you can change which variable you are
   looking at, select a region and click “Dataset” (under View) in
   order to look at the actual numbers, etc. You can also export the
   pictures in several different formats under "File/Export".

   Some users have found that Amrvis does not work properly under X
   with the proprietary Nvidia graphics driver. A fix for this is
   provided in the FAQ (§ :ref:`ch:faq:vis`)—this is due
   to the default behavior of the DAC in mappuing colors.


yt
^^

yt is a great alternative to using Amrvis for visualization,
and understands Castro plotfiles well.  Install yt following
their instructions: `Getting yt <https://yt-project.org/#getyt>`_ .

You should be able to read in your plotfiles using ``yt.load()``
and do any of the plots described in the `yt Cookbook <https://yt-project.org/doc/cookbook/index.html>`_ .


Other Distributed Problem Setups
================================

There are a number of standard problem setups that come with Castro.
These can be used as a starting point toward writing your own setup.
We organize these into subdirectories by broad type (radiation, hydro,
gravity, etc.): The standard categories and *some* of the included
problems are:

* ``gravity_tests``:

   * ``DustCollapse``:

     A pressureless cloud collapse that is a standard test problem for
     gravity. An analytic solution that describes the radius of the
     sphere as a function of time is found in Colgate and
     White :cite:`colgwhite`. This problem is also found
     in the FLASH User’s Guide.

   * ``hydrostatic_adjust``:

     Model a 1-d stellar atmosphere (plane-parallel or
     spherical/self-gravitating) and dump energy in via an analytic
     heat source and watch the atmosphere’s hydrostatic state adjust
     in response. This is the counterpart to the Maestro
     ``test_basestate`` unit test.

* ``hydro_tests``:

   * ``double_bubble``:

     Initialize 1 or 2 bubbles in a stratified atmosphere (isothermal
     or isentropic) and allow for the bubbles to have the same or a
     different :math:`\gamma` from one another / the background
     atmosphere.  This uses the multigamma EOS.

     An analogous problem is implemented in Maestro.

   * ``HCBubble``:

   * ``KH``:

     A Kelvin-Helmholtz shear instability problem.

   * ``oddeven``:

     A grid-aligned shock hitting a very small density perturbation.
     This demonstrates the odd-even decoupling problem discussed in
     :cite:`quirk1997`. This setup serves to test the
     castro.hybrid_riemann option to hydrodynamics.

   * ``reacting_bubble``:

     A reacting bubble in a stratified white dwarf atmosphere. This
     problem was featured in the Maestro reaction
     paper :cite:`maestro:III`.

   * ``RT``:

     A single-model Rayleigh-Taylor instability problem.

   * ``RT_particles``:

   * ``Sedov``:

     The standard Sedov-Taylor blast wave problem. This setup was used
     in the first Castro paper :cite:`castro_I`.

   * ``Sod``:

     A one-dimensional shock tube setup, including the classic Sod
     problem. This setup was used in the original Castro paper.

   * ``Sod_stellar``:

     A version of the Sod shock tube for the general stellar equation
     of state. This setup and the included inputs files was used
     in :cite:`zingalekatz`.

   * ``toy_convect``:

     A simple nova-like convection problem with an external heating
     source. This problem shows how to use the model parser to
     initialize a 1-d atmosphere on the Castro grid, incorporate a
     custom tagging routine, sponge the fluid above the atmosphere,
     and write a custom diagnostics routine.

     A Maestro version of this problem setup also exists.

* ``radiation_tests``:

* ``science``:

* ``unit_tests``:

.. [1]
   Note: previously the radiation
   solver was distributed separately as ``CastroRadiation.git``,
   but this has been merged into the main Castro respository

.. [2]
   Each of these will recognize it as the
   BoxLib format.
