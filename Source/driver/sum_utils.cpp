#include <iomanip>

#include <Castro.H>
#include <Castro_F.H>

#ifdef GRAVITY
#include <Gravity.H>
#include <Gravity_F.H>
#endif

using namespace amrex;

Real Castro::volWgtSum(const std::string& name, Real time, bool local, bool finemask) {
    BL_PROFILE("Castro::volWgtSum()");

    auto mf = derive(name, time, 0);

    BL_ASSERT(mf);

    if (level < parent->finestLevel() && finemask) {
        const MultiFab& mask = getLevel(level + 1).build_fine_mask();
        MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
    }

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const& fab = (*mf).array(mfi);
        auto const& vol = volume.array(mfi);

        const Box& box = mfi.tilebox();

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //

        reduce_op.eval(box, reduce_data,
                       [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept->ReduceTuple {
                           return {fab(i, j, k) * vol(i, j, k)};
                       });
    }

    ReduceTuple hv = reduce_data.value();
    Real sum = amrex::get<0>(hv);

    if (!local) ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real Castro::volWgtSquaredSum(const std::string& name, Real time, bool local) {
    BL_PROFILE("Castro::volWgtSquaredSum()");

    auto mf = derive(name, time, 0);

    BL_ASSERT(mf);

    if (level < parent->finestLevel()) {
        const MultiFab& mask = getLevel(level + 1).build_fine_mask();
        MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
    }

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const& fab = (*mf).array(mfi);
        auto const& vol = volume.array(mfi);

        const Box& box = mfi.tilebox();

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //

        reduce_op.eval(box, reduce_data,
                       [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept->ReduceTuple {
                           return {fab(i, j, k) * fab(i, j, k) * vol(i, j, k)};
                       });
    }

    ReduceTuple hv = reduce_data.value();
    Real sum = amrex::get<0>(hv);

    if (!local) ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real Castro::locWgtSum(const std::string& name, Real time, int idir, bool local) {
    BL_PROFILE("Castro::locWgtSum()");

    auto mf = derive(name, time, 0);

    BL_ASSERT(mf);

    if (level < parent->finestLevel()) {
        const MultiFab& mask = getLevel(level + 1).build_fine_mask();
        MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
    }

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    auto dx = geom.CellSizeArray();
    auto problo = geom.ProbLoArray();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const& fab = (*mf).array(mfi);

        const Box& box = mfi.tilebox();

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //

        reduce_op.eval(box, reduce_data,
                       [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept->ReduceTuple {
                           Real loc[3];

                           loc[0] = problo[0] + (0.5_rt + i) * dx[0];

#if AMREX_SPACEDIM >= 2
                           loc[1] = problo[1] + (0.5_rt + j) * dx[1];
#else
                           loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
                           loc[2] = problo[2] + (0.5_rt + k) * dx[2];
#else
                           loc[2] = 0.0_rt;
#endif

                           Real ds;

                           if (idir == 0) {  // sum(mass * x)
                               ds = fab(i, j, k) * loc[0];
                           } else if (idir == 1) {  // sum(mass * y)
                               ds = fab(i, j, k) * loc[1];
                           } else {  // sum(mass * z)
                               ds = fab(i, j, k) * loc[2];
                           }

                           return {ds};
                       });
    }

    ReduceTuple hv = reduce_data.value();
    Real sum = amrex::get<0>(hv);

    if (!local) ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real Castro::volProductSum(const std::string& name1, const std::string& name2, Real time,
                           bool local) {
    BL_PROFILE("Castro::volProductSum()");

    auto mf1 = derive(name1, time, 0);
    auto mf2 = derive(name2, time, 0);

    BL_ASSERT(mf1);
    BL_ASSERT(mf2);

    if (level < parent->finestLevel()) {
        const MultiFab& mask = getLevel(level + 1).build_fine_mask();
        MultiFab::Multiply(*mf1, mask, 0, 0, 1, 0);
        MultiFab::Multiply(*mf2, mask, 0, 0, 1, 0);
    }

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*mf1, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const& fab1 = (*mf1).array(mfi);
        auto const& fab2 = (*mf2).array(mfi);
        auto const& vol = volume.array(mfi);

        const Box& box = mfi.tilebox();

        reduce_op.eval(box, reduce_data,
                       [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept->ReduceTuple {
                           return {fab1(i, j, k) * fab2(i, j, k) * vol(i, j, k)};
                       });
    }

    ReduceTuple hv = reduce_data.value();
    Real sum = amrex::get<0>(hv);

    if (!local) ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real Castro::locSquaredSum(const std::string& name, Real time, int idir, bool local) {
    BL_PROFILE("Castro::locSquaredSum()");

    auto mf = derive(name, time, 0);

    BL_ASSERT(mf);

    if (level < parent->finestLevel()) {
        const MultiFab& mask = getLevel(level + 1).build_fine_mask();
        MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
    }

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    auto dx = geom.CellSizeArray();
    auto problo = geom.ProbLoArray();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto const& fab = (*mf).array(mfi);

        const Box& box = mfi.tilebox();

        reduce_op.eval(
            box, reduce_data, [=] AMREX_GPU_HOST_DEVICE(int i, int j, int k) noexcept->ReduceTuple {
                Real loc[3];

                loc[0] = problo[0] + (0.5_rt + i) * dx[0];

#if AMREX_SPACEDIM >= 2
                loc[1] = problo[1] + (0.5_rt + j) * dx[1];
#else
            loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
                loc[2] = problo[2] + (0.5_rt + k) * dx[2];
#else
            loc[2] = 0.0_rt;
#endif

                Real ds;

                if (idir == 0) {  // sum(mass * x^2)
                    ds = fab(i, j, k) * loc[0] * loc[0];
                } else if (idir == 1) {  // sum(mass * y^2)
                    ds = fab(i, j, k) * loc[1] * loc[1];
                } else if (idir == 2) {  // sum(mass * z^2)
                    ds = fab(i, j, k) * loc[2] * loc[2];
                } else {  // sum(mass * r^2)
                    ds = fab(i, j, k) * (loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2]);
                }

                return {ds};
            });
    }

    ReduceTuple hv = reduce_data.value();
    Real sum = amrex::get<0>(hv);

    if (!local) ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}
