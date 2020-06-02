#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

using namespace amrex;

void apply_metric(const Box& bx, Array4<Real> const rhs, const Box& rbx, Array4<Real> const ecx,
                  const Box& xbx,
#if AMREX_SPACEDIM >= 2
                  Array4<Real> const ecy, const Box& ybx,
#endif
                  const Real* dx, const int coord_type) {

    // r-z
    if (coord_type == 1) {

        AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
            IntVect idx(D_DECL(i, j, k));

            // at centers
            if (rbx.contains(idx)) {
                Real r = (static_cast<Real>(i) + 0.5_rt) * dx[0];
                rhs(i, j, k) *= r;
            }

            // On x-edges
            if (xbx.contains(idx)) {
                Real r = static_cast<Real>(i) * dx[0];
                ecx(i, j, k) *= r;
            }

#if AMREX_SPACEDIM >= 2
            // On y-edges
            if (ybx.contains(idx)) {
                Real r = (static_cast<Real>(i) + 0.5_rt) * dx[0];
                ecy(i, j, k) *= r;
            }
#endif
        });

#ifndef AMREX_USE_CUDA
    } else {

        amrex::Print() << "Bogus coord_type in apply_metric " << coord_type << std::endl;

        amrex::Error("Error:: MGutils.cpp :: ca_apply_metric");
#endif
    }
}

void do_weight_cc(const Box& bx, Array4<Real> const cc, const Real* dx, const int coord_type) {

    // r-z
    if (coord_type == 1) {

        // At centers
        AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
            Real r = (static_cast<Real>(i) + 0.5_rt) * dx[0];
            cc(i, j, k) *= r;
        });

#ifndef AMREX_USE_CUDA
    } else {

        amrex::Print() << "Bogus coord_type in weight_cc " << coord_type << std::endl;
        amrex::Error("Error:: MGutils.cpp :: ca_weight_cc");
#endif
    }
}

void do_unweight_cc(const Box& bx, Array4<Real> const cc, const Real* dx, const int coord_type) {

    // r-z
    if (coord_type == 1) {

        // At centers
        AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
            Real r = (static_cast<Real>(i) + 0.5_rt) * dx[0];
            cc(i, j, k) /= r;
        });

#ifndef AMREX_USE_CUDA
    } else {

        amrex::Print() << "Bogus coord_type in unweight_cc " << coord_type << std::endl;
        amrex::Error("Error:: MGutils.cpp :: ca_unweight_cc");
#endif
    }
}

void do_unweight_edges(const Box& bx, Array4<Real> const ec, const int idir, const Real* dx,
                       const int coord_type) {

    // r-z
    if (coord_type == 1) {

        if (idir == 0) {

            // On x-edges
            AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
                if (i != 0) {
                    Real r = static_cast<Real>(i) * dx[0];
                    ec(i, j, k) /= r;
                }
            });

        } else {

            // On y-edges
            AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
                Real r = (static_cast<Real>(i) + 0.5_rt) * dx[0];
                ec(i, j, k) /= r;
            });
        }

#ifndef AMREX_USE_CUDA
    } else {

        amrex::Print() << "Bogus coord_type in unweight_edges " << coord_type << std::endl;
        amrex::Error("Error:: MGutils.cpp :: ca_unweight_edges");
#endif
    }
}
