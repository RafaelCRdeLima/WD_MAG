#include <AMReX_REAL.H>

#include <Derive.H>
#include <Castro.H>
#include <Problem_Derive.H>

using namespace amrex;

#ifdef __cplusplus
extern "C"
{
#endif

    void ca_deremag(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                    const FArrayBox& datfab, const Geometry& /*geom*/,
                    Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {
        auto const dat = datfab.array();
        auto const der = derfab.array();

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real Bx_c = 0.5_rt * (dat(i,j,k,0) + dat(i+1,j,k,0));
            Real By_c = 0.5_rt * (dat(i,j,k,1) + dat(i,j+1,k,1));
            Real Bz_c = 0.5_rt * (dat(i,j,k,2) + dat(i,j,k+1,2));

            der(i,j,k,0) = 0.5_rt * (Bx_c*Bx_c + By_c*By_c + Bz_c*Bz_c);
        });
    }

    void ca_deretor(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                    const FArrayBox& datfab, const Geometry& geom,
                    Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {
        auto const dat = datfab.array();
        auto const der = derfab.array();

        auto dx = geom.CellSizeArray();
        auto problo = geom.ProbLoArray();

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real Bx_c = 0.5_rt * (dat(i,j,k,0) + dat(i+1,j,k,0));
            Real By_c = 0.5_rt * (dat(i,j,k,1) + dat(i,j+1,k,1));

            Real x = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);
            Real y = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
            Real varpi = std::sqrt(x*x + y*y);

            Real B_tor = (varpi > 0.0_rt) ? (-y*Bx_c + x*By_c) / varpi : 0.0_rt;

            der(i,j,k,0) = 0.5_rt * B_tor*B_tor;
        });
    }

#ifdef __cplusplus
}
#endif
