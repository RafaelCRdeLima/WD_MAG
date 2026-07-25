#ifndef PROBLEM_DERIVE_H
#define PROBLEM_DERIVE_H

#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

using namespace amrex;

#ifdef __cplusplus
extern "C"
{
#endif

    // Magnetic energy density (erg/cm^3, Gaussian-equivalent: Castro
    // stores B' = B/sqrt(4 pi), and B'^2/2 = B_gauss^2/(8 pi) -- see
    // problem_initialize_mhd_data.H's comment on the same factor).
    void ca_deremag(const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                    const FArrayBox& datfab, const Geometry& geom,
                    Real time, const int* bcrec, int level);

    // Toroidal magnetic energy density -- B_phi^2/2 in a CYLINDRICAL
    // decomposition around the z-axis. This background star is
    // non-rotating (Step 1 requires field-free, non-rotating), so the
    // z-axis is an arbitrary bookkeeping choice, not a dynamically
    // preferred axis -- see docs/teoria.md Sec 6.3/6.4 for the caveat.
    void ca_deretor(const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                    const FArrayBox& datfab, const Geometry& geom,
                    Real time, const int* bcrec, int level);

#ifdef __cplusplus
}
#endif

#endif
