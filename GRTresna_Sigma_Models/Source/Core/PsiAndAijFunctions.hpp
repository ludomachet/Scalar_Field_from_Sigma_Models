/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef PSIANDAIJFUNCTIONS_HPP_
#define PSIANDAIJFUNCTIONS_HPP_

#include "DerivativeOperators.hpp"
#include "FArrayBox.H"
#include "GRParmParse.hpp"
#include "Interval.H"
#include "REAL.H"
#include "RealVect.H"
#include "TensorAlgebra.hpp"
#include "UsingNamespace.H"

class PsiAndAijFunctions
{
  public:
    struct params_t
    {
        Real bh1_bare_mass;
        Real bh2_bare_mass;
        RealVect bh1_spin;
        RealVect bh2_spin;
        RealVect bh1_momentum;
        RealVect bh2_momentum;
        RealVect bh1_offset;
        RealVect bh2_offset;
        bool use_compact_Vi_ansatz;
    };

    static void read_params(GRParmParse &pp, params_t &a_psi_and_Aij_params);

    explicit PsiAndAijFunctions(params_t a_psi_and_Aij_params)
        : m_psi_and_Aij_params(a_psi_and_Aij_params)
    {
    }

    void get_bh_coords(Real &bh_radius, RealVect &loc_bh, const RealVect &loc,
                       const RealVect &bh_offset);

    Real compute_bowenyork_psi(const RealVect &loc);

    Real compute_bowenyork_psiBL(const RealVect &loc);

    void compute_bowenyork_Aij(Tensor<2, Real> &Aij, // const IntVect &iv,
                               const RealVect &loc);

    void compute_ctt_Aij(Tensor<2, Real> &Aij,
                         const FArrayBox &multigrid_vars_box, const IntVect &iv,
                         const RealVect &a_dx, const RealVect &loc) const;

    params_t m_psi_and_Aij_params;
};

#endif /* PSIANDAIJFUNCTIONS_HPP_ */
