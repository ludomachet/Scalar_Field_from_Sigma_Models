/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(AXIDILATON_HPP_)
#error "This file should only be included through AxiDilaton.hpp"
#endif

#ifndef AXIDILATON_IMPL_HPP_
#define AXIDILATON_IMPL_HPP_

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> AxiDilaton<potential_t>::compute_emtensor(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const Tensor<2, data_t> &h_UU, const Tensor<3, data_t> &chris_ULL) const
{
    emtensor_t<data_t> out;

    const data_t chi_regularised = simd_max(1e-6, vars.chi);

    // Copy the field vars into SFObject
    CSFObject<data_t> vars_csf;
    vars_csf.phi_Re = vars.phi_Re;
    vars_csf.phi_Im = vars.phi_Im;
    vars_csf.Pi_Re = vars.Pi_Re;
    vars_csf.Pi_Im = vars.Pi_Im;

    // call the function which computes the em tensor excluding the potential
    emtensor_excl_potential(out, vars, vars_csf, d1.phi_Re, d1.phi_Im, h_UU,
                            chris_ULL);

    // set the default potential values
    data_t V_of_modulus_phi_squared = 0.0;
    data_t dVdmodulus_phi_squared = 0.0;

    // compute potential and add constributions to EM Tensor
    my_potential.compute_potential(V_of_modulus_phi_squared,
                                   dVdmodulus_phi_squared, vars,
                                   m_gamma_squared_coeff);

    out.rho += V_of_modulus_phi_squared;
    out.S += -3. * V_of_modulus_phi_squared; // trace of hij I guess
    FOR2(i, j)
    {
        out.Sij[i][j] +=
            -1. * vars.h[i][j] * V_of_modulus_phi_squared / chi_regularised;
    }

    data_t BH_cutoff = 1.0 / (1.0 + exp(-200.0 * (vars.chi - 0.09)));

    out.rho *= BH_cutoff;
    out.S *= BH_cutoff;
    FOR2(i, j) { out.Sij[i][j] *= BH_cutoff; }

    return out;
}

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
void AxiDilaton<potential_t>::emtensor_excl_potential(
    emtensor_t<data_t> &out, const vars_t<data_t> &vars,
    const CSFObject<data_t> &vars_csf, const Tensor<1, data_t> &d1_phi_Re,
    const Tensor<1, data_t> &d1_phi_Im, const Tensor<2, data_t> &h_UU,
    const Tensor<3, data_t> &chris_ULL) const
{
    // initialise some useful quantities
    data_t modulus_d1_phi_squared = 0.;
    data_t modulus_phi_squared =
        vars_csf.phi_Re * vars_csf.phi_Re + vars_csf.phi_Im * vars_csf.phi_Im;
    data_t modulus_Pi_squared =
        vars_csf.Pi_Re * vars_csf.Pi_Re + vars_csf.Pi_Im * vars_csf.Pi_Im;
    // data_t KAxiDilaton_phi_squared =
    //    0.5 / pow(1- m_gamma_squared_coeff * modulus_phi_squared / 4., 2.0);
    data_t KAxiDilaton_phi_squared =
        1.0 / pow(4.0 - m_gamma_squared_coeff * modulus_phi_squared, 2.0);
    const data_t chi_regularised = simd_max(1e-6, vars.chi);
    FOR2(i, j)
    {
        modulus_d1_phi_squared +=
            chi_regularised * h_UU[i][j] *
            (d1_phi_Re[i] * d1_phi_Re[j] + d1_phi_Im[i] * d1_phi_Im[j]);
    }

    out.rho = 8. * KAxiDilaton_phi_squared *
              (modulus_Pi_squared + modulus_d1_phi_squared);

    // S_i (note lower index) = n^a T_a0
    FOR1(i)
    {
        out.Si[i] =
            -16. * KAxiDilaton_phi_squared *
            (vars_csf.Pi_Re * d1_phi_Re[i] + vars_csf.Pi_Im * d1_phi_Im[i]);
    }

    FOR2(i, j)
    {
        out.Sij[i][j] =
            8. * KAxiDilaton_phi_squared *
            ((2. *
              (d1_phi_Re[i] * d1_phi_Re[j] + d1_phi_Im[i] * d1_phi_Im[j])) +
             (vars.h[i][j] * (modulus_Pi_squared - modulus_d1_phi_squared) /
              chi_regularised));
    }

    out.S = chi_regularised * TensorAlgebra::compute_trace(out.Sij, h_UU);
}

// Adds in the RHS for the matter vars
template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void AxiDilaton<potential_t>::add_matter_rhs(
    rhs_vars_t<data_t> &total_rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec) const
{
    CSFObject<data_t> rhs_csf;
    // advection terms
    CSFObject<data_t> advec_csf;
    advec_csf.phi_Re = advec.phi_Re;
    advec_csf.phi_Im = advec.phi_Im;
    advec_csf.Pi_Re = advec.Pi_Re;
    advec_csf.Pi_Im = advec.Pi_Im;

    // the vars
    CSFObject<data_t> vars_csf;
    vars_csf.phi_Re = vars.phi_Re;
    vars_csf.phi_Im = vars.phi_Im;
    vars_csf.Pi_Re = vars.Pi_Re;
    vars_csf.Pi_Im = vars.Pi_Im;

    const data_t chi_regularised = simd_max(1e-6, vars.chi);

    // call the function for the rhs excluding the potential
    matter_rhs_excl_potential(rhs_csf, vars, vars_csf, d1, d1.phi_Re, d1.phi_Im,
                              d2.phi_Re, d2.phi_Im, advec_csf);

    // set the default potential values
    data_t V_of_modulus_phi_squared = 0.0;
    data_t dVdmodulus_phi_squared = 0.0;

    // compute potential and add constributions to EM Tensor
    my_potential.compute_potential(V_of_modulus_phi_squared,
                                   dVdmodulus_phi_squared, vars,
                                   m_gamma_squared_coeff);

    // adjust RHS for the potential term
    total_rhs.phi_Re = rhs_csf.phi_Re;
    total_rhs.phi_Im = rhs_csf.phi_Im;
    total_rhs.Pi_Re = rhs_csf.Pi_Re - vars.lapse * dVdmodulus_phi_squared * 2. *
                                          vars_csf.phi_Re; // Check factor of 2
    total_rhs.Pi_Im = rhs_csf.Pi_Im - vars.lapse * dVdmodulus_phi_squared * 2. *
                                          vars_csf.phi_Im;

    data_t BH_cutoff = 1.0 / (1.0 + exp(-200.0 * (vars.chi - 0.9)));

    // total_rhs.phi_Re *= BH_cutoff;
    // total_rhs.phi_Im *= BH_cutoff;
    // total_rhs.Pi_Re *= BH_cutoff;
    // total_rhs.Pi_Im *= BH_cutoff;
}

// the RHS excluding the potential terms
template <class potential_t>
template <class data_t, template <typename> class vars_t>
void AxiDilaton<potential_t>::matter_rhs_excl_potential(
    CSFObject<data_t> &rhs_csf, const vars_t<data_t> &vars,
    const CSFObject<data_t> &vars_csf, const vars_t<Tensor<1, data_t>> &d1,
    const Tensor<1, data_t> &d1_phi_Re, const Tensor<1, data_t> &d1_phi_Im,
    const Tensor<2, data_t> &d2_phi_Re, const Tensor<2, data_t> &d2_phi_Im,
    const CSFObject<data_t> &advec_csf) const
{
    using namespace TensorAlgebra;

    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);
    data_t modulus_phi_squared =
        vars_csf.phi_Re * vars_csf.phi_Re + vars_csf.phi_Im * vars_csf.phi_Im;

    const data_t chi_regularised = simd_max(1e-6, vars.chi);

    data_t KAxiDilaton_phi_squared_sqrt =
        1.0 / (4.0 - m_gamma_squared_coeff * modulus_phi_squared);

    // evolution equations for scalar field and (minus) its conjugate momentum
    rhs_csf.phi_Re = advec_csf.phi_Re + vars.lapse * vars_csf.Pi_Re;
    rhs_csf.phi_Im = advec_csf.phi_Im + vars.lapse * vars_csf.Pi_Im;
    rhs_csf.Pi_Re =
        advec_csf.Pi_Re + vars.lapse * vars_csf.Pi_Re * vars.K +
        vars.lapse * 2. * m_gamma_squared_coeff * KAxiDilaton_phi_squared_sqrt *
            (-2. * vars_csf.Pi_Im * vars_csf.Pi_Re * vars_csf.phi_Im +
             vars_csf.Pi_Im * vars_csf.Pi_Im * vars_csf.phi_Re -
             vars_csf.Pi_Re * vars_csf.Pi_Re * vars_csf.phi_Re);
    rhs_csf.Pi_Im =
        advec_csf.Pi_Im + vars.lapse * vars_csf.Pi_Im * vars.K +
        vars.lapse * 2. * m_gamma_squared_coeff * KAxiDilaton_phi_squared_sqrt *
            (-2. * vars_csf.Pi_Im * vars_csf.Pi_Re * vars_csf.phi_Re -
             vars_csf.Pi_Im * vars_csf.Pi_Im * vars_csf.phi_Im +
             vars_csf.Pi_Re * vars_csf.Pi_Re * vars_csf.phi_Im);

    FOR(k)
    {
        rhs_csf.Pi_Re +=
            -vars.lapse * chi_regularised * chris.contracted[k] * d1_phi_Re[k];
        rhs_csf.Pi_Im +=
            -vars.lapse * chi_regularised * chris.contracted[k] * d1_phi_Im[k];
    }

    FOR(k, l)
    {
        rhs_csf.Pi_Re +=
            -h_UU[k][l] * (-chi_regularised * d1.lapse[k] * d1_phi_Re[l] +
                           vars.lapse * (0.5 * d1.chi[k] * d1_phi_Re[l] -
                                         chi_regularised * d2_phi_Re[k][l]));
        rhs_csf.Pi_Im +=
            -h_UU[k][l] * (-chi_regularised * d1.lapse[k] * d1_phi_Im[l] +
                           vars.lapse * (0.5 * d1.chi[k] * d1_phi_Im[l] -
                                         chi_regularised * d2_phi_Im[k][l]));
    }

    FOR(k, l)
    {
        rhs_csf.Pi_Re += 2. * m_gamma_squared_coeff *
                         KAxiDilaton_phi_squared_sqrt * h_UU[k][l] *
                         vars.lapse * vars.chi *
                         (-vars_csf.phi_Re * d1_phi_Im[k] * d1_phi_Im[l] +
                          2. * vars_csf.phi_Im * d1_phi_Im[k] * d1_phi_Re[l] +
                          vars_csf.phi_Re * d1_phi_Re[k] * d1_phi_Re[l]);
        rhs_csf.Pi_Im += 2. * m_gamma_squared_coeff *
                         KAxiDilaton_phi_squared_sqrt * h_UU[k][l] *
                         vars.lapse * vars.chi *
                         (vars_csf.phi_Im * d1_phi_Im[k] * d1_phi_Im[l] +
                          2. * vars_csf.phi_Re * d1_phi_Im[k] * d1_phi_Re[l] -
                          vars_csf.phi_Im * d1_phi_Re[k] * d1_phi_Re[l]);
    }
}

#endif /* AXIDILATON_IMPL_HPP_ */
