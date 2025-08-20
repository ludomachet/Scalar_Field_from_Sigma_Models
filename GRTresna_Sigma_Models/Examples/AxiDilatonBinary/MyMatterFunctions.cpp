/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#include "AxiDilaton.hpp"
#include "SphericalHarmonics.hpp"

Real AxiDilaton::my_potential_function(const Real &phi_Re_here, const Real &phi_Im_here) const
{
    return 0.5 * pow(m_matter_params.scalar_mass, 2.0) * 
            (phi_Re_here * phi_Re_here + phi_Im_here * phi_Im_here) / (1.0 -
            m_matter_params.gamma_squared_coeff * (phi_Re_here * phi_Re_here + phi_Im_here * phi_Im_here) / 4.0);
}

Real AxiDilaton::my_phi_Re_function(const RealVect &loc) const
{
    SphericalHarmonics::Y_lm_t<Real> Y_lm011 =
        SphericalHarmonics::spin_Y_lm(loc[0], loc[1], loc[2], 0, 1, 1);
    SphericalHarmonics::Y_lm_t<Real> Y_lm020 =
        SphericalHarmonics::spin_Y_lm(loc[0], loc[1], loc[2], 0, 2, 0);
    SphericalHarmonics::Y_lm_t<Real> Y_lm022 =
        SphericalHarmonics::spin_Y_lm(loc[0], loc[1], loc[2], 0, 2, 2);
    SphericalHarmonics::Y_lm_t<Real> Y_lm02_2 =
        SphericalHarmonics::spin_Y_lm(loc[0], loc[1], loc[2], 0, 2, -2); 

    Real rr = sqrt(loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2]);

    return m_matter_params.amplitude/  (1. + exp(-(480.-rr)/10.));;
    //return m_matter_params.amplitude * Y_lm011.Real / (1. + rr / 350.);
    //return m_matter_params.amplitude * Y_lm011.Real /(1. + exp(-50. * (rr - .5)));
    //return m_matter_params.amplitude * (1.0-exp(-pow((rr - m_matter_params.r0) / m_matter_params.width, 6.0)))* (Y_lm011.Real);
    //return m_matter_params.amplitude * (exp(-pow((rr - m_matter_params.r0) / m_matter_params.width, 2.0)))* (Y_lm011.Real);
    //return m_matter_params.amplitude * (1.0 + exp(-pow((rr - m_matter_params.r0) / m_matter_params.width, 2.0))) * Y_lm011.Real;
    //return m_matter_params.amplitude * exp(-pow((rr - m_matter_params.r0) / m_matter_params.width, 2.0))*(Y_lm011.Real)/(1.+ exp(-10.*(rr - 2.)));
}

Real AxiDilaton::my_phi_Im_function(const RealVect &loc) const
{
    return 0.;
}

Real AxiDilaton::my_Pi_Re_function(const RealVect &loc) const
{
    return 0.;
}

Real AxiDilaton::my_Pi_Im_function(const RealVect &loc) const
{
    return - m_matter_params.scalar_mass * my_phi_Re_function(loc);
}
