/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef AXIDILATONPOTENTIAL_HPP_
#define AXIDILATONPOTENTIAL_HPP_

#include "simd.hpp"

class AxiDilatonPotential
{
  public:
    struct params_t
    {
        double scalar_mass;
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    AxiDilatonPotential(params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_modulus_phi_squared,
                           data_t &dVdmodulus_phi_squared,
                           const vars_t<data_t> &vars,
                           double gamma_squared_coeff) const
    {
        // First calculate |phi|^2
        data_t modulus_phi_squared =
            vars.phi_Re * vars.phi_Re + vars.phi_Im * vars.phi_Im;
        // The potential value at phi
        V_of_modulus_phi_squared =
            0.5 * pow(m_params.scalar_mass, 2.0) * modulus_phi_squared /
            (1.0 - gamma_squared_coeff * modulus_phi_squared / 4.0);

        // The potential gradient at phi
        // THIS IS WRONG!!
        // dVdmodulus_phi_squared = pow(m_params.scalar_mass, 2.0) /(2.0 *
        // pow(1.0-gamma_squared_coeff * modulus_phi_squared / 4.0,2.0));
        dVdmodulus_phi_squared = pow(m_params.scalar_mass, 2.0) / 2.0;
    }
};

#endif /* AXIDILATONPOTENTIAL_HPP_ */
