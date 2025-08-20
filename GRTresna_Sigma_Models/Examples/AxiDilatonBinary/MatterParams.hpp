/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef MATTERPARAMS_HPP_
#define MATTERPARAMS_HPP_

#include "GRParmParse.hpp"
#include "REAL.H"
#include "RealVect.H"

namespace MatterParams
{

struct params_t
{
    //Real phi_Re_0;
    Real dphi;
    Real dphi_length;
    Real pi_Re_0;
    Real dpi;
    Real dpi_length;
    // Real phi_Im_0;
    // Real dphi_Im;
    // Real dphi_Im_length;
    // Real pi_Im_0;
    // Real dpi_Im;
    // Real dpi_Im_length;
    Real scalar_mass;
    Real gamma_squared_coeff;
    Real amplitude;
    Real width;
    Real r0;
    Real kerr_mass;
    RealVect kerr_spin;
};

inline void read_params(GRParmParse &pp, params_t &matter_params)
{
    // pp.get("phi_Re_0", matter_params.phi_Re_0);
    // pp.get("dphi", matter_params.dphi);
    // pp.get("dphi_length", matter_params.dphi_length);
    // pp.get("pi_Re_0", matter_params.pi_Re_0);
    // pp.get("dpi", matter_params.dpi);
    // pp.get("dpi_length", matter_params.dpi_length);

    // pp.get("phi_Im_0", matter_params.phi_Im_0);
    // pp.get("dphi_Im", matter_params.dphi_Im);
    // pp.get("dphi_Im_length", matter_params.dphi_Im_length);
    // pp.get("pi_Im_0", matter_params.pi_Im_0);
    // pp.get("dpi_Im", matter_params.dpi_Im);
    // pp.get("dpi_Im_length", matter_params.dpi_Im_length);
    std::vector<double> temp_spin1(3);

    pp.get("scalar_mass", matter_params.scalar_mass);
    pp.get("scalar_amplitude", matter_params.amplitude);
    pp.get("gamma_squared_coeff", matter_params.gamma_squared_coeff);
    pp.get("scalar_width", matter_params.width);
    pp.get("scalar_r0", matter_params.r0);
    pp.get("bh1_bare_mass", matter_params.kerr_mass);
    std::vector<double> temp_spin1b(3);
    pp.getarr("bh1_spin", temp_spin1b, 0, 3);

    for (int idir = 0; idir < 3; idir++)
    {
        matter_params.kerr_spin[idir] = temp_spin1b[idir];
    }
}

}; // namespace MatterParams

#endif
