/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#include "AxiDilaton.hpp"
#include "DerivativeOperators.hpp"
#include "EMTensor.hpp"
#include "FArrayBox.H"
#include "GRParmParse.hpp"
#include "Grids.hpp"
#include "IntVect.H"
#include "LevelData.H"
#include "MultigridVariables.hpp" 
#include "PsiAndAijFunctions.hpp"
#include "REAL.H"
#include "RealVect.H"
#include "Tensor.hpp"

void AxiDilaton::initialise_matter_vars(LevelData<FArrayBox> &a_multigrid_vars,
                                         const RealVect &a_dx) const
{
    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    DataIterator dit = a_multigrid_vars.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        // These contain the vars in the boxes, set them all to zero
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];

        // Iterate over the box and set non zero comps
        Box ghosted_box = multigrid_vars_box.box();
        BoxIterator bit(ghosted_box);
        for (bit.begin(); bit.ok(); ++bit)
        {

            // work out location on the grid
            IntVect iv = bit();
            RealVect loc;
            Grids::get_loc(loc, iv, a_dx, center); 

            multigrid_vars_box(iv, c_phi_Re_0) = my_phi_Re_function(loc);
            multigrid_vars_box(iv, c_phi_Im_0) = my_phi_Im_function(loc);
            multigrid_vars_box(iv, c_Pi_Re_0) = my_Pi_Re_function(loc);
            multigrid_vars_box(iv, c_Pi_Im_0) = my_Pi_Im_function(loc);
        }
    }
}

// template <class data_t>
emtensor_t AxiDilaton::compute_emtensor(const IntVect a_iv,
                                         const RealVect &a_dx,
                                         FArrayBox &a_multigrid_vars_box) const
{
    emtensor_t out;

    DerivativeOperators derivs(a_dx);
    RealVect loc;
    Grids::get_loc(loc, a_iv, a_dx, center);

    Real rr = 0.;
    FOR(i){rr += loc[i]*loc[i];}
    rr = sqrt(rr);

    Real BH_cutoff =  1.0 / (1.0 + exp(- 200.0 * (rr - 0.7))) ;

    Real psi_reg = a_multigrid_vars_box(a_iv, c_psi_reg);
    Real psi_bh = psi_and_Aij_functions->compute_bowenyork_psi(loc);
    Real psi_0 = psi_reg + psi_bh;

    Real Pi_Re_0 = a_multigrid_vars_box(a_iv, c_Pi_Re_0) * BH_cutoff;
    Real phi_Re_0 = a_multigrid_vars_box(a_iv, c_phi_Re_0) * BH_cutoff;
    Real Pi_Im_0 = a_multigrid_vars_box(a_iv, c_Pi_Im_0) * BH_cutoff;
    Real phi_Im_0 = a_multigrid_vars_box(a_iv, c_phi_Im_0) * BH_cutoff;

    Real modulus_phi_squared = phi_Re_0 * phi_Re_0 + phi_Im_0 * phi_Im_0;
    Real modulus_Pi_squared = Pi_Re_0 * Pi_Re_0 + Pi_Im_0 * Pi_Im_0;
    Real KAxiDilaton_phi_squared = 1. / pow(4.0- m_matter_params.gamma_squared_coeff * modulus_phi_squared , 2.0);

    Tensor<1, Real, SpaceDim> d1_phi_Re;
    Tensor<1, Real, SpaceDim> d1_phi_Im;
    derivs.get_d1(d1_phi_Re, a_iv, a_multigrid_vars_box, c_phi_Re_0);
    derivs.get_d1(d1_phi_Im, a_iv, a_multigrid_vars_box, c_phi_Im_0);
    Real d1_phi_squared = 0.;
    FOR(i) { d1_phi_squared += (d1_phi_Re[i] * d1_phi_Re[i] + d1_phi_Im[i] * d1_phi_Im[i]) * BH_cutoff; }

    Real V_of_phi = my_potential_function(phi_Re_0, phi_Im_0);

    out.rho =
        8. * KAxiDilaton_phi_squared * (pow(psi_0, -4.0) * d1_phi_squared + modulus_Pi_squared ) + V_of_phi;
    FOR(i) { out.Si[i] = - 16.* KAxiDilaton_phi_squared * (Pi_Re_0 * d1_phi_Re[i] + Pi_Im_0 * d1_phi_Im[i]); }

    out.rho *= BH_cutoff;
    FOR(i) { out.Si[i]  *= BH_cutoff; }

    return out;
}
