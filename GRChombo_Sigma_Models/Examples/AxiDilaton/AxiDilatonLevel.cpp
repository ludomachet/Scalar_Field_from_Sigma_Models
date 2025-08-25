/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "AxiDilatonLevel.hpp"
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "CustomExtraction.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "SixthOrderDerivatives.hpp"
#include "TraceARemoval.hpp"

// For RHS update
#include "MatterCCZ4RHS.hpp"

// For constraints calculation
#include "NewMatterConstraints.hpp"

// For tag cells
#include "FixedGridsTaggingCriterion.hpp"
// #include "ChiAndPhiTaggingCriterion.hpp"
// #include "FixedGridsTaggingCriterion.hpp"

// Problem specific includes
#include "AxiDilaton.hpp"
#include "AxiDilatonPotential.hpp"
#include "ComputePack.hpp"
#include "GammaCalculator.hpp"
#include "InitialAxiDilatonData.hpp"
#include "KerrBH.hpp"
#include "SetValue.hpp"

// ADM quantities
#include "ADMQuantities.hpp"
#include "ADMQuantitiesExtraction.hpp"

// Things to do at each advance step, after the RK4 is calculated
void AxiDilatonLevel::specificAdvance()
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(),
                          PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(
            NanCheck(m_dx, m_p.center, "NaNCheck in specific Advance"),
            m_state_new, m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

// Initial data for field and metric variables
void AxiDilatonLevel::initialData()
{
    CH_TIME("AxiDilatonLevel::initialData");
    if (m_verbosity)
        pout() << "AxiDilatonLevel::initialData " << m_level << endl;

    // First set everything to zero then initial conditions for scalar field -
    // here a Kerr BH and a scalar field profile
    BoxLoops::loop(
        make_compute_pack(SetValue(0.), KerrBH(m_p.kerr_params, m_dx),
                          InitialAxiDilatonData(m_p.initial_params, m_dx)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    fillAllGhosts();
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);
}

// Things to do when restarting from a checkpoint, including
// restart from the initial condition solver output
void AxiDilatonLevel::postRestart()
{
    // On restart calculate the constraints on every level
    fillAllGhosts();
    AxiDilatonPotential potential(m_p.potential_params);
    AxiDilatonWithPotential axi_dilaton_field(potential,
                                              m_p.gamma_squared_coeff);
    BoxLoops::loop(MatterConstraints<AxiDilatonWithPotential>(
                       axi_dilaton_field, m_dx, m_p.G_Newton, c_Ham,
                       Interval(c_Mom, c_Mom)),
                   m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

    // Use AMR Interpolator and do lineout data extraction
    // pass the boundary params so that we can use symmetries
    AMRInterpolator<Lagrange<2>> interpolator(
        m_gr_amr, m_p.origin, m_p.dx, m_p.boundary_params, m_p.verbosity);

    // this should fill all ghosts including the boundary ones according
    // to the conditions set in params.txt
    interpolator.refresh();

    // restart works from level 0 to highest level, so want this to happen last
    // on finest level
    int write_out_level = m_p.max_level;
    if (m_level == write_out_level)
    {
        // AMRReductions for diagnostic variables
        AMRReductions<VariableType::diagnostic> amr_reductions_diagnostic(
            m_gr_amr);
        double L2_Ham = amr_reductions_diagnostic.norm(c_Ham);
        double L2_Mom = amr_reductions_diagnostic.norm(c_Mom);

        // only on rank zero write out the result
        if (procID() == 0)
        {
            pout() << "The initial norm of the constraint vars on restart is "
                   << L2_Ham << " for the Hamiltonian constraint and " << L2_Mom
                   << " for the momentum constraints" << endl;
        }

        // set up the query and execute it
        int num_points = 3 * m_p.ivN[0];
        CustomExtraction constraint_extraction(c_Ham, c_Mom, num_points, m_p.L,
                                               m_p.center, m_dt, m_time);
        constraint_extraction.execute_query(
            &interpolator, m_p.data_path + "constraint_lineout");
    }
}

#ifdef CH_USE_HDF5
// Things to do before outputting a checkpoint file
void AxiDilatonLevel::prePlotLevel()
{
    fillAllGhosts();
    AxiDilatonPotential potential(m_p.potential_params);
    AxiDilatonWithPotential axi_dilaton_field(potential,
                                              m_p.gamma_squared_coeff);
    BoxLoops::loop(MatterConstraints<AxiDilatonWithPotential>(
                       axi_dilaton_field, m_dx, m_p.G_Newton, c_Ham,
                       Interval(c_Mom, c_Mom)),
                   m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}
#endif

// Things to do in RHS update, at each RK4 step
void AxiDilatonLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                      const double a_time)
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(),
                          PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
    AxiDilatonPotential potential(m_p.potential_params);
    AxiDilatonWithPotential axi_dilaton_field(potential,
                                              m_p.gamma_squared_coeff);
    if (m_p.max_spatial_derivative_order == 4)
    {
        MatterCCZ4RHS<AxiDilatonWithPotential, MovingPunctureGauge,
                      FourthOrderDerivatives>
            my_ccz4_matter(axi_dilaton_field, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.rescale_sigma, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else if (m_p.max_spatial_derivative_order == 6)
    {
        MatterCCZ4RHS<AxiDilatonWithPotential, MovingPunctureGauge,
                      SixthOrderDerivatives>
            my_ccz4_matter(axi_dilaton_field, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.rescale_sigma, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
}

// Things to do at ODE update, after soln + rhs
void AxiDilatonLevel::specificUpdateODE(GRLevelData &a_soln,
                                        const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void AxiDilatonLevel::preTagCells()
{
    fillAllGhosts(VariableType::evolution, Interval(c_chi, c_chi));
    fillAllGhosts(VariableType::evolution, Interval(c_phi_Re, c_phi_Re));
    fillAllGhosts(VariableType::evolution, Interval(c_phi_Im, c_phi_Im));
}

void AxiDilatonLevel::computeTaggingCriterion(
    FArrayBox &tagging_criterion, const FArrayBox &current_state,
    const FArrayBox &current_state_diagnostics)
{
    BoxLoops::loop(FixedGridsTaggingCriterion(m_dx, m_level, m_p.L, m_p.center),
                   current_state, tagging_criterion);
}

// void AxiDilatonLevel::computeTaggingCriterion(
//     FArrayBox &tagging_criterion, const FArrayBox &current_state,
//     const FArrayBox &current_state_diagnostics)
// {
//     BoxLoops::loop(ChiAndPhiTaggingCriterion(m_dx, m_p.threshold_chi,
//     m_p.threshold_phi),
//                    current_state, tagging_criterion);
// }

void AxiDilatonLevel::specificPostTimeStep()
{
#ifdef USE_AHFINDER
    if (m_p.AH_activate && m_level == m_p.AH_params.level_to_run)
        m_bh_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);
#endif

    CH_TIME("AxiDilatonLevel::specificPostTimeStep");
    // Do the extraction on the min extraction level
    if (m_p.activate_extraction == 1)
    {
        int min_level = m_p.extraction_params.min_extraction_level();
        bool calculate_adm = at_level_timestep_multiple(min_level);
        if (calculate_adm)
        {
            // Populate the ADM Mass and Spin values on the grid
            fillAllGhosts();
            BoxLoops::loop(ADMQuantities(m_p.extraction_params.center, m_dx,
                                         c_Madm, c_Jadm),
                           m_state_new, m_state_diagnostics,
                           EXCLUDE_GHOST_CELLS);

            if (m_level == min_level)
            {
                CH_TIME("ADMExtraction");
                // Now refresh the interpolator and do the interpolation
                m_gr_amr.m_interpolator->refresh();
                ADMQuantitiesExtraction my_extraction(
                    m_p.extraction_params, m_dt, m_time, m_restart_time, c_Madm,
                    c_Jadm);
                my_extraction.execute_query(m_gr_amr.m_interpolator);
            }
        }
    }
}
