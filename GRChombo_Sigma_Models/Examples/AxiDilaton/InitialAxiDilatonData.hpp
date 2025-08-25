/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALAXIDILATONDATA_HPP_
#define INITIALAXIDILATONDATA_HPP_

#include "AxiDilaton.hpp"
#include "Cell.hpp"
#include "Complex.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4RHS.hpp"
#include "SphericalHarmonics.hpp"
#include "SuperradiantInitialCondition.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which sets the initial scalar field matter config
class InitialAxiDilatonData
{
  public:
    //! A structure for the input params for scalar field properties and initial
    //! conditions
    struct params_t
    {
        double amplitude; //!< Amplitude of bump in initial SF bubble
        std::array<double, CH_SPACEDIM>
            center;         //!< Centre of perturbation in initial SF bubble
        double width;       //!< Width of bump in initial SF bubble
        double scalar_mass; // Mass of the SF
        double r0;          // location of gaussian pulse
        double kerr_mass;
        double kerr_spin;
    };

    //! The constructor
    InitialAxiDilatonData(params_t a_params, double a_dx)
        : m_dx(a_dx), m_params(a_params)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        // Coordinates<double> coords(current_cell, m_dx, m_params.center);
        data_t rr = coords.get_radius();

        // double my_x = current_cell[0] * m_dx - m_params.center;

        // cout << m_params.kerr_spin << endl;

        double kerrA = m_params.kerr_spin;
        double kerrM = m_params.kerr_mass;
        double scalarM = m_params.scalar_mass;

        // SuperradiantInitCond::Superradiant Initial =
        // SuperradiantInitCond::InitialConditionFunction(coords.x, coords.y,
        // coords.z, m_params.scalar_mass, m_params.kerr_spin,
        // m_params.kerr_mass);
        //  SuperradiantInitCond::InitialConditionFunction(coords.x, coords.y,
        //  coords.z, scalarM, kerrA, kerrM);

        // cout << "Initial condition = " << Initial.Real << "+I " << Initial.Im
        // << endl;

        // data_t phi_Re = Initial.Real;
        // data_t phi_Im = Initial.Im;

        // current_cell.store_vars(phi_Re, c_phi_Re);
        // current_cell.store_vars(phi_Im, c_phi_Im);
        // current_cell.store_vars(- m_params.scalar_mass * phi_Im, c_Pi_Re);
        // current_cell.store_vars(- m_params.scalar_mass * phi_Re, c_Pi_Im);

        // calculate the field value WITH SPHERICAL HARMONICS

        SphericalHarmonics::Y_lm_t<data_t> Y_lm011 =
            SphericalHarmonics::spin_Y_lm(coords.x, coords.y, coords.z, 0, 1,
                                          1);
        SphericalHarmonics::Y_lm_t<data_t> Y_lm020 =
            SphericalHarmonics::spin_Y_lm(coords.x, coords.y, coords.z, 0, 2,
                                          0);
        SphericalHarmonics::Y_lm_t<data_t> Y_lm022 =
            SphericalHarmonics::spin_Y_lm(coords.x, coords.y, coords.z, 0, 2,
                                          2);
        SphericalHarmonics::Y_lm_t<data_t> Y_lm02_2 =
            SphericalHarmonics::spin_Y_lm(coords.x, coords.y, coords.z, 0, 2,
                                          -2);

        // data_t phi = m_params.amplitude / (1+1e5/(1e5-rr));
        // data_t phi = m_params.amplitude / (1. + exp(-(480.-rr)/10.));
        // data_t phi = m_params.amplitude / (1. + exp(-(480.-rr)/10.));
        data_t phi = m_params.amplitude / (1. + rr / 1024.);

        current_cell.store_vars(phi, c_phi_Re);
        current_cell.store_vars(0.0, c_phi_Im);
        current_cell.store_vars(0.0, c_Pi_Re);
        current_cell.store_vars(m_params.scalar_mass * phi, c_Pi_Im);
    }

  protected:
    double m_dx;
    const params_t m_params; //!< The matter initial condition params
};

#endif /* INITIALAXIDILATONDATA_HPP_ */

// class InitialAxiDilatonData
// {
//   public:
//     //! A structure for the input params for scalar field properties and
//     initial
//     //! conditions
//     struct params_t
//     {
//         double amplitude; //!< Amplitude of bump in initial SF bubble
//         std::array<double, CH_SPACEDIM>
//             center;   //!< Centre of perturbation in initial SF bubble
//         double width; //!< Width of bump in initial SF bubble
//         double scalar_mass; // Mass of the SF
//         double r0; // location of gaussian pulse
//         double kerr_mass;
//         double kerr_spin;
//     };

//     //! The constructor
//     InitialAxiDilatonData(params_t a_params, double a_dx)
//         : m_dx(a_dx), m_params(a_params)
//     {
//     }

//     //! Function to compute the value of all the initial vars on the grid
//     template <class data_t> void compute(Cell<data_t> current_cell) const
//     {
//         // where am i?
//         Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

//         cout << m_params.kerr_spin << endl;

//         data_t kerrA = m_params.kerr_spin;
//         data_t kerrM = m_params.kerr_mass;
//         data_t scalarM = m_params.scalar_mass;

//         SuperradiantInitCond::Superradiant<data_t> Initial =
//             //SuperradiantInitCond::InitialConditionFunction(coords.x,
//             coords.y, coords.z, m_params.scalar_mass, m_params.kerr_spin,
//             m_params.kerr_mass);
//             SuperradiantInitCond::InitialConditionFunction(coords.x,
//             coords.y, coords.z, scalarM, kerrA, kerrM);

//         data_t phi_Re = Initial.Real;
//         data_t phi_Im = Initial.Im;

//         current_cell.store_vars(phi_Re, c_phi_Re);
//         current_cell.store_vars(phi_Im, c_phi_Im);
//         current_cell.store_vars(- m_params.scalar_mass * phi_Im, c_Pi_Re);
//         current_cell.store_vars(- m_params.scalar_mass * phi_Re, c_Pi_Im);

//         // calculate the field value WITH SPHERICAL HARMONICS

//         // SphericalHarmonics::Y_lm_t<data_t> Y_lm011 =
//         //   SphericalHarmonics::spin_Y_lm(coords.x, coords.y, coords.z, 0,
//         1, 1);
//         // SphericalHarmonics::Y_lm_t<data_t> Y_lm020 =
//         //   SphericalHarmonics::spin_Y_lm(coords.x, coords.y, coords.z, 0,
//         2, 0);
//         // SphericalHarmonics::Y_lm_t<data_t> Y_lm022 =
//         //    SphericalHarmonics::spin_Y_lm(coords.x, coords.y, coords.z, 0,
//         2, 2);
//         // SphericalHarmonics::Y_lm_t<data_t> Y_lm02_2 =
//         //    SphericalHarmonics::spin_Y_lm(coords.x, coords.y, coords.z, 0,
//         2, -2);

//         //data_t phi = m_params.amplitude / (1. + rr / 350.);
//         // data_t phi = m_params.amplitude * Y_lm011.Real / (1. + rr / 400.);

//         // current_cell.store_vars(phi, c_phi_Re);
//         // current_cell.store_vars(0.0, c_phi_Im);
//         // current_cell.store_vars(0.0, c_Pi_Re);
//         // current_cell.store_vars(- m_params.scalar_mass * phi, c_Pi_Im);
//     }

//   protected:
//     double m_dx;
//     const params_t m_params; //!< The matter initial condition params
// };

// #endif /* INITIALAXIDILATONDATA_HPP_ */
