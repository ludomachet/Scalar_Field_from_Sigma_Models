/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALAXIDILATONDATA_HPP_
#define INITIALAXIDILATONDATA_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4RHS.hpp"
#include "AxiDilaton.hpp"
#include "Tensor.hpp"
#include "Complex.hpp"
#include "SphericalHarmonics.hpp"
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
            center;   //!< Centre of perturbation in initial SF bubble
        double width; //!< Width of bump in initial SF bubble
        double scalar_mass; // Mass of the SF 
        double r0; // location of gaussian pulse
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
        data_t rr = coords.get_radius();
        data_t rr2 = rr * rr;

        // calculate the field value WITH SPHERICAL HARMONICS
     
        SphericalHarmonics::Y_lm_t<data_t> Y_lm011 =
          SphericalHarmonics::spin_Y_lm(coords.x, coords.y, coords.z, 0, 1, 1);
        SphericalHarmonics::Y_lm_t<data_t> Y_lm020 =
          SphericalHarmonics::spin_Y_lm(coords.x, coords.y, coords.z, 0, 2, 0);
        SphericalHarmonics::Y_lm_t<data_t> Y_lm022 =
           SphericalHarmonics::spin_Y_lm(coords.x, coords.y, coords.z, 0, 2, 2);
        SphericalHarmonics::Y_lm_t<data_t> Y_lm02_2 =
           SphericalHarmonics::spin_Y_lm(coords.x, coords.y, coords.z, 0, 2, -2);

	data_t phi = m_params.amplitude / (1. + exp(-(480.-rr)/10.));
	//data_t phi = m_params.amplitude * Y_lm011.Real / (1. + rr / 400.);

        current_cell.store_vars(phi, c_phi_Re);
        current_cell.store_vars(0.0, c_phi_Im);
        current_cell.store_vars(0.0, c_Pi_Re);
        current_cell.store_vars(- m_params.scalar_mass * phi, c_Pi_Im);
    }

  protected:
    double m_dx;
    const params_t m_params; //!< The matter initial condition params
};

#endif /* INITIALAXIDILATONDATA_HPP_ */
