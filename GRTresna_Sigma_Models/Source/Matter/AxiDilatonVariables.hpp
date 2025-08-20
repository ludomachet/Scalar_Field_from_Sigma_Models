/* GRTresna
 * Copyright 2024 The GRTL Collaboration.
 * Please refer to LICENSE in GRTresna's root directory.
 */

#ifndef AXIDILATONVARIABLES_HPP_
#define AXIDILATONVARIABLES_HPP_

#include "MetricVariables.hpp"
#include "ParityDefinitions.hpp"

// Matter Vars
enum
{
    c_phi_Re_0 = NUM_METRIC_VARS,
    c_phi_Im_0,
    c_Pi_Re_0,
    c_Pi_Im_0,
    NUM_MULTIGRID_VARS
};

namespace MatterVariables
{

static const std::array<std::string, NUM_MULTIGRID_VARS - NUM_METRIC_VARS>
    variable_names = {"phi_Re_0", "phi_Im_0", "Pi_Re_0", "Pi_Im_0"};

static constexpr std::array<int, NUM_MULTIGRID_VARS - NUM_METRIC_VARS> const
    vars_parity = {EVEN, EVEN, EVEN, EVEN};

} // namespace MatterVariables

#endif // AXIDILATONDVARIABLES_HPP_