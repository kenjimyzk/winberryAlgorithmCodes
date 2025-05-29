// Dynare shell which declares model and solves for aggregate dynamics using
// first order approximation (when approximating conditional expectation with 
// polynomials) - Updated for Dynare 5.x/6.x compatibility
//
// Thomas Winberry, July 26th, 2016
// Updated for Dynare 5.x/6.x compatibility

//----------------------------------------------------------------
// Load parameters
//----------------------------------------------------------------

@#include "parameters_polynomials.mod"

//----------------------------------------------------------------
// Define variables
//----------------------------------------------------------------

@#include "variables_polynomials.mod"

//----------------------------------------------------------------
// Model equations
//----------------------------------------------------------------

model;

@#include "equations_polynomials.mod"

end;

//----------------------------------------------------------------
// Computation
//----------------------------------------------------------------

// Specify shock process
shocks;
    var aggregateTFPShock = 1;
end;

// Set Dynare options for compatibility with version 5.x/6.x
options_.solve_algo = 4;                    // Algorithm for solving the steady state
options_.steadystate.nocheck = 1;          // Don't check steady state residuals
options_.qz_criterium = 1.000001;          // Criterion for stable eigenvalues
options_.lyapunov_fp = 1;                  // Use fixed point iteration for Lyapunov equation
options_.sylvester_fp = 1;                 // Use fixed point iteration for Sylvester equation

// Alternative solver options that can be tried if the default fails:
// options_.solve_algo = 9;                 // Trust region algorithm
// options_.solve_algo = 10;                // Levenberg-Marquardt mixed complementarity problem
// options_.maxit_ = 1000;                  // Maximum number of iterations
// options_.tolf = 1e-5;                    // Tolerance on function values

// Compute steady state
steady(nocheck);

// Check regularity conditions (uncomment to check)
// check;
// model_diagnostics;
// model_info;

// Simulate the model
stoch_simul(order=1,hp_filter=100,irf=40) aggregateTFP logAggregateOutput 
	logAggregateConsumption logAggregateInvestment logWage r;