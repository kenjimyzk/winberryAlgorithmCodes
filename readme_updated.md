# Krusell-Smith Model: Dynare 5.x/6.x Compatible Version

This repository contains an updated version of the Krusell-Smith (1998) heterogeneous agent macroeconomic model toolkit, made compatible with **Dynare 5.x and 6.x**.

## Original Work
Based on the original toolkit by **Thomas Winberry** (2016): "A Toolbox for Solving and Estimating Heterogeneous Agent Macro Models"

## Updates for Dynare 5.x/6.x Compatibility

### Key Changes Made:
1. **Function signatures**: Updated `_steadyState` to `_steadystate` (lowercase)
2. **Options structure**: Changed `options` to `options_` throughout
3. **Parameter handling**: Added `deblank()` for parameter name processing
4. **Array indexing**: Fixed `ys(ii)` to `ys(ii,1)` for proper matrix indexing
5. **Optimization options**: Updated `optimoptions` syntax for modern MATLAB versions
6. **Error handling**: Enhanced error reporting and debugging capabilities
7. **Dynare calling**: Added support for both old and new Dynare calling conventions

## Quick Start

### Prerequisites
- MATLAB R2015b or later
- Dynare 5.0 or later (tested with Dynare 5.4 and 6.1)

### Basic Usage

1. **Steady State Analysis**:
```matlab
steadyState
```

2. **Dynamic Analysis**:
```matlab
dynamics
```

### File Structure

#### Main Files:
- `steadyState.m` - Compute and analyze steady state equilibrium
- `dynamics.m` - Solve for aggregate dynamics with TFP shocks

#### Dynare Model Files:
- `firstOrderDynamics_polynomials.mod` - Model using polynomial approximation
- `firstOrderDynamics_splines.mod` - Model using spline approximation
- `parameters_polynomials.mod` / `parameters_splines.mod` - Parameter declarations
- `variables_polynomials.mod` / `variables_splines.mod` - Variable declarations
- `equations_polynomials.mod` / `equations_splines.mod` - Model equations

#### Steady State Functions:
- `firstOrderDynamics_polynomials_steadystate.m` - Steady state for polynomial version
- `firstOrderDynamics_splines_steadystate.m` - Steady state for spline version

#### Core Functions (in Auxiliary Functions/):
- `setParameters.m` - Parameter configuration
- `coreSteadyState.m` - Core steady state computation
- `computeMCResidualPolynomials.m` - Market clearing with parametric distribution
- `computeMCResidualHistogram.m` - Market clearing with histogram method
- `updateCoefficients_*.m` - Policy function updates

## Model Description

### Economic Environment:
- **Households**: Face idiosyncratic employment/unemployment risk
- **Technology**: Cobb-Douglas production with capital and labor
- **Government**: Provides unemployment insurance
- **Equilibrium**: Competitive equilibrium with aggregate uncertainty

### Solution Method:
1. **Individual Problem**: Solved using either:
   - Chebyshev polynomial approximation of conditional expectation
   - Linear spline approximation of savings policy
2. **Distribution**: Approximated using parametric exponential polynomial family
3. **Dynamics**: First-order approximation around steady state

## New Features in Updated Version

### Enhanced Error Handling:
```matlab
try
    % Computation code
    coreSteadyState;
catch ME
    fprintf('Error: %s\n', ME.message);
    % Debugging information provided
end
```

### Version Compatibility Check:
```matlab
% Automatic detection of Dynare version
try
    dyn_ver = dynare_version();
    fprintf('Dynare version: %s.%s\n', dyn_ver.major, dyn_ver.minor);
catch
    warning('Could not determine Dynare version');
end
```

### Improved Plotting:
- Automatic figure saving to `/Figures` directory
- Enhanced plot formatting and styling
- Better axis labels and legends

## Parameters

### Key Economic Parameters:
- `bbeta = 0.96` - Discount factor
- `ssigma = 2` - Risk aversion
- `aalpha = 0.36` - Capital share
- `ddelta = 0.1` - Depreciation rate
- `aggEmployment = 0.93` - Employment rate
- `mmu = 0.15` - Unemployment benefit replacement rate

### Numerical Parameters:
- `nAssets = 25` - Grid points for assets
- `nMeasure = 3` - Moments of distribution
- `splineOpt = 0` - Approximation method (0=polynomials, 1=splines)

## Troubleshooting

### Common Issues:

1. **"Blanks in variable names"**
   - Solution: Updated with `deblank()` function calls

2. **"Steady state not found"**
   - Try changing: `options_.solve_algo = 9` or `options_.solve_algo = 10`

3. **"Undefined function or variable"**
   - Check MATLAB path includes all subdirectories
   - Verify all global variables are properly declared

4. **Convergence Issues**:
   - Adjust `tolerance` and `maxIterations` in `setParameters.m`
   - Try different dampening factors

### Performance Tips:
- Use polynomial approximation (`splineOpt = 0`) for faster computation
- Reduce grid size for initial testing
- Increase tolerance for quicker but less accurate results

## Output

### Steady State Results:
- Aggregate capital stock
- Interest rate and wage
- Policy functions (savings and consumption)
- Wealth distribution (histogram vs parametric approximation)

### Dynamic Results:
- Impulse response functions to TFP shocks
- Variance decomposition
- Theoretical moments

## Validation

The updated code produces results consistent with the original Winberry (2016) toolkit while being compatible with modern Dynare versions.

## Support

For issues specific to this updated version, please check:
1. Dynare version compatibility (requires 5.0+)
2. MATLAB version compatibility (requires R2015b+)
3. All auxiliary functions are in the correct path

## Citation

If you use this code, please cite:
- Winberry, Thomas (2016): "A Toolbox for Solving and Estimating Heterogeneous Agent Macro Models"
- Original Krusell-Smith model: Krusell, Per, and Anthony A. Smith Jr. (1998): "Income and wealth heterogeneity in the macroeconomy." Journal of Political Economy 106.5: 867-896.

## License

This code is provided for academic and research purposes. Please respect the original author's licensing terms.