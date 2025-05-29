# Installation and Setup Guide

## System Requirements

### Minimum Requirements:
- **MATLAB**: R2015b or later
- **Dynare**: 5.0 or later 
- **Operating System**: Windows, macOS, or Linux
- **RAM**: 4GB minimum, 8GB recommended
- **Storage**: 500MB free space

### Recommended Setup:
- **MATLAB**: R2020a or later
- **Dynare**: 6.1 or later
- **RAM**: 16GB for large-scale computations

## Installation Steps

### 1. Install MATLAB
Download and install MATLAB from MathWorks. Ensure you have the following toolboxes:
- Optimization Toolbox
- Statistics and Machine Learning Toolbox (optional but recommended)

### 2. Install Dynare

#### Option A: Download from Official Site
1. Visit: https://www.dynare.org/download/
2. Download Dynare 6.1 or later
3. Follow installation instructions for your OS

#### Option B: MATLAB Add-On (Recommended)
```matlab
% In MATLAB Command Window
% Go to Home > Add-Ons > Get Add-Ons
% Search for "Dynare" and install
```

### 3. Verify Installation
```matlab
% Test Dynare installation
dynare_version()

% Should display something like:
% This is Dynare version 6.1
```

### 4. Download the Updated Krusell-Smith Code
1. Download all the updated files to a folder (e.g., `KS_Model/`)
2. Create the following directory structure:
```
KS_Model/
├── steadyState.m
├── dynamics.m
├── firstOrderDynamics_polynomials.mod
├── firstOrderDynamics_splines.mod
├── Auxiliary Functions/
│   ├── setParameters.m
│   ├── coreSteadyState.m
│   ├── computeGrids.m
│   ├── computePolynomials.m
│   ├── firstOrderDynamics_polynomials_steadystate.m
│   ├── firstOrderDynamics_splines_steadystate.m
│   ├── parameters_polynomials.mod
│   ├── parameters_splines.mod
│   ├── variables_polynomials.mod
│   ├── variables_splines.mod
│   ├── equations_polynomials.mod
│   ├── equations_splines.mod
│   └── [other auxiliary functions...]
└── Figures/ (will be created automatically)
```

## Setup and Configuration

### 1. Set MATLAB Path
```matlab
% Navigate to your KS_Model directory
cd('path/to/KS_Model')

% Add current directory and subdirectories to path
addpath(genpath(pwd))

% Save the path for future sessions
savepath
```

### 2. Test Installation
```matlab
% Quick test - this should run without errors
cd('path/to/KS_Model')
setParameters
fprintf('Installation successful!\n')
```

### 3. Configure Dynare (if needed)
```matlab
% If Dynare is not automatically detected, add it manually:
addpath('path/to/dynare/matlab')
```

## Running Your First Simulation

### Test 1: Steady State Computation
```matlab
% Navigate to model directory
cd('path/to/KS_Model')

% Run steady state computation (takes 2-5 minutes)
steadyState

% Expected output:
% - Text output showing convergence
% - 4 figures showing policy functions and distributions
% - Figures saved to ./Figures/ directory
```

### Test 2: Dynamic Analysis
```matlab
% Run dynamic analysis (takes 5-10 minutes)
dynamics

% Expected output:
% - Dynare output showing model solution
% - Impulse response functions
% - Theoretical moments
```

## Troubleshooting Common Installation Issues

### Issue 1: "Dynare not found"
**Solution:**
```matlab
% Check if Dynare is in path
which dynare

% If empty, add Dynare to path:
addpath('/Applications/Dynare/6.1/matlab')  % macOS example
addpath('C:\dynare\6.1\matlab')             % Windows example
savepath
```

### Issue 2: "optimoptions not found"
**Cause:** Using very old MATLAB version
**Solution:** Update to MATLAB R2013a or later, or modify code to use `optimset`

### Issue 3: "Permission denied" errors
**Solution:**
```matlab
% Make sure you have write permissions in the directory
% Try running MATLAB as administrator (Windows) or with sudo (Linux/macOS)
```

### Issue 4: "Out of memory" errors
**Solutions:**
1. **Reduce grid size:**
```matlab
% Edit setParameters.m
nAssets = 15;        % Reduce from 25
nAssetsFine = 50;    % Reduce from 100
```

2. **Increase MATLAB memory:**
```matlab
% In MATLAB startup file
java.lang.System.setProperty('java.awt.headless','true')
```

### Issue 5: Convergence problems
**Solutions:**
1. **Adjust tolerance:**
```matlab
% In setParameters.m
tolerance = 1e-4;    % Increase from 1e-5
```

2. **Change solver:**
```matlab
% In .mod files, try different solver:
options_.solve_algo = 9;  % or 10
```

## Performance Optimization

### For Faster Computation:
```matlab
% In setParameters.m
splineOpt = 0;           % Use polynomials (faster)
nAssets = 20;            % Reduce grid size
maxIterations = 1e4;     % Reduce max iterations
tolerance = 1e-4;        % Increase tolerance
```

### For Higher Accuracy:
```matlab
% In setParameters.m
nAssets = 30;            % Increase grid size
nAssetsFine = 200;       % Finer grid for plots
tolerance = 1e-6;        % Tighter tolerance
maxIterations = 5e4;     % More iterations
```

## Validation Tests

Run these tests to ensure everything is working correctly:

### Test 1: Parameter Loading
```matlab
setParameters
fprintf('Beta: %.3f, Alpha: %.3f\n', bbeta, aalpha)
% Expected: Beta: 0.960, Alpha: 0.360
```

### Test 2: Grid Computation
```matlab
setParameters
computeGrids
fprintf('Assets grid size: %d\n', length(vAssetsGrid))
% Expected: Assets grid size: 25
```

### Test 3: Steady State (Quick)
```matlab
% Set fast parameters for testing
splineOpt = 0;
nAssets = 10;
tolerance = 1e-3;
steadyState
% Should complete in under 1 minute
```

## Getting Help

### Built-in Diagnostics:
The updated code includes diagnostic information:
```matlab
% Run dynamics with diagnostics
dynamics  % Will show version info and file checks
```

### Common Commands for Debugging:
```matlab
% Check Dynare version
dynare_version()

% Check current directory
pwd

% List files in current directory
ls

% Check if specific functions exist
which setParameters
which coreSteadyState
```

### Performance Monitoring:
```matlab
% Time your computations
tic
steadyState
elapsed_time = toc
fprintf('Computation took %.2f seconds\n', elapsed_time)
```

Now you're ready to explore the Krusell-Smith heterogeneous agent model with the updated, Dynare-compatible code!