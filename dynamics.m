% Shell which declares parameters and calls Dynare to solve model using
% first order approximation of aggregate dynamics
% Updated for Dynare 5.x/6.x compatibility
%
% Thomas Winberry, July 26th, 2016
% Updated for Dynare 5.x/6.x compatibility

clear all
close all
clc

% Store current directory and change to Auxiliary Functions
oldFolder = pwd;
cd('./Auxiliary Functions');

try
    %----------------------------------------------------------------
    % Set parameters
    %----------------------------------------------------------------
    setParameters;
    
    %----------------------------------------------------------------
    % Compute approximation tools
    %----------------------------------------------------------------
    
    % Grids
    computeGrids;
    
    % Polynomials over grids (only if using polynomials to approximate conditional expectation)
    if splineOpt == 0
        computePolynomials;
    end
    
    %----------------------------------------------------------------
    % Save parameters in .mat files to import into Dynare 
    %----------------------------------------------------------------
    
    if splineOpt == 0	% if using polynomials to approximate individual decisions
        
        % Economic parameters
        save economicParameters.mat bbeta ssigma aaBar aalpha ddelta vEpsilonGrid aggEmployment ...
            mmu ttau rrhoTFP ssigmaTFP mEpsilonTransition
            
        % Approximation parameters
        save approximationParameters.mat nEpsilon nAssets nState assetsMin assetsMax nAssetsFine nStateFine nAssetsQuadrature nStateQuadrature ...
            nMeasure nMeasureCoefficients kRepSS maxIterations tolerance dampening
            
        % Grids
        save grids.mat vAssetsGridZeros vAssetsGrid mEpsilonGrid mAssetsGrid mEpsilonPrimeGrid vAssetsGridFine ...
            vAssetsGridFineZeros mEpsilonGridFine mAssetsGridFine mEpsilonPrimeGridFine vQuadratureWeights ...
            vAssetsGridQuadratureZeros vAssetsGridQuadrature mEpsilonGridQuadrature mAssetsGridQuadrature
            
        % Polynomials
        save polynomials.mat vAssetsPoly vAssetsPolySquared vAssetsPolyFine vAssetsPolyQuadrature vAssetsPolyBC
        
    else	% if using splines to approximate individual decisions
        
        % Economic parameters
        save economicParameters.mat bbeta ssigma aaBar aalpha ddelta vEpsilonGrid aggEmployment ...
            mmu ttau rrhoTFP ssigmaTFP mEpsilonTransition
            
        % Approximation parameters
        save approximationParameters.mat nEpsilon nAssets nState assetsMin assetsMax nAssetsFine nStateFine nAssetsQuadrature nStateQuadrature ...
            nMeasure nMeasureCoefficients kRepSS maxIterations tolerance dampening
            
        % Grids
        save grids.mat vAssetsGrid mEpsilonGrid mAssetsGrid mEpsilonPrimeGrid vAssetsGridFine ...
            mEpsilonGridFine mAssetsGridFine mEpsilonPrimeGridFine vQuadratureWeights ...
            vAssetsGridQuadrature mEpsilonGridQuadrature mAssetsGridQuadrature
        
    end
    
    %----------------------------------------------------------------
    % Run Dynare
    %----------------------------------------------------------------
    
    fprintf('\n=== Running Dynare ===\n');
    
    if splineOpt == 0	% if using polynomials to approximate individual decisions
        
        fprintf('Using polynomial approximation...\n');
        
        % Different calling conventions for different environments
        if exist('OCTAVE_VERSION', 'builtin')
            % Octave
            dynare('firstOrderDynamics_polynomials.mod');
        else
            % MATLAB - try modern syntax first, fall back to old syntax
            try
                dynare('firstOrderDynamics_polynomials');
            catch ME1
                try
                    dynare firstOrderDynamics_polynomials;
                catch ME2
                    fprintf('Error with both Dynare calling methods:\n');
                    fprintf('Method 1 error: %s\n', ME1.message);
                    fprintf('Method 2 error: %s\n', ME2.message);
                    rethrow(ME2);
                end
            end
        end
        
    else	% if using splines to approximate individual decisions
        
        fprintf('Using spline approximation...\n');
        
        % Different calling conventions for different environments
        if exist('OCTAVE_VERSION', 'builtin')
            % Octave
            dynare('firstOrderDynamics_splines.mod');
        else
            % MATLAB - try modern syntax first, fall back to old syntax
            try
                dynare('firstOrderDynamics_splines');
            catch ME1
                try
                    dynare firstOrderDynamics_splines;
                catch ME2
                    fprintf('Error with both Dynare calling methods:\n');
                    fprintf('Method 1 error: %s\n', ME1.message);
                    fprintf('Method 2 error: %s\n', ME2.message);
                    rethrow(ME2);
                end
            end
        end
        
    end
    
    fprintf('=== Dynare completed successfully ===\n');
    
catch ME
    fprintf('\n=== Error occurred ===\n');
    fprintf('Error in %s at line %d\n', ME.stack(1).name, ME.stack(1).line);
    fprintf('Error message: %s\n', ME.message);
    
    % Provide helpful debugging information
    fprintf('\n=== Debugging Information ===\n');
    fprintf('Current splineOpt: %d\n', splineOpt);
    fprintf('Current directory: %s\n', pwd);
    
    % Check if required files exist
    if splineOpt == 0
        fprintf('Looking for: firstOrderDynamics_polynomials.mod\n');
        if exist('firstOrderDynamics_polynomials.mod', 'file')
            fprintf('✓ File found\n');
        else
            fprintf('✗ File not found\n');
        end
    else
        fprintf('Looking for: firstOrderDynamics_splines.mod\n');
        if exist('firstOrderDynamics_splines.mod', 'file')
            fprintf('✓ File found\n');
        else
            fprintf('✗ File not found\n');
        end
    end
    
    % Return to original directory before rethrowing
    cd(oldFolder);
    rethrow(ME);
end

% Return to original directory
cd(oldFolder);

fprintf('\nDynare simulation completed. Check the workspace for results.\n');