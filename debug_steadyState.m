% Debug script for steadyState.m
% Use this to identify the exact issue at line 157
%
% Updated for Dynare 5.x/6.x compatibility

clear all
close all
clc

% Store current directory and change to Auxiliary Functions
oldFolder = pwd;
cd('./Auxiliary Functions');

try
    %----------------------------------------------------------------
    % Step 1: Set parameters
    %----------------------------------------------------------------
    fprintf('=== Step 1: Setting Parameters ===\n');
    setParameters;
    fprintf('✓ Parameters set successfully\n');
    
    % Display key parameters
    fprintf('  nEpsilon = %d\n', nEpsilon);
    fprintf('  nAssets = %d\n', nAssets);
    fprintf('  nMeasure = %d\n', nMeasure);
    fprintf('  nAssetsFine = %d\n', nAssetsFine);
    fprintf('  splineOpt = %d\n', splineOpt);
    
    % Initialize global price variables
    global r w
    kRepSS = ((aalpha * (aggEmployment ^ (1 - aalpha))) / ((1 / bbeta) - (1 - ddelta))) ^ (1 / (1 - aalpha));
    r = aalpha * (kRepSS ^ (aalpha - 1)) * (aggEmployment ^ (1 - aalpha)) - ddelta;
    w = (kRepSS ^ aalpha) * (1 - aalpha) * (aggEmployment ^ (-aalpha));
    fprintf('  Initial price guess: r = %.4f, w = %.4f\n', r, w);
    
    %----------------------------------------------------------------
    % Step 2: Compute grids
    %----------------------------------------------------------------
    fprintf('\n=== Step 2: Computing Grids ===\n');
    computeGrids;
    fprintf('✓ Grids computed successfully\n');
    
    % Check grid variables
    fprintf('  vAssetsGridFine size: %d\n', length(vAssetsGridFine));
    fprintf('  assetsMin = %.4f, assetsMax = %.4f\n', assetsMin, assetsMax);
    
    %----------------------------------------------------------------
    % Step 3: Compute polynomials (if needed)
    %----------------------------------------------------------------
    if splineOpt == 0
        fprintf('\n=== Step 3: Computing Polynomials ===\n');
        computePolynomials;
        fprintf('✓ Polynomials computed successfully\n');
    end
    
    %----------------------------------------------------------------
    % Step 4: Test coreSteadyState step by step
    %----------------------------------------------------------------
    fprintf('\n=== Step 4: Testing coreSteadyState Components ===\n');
    
    % Test histogram computation first
    fprintf('Testing histogram computation...\n');
    kRepSS = ((aalpha * (aggEmployment ^ (1 - aalpha))) / ((1 / bbeta) - (1 - ddelta))) ^ (1 / (1 - aalpha));
    fprintf('  Representative agent capital: %.4f\n', kRepSS);
    
    % Test market clearing residual function
    try
        residual = computeMCResidualHistogram(1.01 * kRepSS);
        fprintf('✓ Market clearing residual computed: %.6f\n', residual);
    catch ME
        fprintf('✗ Error in computeMCResidualHistogram: %s\n', ME.message);
        cd(oldFolder);
        return;
    end
    
    %----------------------------------------------------------------
    % Step 5: Run full coreSteadyState with error checking
    %----------------------------------------------------------------
    fprintf('\n=== Step 5: Running coreSteadyState ===\n');
    
    % Set display option for debugging
    displayOpt = 'iter-detailed';
    
    try
        coreSteadyState;
        
        if exist('check', 'var') && check == 1
            fprintf('✗ coreSteadyState failed (check = 1)\n');
            cd(oldFolder);
            return;
        end
        
        fprintf('✓ coreSteadyState completed\n');
        
    catch ME
        fprintf('✗ Error in coreSteadyState: %s\n', ME.message);
        fprintf('  Error occurred in: %s at line %d\n', ME.stack(1).name, ME.stack(1).line);
        cd(oldFolder);
        return;
    end
    
    %----------------------------------------------------------------
    % Step 6: Check variables created by coreSteadyState
    %----------------------------------------------------------------
    fprintf('\n=== Step 6: Checking Variables ===\n');
    
    variables_to_check = {'aggregateCapital', 'mMoments', 'mParameters', 'mHat', 'mCoefficients'};
    
    for i = 1:length(variables_to_check)
        var_name = variables_to_check{i};
        if exist(var_name, 'var')
            eval(['var_size = size(' var_name ');']);
            fprintf('✓ %s exists with size [%s]\n', var_name, num2str(var_size));
        else
            fprintf('✗ %s does not exist\n', var_name);
        end
    end
    
    %----------------------------------------------------------------
    % Step 7: Test the specific computation at line 157
    %----------------------------------------------------------------
    fprintf('\n=== Step 7: Testing Distribution Computation (around line 157) ===\n');
    
    if exist('mMoments', 'var') && exist('mParameters', 'var')
        
        % Check dimensions
        fprintf('  mMoments size: [%d, %d]\n', size(mMoments,1), size(mMoments,2));
        fprintf('  mParameters size: [%d, %d]\n', size(mParameters,1), size(mParameters,2));
        fprintf('  vAssetsGridFine size: [%d, %d]\n', size(vAssetsGridFine,1), size(vAssetsGridFine,2));
        
        try
            % Test the computation that likely fails at line 157
            mDistributionFine = zeros(nEpsilon,nAssetsFine);
            
            for iEpsilon = 1 : nEpsilon
                fprintf('  Processing epsilon state %d...\n', iEpsilon);
                
                % First moment (uncentered)
                mGridMoments = zeros(nAssetsFine,nMeasure);
                mGridMoments(:,1) = (vAssetsGridFine - mMoments(iEpsilon,1));
                
                % Higher order moments (centered)
                for iMoment = 2 : nMeasure
                    mGridMoments(:,iMoment) = (vAssetsGridFine - mMoments(iEpsilon,1)) .^ iMoment - ...
                        mMoments(iEpsilon,iMoment);
                end
                
                % Compute density away from borrowing constraint
                mDistributionFine(iEpsilon,:) = mParameters(iEpsilon,1) * exp(mGridMoments * ...
                    mParameters(iEpsilon,2:nMeasure+1)');
                
                fprintf('    ✓ Epsilon state %d completed\n', iEpsilon);
            end
            
            fprintf('✓ Distribution computation successful\n');
            
        catch ME
            fprintf('✗ Error in distribution computation: %s\n', ME.message);
            fprintf('  This is likely the error at line 157\n');
            
            % More detailed error information
            if exist('iEpsilon', 'var')
                fprintf('  Error occurred at iEpsilon = %d\n', iEpsilon);
            end
            
            cd(oldFolder);
            return;
        end
        
    else
        fprintf('✗ Required variables (mMoments, mParameters) not available\n');
    end
    
    %----------------------------------------------------------------
    % Success message
    %----------------------------------------------------------------
    fprintf('\n=== DEBUG COMPLETED SUCCESSFULLY ===\n');
    fprintf('All components are working. You should be able to run steadyState.m now.\n');
    
catch ME
    fprintf('\n=== DEBUG FAILED ===\n');
    fprintf('Error in %s at line %d\n', ME.stack(1).name, ME.stack(1).line);
    fprintf('Error message: %s\n', ME.message);
end

% Return to original directory
cd(oldFolder);