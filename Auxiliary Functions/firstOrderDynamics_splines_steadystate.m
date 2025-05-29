function [ys,params,check] = firstOrderDynamics_splines_steadystate(ys,exo,M_,options_)
% Computes stationary equilibrium of the model for Dynare; format is required
% to be called by Dynare (follows example of NK_baseline.mod in Dynare examples)
% Updated for Dynare 5.x/6.x compatibility
%
% Thomas Winberry, July 26th, 2016
% Updated for Dynare 5.x/6.x compatibility

tStart = tic;
fprintf('\nComputing steady state...\n')

%----------------------------------------------------------------
% Call parameters (the next set of commands will overwrite some)
%----------------------------------------------------------------
setParameters;

%----------------------------------------------------------------
% Read in parameters from Dynare declaration
%----------------------------------------------------------------

% Initialize indicator
check = 0;

% Read out parameters to access them with their name
for iParameter = 1:M_.param_nbr
    paramname = deblank(M_.param_names{iParameter,:});  % Added deblank for Dynare 5.x/6.x
    eval(['global ' paramname]);
    eval([paramname ' = M_.params(' int2str(iParameter) ');']);
end

%----------------------------------------------------------------
% Steady State Computation
%----------------------------------------------------------------
try
    displayOpt = 'off';       % 'iter-detailed' or 'off'
    coreSteadyState;
    
    % Check if steady state computation was successful
    if check == 1
        fprintf('Steady state computation failed in coreSteadyState\n');
        return;
    end
    
    % Prices
    r = aalpha * (aggregateCapital ^ (aalpha - 1)) * (aggEmployment ^ (1 - aalpha)) - ddelta;
    w = (aggregateCapital ^ aalpha) * (1 - aalpha) * (aggEmployment ^ (-aalpha));
    
catch ME
    fprintf('Error in steady state computation: %s\n', ME.message);
    check = 1;
    ys = NaN(M_.orig_endo_nbr,1);
    params = NaN(M_.param_nbr,1);
    return;
end

%----------------------------------------------------------------
% Save values of steady state variables for Dynare (must be exactly
% as declared in Dynare)
%----------------------------------------------------------------

try
    % Coefficients on spline approximation (different from polynomial version)
    for iEpsilon = 1 : nEpsilon
        for iAsset = 1 : nAssets
            eval(sprintf('coefficient_%d_%d = mCoefficients(iEpsilon,iAsset);',...
                iEpsilon,iAsset));
        end
    end
    
    % Moments and parameters of density away from borrowing constraint
    for iEpsilon = 1 : nEpsilon
        for iMoment = 1 : nMeasure
            eval(sprintf('moment_%d_%d = mMoments(iEpsilon,iMoment);',iEpsilon,iMoment));
            eval(sprintf('measureCoefficient_%d_%d = mParameters(iEpsilon,iMoment+1);',iEpsilon,iMoment));
        end
    end
    
    % Mass at borrowing constraint
    for iEpsilon = 1 : nEpsilon
        eval(sprintf('mHat_%d = mHat(iEpsilon);',iEpsilon));
    end
    
    % Other variables
    aggregateCapital = (1 - mHat(1,1)) * (1 - aggEmployment) * mMoments(1,1) + (1 - mHat(2,1)) * aggEmployment * mMoments(2,1);
    aggregateTFP = 0;
    logAggregateOutput = log(exp(aggregateTFP) * (aggregateCapital ^ aalpha) * (aggEmployment ^ (1 - aalpha)));
    logAggregateInvestment = log(ddelta * aggregateCapital);
    logAggregateConsumption = log(exp(logAggregateOutput) - exp(logAggregateInvestment));
    logWage = log(w);
    
catch ME
    fprintf('Error in setting steady state variables: %s\n', ME.message);
    check = 1;
    ys = NaN(M_.orig_endo_nbr,1);
    params = NaN(M_.param_nbr,1);
    return;
end

%----------------------------------------------------------------
% Prepare output for Dynare
%----------------------------------------------------------------

% Initialize parameter vector
params = NaN(M_.param_nbr,1);
for iter = 1:length(M_.params) % update parameters set in the file
    eval(['params(' num2str(iter) ',1) = ' deblank(M_.param_names{iter}) ';' ])
end

% Initialize endogenous variables vector
ys = NaN(M_.orig_endo_nbr,1);

% Save endogenous variables back into ys
for ii = 1 : M_.orig_endo_nbr
    varname = deblank(M_.endo_names{ii,:});  % Added deblank for Dynare 5.x/6.x
    try
        eval(['ys(' int2str(ii) ',1) = ' varname ';']);  % Fixed indexing for Dynare 5.x/6.x
    catch
        fprintf('Warning: Could not set variable %s\n', varname);
        ys(ii,1) = 0;  % Set to zero if variable not found
    end
end

% Final check for NaN values
if any(isnan(ys))
    fprintf('Warning: Some steady state values are NaN\n');
    nan_indices = find(isnan(ys));
    for i = 1:length(nan_indices)
        idx = nan_indices(i);
        varname = deblank(M_.endo_names{idx,:});
        fprintf('  Variable %s (index %d) is NaN\n', varname, idx);
    end
    check = 1;
end

if any(isnan(params))
    fprintf('Warning: Some parameter values are NaN\n');
    check = 1;
end

fprintf('... Done!  Elapsed time: %2.2f seconds \n\n',toc(tStart))