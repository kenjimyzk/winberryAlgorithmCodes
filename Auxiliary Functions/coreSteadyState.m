% Computes market clearing capital stock and associated distribution and 
% decision rules in steady state
% Updated for Dynare 5.x/6.x compatibility
%
% Thomas Winberry, July 26th, 2016
% Updated for Dynare 5.x/6.x compatibility

%----------------------------------------------------------------
% Initialize check variable for error handling
%----------------------------------------------------------------
check = 0;

%----------------------------------------------------------------
% Initialize global price variables
%----------------------------------------------------------------

% Compute initial price guess based on representative agent steady state
global r w

% Representative agent steady state capital
kRepSS = ((aalpha * (aggEmployment ^ (1 - aalpha))) / ((1 / bbeta) - (1 - ddelta))) ^ (1 / (1 - aalpha));

% Initial price guess
r = aalpha * (kRepSS ^ (aalpha - 1)) * (aggEmployment ^ (1 - aalpha)) - ddelta;
w = (kRepSS ^ aalpha) * (1 - aalpha) * (aggEmployment ^ (-aalpha));

fprintf('Initial price guess: r = %.4f, w = %.4f\n', r, w);

%----------------------------------------------------------------
% Compute approximation tools
%----------------------------------------------------------------

try
    % Grids
    computeGrids;
    
    % Polynomials over grids (only if using polynomials to approximate conditional expectation)
    if splineOpt == 0
        computePolynomials;
    end
    
catch ME
    fprintf('Error in computing grids/polynomials: %s\n', ME.message);
    check = 1;
    return;
end

%----------------------------------------------------------------
% Compute initial guess of market-clearing capital stock using
% histogram approximation of distribution, from Young (2010)
%----------------------------------------------------------------

t0 = tic;
fprintf('Computing initial guess from histogram...\n')

try
    % Solve for market clearing capital stock
    f = @(capital) computeMCResidualHistogram(capital);
    
    % Set solver options based on MATLAB version
    if exist('optimoptions', 'file') == 2
        % Modern MATLAB
        options = optimoptions('fsolve','Display',displayOpt,'FunctionTolerance',1e-2,...
            'MaxIterations',1000,'MaxFunctionEvaluations',1000);
    else
        % Older MATLAB
        options = optimset('Display',displayOpt,'TolFun',1e-2,...
            'MaxIter',1000,'MaxFunEvals',1000);
    end
    
    % Compute representative agent steady state capital as initial guess
    kRepSS = ((aalpha * (aggEmployment ^ (1 - aalpha))) / ((1 / bbeta) - (1 - ddelta))) ^ (1 / (1 - aalpha));
    
    [aggregateCapitalInit,err,exitflag] = fsolve(f,1.01*kRepSS,options);
    
    % Return exitflag if market clearing not solved
    if exitflag < 1
        fprintf('Error: Could not solve for initial capital stock. Exit flag: %d\n', exitflag);
        check = 1;
        return; 
    end	
    
    aggregateCapital = aggregateCapitalInit;
    if strcmp(displayOpt,'iter-detailed') == 1
        fprintf('Done! Time to compute: %2.2f seconds \n\n',toc(t0))
    end
    
catch ME
    fprintf('Error in histogram computation: %s\n', ME.message);
    check = 1;
    return;
end

%----------------------------------------------------------------
% Compute moments of histogram to use as initial guess for parametric family
%----------------------------------------------------------------

try
    % Compute histogram
    [~, mHistogram] = computeMCResidualHistogram(aggregateCapital);
    
    % Compute moments from histogram
    mMomentsHistogram = zeros(nEpsilon,nMeasure);
    aGridMoments = zeros(nEpsilon,nAssetsQuadrature,nMeasure); % grid for computing PDF
    
    for iEpsilon = 1 : nEpsilon
        
        % First moment (uncentered)
        mMomentsHistogram(iEpsilon,1) = sum(vAssetsGridFine' .* (mHistogram(iEpsilon,:) ./ ...
            sum(mHistogram(iEpsilon,:))));
        aGridMoments(iEpsilon,:,1) = vAssetsGridQuadrature - mMomentsHistogram(iEpsilon,1);
        
        % Higher order moments (centered)
        for iMoment = 2 : nMeasure
            mMomentsHistogram(iEpsilon,iMoment) = sum(((vAssetsGridFine' - mMomentsHistogram(iEpsilon,1)) .^ iMoment) .* ...
                (mHistogram(iEpsilon,:) ./ sum(mHistogram(iEpsilon,:))));
            aGridMoments(iEpsilon,:,iMoment) = (vAssetsGridQuadrature' - mMomentsHistogram(iEpsilon,1)) .^ ...
                iMoment - mMomentsHistogram(iEpsilon,iMoment);
        end	
        
    end
    
    % Mass at borrowing constraint
    mHatHistogram = [mHistogram(1,1) / sum(mHistogram(1,:));mHistogram(2,1) / sum(mHistogram(2,:))];
    
catch ME
    fprintf('Error in computing moments from histogram: %s\n', ME.message);
    check = 1;
    return;
end

%----------------------------------------------------------------
% Compute market-clearing capital stock from parametric family
%----------------------------------------------------------------

t0 = tic; 
fprintf('Compute steady state from parametric family...\n')

try
    % Solve for market clearing capital stock
    f = @(capital) computeMCResidualPolynomials(capital,mMomentsHistogram,aGridMoments,mHatHistogram);
    
    % Set solver options
    if exist('optimoptions', 'file') == 2
        options = optimoptions('fsolve','Display',displayOpt,'FunctionTolerance',1e-2,...
            'MaxIterations',1000,'MaxFunctionEvaluations',1000);
    else
        options = optimset('Display',displayOpt,'TolFun',1e-2,...
            'MaxIter',1000,'MaxFunEvals',1000);
    end
    
    if abs(f(aggregateCapitalInit)) > 1e-4
        [aggregateCapital,err,exitflag] = fsolve(f,aggregateCapitalInit,options);
    end
    
    % Return error if market clearing not solved
    if exitflag < 1
        fprintf('Error: Could not solve for final capital stock. Exit flag: %d\n', exitflag);
        check = 1;
        return; 
    end	
    
    if strcmp(displayOpt,'iter-detailed') == 1
        fprintf('Done! Time to compute: %2.2f seconds \n\n',toc(t0))
    end
    
catch ME
    fprintf('Error in parametric family computation: %s\n', ME.message);
    check = 1;
    return;
end

%----------------------------------------------------------------
% Compute other objects from steady state
%----------------------------------------------------------------

try
    [~,mCoefficients,mParameters,mMoments,mHat] = ...
        computeMCResidualPolynomials(aggregateCapital,mMomentsHistogram,aGridMoments,mHatHistogram);
    
    % Verify outputs
    if isempty(mCoefficients) || isempty(mParameters) || isempty(mMoments) || isempty(mHat)
        fprintf('Error: Some output variables are empty\n');
        check = 1;
        return;
    end
    
    % Check dimensions
    [n1,n2] = size(mMoments);
    if n1 ~= nEpsilon || n2 ~= nMeasure
        fprintf('Warning: mMoments has unexpected dimensions [%d,%d], expected [%d,%d]\n', n1, n2, nEpsilon, nMeasure);
    end
    
    [n1,n2] = size(mParameters);  
    if n1 ~= nEpsilon || n2 ~= nMeasure+1
        fprintf('Warning: mParameters has unexpected dimensions [%d,%d], expected [%d,%d]\n', n1, n2, nEpsilon, nMeasure+1);
    end
    
    if length(mHat) ~= nEpsilon
        fprintf('Warning: mHat has unexpected length %d, expected %d\n', length(mHat), nEpsilon);
    end
    
    fprintf('Steady state computation completed successfully.\n');
    
catch ME
    fprintf('Error in computing final steady state objects: %s\n', ME.message);
    check = 1;
    return;
end