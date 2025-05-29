% Computes and analyzes steady state with no aggregate shocks
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
    % Compute Steady State
    %----------------------------------------------------------------
    
    fprintf('=== Computing Steady State ===\n');
    
    % Solve for steady state capital stock, distribution, and decision rules
    coreSteadyState;
    
    % Check if computation was successful
    if exist('check', 'var') && check == 1
        error('Steady state computation failed');
    end
    
    % Compute decision rules along fine grid for analysis
    [~,mHistogram,mAssetsPrime,mConsumption] = computeMCResidualHistogram(aggregateCapital);
    
    % Compute density along fine grid
    mDistributionFine = zeros(nEpsilon,nAssetsFine);
    for iEpsilon = 1 : nEpsilon
        
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
        
        % Mass at borrowing constraint
        %mDistributionFine(iEpsilon,1) = mHat(iEpsilon,1); % Commented out for now; need fine quadrature grid to capture correctly
        
    end
    
    fprintf('=== Steady State Computation Completed ===\n');
    
    %----------------------------------------------------------------
    % Plot results 
    %----------------------------------------------------------------
    
    fprintf('=== Generating Plots ===\n');
    
    % Create figure directory if it doesn't exist
    if ~exist('../Figures', 'dir')
        mkdir('../Figures');
    end
    
    % Savings function
    figure('Position', [100, 100, 800, 600]);
    hold on
    plot(vAssetsGridFine,mAssetsPrime(1,:),'linewidth',2,'color',[0.7,0.1,0.1],'DisplayName','Unemployed')
    plot(vAssetsGridFine,mAssetsPrime(2,:),'linewidth',2,'color',[0.0,0.2,0.5],'DisplayName','Employed')
    plot(vAssetsGridFine,vAssetsGridFine,'k--','linewidth',1,'DisplayName','45° line')
    xlabel('Assets, $a$','interpreter','latex','FontSize',14)
    ylabel('Savings, $s(\varepsilon,a)$','interpreter','latex','FontSize',14)
    xlim([aaBar 0.9*assetsMax])
    title('Savings Decision Rule','FontSize',16)
    legend('Location','southeast','FontSize',12)
    grid on
    set(gcf,'color','w')
    hold off
    
    % Save figure
    saveas(gcf, '../Figures/savings_policy.png');
    fprintf('Saved: ../Figures/savings_policy.png\n');
    
    % Consumption function
    figure('Position', [150, 150, 800, 600]);
    hold on
    plot(vAssetsGridFine,mConsumption(1,:),'linewidth',2,'color',[0.7,0.1,0.1],'DisplayName','Unemployed')
    plot(vAssetsGridFine,mConsumption(2,:),'linewidth',2,'color',[0.0,0.2,0.5],'DisplayName','Employed')
    xlabel('Assets, $a$','interpreter','latex','FontSize',14)
    ylabel('Consumption, $c(\varepsilon,a)$','interpreter','latex','FontSize',14)
    xlim([aaBar 0.9*assetsMax])
    title('Consumption Decision Rule','FontSize',16)
    legend('Location','southeast','FontSize',12)
    grid on
    set(gcf,'color','w')
    hold off
    
    % Save figure
    saveas(gcf, '../Figures/consumption_policy.png');
    fprintf('Saved: ../Figures/consumption_policy.png\n');
    
    % Distribution of unemployed households
    figure('Position', [200, 200, 800, 600]);
    hold on
    plot(vAssetsGridFine,mHistogram(1,:) / sum(mHistogram(1,:)),...
        'linewidth',2,'color',[0.0,0.2,0.5],'DisplayName','Histogram')
    plot(vAssetsGridFine,mDistributionFine(1,:) ./ sum(mDistributionFine(1,:)),...
        'linewidth',2,'color',[0.7,0.1,0.1],'linestyle','--','DisplayName','Parametric Family')
    xlabel('Assets, $a$','interpreter','latex','FontSize',14)
    ylabel('Mass of households, $g(\varepsilon,a)$','interpreter','latex','FontSize',14)
    xlim([aaBar 0.9*assetsMax])
    title('Invariant Distribution of Households (Unemployed)','FontSize',16)
    legend('Location','northeast','FontSize',12)
    grid on
    set(gcf,'color','w')
    hold off
    
    % Save figure
    saveas(gcf, '../Figures/distribution_unemployed.png');
    fprintf('Saved: ../Figures/distribution_unemployed.png\n');
    
    % Distribution of employed households
    figure('Position', [250, 250, 800, 600]);
    hold on
    plot(vAssetsGridFine,mHistogram(2,:) / sum(mHistogram(2,:)),...
        'linewidth',2,'color',[0.0,0.2,0.5],'DisplayName','Histogram')
    plot(vAssetsGridFine,mDistributionFine(2,:) ./ sum(mDistributionFine(2,:)),...
        'linewidth',2,'color',[0.7,0.1,0.1],'linestyle','--','DisplayName','Parametric Family')
    xlabel('Assets, $a$','interpreter','latex','FontSize',14)
    ylabel('Mass of households, $g(\varepsilon,a)$','interpreter','latex','FontSize',14)
    xlim([aaBar 0.9*assetsMax])
    title('Invariant Distribution of Households (Employed)','FontSize',16)
    legend('Location','northeast','FontSize',12)
    grid on
    set(gcf,'color','w')
    hold off
    
    % Save figure
    saveas(gcf, '../Figures/distribution_employed.png');
    fprintf('Saved: ../Figures/distribution_employed.png\n');
    
    %----------------------------------------------------------------
    % Display key results
    %----------------------------------------------------------------
    
    fprintf('\n=== KEY RESULTS ===\n');
    fprintf('Aggregate Capital: %.4f\n', aggregateCapital);
    fprintf('Interest Rate: %.4f (%.2f%%)\n', r, r*100);
    fprintf('Wage: %.4f\n', w);
    
    % Compute some distributional statistics
    meanAssets = [sum(vAssetsGridFine .* (mHistogram(1,:) / sum(mHistogram(1,:)))), ...
                  sum(vAssetsGridFine .* (mHistogram(2,:) / sum(mHistogram(2,:))))];
    
    fprintf('\nDistributional Statistics:\n');
    fprintf('Mean assets (unemployed): %.4f\n', meanAssets(1));
    fprintf('Mean assets (employed): %.4f\n', meanAssets(2));
    fprintf('Mass at borrowing constraint (unemployed): %.4f\n', mHat(1));
    fprintf('Mass at borrowing constraint (employed): %.4f\n', mHat(2));
    
    fprintf('\n=== All computations completed successfully ===\n');
    
catch ME
    fprintf('\n=== Error occurred ===\n');
    fprintf('Error in %s at line %d\n', ME.stack(1).name, ME.stack(1).line);
    fprintf('Error message: %s\n', ME.message);
    
    % Return to original directory before rethrowing
    cd(oldFolder);
    rethrow(ME);
end

% Return to original directory
cd(oldFolder);