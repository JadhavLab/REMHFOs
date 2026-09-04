function jds_rippleRateSleepState(animalprefixlist)
% Compute ripple rates in NREM and REM sleep and compare them statistically.
%
% Inputs:
%   animalprefixlist - cell array of animal prefix strings
%
% Outputs:
%   Creates a figure showing ripple rate comparison between NREM and REM
%   Returns p-value from statistical comparison
%
% Example:
%   p = jds_rippleRateSleepState({'ER1', 'JS14'})
%
% Dependencies:
%   config.m, utils.m

% Validate inputs
config.validateInputs(animalprefixlist);

% Initialize variables
day = config.DEFAULT_DAY;
nremRates = [];
remRates = [];
animalInfo = [];

try
    % Process each animal
    for animalIdx = 1:length(animalprefixlist)
        animalPrefix = animalprefixlist{animalIdx};
        
        % Load sleep and ripple data
        [remData, nremData] = utils.loadSleepData(animalPrefix, day);
        [remRipples, nremRipples] = utils.loadRippleData(animalPrefix, day);
        epochs = utils.loadEpochs(animalPrefix, day);
        
        % Process each epoch
        for epochIdx = 1:length(epochs)
            epoch = epochs(epochIdx);
            
            % Check if REM data exists for this epoch
            if ~isfield(remData{day}{epoch}, 'total_duration')
                continue;
            end
            
            % Get durations
            nremDuration = nremData{day}{epoch}.total_duration;
            remDuration = remData{day}{epoch}.total_duration;
            
            % Count ripples
            nremRippleCount = length(nremRipples{day}{epoch}.starttime);
            remRippleCount = length(remRipples{day}{epoch}.starttime);
            
            % Calculate rates
            nremRate = utils.calculateEventRate(nremRippleCount, nremDuration);
            remRate = utils.calculateEventRate(remRippleCount, remDuration);
            
            % Store results
            nremRates = [nremRates; nremRate];
            remRates = [remRates; remRate];
            animalInfo = [animalInfo; animalIdx, epoch];
        end
    end
    
    % Perform statistical comparison
    [pValue, stats] = utils.compareRates(nremRates, remRates, 'ranksum');
    
    % Create visualization
    utils.setupFigure('Ripple Rate Comparison');
    
    % Prepare data for boxplot
    allRates = [nremRates; remRates];
    groupLabels = [repmat({'NREM ripple rate'}, length(nremRates), 1); ...
                   repmat({'REM ripple rate'}, length(remRates), 1)];
    
    % Create boxplot
    boxplot(allRates, groupLabels, 'OutlierSize', 7, 'Symbol', 'k+');
    
    % Format plot
    utils.formatPlot('', 'Ripple rate (Hz)', [0.5, 2.5], [-0.2, 1.6]);
    yticks(0:0.4:1.6);
    
    % Add title with p-value
    title(sprintf('Ripple rate (NREM vs REM) - p = %.4f', pValue));
    
    % Display summary statistics
    fprintf('NREM ripple rate: %.3f ± %.3f Hz (n = %d)\n', ...
            mean(nremRates), std(nremRates), length(nremRates));
    fprintf('REM ripple rate: %.3f ± %.3f Hz (n = %d)\n', ...
            mean(remRates), std(remRates), length(remRates));
    fprintf('Statistical test: Wilcoxon rank-sum, p = %.4f\n', pValue);
    
catch ME
    error('Error in jds_rippleRateSleepState: %s', ME.message);
end

end
