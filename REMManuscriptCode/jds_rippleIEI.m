function jds_rippleIEI(animalprefixlist)
% Calculate and plot distributions of ripple inter-event-intervals.
% Only plots IEIs up to 1 second for visualization purposes.
%
% Inputs:
%   animalprefixlist - cell array of animal prefix strings
%
% Outputs:
%   Creates a figure showing IEI distributions for REM and NREM states
%
% Example:
%   jds_rippleIEI({'ER1', 'JS14'})
%
% Dependencies:
%   config.m, utils.m

% Validate inputs
config.validateInputs(animalprefixlist);

% Initialize variables
day = config.DEFAULT_DAY;
nremIEI = [];
remIEI = [];

% Define sleep states to process
sleepStates = {'REM', 'NREM'};

try
    % Process each sleep state
    for stateIdx = 1:length(sleepStates)
        state = sleepStates{stateIdx};
        
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
                
                if strcmp(state, 'REM')
                    % Process REM data
                    if ~isfield(remData{day}{epoch}, 'starttime')
                        continue;
                    end
                    
                    remList = [remData{day}{epoch}.starttime, remData{day}{epoch}.endtime];
                    rippleTimes = remRipples{day}{epoch}.starttime;
                    
                    % Calculate IEIs for each REM period
                    for remIdx = 1:size(remList, 1)
                        remWindow = remList(remIdx, :);
                        remRippleTimes = rippleTimes(logical(isExcluded(rippleTimes(:, 1), remWindow)), 1);
                        if length(remRippleTimes) > 1
                            iei = utils.calculateIEI(remRippleTimes, config.RIPPLE_MAX_IEI);
                            remIEI = [remIEI; iei];
                        end
                    end
                    
                else
                    % Process NREM data
                    if ~isfield(nremData{day}{epoch}, 'starttime')
                        continue;
                    end
                    
                    nremList = [nremData{day}{epoch}.starttime, nremData{day}{epoch}.endtime];
                    rippleTimes = nremRipples{day}{epoch}.starttime;
                    
                    % Calculate IEIs for each NREM period
                    for nremIdx = 1:size(nremList, 1)
                        nremWindow = nremList(nremIdx, :);
                        nremRippleTimes = rippleTimes(logical(isExcluded(rippleTimes(:, 1), nremWindow)), 1);
                        if length(nremRippleTimes) > 1
                            iei = utils.calculateIEI(nremRippleTimes, config.RIPPLE_MAX_IEI);
                            nremIEI = [nremIEI; iei];
                        end
                    end
                end
            end
        end
    end
    
    % Create visualization
    utils.setupFigure('Ripple Inter-Event Intervals');
    
    % Plot histograms
    hold on;
    histogram(remIEI, 50, 'DisplayStyle', 'stairs', 'Normalization', 'probability', ...
             'DisplayName', 'REM', 'EdgeColor', 'blue');
    histogram(nremIEI, 50, 'DisplayStyle', 'stairs', 'Normalization', 'probability', ...
             'DisplayName', 'NREM', 'EdgeColor', 'red');
    
    % Format plot
    utils.formatPlot('Inter-event interval (s)', 'Probability', [-0.05, config.RIPPLE_MAX_IEI]);
    legend('Location', 'best');
    
    % Add statistics
    [pValue, ~] = utils.compareRates(remIEI, nremIEI, 'ranksum');
    title(sprintf('Ripple IEI Distribution (p = %.4f)', pValue));
    
catch ME
    error('Error in jds_rippleIEI: %s', ME.message);
end

end