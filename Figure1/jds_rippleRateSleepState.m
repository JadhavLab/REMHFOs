function jds_rippleRateSleepState(animalprefixlist)
%JDS_RIPPLERATESLEEPSTATE Compare cortical ripple rates in NREM vs REM.
%   jds_rippleRateSleepState(animalprefixlist) computes the cortical ripple
%   rate (events/s) within NREM and REM for each sleep epoch, compares the
%   two states (rank-sum), and plots a boxplot.
%
%   animalprefixlist - cell array of animal prefix strings

day = 1;

nremRates = [];
remRates = [];

for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);

    remrips = load(sprintf('%s%sctxrippletime_REM0%d.mat',dir,animalprefix,day));
    nremrips = load(sprintf('%s%sctxrippletime_SWS0%d.mat',dir,animalprefix,day));
    rem = load(sprintf('%s%srem0%d.mat',dir,animalprefix,day));
    rem = rem.rem;
    load(sprintf('%s%sswsALL0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%sremeps0%d.mat',dir,animalprefix,day));
    epochs = remeps;

    for e = 1:length(epochs)
        epoch = epochs(e);

        if isempty(rem{day}{epoch}.starttime)
            continue
        end

        nremDur = sws{day}{epoch}.total_duration;
        remDur = rem{day}{epoch}.total_duration;

        nremCount = length(nremrips.ctxripple{day}{epoch}.starttime);
        remCount = length(remrips.ctxripple{day}{epoch}.starttime);

        nremRates = [nremRates; nremCount/nremDur];
        remRates = [remRates; remCount/remDur];
    end
end

[p h] = ranksum(nremRates,remRates)

datacombinedRates = [nremRates; remRates];
g1 = repmat({'NREM ripple rate'},length(nremRates),1);
g2 = repmat({'REM ripple rate'},length(remRates),1);
g = [g1;g2];

figure;
h = boxplot(datacombinedRates,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
xlim([0.5 2.5])
ylim([-0.2 1.6])
yticks([0:0.4:1.6])
ylabel('Ripple rate (Hz)')
title(['Ripple rate (NREM vs REM) - p = ' num2str(p)])
set(gcf, 'renderer', 'painters')
