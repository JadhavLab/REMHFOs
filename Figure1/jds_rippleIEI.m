function jds_rippleIEI(animalprefixlist)
%JDS_RIPPLEIEI Distributions of cortical ripple inter-event intervals.
%   jds_rippleIEI(animalprefixlist) computes inter-event intervals between
%   successive cortical ripple start times within REM bouts and within NREM
%   bouts, compares the two distributions (rank-sum), and plots probability
%   histograms. Only IEIs up to 1 s are kept.
%
%   animalprefixlist - cell array of animal prefix strings

day = 1;
maxIEI = 1; %seconds

nremIEI = [];
remIEI = [];

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

        %REM bouts
        if ~isempty(rem{day}{epoch}.starttime)
            remList = [rem{day}{epoch}.starttime rem{day}{epoch}.endtime];
            ripTimes = remrips.ctxripple{day}{epoch}.starttime;
            for r = 1:size(remList,1)
                boutRips = ripTimes(logical(isExcluded(ripTimes, remList(r,:))));
                if length(boutRips) > 1
                    iei = diff(boutRips);
                    remIEI = [remIEI; iei(iei <= maxIEI)];
                end
            end
        end

        %NREM bouts
        if ~isempty(sws{day}{epoch}.starttime)
            nremList = [sws{day}{epoch}.starttime sws{day}{epoch}.endtime];
            ripTimes = nremrips.ctxripple{day}{epoch}.starttime;
            for r = 1:size(nremList,1)
                boutRips = ripTimes(logical(isExcluded(ripTimes, nremList(r,:))));
                if length(boutRips) > 1
                    iei = diff(boutRips);
                    nremIEI = [nremIEI; iei(iei <= maxIEI)];
                end
            end
        end
    end
end

[p h] = ranksum(remIEI,nremIEI)

figure; hold on
histogram(remIEI,50,'DisplayStyle','stairs','Normalization','probability');
histogram(nremIEI,50,'DisplayStyle','stairs','Normalization','probability');
legend({'REM','NREM'})
xlim([-0.05 maxIEI])
xlabel('Inter-event interval (s)')
ylabel('Probability')
title(['Ripple IEI Distribution (p = ' num2str(p) ')'])
set(gcf, 'renderer', 'painters')
