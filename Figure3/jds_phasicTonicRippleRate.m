function jds_phasicTonicRippleRate(animalprefixlist)
%JDS_PHASICTONICRIPPLERATE Cortical ripple rate in phasic vs tonic REM.
%   jds_phasicTonicRippleRate(animalprefixlist) compares the rate of
%   cortical REM ripples within phasic (high-theta) REM bouts vs the
%   remaining (tonic) REM time (sign-rank, boxplot).
%
%   animalprefixlist - cell array of animal prefix strings

day = 1;

phasicRate = [];
tonicRate = [];

phasicTime = [];
tonicTime = [];
totalDur = [];
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);

    load(sprintf('%s%sctxrippletime_REM0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%sphasicrembouts0%d.mat',dir,animalprefix,day)); %high-theta REM bouts

    rem = load(sprintf('%s%srem0%d.mat',dir,animalprefix,day));
    rem = rem.rem;
    load(sprintf('%s%sremeps0%d.mat',dir,animalprefix,day));
    epochs = remeps;

    for ep = 1:length(epochs)
        epoch = epochs(ep);

        ripmidtimes = (ctxripple{day}{epoch}.endtime - ctxripple{day}{epoch}.starttime)/2 +...
            ctxripple{day}{epoch}.starttime;

        phasicList = [phasicrem{day}{epoch}.starttime phasicrem{day}{epoch}.endtime];

        inPhasic = isExcluded(ripmidtimes, phasicList); %in phasic
        inTonic = ~isExcluded(ripmidtimes, phasicList); %not in phasic

        remdur = rem{day}{epoch}.total_duration;
        phasicDur = phasicrem{day}{epoch}.total_duration;
        tonicDur = remdur - phasicDur;
        totalDur = [totalDur; remdur];

        phasicRate = [phasicRate; sum(inPhasic)/phasicDur];
        tonicRate = [tonicRate; sum(inTonic)/tonicDur];
        phasicTime = [phasicTime; phasicDur];
        tonicTime = [tonicTime; tonicDur];
    end
end

p = signrank(phasicRate,tonicRate)

datacombinedRate = [phasicRate; tonicRate];
g1 = repmat({'High theta REM'},length(phasicRate),1);
g2 = repmat({'Low theta REM'},length(tonicRate),1);
g = [g1;g2];

figure
h = boxplot(datacombinedRate,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
title(['Ripple rate High Theta vs Low Theta REM-p = ' num2str(p)])

