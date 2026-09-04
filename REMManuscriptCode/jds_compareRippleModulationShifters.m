clear; savedir = '/Volumes/JUSTIN/SingleDay/ProcessedDataREM/';
% load([savedir 'Allanim_250ctxnremnoncoordripplemod400_0mscrit_sleep_CA1_alldata_largewin_sepeps_gather_X6.mat'])
% load([savedir 'CA1nrempfcindripmodsigdata.mat'])
load([savedir 'CA1nremca1allripmodsigdata.mat'])
% load([savedir 'CA1nremca1indripmodsigdata.mat'])

% 
% 
% allripplemod_idx = [];
% allripmodhists = [];
% allripmodMI = [];
% for i=1:length(allripplemod)
%     allripplemod_idx=[allripplemod_idx;allripplemod(i).index];
%     allripmodhists=[allripmodhists; zscore(filtfilt(b,1,mean(rast2mat_lrg(allripplemod(i).raster))))];
%     allripmodMI = [allripmodMI; allripplemod(i).Dm];
% end
% 
% clearvars -except allripmodMI allripmodhists allripplemod_idx

b=gaussian(20,61);
day = 1;

animalprefixlist = {'KL8','JS17','JS15','JS14','JS12','JS13','JS34','BG1','JS21','ZT2'};

shiftSupp = [];
nonShiftSupp = [];

shiftHist = [];
nonShiftHist = [];

for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
    idxs = find(allripplemod_idx(:,1) == a);
    animIdx = allripplemod_idx(idxs,2:end);
    animhists = allripmodhists(idxs,:);
    animmod = allripmodMI(idxs);

    load(sprintf('%s%sremeps0%d.mat',dir,animalprefix,day));% get sws time
    epochs = remeps;

    load(sprintf('%s%sCA1_RTimeStrengthSleepNewSpk_20_%02d.mat',dir,animalprefix,day));
    CA1_R = RtimeStrength; clear RtimeStrength
    load(sprintf('%s%sPFC_RTimeStrengthSleepNewSpk_20_%02d.mat',dir,animalprefix,day));
    PFC_R = RtimeStrength; clear RtimeStrength
    load(sprintf('%s%sremeps%02d.mat',dir,animalprefix,day));
    load(sprintf('%s%sswsALL%02d.mat',dir,animalprefix,day));
    rem = load(sprintf('%s%srem%02d.mat',dir,animalprefix,day));
    rem = rem.rem;

    %Load ripples
    indrips = load(sprintf('%s%sctxrippletime_chainREM%02d.mat',dir,animalprefix,day));
    load(sprintf('%s/%sCA1remshiftershighthresh%02d.mat',dir,animalprefix,day));

    epochs = remeps;

    for e = 1:length(epochs)

        ep = epochs(e);
        if ep == 1
            continue
        end

        shifters = shiftList{day}{ep}.cellidx;

        epidx = find(animIdx(:,2) == ep);
        animIdxEp = animIdx(epidx,3:4);
        animhistsEp = animhists(epidx,:);
        animmodEp = animmod(epidx);
        
        rips = [indrips.ctxripple{day}{ep}.starttimeC indrips.ctxripple{day}{ep}.endtimeC];
        if isempty(rips)
            continue
        end
        if length(rips(:,1)) > 10
            load(sprintf('%s%sspikes%02d.mat',dir,animalprefix,day));
            numncrips = length(rips(:,1));

            CA1idx = animIdxEp;

            for c = 1:length(CA1idx(:,1))
                tmp = [];
                idx1 = find(CA1idx(c,1) == shifters(:,1));
                idx2 = find(CA1idx(c,2) == shifters(:,2));
                idx3 = intersect(idx1,idx2);
                if isempty(idx3)
                    continue
                end

                isShift = shifters(idx3,3);

                if isShift == 1
                    shiftSupp = [shiftSupp; animmodEp(c)];
                    shiftHist = [shiftHist; animhistsEp(c,:)];
                else
                    nonShiftSupp = [nonShiftSupp; animmodEp(c)];
                    nonShiftHist = [nonShiftHist; animhistsEp(c,:)];
                end
            end
        end
    end
end

[p h] = ranksum(shiftSupp,nonShiftSupp)
datacombinedShift = [nonShiftSupp; shiftSupp];
g1 = repmat({'Non shift'},length(nonShiftSupp),1);
g2 = repmat({'Shift'},length(shiftSupp),1);
g = [g1;g2];

figure; hold on
h = boxplot(datacombinedShift,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-0.5 0.3])
xlim([0.5 2.5])
title(['Shifters NREM ripple suppression-p = ' num2str(p)])
ylabel('CA1 suppression')
set(gcf, 'renderer', 'painters')

figure; hold on
xaxis=-1049:1050;
pl1 = plot(xaxis,mean(shiftHist,1),'-r')
boundedline(xaxis,mean(shiftHist,1),std(shiftHist)./sqrt(size(shiftHist,1)),'-r');
pl2 = plot(xaxis,mean(nonShiftHist,1),'-b')
boundedline(xaxis,mean(nonShiftHist,1),std(nonShiftHist)./sqrt(size(nonShiftHist,1)),'-b');
xlim([-250 250])
xticks([-250:250:250])
x = [0 0];
y = [-1.2 0.2];
plot(x,y,'--k')
ylabel('Mean z-scored psth')
xlabel('Time from noncoordinated NREM ripple (ms)')
legend([pl1 pl2],{'Shift','Non shift'});
title('CA1 Modulation - NREM All SWRs')
set(gcf, 'renderer', 'painters')

keyboard

