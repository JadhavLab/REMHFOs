clear all
savedir = '/Volumes/JUSTIN/SingleDay/ProcessedDataREM/';
% load([savedir 'Allanim_250ctxnremnoncoordripplemod400_0mscrit_sleep_CA1_alldata_largewin_sepeps_gather_X6.mat'])
% load([savedir 'Allanim_250swrnremcoordripplemodbwin200_0mscrit_sleep_CA1_alldata_largewin_sepeps_gather_X6'])
load([savedir 'CA1nrempfcindripmodsigdata.mat'])
% load([savedir 'CA1nremca1allripmodsigdata.mat'])
% load([savedir 'CA1nremca1coordripmodsigdata200bwin.mat'])
%
% allripplemod_idx = [];
% allripmodhists = [];
% allripmodMI = [];
% allripmodPeak = [];
% for i=1:length(allripplemod)
%     if allripplemod(i).rasterShufP2<0.05
%         allripplemod_idx=[allripplemod_idx;allripplemod(i).index];
%         allripmodhists=[allripmodhists; zscore(filtfilt(b,1,mean(rast2mat_lrg(allripplemod(i).raster))))];
%         allripmodMI = [allripmodMI; allripplemod(i).Dm];
%         allripmodPeak = [allripmodPeak; allripplemod(i).modln_div];
%     end
% end
% 
% clearvars -except allripmodMI allripmodhists allripplemod_idx allripmodPeak

b=gaussian(20,61);
day = 1;

animalprefixlist = {'KL8','JS17','JS15','JS14','JS12','JS13','JS34','BG1','JS21','ZT2'};

cofiringMod = [];
cofiringModHists = [];
rateChangeMod = [];

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

    epochs = remeps;

    for e = 1:length(epochs)

        ep = epochs(e);
        if ep == 1
            continue
        end

        epidx = find(animIdx(:,2) == ep);
        animIdxEp = animIdx(epidx,3:4);
        animhistsEp = animhists(epidx,:);
        animmodEp = animmod(epidx);

        rips = [indrips.ctxripple{day}{ep}.starttimeC indrips.ctxripple{day}{ep}.endtimeC];
        if isempty(rips)
            continue
        end
        firstsws = [sws{day}{ep}.starttime(1) sws{day}{ep}.endtime(1)];
        lastsws = [sws{day}{ep}.starttime(end) sws{day}{ep}.endtime(end)];
        allrem = [rem{day}{ep}.starttime rem{day}{ep}.endtime];
        remdur = rem{day}{ep}.total_duration;
        %         if ((firstsws(2)-firstsws(1) < 30)) || ((lastsws(2)-lastsws(1) < 30))
        %             continue
        %         end
        if length(rips(:,1)) > 10
            load(sprintf('%s%sspikes%02d.mat',dir,animalprefix,day));
            numncrips = length(rips(:,1));

            CA1idx2 = CA1_R{ep}.cellidx;
            CA1idx = animIdxEp;
            PFCidx = PFC_R{ep}.cellidx;

            for c = 1:length(CA1idx(:,1))
                tmp = [];
                idx1 = find(CA1idx(c,1) == CA1idx2(:,1));
                idx2 = find(CA1idx(c,2) == CA1idx2(:,2));
                idx3 = intersect(idx1,idx2);
                if isempty(idx3)
                    continue
                end
                cell1spks = spikes{day}{ep}{CA1idx(c,1)}{CA1idx(c,2)}.data(:,1);
                firstnremFR = (sum(isExcluded(cell1spks, firstsws)))/(firstsws(2)-firstsws(1));
                lastnremFR = (sum(isExcluded(cell1spks, lastsws)))/(lastsws(2)-lastsws(1));
                remFR = (sum(isExcluded(cell1spks, allrem)))/remdur;
                FRchange = lastnremFR - firstnremFR;
                for pp = 1:length(PFCidx(:,1))
                    if (~isempty(spikes{day}{ep}{CA1idx(c,1)}{CA1idx(c,2)})) &&...
                            (~isempty(spikes{day}{ep}{PFCidx(pp,1)}{PFCidx(pp,2)}))
                        cell1spks = spikes{day}{ep}{CA1idx(c,1)}{CA1idx(c,2)}.data(:,1);
                        cell2spks = spikes{day}{ep}{PFCidx(pp,1)}{PFCidx(pp,2)}.data(:,1);

                        spkbins1_c = periodAssign(cell1spks, rips);
                        spkbins2_p = periodAssign(cell2spks, rips);

                        ripnum1_c = unique(spkbins1_c);
                        activeinrip1_c = ripnum1_c(find(ripnum1_c ~= 0));

                        ripnum2_p = unique(spkbins2_p);
                        activeinrip2_p = ripnum2_p(find(ripnum2_p ~= 0));

                        common_hprip = length(find(ismember(activeinrip1_c, activeinrip2_p)));

                        %calculate zscored ripple coactivity

                        nAB_hp = common_hprip;
                        nA_hp = length(activeinrip1_c);
                        nB_hp = length(activeinrip2_p);

                        coact_hprip = (nAB_hp - (nA_hp*nB_hp/numncrips))/...
                            sqrt(nA_hp*nB_hp*(numncrips - nA_hp)*(numncrips - nB_hp)/...
                            (numncrips^2*(numncrips-1)));

                        tmp = [tmp; coact_hprip];
                    end
                end
                cofiringMod = [cofiringMod; [nanmean(tmp) animmodEp(c)]];
                cofiringModHists = [cofiringModHists; animhistsEp(c,:)];
                rateChangeMod = [rateChangeMod; [FRchange animmodEp(c)]];
            end
        end
    end
end

cofiringMod2 = cofiringMod(~isnan(cofiringMod(:,1)),:);
[r pcorr] = corrcoef(cofiringMod2)
% figure
% scatter(cofiringMod2(:,1),cofiringMod2(:,2));
% xlabel('Cofiring')
% ylabel('NREM PFC ripple suppression')
% hold on; lsline
% title(['Cofiring CA1 suppression corr r = ' num2str(r(1,2)) ' p = ' num2str(p(1,2))])
% 
% low = cofiringMod(find(cofiringMod(:,1) < nanmean(cofiringMod(:,1))),2);
% high = cofiringMod(find(cofiringMod(:,1) > nanmean(cofiringMod(:,1))),2);
low = cofiringMod(find(cofiringMod(:,1) < 0),2);
high = cofiringMod(find(cofiringMod(:,1) > 0),2);
% lowHist = cofiringModHists(find(cofiringMod(:,1) < nanmean(cofiringMod(:,1))),:);
% highHist = cofiringModHists(find(cofiringMod(:,1) > nanmean(cofiringMod(:,1))),:);
% nanHist = cofiringModHists(find(cofiringMod(:,1) < nanmean(cofiringMod(:,1))),:);

[p h] = ranksum(low,high)
datacombinedSuppression = [low; high];
g1 = repmat({'Low cofiring'},length(low),1);
g2 = repmat({'High cofiring'},length(high),1);
g = [g1;g2];

figure; hold on
h = boxplot(datacombinedSuppression,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
xlim([0.5 2.5])
% ylim([-0.02 0.2])
% yticks([0.11:0.01:0.14])
title(['Suppression-p = ' num2str(p)])
ylabel('Suppression')
set(gcf, 'renderer', 'painters')

lowHist = cofiringModHists(find(cofiringMod(:,1) < 0),:);
highHist = cofiringModHists(find(cofiringMod(:,1) > 0),:);
nanHist = cofiringModHists(find(isnan(cofiringMod(:,1))),:);
lowHist = [lowHist; nanHist];

figure; hold on
xaxis=-1049:1050;
pl1 = plot(xaxis,mean(lowHist,1),'-r')
boundedline(xaxis,mean(lowHist,1),std(lowHist)./sqrt(size(lowHist,1)),'-r');
pl2 = plot(xaxis,mean(highHist,1),'-b')
boundedline(xaxis,mean(highHist,1),std(highHist)./sqrt(size(highHist,1)),'-b');
xlim([-250 250])
xticks([-250:250:250])
x = [0 0];
y = [-1.2 0.2];
plot(x,y,'--k')
ylabel('Mean z-scored psth')
xlabel('Time from noncoordinated NREM ripple (ms)')
legend([pl1 pl2],{'Other','High cofiring'});
title('CA1 suppression')
set(gcf, 'renderer', 'painters')

quartile_sep = floor(length(cofiringMod2(:,1))/4);
spklatquar = sortrows(cofiringMod2,1);
vals = [];
cnt = 1;
for s = 1:4
    if s < 4
        tmp = spklatquar(cnt:quartile_sep*s,1:2);
        tmp(:,3) = s;
    else
        tmp = spklatquar(cnt:end,1:2);
        tmp(:,3) = s;
    end
    vals{s} = tmp;
    cnt = cnt + quartile_sep;
    clear tmp
end

v = cellfun(@mean,vals,'UniformOutput',false);

v2 = cellfun((@(x) std(x)./sqrt(length(x(:,1)))),vals,'UniformOutput',false);

data_sems = vertcat(v2{:});

data_means = vertcat(v{:});

X = [1:4];

[p2 h2] = ranksum(vals{1}(:,2),vals{4}(:,2))

figure
errorbar(X, data_means(:,2), data_sems(:,2),'b','LineWidth',3,'Capsize',0);
hold on
xlim([0.5 4.5])
xticks([1:4])
ylabel('CA1 Coord SWR mod')
xlabel(['Cofiring quartile - Q1vsQ4 p = ' num2str(p2)])
title(['Cofiring CA1 suppression corr r = ' num2str(r(1,2)) ' p = ' num2str(pcorr(1,2))])
set(gcf, 'renderer', 'painters')

keyboard

