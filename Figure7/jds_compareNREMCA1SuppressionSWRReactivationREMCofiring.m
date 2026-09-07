%jds_compareNREMCA1SuppressionSWRReactivationREMCofiring
%For CA1 cells with both a noncoordinated-ripple suppression index and a
%coordinated-SWR reactivation index in NREM, relates each cell's REM
%chain-ripple cofiring with PFC to those NREM modulation indices
%(high/low cofiring rank-sum, quartile plot).

clear all
close all
savedir = '/Volumes/JUSTIN/SingleDay/ProcessedDataREM/';

load([savedir 'CA1nrempfcindripmodsigdata.mat'])
allripmodMI_C = allripmodMI;
allripplemod_idx_C = allripplemod_idx;

load([savedir 'CA1nremca1coordripmodsigdata.mat'])
allripmodMI_H = allripmodMI;
allripplemod_idx_H = allripplemod_idx;

animalprefixlist = {'KL8','JS17','JS15','JS14','JS12','JS13','JS34','BG1','JS21','ZT2'};
day = 1;
cofiringMod = [];
cofiringModHists = [];
rateChangeMod = [];

for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);

    idxsH = find(allripplemod_idx_H(:,1) == a);
    animIdxH = allripplemod_idx_H(idxsH,2:end);
    animmodH = allripmodMI_H(idxsH);

    idxsC = find(allripplemod_idx_C(:,1) == a);
    animIdxC = allripplemod_idx_C(idxsC,2:end);
    animmodC = allripmodMI_C(idxsC);

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

        rips = [indrips.ctxripple{day}{ep}.starttimeC indrips.ctxripple{day}{ep}.endtimeC];
        if isempty(rips)
            continue
        end

        if length(rips(:,1)) > 10
            load(sprintf('%s%sspikes%02d.mat',dir,animalprefix,day));
            numncrips = length(rips(:,1));

            CA1idx2 = CA1_R{ep}.cellidx;
            PFCidx = PFC_R{ep}.cellidx;

            epidxH = find(animIdxH(:,2) == ep);
            animIdxEpH = animIdxH(epidxH,3:4);
            animmodEpH = animmodH(epidxH);

            epidxC = find(animIdxC(:,2) == ep);
            animIdxEpC = animIdxC(epidxC,3:4);
            animmodEpC = animmodC(epidxC);

            CA1idx = animIdxEpH;
            CA1idx3 = animIdxEpC;

            for c = 1:length(CA1idx(:,1))
                tmp = [];
                idx1 = find(CA1idx(c,1) == CA1idx2(:,1));
                idx2 = find(CA1idx(c,2) == CA1idx2(:,2));
                idx3 = intersect(idx1,idx2);
                if isempty(idx3)
                    continue
                end
                idx1 = find(CA1idx(c,1) == CA1idx3(:,1));
                idx2 = find(CA1idx(c,2) == CA1idx3(:,2));
                idx3 = intersect(idx1,idx2);
                if isempty(idx3)
                    continue
                end
                Hmod = animmodEpH(c);
                Cmod = animmodEpC(idx3);
                cell1spks = spikes{day}{ep}{CA1idx(c,1)}{CA1idx(c,2)}.data(:,1);
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
                cofiringMod = [cofiringMod; [nanmean(tmp) Cmod Hmod]];
            end
        end
    end
end

cofiringMod2 = cofiringMod(~isnan(cofiringMod(:,1)),:);
[r pcorr] = corrcoef(cofiringMod2)

low = cofiringMod(find(cofiringMod(:,1) < 0),2);
high = cofiringMod(find(cofiringMod(:,1) > 0),2);

[p h] = ranksum(low,high)
datacombinedSuppression = [low; high];
g1 = repmat({'Low cofiring'},length(low),1);
g2 = repmat({'High cofiring'},length(high),1);
g = [g1;g2];

figure; hold on
h = boxplot(datacombinedSuppression,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
xlim([0.5 2.5])
title(['Suppression-p = ' num2str(p)])
ylabel('Suppression')
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
ylabel('CA1 suppression')
xlabel(['Cofiring - Q1vsQ4 p = ' num2str(p2)])
title(['Cofiring CA1 suppression corr r = ' num2str(r(1,2)) ' p = ' num2str(pcorr(1,2))])
set(gcf, 'renderer', 'painters')


