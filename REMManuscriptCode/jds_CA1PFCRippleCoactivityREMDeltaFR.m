function jds_CA1PFCRippleCoactivityREMDeltaFR(animalprefixlist)

%To do: plot the peak strength within specified ripple events. compare
%ripple types for different assemblies - might need to change binsize for
%raster generation
day = 1;

ncRippleCo = [];
cRippleCo = [];
firingRates = [];
firingRatesWithRem = [];
meanFRLow = [];
meanFRHigh = [];
upFR = [];
downFR = [];

for a = 1:length(animalprefixlist)

    animalprefix = char(animalprefixlist(a));
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    %Load reactivation strength file for all assemblies and epochs
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
        firstsws = [sws{day}{ep}.starttime(1) sws{day}{ep}.endtime(1)];
        lastsws = [sws{day}{ep}.starttime(end) sws{day}{ep}.endtime(end)];
        allrem = [rem{day}{ep}.starttime rem{day}{ep}.endtime];
        remdur = rem{day}{ep}.total_duration;
        if ((firstsws(2)-firstsws(1) < 30)) || ((lastsws(2)-lastsws(1) < 30))
            continue
        end
        [ctxidx, hpidx] = jds_getallepcells(dir, animalprefix, day, ep, []);
        CA1assemblytmp = CA1_R{ep}.reactivationStrength;
        CA1num = length(CA1assemblytmp);
        PFCassemblytmp = PFC_R{ep}.reactivationStrength;
        PFCnum = length(PFCassemblytmp);

%         ctxrip = indrips.ctxripple{day}{ep}.C_sep;
%         tmp = [];
%         for c = 1:length(ctxrip)
%             tmp = [tmp; [ctxrip{c}(2,1) ctxrip{c}(2,2)]]; 
% %             sub = ctxrip{c};
% %             sub(2,:) = [];
% %             tmp = [tmp; [ctxrip{c}(2:end,1) ctxrip{c}(2:end,2)]]; 
% %             tmp = [tmp; sub]; 
%         end

        ncrips = [indrips.ctxripple{day}{ep}.starttimeC indrips.ctxripple{day}{ep}.endtimeC];
%         ncrips = tmp;
        if isempty(ncrips)
            continue
        end
        if length(ncrips(:,1)) > 10
            load(sprintf('%s%sspikes%02d.mat',dir,animalprefix,day));
            numncrips = length(ncrips(:,1));

            CA1idx = CA1_R{ep}.cellidx;
            PFCidx = PFC_R{ep}.cellidx;

            Chigh = CA1idx;
%             Chigh = PFCidx;
            Phigh = PFCidx;
%             Phigh = CA1idx;
            for c = 1:length(Chigh(:,1))
                tmp = [];
                cell1spks = spikes{day}{ep}{Chigh(c,1)}{Chigh(c,2)}.data(:,1);
                firstnremFR = (sum(isExcluded(cell1spks, firstsws)))/(firstsws(2)-firstsws(1));
                lastnremFR = (sum(isExcluded(cell1spks, lastsws)))/(lastsws(2)-lastsws(1));
                remFR = (sum(isExcluded(cell1spks, allrem)))/remdur;
                FRchange = lastnremFR - firstnremFR;
                meanRt = spikes{day}{ep}{Chigh(c,1)}{Chigh(c,2)}.meanrate;
                if FRchange < -5
                    continue
                end
                firingRates = [firingRates; [firstnremFR lastnremFR]];
                firingRatesWithRem = [firingRatesWithRem; [firstnremFR remFR lastnremFR]];
                for pp = 1:length(Phigh(:,1))
                    if (~isempty(spikes{day}{ep}{Chigh(c,1)}{Chigh(c,2)})) &&...
                            (~isempty(spikes{day}{ep}{Phigh(pp,1)}{Phigh(pp,2)}))
                        cell1spks = spikes{day}{ep}{Chigh(c,1)}{Chigh(c,2)}.data(:,1);
                        cell2spks = spikes{day}{ep}{Phigh(pp,1)}{Phigh(pp,2)}.data(:,1);

                        spkbins1_c = periodAssign(cell1spks, ncrips);
                        spkbins2_p = periodAssign(cell2spks, ncrips);

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
                ncRippleCo = [ncRippleCo; [nanmean(tmp) FRchange meanRt remFR]];
                if FRchange > -0.0447% mean = -0.0363 %median = -0.0447
                    upFR = [upFR; nanmean(tmp)];
                elseif FRchange < -0.0447
                    downFR = [downFR; nanmean(tmp)];
                end
            end
        end
    end
end
p = friedman(firingRatesWithRem,1,'on')
[p, tbl, stats] = friedman(firingRatesWithRem, 1, 'off');
[COMPARISON,MEANS,H,GNAMES] = multcompare(stats, 'CType', 'bonferroni');
meanFrs = mean(firingRatesWithRem);
semFrs = std(firingRatesWithRem)./sqrt(size(firingRatesWithRem,1));
figure
errorbar(1:3,meanFrs,semFrs,'-k')
xlim([0.5 3.5])
xticks([1:3])
xticklabels({'First NREM','REM','Last NREM'})
set(gcf, 'renderer', 'painters')

datacombinedEpFrs = [firingRatesWithRem(:,1); firingRatesWithRem(:,2); firingRatesWithRem(:,3)];
g1 = repmat({'First NREM'},length(firingRatesWithRem(:,1)),1);
g2 = repmat({['REM p=' num2str(COMPARISON(1,6))]},length(firingRatesWithRem(:,2)),1);
g3 = repmat({['Last NREM p=' num2str(COMPARISON(2,6))]},length(firingRatesWithRem(:,3)),1);
g = [g1;g2;g3];

figure; 
h = boxplot(datacombinedEpFrs,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-1 1])
title(['CA1 FR eps-p = ' num2str(COMPARISON(3,6)) ' compare 2-3'])
ylim([-0.4 1.6])
yticks([-0.4:0.4:1.6])
ylabel('Firing Rate (Hz)')
set(gcf, 'renderer', 'painters')

meanCo = nanmean(ncRippleCo(:,1));
meanCo = 0;

highCo = ncRippleCo(find(ncRippleCo(:,1)>meanCo),:);
lowCo = ncRippleCo(find(ncRippleCo(:,1)<meanCo),:);

[p1 h1] = ranksum(highCo(:,2),lowCo(:,2))

datacombinedFRchange = [highCo(:,2); lowCo(:,2)];
g1 = repmat({'HighCofiring'},length(highCo(:,2)),1);
g2 = repmat({'LowCofiring'},length(lowCo(:,2)),1);
g = [g1;g2];

figure; 
h = boxplot(datacombinedFRchange,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-1 1])
title(['CA1 FR change-p = ' num2str(p1)])
set(gcf, 'renderer', 'painters')

[pFR hFR] = ranksum(highCo(:,3),lowCo(:,3))

datacombinedFR = [highCo(:,3); lowCo(:,3)];
g1 = repmat({'HighCofiring'},length(highCo(:,3)),1);
g2 = repmat({'LowCofiring'},length(lowCo(:,3)),1);
g = [g1;g2];

[p1 h1] = ranksum(highCo(:,2),lowCo(:,2))

datacombinedFRchange = [highCo(:,2); lowCo(:,2)];
g1 = repmat({'HighCofiring'},length(highCo(:,2)),1);
g2 = repmat({'LowCofiring'},length(lowCo(:,2)),1);
g = [g1;g2];

figure; 
h = boxplot(datacombinedFR,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
ylim([-1 1])
title(['CA1 mean FR-p = ' num2str(pFR)])
set(gcf, 'renderer', 'painters')

figure;
bar([mean(highCo(:,2)) mean(lowCo(:,2))],'k');
hold on
errorbar([1 2], [mean(highCo(:,2)) mean(lowCo(:,2))],...
    [std(highCo(:,2))./sqrt(length(highCo(:,2))) std(lowCo(:,2))./sqrt(length(lowCo(:,2)))],...
    'LineStyle','none')
xticklabels({'High cofiring','Low'})
title(['CA1 FR change-p = ' num2str(p1)])
ylabel('FR change first-last NREM')
set(gcf, 'renderer', 'painters')

first = firingRates(:,1);
last = firingRates(:,2);

[p2 h2] = signrank(firingRates(:,1),firingRates(:,2))

datacombinedFRchange = [first; last];
g1 = repmat({'first'},length(first),1);
g2 = repmat({'last'},length(last),1);
g = [g1;g2];

figure; hold on
h = boxplot(datacombinedFRchange,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
xlim([0.5 2.5])
set(gca, 'YScale', 'log')
% ylim([0.05 2])
title(['CA1 FR -p = ' num2str(p2)])
set(gcf, 'renderer', 'painters')

[p3 h3] = ranksum(highCo(:,4),lowCo(:,4))

datacombinedREMFR = [highCo(:,4); lowCo(:,4)];
g1 = repmat({'HighCofiring'},length(highCo(:,4)),1);
g2 = repmat({'LowCofiring'},length(lowCo(:,4)),1);
g = [g1;g2];

figure; 
h = boxplot(datacombinedREMFR,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
ylim([-0.25 2.5])
title(['CA1 REM FR -p = ' num2str(p3)])
set(gcf, 'renderer', 'painters')

[p4 h4] = ranksum(downFR,upFR)

datacombinedFRchangeCofiring = [downFR; upFR];
g1 = repmat({'FR change lower than mean'},length(downFR),1);
g2 = repmat({'FR change higher than mean'},length(upFR),1);
g = [g1;g2];

figure; 
h = boxplot(datacombinedFRchangeCofiring,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-0.25 2.5])
ylabel('Cofiring with PFC')
title(['CA1 Cofiring -p = ' num2str(p4)])
set(gcf, 'renderer', 'painters')

% figure; hold on
set(gca, 'YScale', 'log')
xlim([0.5 2.5])
for c = 1:length(first)
    x = [1 2];
    tmp = firingRates(c,:);
    if tmp(2) > tmp(1)
    plot(x,tmp,'-r')
    elseif tmp(2) < tmp(1)
        plot(x,tmp,'-k')
    elseif tmp(2) == tmp(1)
        plot(x,tmp,'-m')
    end
end
ylim([0.003 10])

meanFr = mean(firingRates);
semFr = std(firingRates)./sqrt(size(firingRates,1))
figure
errorbar([1:2], meanFr, semFr)
xlim([0.5 2.5])
ylim([0.46 0.6])
title(['CA1 FR -p = ' num2str(p2)])

figure
scatter(ncRippleCo(:,1),ncRippleCo(:,2))
mdlr = fitlm(ncRippleCo(:,1),ncRippleCo(:,2),'RobustOpts','on')
pVal = mdlr.Coefficients.pValue(2);
r = sqrt(mdlr.Rsquared.Ordinary); 
hold on
lsline
% [r p] = corrcoef(ncRippleCo(:,1),ncRippleCo(:,2),'rows','complete')
% [r p3] = corr(ncRippleCo(:,1),ncRippleCo(:,2),'Type','Spearman','rows','complete')
title(['CA1 FR change vs cofiring - p=' num2str(pVal) ' r=' num2str(r) ' Robust LR'])
set(gcf, 'renderer', 'painters')

outlier = find(isoutlier(mdlr.Residuals.Raw))
noOutiers = ncRippleCo;
noOutiers(outlier,:) = [];

figure
%% no outliers
% idx = find(~isnan(noOutiers(:,1)));
% tt = noOutiers(idx,:);
%% outliers included 
idx = find(~isnan(ncRippleCo(:,1)));
tt = ncRippleCo(idx,:);
%%
quartile_sep = floor(length(tt(:,1))/4);
spklatquar = sortrows(tt,1);
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

errorbar(X, data_means(:,2), data_sems(:,2),'b','LineWidth',3);
hold on
xlim([0.5 4.5])
ylabel('Firing Rate change')
xlabel('Cofiring Quartile')
title(['CA1 FR change vs cofiring - p=' num2str(pVal) ' r=' num2str(r) ' Robust LR'])
set(gcf, 'renderer', 'painters')

keyboard
