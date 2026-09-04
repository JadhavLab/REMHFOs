function jds_rippleTriggeredAssemblyStrengthREM(animalprefixlist,area,state)

%To do: plot the peak strength within specified ripple events. compare
%ripple types for different assemblies - might need to change binsize for
%raster generation
day = 1;

% bins = 400; %for 5 ms
% bins = 20; %for 100 ms
bins = 50; %for 20ms
% peakbins = find((-bins:bins)<=20 & (-bins:bins)>=0);
peakbins = find(abs(-bins:bins)<=25);
% peakbins = find(abs(-bins:bins)<=13);
allbins = -bins:bins;


allevents_ctxriptrig = [];
allevents_ca1riptrig = [];
allevents_ctxriptrigNorm = [];
allevents_ca1riptrigNorm = [];
allevents_ctxripstrength = [];
allevents_ca1ripstrength = [];
allevents_ctxriptrigMeanReact = [];
allevents_ca1riptrigMeanReact = [];
shuf = 0;
g1 = gaussian(3, 10);
firstripinchain = 1;

for a = 1:length(animalprefixlist)

    animalprefix = char(animalprefixlist(a));
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    %Load reactivation strength file for all assemblies and epochs
    load(sprintf('%s%s%s_RTimeStrength%sNewSpk_20_%02d.mat',dir,animalprefix,area,state,day));
    %Load ripples
    load(sprintf('%s%sctxrippletime_chainREM%02d.mat',dir,animalprefix,day));
    ripple = ctxripple;
    load(sprintf('%s%sctxrippletime_chainSWS%02d.mat',dir,animalprefix,day));
    load(sprintf('%s%sremeps%02d.mat',dir,animalprefix,day));

    %Which epochs to analyze
    epochs = remeps;
    epochs = epochs(find(epochs>1));

    for e = 1:length(epochs)
        ep = epochs(e);
        assemblytmp = RtimeStrength{ep}.reactivationStrength;
%         ctxripstarts = ctxripple{day}{ep}.starttime;
        if firstripinchain == 1
            tmp = ripple{day}{ep}.C_sep;
            tmp2 = ctxripple{day}{ep}.C_sep;
            ripstarts = [];
            ctxripstarts = [];
            for c = 1:length(tmp)
                ripstarts = [ripstarts; tmp{c}(1,1)];
%                 ripstarts = [ripstarts; tmp{c}(1,1)+(tmp{c}(end,2)-tmp{c}(1,1))/2]; %middle
            end
            for cc = 1:length(tmp2)
                ctxripstarts = [ctxripstarts; tmp2{cc}(1,1)];
            end
        else
            ripstarts = ripple{day}{ep}.starttimeC;
        end

        %Do PFC ripples
        if (length(ctxripstarts) > 10) && (length(ripstarts) > 10)
            if ~isempty(assemblytmp)
                for ii = 1:length(assemblytmp)
                    react_idx = [];
                    for t = 1:length(ctxripstarts)
                        idxtmp = lookup(ctxripstarts(t), assemblytmp{ii}(:,1));
                        react_idx = [react_idx; idxtmp];
                    end
                    atmp = [];
                    atmpMean = [];
                    strengthstmp = assemblytmp{ii}(:,2);
                    if shuf == 1
                        strengthstmp = strengthstmp(randperm(length(strengthstmp)));
                    end
                    for r = 1:length(react_idx)
                        if ((react_idx(r) + bins) < length(strengthstmp)) && ((react_idx(r) - bins) > 1)
                            tmp = strengthstmp((react_idx(r) - bins):(react_idx(r) + bins)); %get vector of reactivation strenths for specified time period
                            atmp = [atmp; tmp'];
                            atmpMean = [atmpMean; max(tmp)];
                        end
                    end
                    allevents_ctxriptrigMeanReact = [allevents_ctxriptrigMeanReact; mean(atmpMean)];
                    allevents_ctxriptrig = [allevents_ctxriptrig; zscore(smoothvect(mean(atmp,1),g1))];
                    allevents_ctxriptrigNorm = [allevents_ctxriptrigNorm; normalize(smoothvect(mean(atmp,1),g1),'range')];
                end
                %Do CA1 ripples
                for ii = 1:length(assemblytmp)
                    react_idx = [];
                    for t = 1:length(ripstarts)
                        idxtmp = lookup(ripstarts(t), assemblytmp{ii}(:,1));
                        react_idx = [react_idx; idxtmp];
                    end
                    atmp = [];
                    atmpMean = [];
                    strengthstmp = assemblytmp{ii}(:,2);
                    if shuf == 1
                        strengthstmp = strengthstmp(randperm(length(strengthstmp)));
                    end
                    for r = 1:length(react_idx)
                        if ((react_idx(r) + bins) < length(strengthstmp)) && ((react_idx(r) - bins) > 1)
                            tmp = strengthstmp((react_idx(r) - bins):(react_idx(r) + bins));
                            atmp = [atmp; tmp'];
                            atmpMean = [atmpMean; max(tmp)];
                        end
                    end
                    allevents_ca1riptrigMeanReact = [allevents_ca1riptrigMeanReact; mean(atmpMean)];
                    allevents_ca1riptrig = [allevents_ca1riptrig; zscore(smoothvect(mean(atmp,1),g1))];
                    allevents_ca1riptrigNorm = [allevents_ca1riptrigNorm; normalize(smoothvect(mean(atmp,1),g1),'range')];
                end
            end
        end
    end
end

[p h] = ranksum(allevents_ctxriptrigMeanReact, allevents_ca1riptrigMeanReact);
datacombinedReactivation = [allevents_ca1riptrigMeanReact; allevents_ctxriptrigMeanReact];
g1 = repmat({'REM'},length(allevents_ca1riptrigMeanReact),1);
g2 = repmat({'NREM'},length(allevents_ctxriptrigMeanReact),1);
g = [g1;g2];

figure;
h = boxplot(datacombinedReactivation,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-0.02 0.2])
title(['Reactivation strength-p = ' num2str(p)])
ylabel('Reactivation strength (AU)')
set(gcf, 'renderer', 'painters')

allevents_ctxriptrig_z_mn = mean(allevents_ctxriptrig);
allevents_ca1riptrig_z_mn = mean(allevents_ca1riptrig);

allevents_ctxriptrig_z_sem = (std(allevents_ctxriptrig)./sqrt(length(allevents_ctxriptrig(:,1))));
allevents_ca1riptrig_z_sem = (std(allevents_ca1riptrig)./sqrt(length(allevents_ca1riptrig(:,1))));

mnPeakCtxRip = mean(allevents_ctxriptrig(:,peakbins),2);
mnPeakCa1Rip = mean(allevents_ca1riptrig(:,peakbins),2);
%
% rDiff = mnPeakCtxRip-mnPeakCa1Rip;
% stem(rDiff); view(90,90)
% g1 = gaussian(3, 10);

%% Example assemblies
%CA1
% subplot(2,1,1)
% hold on
% plot(smoothvect(allevents_ca1riptrig(102,:),g1))
% xlim([50 150])
% xticks([50 100 150])
% subplot(2,1,2)
% hold on
% plot(smoothvect(allevents_ca1riptrig(103,:),g1))
% xlim([50 150])
% xticks([50 100 150])
% subplot(2,1,1)
% plot(smoothvect(allevents_ctxriptrig(102,:),g1))
% subplot(2,1,2)
% plot(smoothvect(allevents_ctxriptrig(103,:),g1))

%PFC
% subplot(2,1,1)
% hold on
% plot(smoothvect(allevents_ca1riptrig(9,:),g1))
% xlim([50 150])
% xticks([50 100 150])
% subplot(2,1,2)
% hold on
% plot(smoothvect(allevents_ca1riptrig(14,:),g1))
% xlim([50 150])
% xticks([50 100 150])
% subplot(2,1,1)
% plot(smoothvect(allevents_ctxriptrig(9,:),g1))
% subplot(2,1,2)
% plot(smoothvect(allevents_ctxriptrig(14,:),g1))
% set(gcf, 'renderer', 'painters')
%%

figure; hold on
ax1 = gca;
ax1.FontSize = 14;
pl1 = plot([-bins:bins],allevents_ca1riptrig_z_mn,'-k','LineWidth',1)
boundedline([-bins:bins],allevents_ca1riptrig_z_mn,allevents_ca1riptrig_z_sem,'-k');
pl2 = plot([-bins:bins],allevents_ctxriptrig_z_mn,'-r','LineWidth',1)
boundedline([-bins:bins],allevents_ctxriptrig_z_mn,allevents_ctxriptrig_z_sem,'-r');

title(sprintf('NC Ripple Triggered %s Reactivation Strength',area))

ylabel('Reactivation Strength')
xlabel('Time from Ripple Onset (s)')
xlim([(-25) 25]);
xticks([-25 0 25])
% xticklabels({'-2','-1.5','-1','-0.5','0','0.5','1','1.5','2'});
xticklabels({'-0.5','0','0.5'});

xx = [0 0];
plot(xx,ax1.YLim,'--k')
ylim = (ax1.YLim);
legend([pl1 pl2],{'PFC Ripples - REM','PFC Ripples - SWS'})
set(gcf, 'renderer', 'painters')

%sort by SWS
[S I] = sort(mnPeakCtxRip,'descend');

figure;
imagesc(allevents_ctxriptrig(I,:))
title('SWS PFC Ripple Triggered Reactivation Strength')
colorbar
set(gcf, 'renderer', 'painters')
clim([-2 4])
xlim([25 75])
xticks([25 50 75])
xticklabels({'-0.5','0','0.5'})
xlabel('Time from Ripple onset')

figure
imagesc(allevents_ca1riptrig(I,:))
title('REM PFC Ripple Triggered Reactivation Strength')
colorbar
set(gcf, 'renderer', 'painters')
clim([-2 4])
xlim([25 75])
xticks([25 50 75])
xticklabels({'-0.5','0','0.5'})
xlabel('Time from Ripple onset')

%sort by rem
[M I] = max(allevents_ca1riptrig(:,peakbins)');
[S I] = sort(I,'ascend');
srtCtx = allevents_ctxriptrig(I,:);
srtCa1 = allevents_ca1riptrig(I,:);
figure
imagesc(allevents_ca1riptrig(I,:))
clim([-2 4])
xlim([25 75])
xticks([25 50 75])
xticklabels({'-0.5','0','0.5'})
figure
imagesc(allevents_ctxriptrig(I,:))
clim([-2 4])
xlim([25 75])
xticks([25 50 75])
xticklabels({'-0.5','0','0.5'})

[M I] = max(srtCtx');
ctxpkloc = allbins(I)*20;
figure; plot(ctxpkloc,flip(1:size(allevents_ctxriptrig,1)),'ko')
xlim([-550 550])

[M I] = max(srtCa1');
ca1pkloc = allbins(I)*20;
figure; plot(ca1pkloc,flip(1:size(allevents_ca1riptrig,1)),'ro')
xlim([-550 550])

idx = find(M > 1);

%sort by rem
[M I] = max(allevents_ctxriptrig(:,peakbins)');
[S I] = sort(I,'ascend');

figure
imagesc(allevents_ctxriptrig(I,:))
clim([-2 4])
xlim([25 75])
xticks([25 50 75])
xticklabels({'-0.5','0','0.5'})
figure
imagesc(allevents_ca1riptrig(I,:))
clim([-2 4])
xlim([25 75])
xticks([25 50 75])
xticklabels({'-0.5','0','0.5'})

%%
keyboard
