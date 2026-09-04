function jds_detectPopulationSuppressionREM(animalprefixlist,area)

% --------------- Parameters ---------------
pret=200; postt=200;

%% --------------------------------------------------
%  %Align Single Cell Firing Rate to event
% -------------------------------------------------
day = 1;
savedata = 1;
binsize = 0.005;
g1 = gaussian(3,3);
g2 = gaussian(3,24);
thresh = [0];
lengthRange = [0.3 0.5];
allEvsThresh = [];
allAnimEvs = [];
allWavs = [];
allMinWavs = [];
allHist = [];
shuf = 1;
shufnum = 1000;
for a = 1:length(animalprefixlist)
    suppressiontimes = [];
    anim_data_tmp = [];
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

%     load(sprintf('%s%sctxrippletime_REM0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%sctxrippletime_SWS0%d.mat',dir,animalprefix,day));

    load(sprintf('%s%s_spikematrix_ev_allepochallcell5_%02d.mat',dir,animalprefix,day));
%     rem = load(sprintf('%s%srem0%d.mat',dir,animalprefix,day));
        load(sprintf('%s%sswsALL0%d.mat',dir,animalprefix,day));
%     rem = rem.rem;
        rem = sws;
    load(sprintf('%s%sremeps0%d.mat',dir,animalprefix,day));
    epochs = remeps;
    for e = 1:length(epochs)
        ep_tmp = [];
        ep = epochs(e);

%         currRips = ctxripple{day}{ep}.starttime;
                currRips = ctxripple{day}{ep}.starttime;

        if strcmp(area,'CA1')
            datamat = observation_matrix{ep}.hpdata;
        else
            datamat = observation_matrix{ep}.ctxdata;
        end

        timevect = observation_matrix{ep}.timeeeg(1)*1000:5:observation_matrix{ep}.timeeeg(end)*1000;
        if ~isempty(rem{day}{ep}.starttime)
            remlist = [rem{day}{ep}.starttime rem{day}{ep}.endtime];
            remlength = rem{day}{ep}.total_duration; %in seconds
            [~,remvec] = wb_list2vec(remlist,timevect/1000);
            remidx = find(remvec == 1);

            meanpopFR = sum(sum(datamat(:,remidx)))/remlength; %mean population rate during rem

            activitysummed = sum(datamat);
            activitysummed = smoothvect(activitysummed, g1);
            binnedFR = activitysummed./0.005;
            normFR = binnedFR./meanpopFR; %normalized population FR normalized by REM avg
            for t = 1:length(thresh)
                minVec = movmin(normFR,5); %100ms moving minimum
                %                 smMin = conv(minVec,g2,'same');
                smMin = conv(normFR,g2,'same');
                zSmMin = zscore(smMin);
                %                 zSmMin = zscore(normFR);
                allEvs = [];
                allMinEvs = [];
                supprTroughTimes = [];
                tmpThresh = thresh(t);
                %                 tmpThresh = min(zSmMin)/2;
                belowThresh = zSmMin < tmpThresh;
                supprList = vec2list(belowThresh,timevect/1000);
                for i = 1:size(supprList,1)
                    tmp = supprList(i,:);
                    if sum(isExcluded(tmp',remlist)) == 2 %if both start and end are within rem sleep
                        if ((tmp(2)-tmp(1)) > lengthRange(1)) && ...
                                ((tmp(2)-tmp(1)) < lengthRange(2)) %threshold for event length
                            stidx = lookup(tmp(1),timevect/1000);
                            endidx = lookup(tmp(2),timevect/1000);
                            tmpVec = zSmMin(stidx:endidx);
                            [M minIdx] = min(tmpVec);
                            minTime = tmp(1) + (binsize*minIdx); %add time to starttime of event to get point of largest suppression
                            midIdx = stidx+minIdx; %get the min idx
                            if (midIdx+postt) < length(normFR)
                                supprWav = conv(normFR(midIdx-pret:midIdx+postt),g2,'same'); %1 second around suppression event
                                supprMin = zSmMin(midIdx-pret:midIdx+postt);
                                allEvs = [allEvs; supprWav]; %compile for plotting
                                allMinEvs = [allMinEvs; supprMin];
                                supprTroughTimes = [supprTroughTimes; minTime];
                            end
                        end
                    end
                end
                for r = 1:size(supprTroughTimes,1)
                    currSup = supprTroughTimes(r,1);
                    currRipsTmp =  currRips(find( (currRips>=(currSup-5))...
                        & (currRips<=(currSup+5)) ));
                    currRipsTmp = currRipsTmp-(currSup);
                    histRips = histc(currRipsTmp,[-5:0.1:5]);
                    allHist = [allHist; histRips(:).'];
                end
                anim_data_tmp = [anim_data_tmp; allEvs];
                allWavs = [allWavs; allEvs];
                allMinWavs = [allMinWavs; allMinEvs];
                allEvsThresh{a}{t} = allEvs;
                suppressiontimes{day}{ep}.troughtime{t} = supprTroughTimes;
                suppressiontimes{day}{ep}.muawavs{t} = allEvs;
                %                 suppressiontimes{day}{ep}.minwavs{t} = allMinEvs;
                suppressiontimes{day}{ep}.descrip = 'population suppression events (z)';
                suppressiontimes{day}{ep}.threshold = thresh; 
                suppressiontimes{day}{ep}.lengthrange = lengthRange;
                if shuf == 1
                    for s = 1:shufnum
                        tmpMat = [];
                        for c = 1:size(datamat,1)
                            randcirc = randi(size(datamat,2),1);
                            tmpCell = circshift(datamat(c,:),randcirc,2);
                            tmpMat = [tmpMat; tmpCell];
                        end
                        datamat_s = tmpMat;
                        activitysummed_s = sum(datamat_s);
                        activitysummed_s = smoothvect(activitysummed_s, g1);
                        binnedFR_s = activitysummed_s./0.005;
                        normFR_s = binnedFR_s./meanpopFR; %normalized population FR normalized by REM avg
                        smMin_s = conv(normFR_s,g2,'same');
                        zSmMin_s = zscore(smMin_s);
                        allEvs_s = [];
                        allMinEvs_s = [];
                        tmpThresh = thresh(t);
                        belowThresh_s = zSmMin_s < tmpThresh;
                        supprList_s = vec2list(belowThresh_s,timevect/1000);
                        supprTroughTimes_s = [];
                        for i = 1:size(supprList_s,1)
                            tmp = supprList_s(i,:);
                            if sum(isExcluded(tmp',remlist)) == 2 %if both start and end are within rem sleep
                                if ((tmp(2)-tmp(1)) > lengthRange(1)) && ...
                                        ((tmp(2)-tmp(1)) < lengthRange(2)) %threshold for event length
                                    stidx = lookup(tmp(1),timevect/1000);
                                    endidx = lookup(tmp(2),timevect/1000);
                                    tmpVec = zSmMin_s(stidx:endidx);
                                    [M minIdx] = min(tmpVec);
                                    minTime = tmp(1) + (binsize*minIdx); %add time to starttime of event to get point of largest suppression
                                    midIdx = stidx+minIdx; %get the min idx
                                    if (midIdx+postt) < length(normFR_s)
                                        supprWav = conv(normFR_s(midIdx-pret:midIdx+postt),g2,'same'); %1 second around suppression event
                                        supprMin = zSmMin_s(midIdx-pret:midIdx+postt);
                                        allEvs_s = [allEvs_s; supprWav]; %compile for plotting
                                        allMinEvs_s = [allMinEvs_s; supprMin];
                                        supprTroughTimes_s = [supprTroughTimes_s; minTime];
                                    end
                                end
                            end
                        end
                        suppressiontimes{day}{ep}.shuf{t}{s} = supprTroughTimes_s;
                        suppressiontimes{day}{ep}.muawavs_s{t}{s} = allEvs_s;
                    end
                end
            end
        end
    end
    allAnimEvs{a} = anim_data_tmp;
    if savedata == 1
        save(sprintf('%s%sPFCNREMmuasuppression%02d.mat', dir,animalprefix,day), 'suppressiontimes', '-v7.3');
    end
end
keyboard
figure; hold on %plot bounded line
ax1 = gca;
ax1.FontSize = 14;
pl1 = plot(-pret*binsize:binsize:postt*binsize,nanmean(allWavs),'-k','LineWidth',1)
boundedline(-pret*binsize:binsize:postt*binsize,nanmean(allWavs),...
    nanstd(allWavs)./sqrt(length(allWavs(:,1))),'-k');
xlim([-0.5 0.5])
set(gcf, 'renderer', 'painters')
ylabel('MUA suppression events')
xlabel('Time from event trough (ms)')
%
% figure; hold on %plot bounded line
% ax1 = gca;
% ax1.FontSize = 14;
% pl1 = plot(-pret*binsize:binsize:postt*binsize,nanmean(allMinWavs),'-r','LineWidth',1)
% boundedline(-pret*binsize:binsize:postt*binsize,nanmean(allMinWavs),...
%     nanstd(allMinWavs)./sqrt(length(allMinWavs(:,1))),'-r');
% xlim([-0.5 0.5])
% set(gcf, 'renderer', 'painters')
% ylabel('Moving minimum')
% xlabel('Time from event minimum (ms)')

figure; hold on
imagesc(-200:200, 1:size(allWavs,1), allWavs(randperm(length(allWavs)),:))
ylabel('Event #')
title('Suppression events')
set(gcf, 'renderer', 'painters')
xlim([-100 100])
ylim([1 size(allWavs,1)])
yticks(size(allWavs,1))
xticks([-100:100:100])
xticklabels({'-500','0','500'})
colormap(magma)

figure; hold on
bar(-50:50, sum(allHist)./size(allHist,1))
y = [0 0.15];
x = [0 0];
plot(x,y,'--k')
% plot(-50:50, conv(sum(allHist)./size(allHist,1),g1,'same'))
xlim([-30 30])
xticks([-30:10:30])
xticklabels({'-3','-2','-1','0','1','2','3'})
xlabel('Time from suppression trough (s)')
ylabel('PFC ripple probability')

keyboard