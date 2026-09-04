%Get firing rate bias and percentage of ripple during which the cell is
%active (control for the cell being active in a different number of PFC vs
%CA1 ripples)
clear all;
close all;
%%
animalprefixlist = {'ZT2','JS17','JS15','JS14','JS12','JS13','JS34','BG1','JS21','KL8'};
day = 1;

%%
blen = 6/1000;

shift_SWRburstingProb = []; %frac of rips with bursting
non_SWRburstingProb = [];
shift_SWRburstingProp = []; %frac of spks in burst over all rips
non_SWRburstingProp = [];
shift_SWRspkLat = [];
non_SWRspkLat = [];
numShift = 0;
numNonshift = 0;

for a = 1:length(animalprefixlist)
    animalprefix = char(animalprefixlist{a});
    animdir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    
    load(sprintf('%s%sctxrippletime_chainREM0%d.mat',animdir,animalprefix,day));% get ripple time
    load(sprintf('%s%sspikes0%d.mat',animdir,animalprefix,day));
    rem = load(sprintf('%s%srem0%d.mat',animdir,animalprefix,day));
    rem = rem.rem;
    load(sprintf('%s%sremeps0%d.mat',animdir,animalprefix,day));
    load(sprintf('%s%sCA1remshiftershighthresh%02d.mat',animdir,animalprefix,day));
    epochs = remeps;

    for e = 1:length(epochs)
        epoch = epochs(e);
                
        if epoch <10
            epochstring = ['0',num2str(epoch)];
        else
            epochstring = num2str(epoch);
        end
        
        riptimes = [ctxripple{day}{epoch}.starttimeNC ctxripple{day}{epoch}.endtimeNC];
        if (isempty(riptimes)) || (epoch == 1)
            continue
        end
        remlist = [rem{day}{epoch}.starttime rem{day}{epoch}.endtime];
        remdur = rem{day}{epoch}.total_duration;
        curreegfile = [animdir,'/EEG/',animalprefix,'eeg', '01','-',epochstring,'-','02']; %use any tetrode
        load(curreegfile);
        time1 = geteegtimes(eeg{day}{epoch}{2}) ; % construct time array
        
        [~,swsvec] = wb_list2vec(remlist,time1);
        
        [~,ripvec] = wb_list2vec(riptimes,time1);
        allmodcells = shiftList{day}{epoch}.cellidx;
        numShift = numShift + sum(allmodcells(:,3) == 1);
        numNonshift = numNonshift + sum(allmodcells(:,3) == 0);
        %%
        for cellcount = 1:length(allmodcells(:,1))
            cellshift = allmodcells(cellcount,3);
            if ~isempty(riptimes)
                index = [day,epoch,allmodcells(cellcount,[1:2])];
                
                if (length(riptimes(:,1)) > 1) && (length(remlist(:,1)) > 1)
                    if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)})
                        spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
                        fRate = spikes{index(1)}{index(2)}{index(3)}{index(4)}.meanrate;
                        
                        ripSpkCnt = 0;
                        burstSpks = [];
                        spkLat = [];
                        inRipSpks = [];
                        for r = 1:length(riptimes(:,1))
                            spkRipTmp = find((spiketimes > riptimes(r,1)) & (spiketimes < riptimes(r,2)));
                            spkTimeTmp = spiketimes(spkRipTmp);
                            ripSpkCnt = ripSpkCnt + length(spkTimeTmp);
                            if (~isempty(spkTimeTmp)) %Latency to first spk
                                spkLatTmp = spkTimeTmp(1) - riptimes(r,1);
                                if isnan(spkLatTmp)
                                    keyboard
                                end
                                spkLat = [spkLat; spkLatTmp];
                            end
                            if (~isempty(spkTimeTmp)) && (length(spkTimeTmp) > 2)
                                inRipSpks = [inRipSpks; spkTimeTmp];
                                tmpisi = diff(spkTimeTmp);
                                bspikes = find(tmpisi < blen);
                                if ~isempty(bspikes)
                                    burstSpks = [burstSpks; 1];
                                else
                                    burstSpks = [burstSpks; 0];
                                end
                            end
                        end
                        tmpisiAll = diff(inRipSpks);
                        
                        % find the intervals less than the burst length
                        bspikesAll = find(tmpisiAll < blen);
                        
                        nadjacent = 0;
                        if (length(bspikesAll) > 2)
                            nadjacent = length(find(diff(bspikesAll) == 1));
                        end
                        propb = (length(bspikesAll) * 2 - nadjacent) / length(inRipSpks);
                        if length(inRipSpks) > 1
                            if cellshift == 1
                                shift_SWRburstingProb = [shift_SWRburstingProb; [sum(burstSpks)/length(riptimes(:,1)) allmodcells(cellcount,6)]];
                                shift_SWRspkLat = [shift_SWRspkLat; [mean(spkLat) allmodcells(cellcount,6)]];
                                shift_SWRburstingProp = [shift_SWRburstingProp; [propb allmodcells(cellcount,6)]];
                            elseif cellshift == 0
                                non_SWRburstingProb = [non_SWRburstingProb; [sum(burstSpks)/length(riptimes(:,1)) allmodcells(cellcount,6)]];
                                non_SWRspkLat = [non_SWRspkLat; [mean(spkLat) allmodcells(cellcount,6)]];
                                non_SWRburstingProp = [non_SWRburstingProp; [propb allmodcells(cellcount,6)]];
                            end
                        end
                    end
                end
            end
        end
    end
end
datameans = [mean(shift_SWRburstingProb(:,1)) mean(non_SWRburstingProb(:,1))]
datasems = [(std(shift_SWRburstingProb(:,1))/sqrt(length(shift_SWRburstingProb(:,1))))...
    (std(non_SWRburstingProb(:,1))/sqrt(length(non_SWRburstingProb(:,1))))]
datamedians = [median(shift_SWRburstingProb(:,1)) median(non_SWRburstingProb(:,1))]

bar([1:2],datameans,'k')
hold on
er = errorbar([1:2],datameans,datasems);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('Bursting Probability')
title('Bursting Probability during NC SWRs')
xticklabels({'Shift','Nonshift'}); xtickangle(45)

[p h] = ranksum(non_SWRburstingProb(:,1),shift_SWRburstingProb(:,1))

datacombinedBurst = [non_SWRburstingProb(:,1); shift_SWRburstingProb(:,1)];
g1 = repmat({'Nonshift'},length(non_SWRburstingProb(:,1)),1);
g2 = repmat({'Shift'},length(shift_SWRburstingProb(:,1)),1);
g = [g1;g2];

% boxplot(datacombinedBurst,g);
figure; 
h = boxplot(datacombinedBurst,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
ylim([-0.02 0.2])
title('CA1 SWR Bursting')
title(['CA1 SWR Bursting-p = ' num2str(p)])
%%
datameansLat = [mean(shift_SWRspkLat(:,1)) mean(non_SWRspkLat(:,1))]
datasemsLat = [(std(shift_SWRspkLat(:,1))/sqrt(length(shift_SWRspkLat(:,1))))...
    (std(non_SWRspkLat(:,1))/sqrt(length(non_SWRspkLat(:,1))))]
datamediansLat = [median(shift_SWRspkLat(:,1)) median(non_SWRspkLat(:,1))]

figure
bar([1:2],datameansLat,'k')
hold on
er = errorbar([1:2],datameansLat,datasemsLat);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('Latency to first spk')
title('Latency to first spk during NC SWRs')
xticklabels({'Shift','Nonshift'}); xtickangle(45)

[p1 h1] = ranksum(shift_SWRspkLat(:,1),non_SWRspkLat(:,1))

datacombinedLatency = [shift_SWRspkLat(:,1); non_SWRspkLat(:,1)];
g1 = repmat({'Shift'},length(shift_SWRspkLat(:,1)),1);
g2 = repmat({'Nonshift'},length(non_SWRspkLat(:,1)),1);
g = [g1;g2];

% boxplot(datacombinedBurst,g);
figure; 
h = boxplot(datacombinedLatency,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
ylim([-0.02 0.2])
title('CA1 SWR Spike Latencuy')
title(['CA1 SWR Spike Latency-p = ' num2str(p1)])

[rburstIn pburstIn] = corrcoef(non_SWRburstingProb)

quartile_sep = floor(length(non_SWRburstingProb(:,1))/4);
spkburstquar = sortrows(non_SWRburstingProb,1);
vals = [];
cnt = 1;
for s = 1:4
    if s < 4
        tmp = spkburstquar(cnt:quartile_sep*s,1:2);
        tmp(:,3) = s;
    else
        tmp = spkburstquar(cnt:end,1:2);
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

figure; hold on; errorbar(X, data_means(:,2), data_sems(:,2),'k','LineWidth',3);
xlim([0.5 4.5])
ylabel('Shift (deg.)')
xlabel('Burst Probability Quartile')

[rburstEx pburstEx] = corrcoef(shift_SWRburstingProb)

quartile_sep = floor(length(shift_SWRburstingProb(:,1))/4);
spkburstquar = sortrows(shift_SWRburstingProb,1);
vals = [];
cnt = 1;
for s = 1:4
    if s < 4
        tmp = spkburstquar(cnt:quartile_sep*s,1:2);
        tmp(:,3) = s;
    else
        tmp = spkburstquar(cnt:end,1:2);
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

errorbar(X, data_means(:,2), data_sems(:,2),'r','LineWidth',3);
xlim([0.5 4.5])
ylabel('Shift (deg.)')
xlabel('Burst Probability Quartile')

legend({['INH-p = ' num2str(pburstIn(1,2)) 'r=' num2str(rburstIn(1,2))] ['EXC-p = ' num2str(pburstEx(1,2)) 'r=' num2str(rburstEx(1,2))]})
title('CA1 SWR Bursting - Quartile')

%%
%Latency Quartile

[rlatIn platIn] = corrcoef(non_SWRspkLat)

quartile_sep = floor(length(non_SWRspkLat(:,1))/4);
spklatquar = sortrows(non_SWRspkLat,1);
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

figure; hold on; errorbar(X, data_means(:,2), data_sems(:,2),'k','LineWidth',3);
xlim([0.5 4.5])
ylabel('Modulation Index')
xlabel('Spike Latency Quartile')

[rlatEx platEx] = corrcoef(shift_SWRspkLat)

quartile_sep = floor(length(shift_SWRspkLat(:,1))/4);
spklatquar = sortrows(shift_SWRspkLat,1);
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

errorbar(X, data_means(:,2), data_sems(:,2),'r','LineWidth',3);
xlim([0.5 4.5])
ylabel('Modulation Index')
xlabel('Spike Latency Quartile')

legend({['INH-p = ' num2str(platIn(1,2)) 'r=' num2str(rlatIn(1,2))] ['EXC-p = ' num2str(platEx(1,2)) 'r=' num2str(rlatEx(1,2))]})
title('CA1 SWR Spike Latency - Quartile')
ylim([-0.7 0.9])

keyboard;