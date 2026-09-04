function jds_getSpindleTrainsNREM(animalprefixlist, day, area)

numNonChainSpins = 0;
numChainSpins = 0;
chainProp = [];
chainRate = [];
totalNumRips = 0;
numsInChain = [];
chainSpinIEI = [];
chainSpinLengths = [];
isoSpinLengths = [];
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    savedata = 0;
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    %     dir = sprintf('/Volumes/JUSTIN/NovelFamiliarNovel/%s_direct/',animalprefix);
    load(sprintf('%s%sremeps0%d.mat',dir,animalprefix,day));

    load(sprintf('%s%sswsALL0%d.mat',dir,animalprefix,day));

    % get ripple time
    load(sprintf('%s%sctxspindletimeG_SWS0%d.mat',dir,animalprefix,day));
    spin = ctxspindle; clear ctxspindle

    eps = remeps;
    for ep = eps

        spin2 = spin{day}{ep};

        spintimes = [spin2.starttime spin2.endtime];
        lengthidx = (spintimes(:,2) - spintimes(:,1) > 0.2);
        spintimes = spintimes(lengthidx,:);
        dur = sws{day}{ep}.total_duration;

        if length(spintimes) > 1
            tdiff = diff(spintimes(:,1));
            trainidx = find(tdiff<=2.78); %ganguly lab paper
            rvec = 1:length(spintimes(:,1))-1;
            spintimenonchain = [];
            spintimechain = [];
            spintimechaincombined = [];
            done = [];
            chainNum = 0;
            for c = 1:length(trainidx)
                if ~isempty(find(done == c))
                    continue
                end
                chainNum = chainNum + 1;
                tmp = [];
                idx = trainidx(c);
                idx2 = trainidx(c) + 1;
                t1 = spintimes(idx,:);
                t2 = spintimes(idx2,:);
                tmp = [tmp; [t1; t2]];
                rvec(idx+1) = idx;
                diffcnt = 1;
                for cc = 1:length(trainidx)
                    if cc > c
                        if trainidx(cc) - trainidx(c) == diffcnt
                            tmp = [tmp; spintimes(trainidx(cc)+1,:)];
                            diffcnt = diffcnt + 1;
                            rvec(cc+1) = cc;
                            done = [done; cc];
                        end
                    end
                end
                chainSpinIEI = [chainSpinIEI; diff(tmp(:,1))];
                chainSpinLengths = [chainSpinLengths; tmp(:,2)-tmp(:,1)];
                spintimechaincombined = [spintimechaincombined; tmp];
                spintimechain{chainNum} = tmp;
                numsInChain = [numsInChain; length(tmp)];
            end
            totalNumRips = totalNumRips + length(spintimechain);
            
            nonchainidx = ~ismember(rvec,trainidx);
            spintimenonchain = spintimes(nonchainidx,:);
%             numsInChain = [numsInChain; ones(size(riptimenonchain,1),1)];
            numNonChainSpins = numNonChainSpins + size(spintimenonchain,1);
            numChainSpins = numChainSpins + size(spintimechaincombined,1);
            propTmp = size(spintimechaincombined,1)/(size(spintimechaincombined,1) + ...
                size(spintimenonchain,1));
            chainProp = [chainProp; propTmp];
            chainRate = [chainRate; size(spintimechaincombined,1)/dur];
            isoSpinLengths = [isoSpinLengths; spintimenonchain(:,2)-spintimenonchain(:,1)];

            if isequal(area,'PFC')
                if ~isempty(spintimechaincombined)
                    ctxspindlenew{day}{ep}.starttimeC = spintimechaincombined(:,1);
                    ctxspindlenew{day}{ep}.endtimeC = spintimechaincombined(:,2);
                else
                    ctxspindlenew{day}{ep}.starttimeC = [];
                    ctxspindlenew{day}{ep}.endtimeC = [];
                end
                ctxspindlenew{day}{ep}.starttimeNC = spintimenonchain(:,1);
                ctxspindlenew{day}{ep}.endtimeNC = spintimenonchain(:,2);
                ctxspindlenew{day}{ep}.C_sep = spintimechain;
            end
        else
            if isequal(area,'PFC')
                ctxspindlenew{day}{ep}.starttimeC = [];
                ctxspindlenew{day}{ep}.endtimeC = [];
                ctxspindlenew{day}{ep}.starttimeNC = [];
                ctxspindlenew{day}{ep}.endtimeNC = [];
                ctxspindlenew{day}{ep}.C_sep = [];
            end

        end
        clear spintimes
    end

    if savedata == 1
        if isequal(area,'PFC')
            ctxspindle = ctxspindlenew;
            clear ctxspindlenew;
            save(sprintf('%s/%sctxspindletime_chainSWS%02d.mat', dir, animalprefix, day), 'ctxspindle');
        end
    end
end
% p = signrank(remrate,swsrate)
% 
% datacombinedRate = [remrate; swsrate];
% g1 = repmat({'REM'},length(remrate),1);
% g2 = repmat({'SWS'},length(swsrate),1);
% g = [g1;g2];
% 
% % boxplot(datacombinedCoact,g);
% figure
% h = boxplot(datacombinedRate,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% % ylim([-0.2 1.2])
% title(['Chain ripple rate-p = ' num2str(p)])
counts = histcounts(numsInChain);
prop = counts./totalNumRips;
bar(prop, 'BarWidth', 1);
xlim([0.5 10.5])
ylabel('Proportion of events')
xlabel('Number of spindles in event')
title('Spindle Chains')
set(gcf, 'renderer', 'painters')

[p1 h1] = ranksum(isoSpinLengths(:,1),chainSpinLengths(:,1))

datacombinedSpindles = [isoSpinLengths(:,1); chainSpinLengths(:,1)];
g1 = repmat({'Isolated'},length(isoSpinLengths(:,1)),1);
g2 = repmat({'Chain'},length(chainSpinLengths(:,1)),1);
g = [g1;g2];

% boxplot(datacombinedBurst,g);
figure; 
h = boxplot(datacombinedSpindles,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
ylim([0.05 1.5])
title(['Spindle Durations-p = ' num2str(p1)])
set(gcf, 'renderer', 'painters')
keyboard

