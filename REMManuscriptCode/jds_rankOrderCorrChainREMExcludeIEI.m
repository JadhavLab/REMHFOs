function jds_rankOrderCorrChainREMExcludeIEI(animalprefixlist,day,cellcountthresh)


remRho = [];
nremRho = [];
remRhoEp = [];
nremRhoEp = [];
plotranks = 0;
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    %%
    %-----match neurons across epochs-----%

    %%
    %-----create the event matrix during SWRs-----%
    load(sprintf('%s/%sctxgammatime_chainREM%02d.mat', dir, animalprefix, day));
    rem_rip = ctxgamma; clear ctxripple
    load(sprintf('%s/%sctxrippletime_chainSWS%02d.mat', dir, animalprefix, day));
    load(sprintf('%s/%sspikes%02d.mat', dir, animalprefix, day));
    nrem_rip = ctxripple; clear ctxripple
    load(sprintf('%s/%sremeps%02d.mat', dir, animalprefix, day));
    eps = remeps;
    for ep = eps
        [ctxidx, hpidx] = jds_getallepcells_includeall(dir, animalprefix, day, ep, []);
        hpnum = length(hpidx(:,1));
        ctxnum = length(ctxidx(:,1));

        nonriptimes = [];
        rem_riptimes = [];
        
        for r = 1:length(rem_rip{day}{ep}.C_sep)
            riptmp = rem_rip{day}{ep}.C_sep{r};
            rem_riptimes = [rem_riptimes; [riptmp(1,1) riptmp(end,2) 1]];
            for rr = 1:size(riptmp,1)-1
                tmpEnd = riptmp(rr,2);
                tmpSt = riptmp(rr+1,1);
                nonriptimes = [nonriptimes; [tmpEnd tmpSt]]; %compiling between rip times
            end
        end

        nrem_riptimes = [];

        for r = 1:length(nrem_rip{day}{ep}.C_sep)
            riptmp = nrem_rip{day}{ep}.C_sep{r};
            nrem_riptimes = [nrem_riptimes; [riptmp(1,1) riptmp(end,2) 2]];
            for rr = 1:size(riptmp,1)-1
                tmpEnd = riptmp(rr,2);
                tmpSt = riptmp(rr+1,1);
                nonriptimes = [nonriptimes; [tmpEnd tmpSt]]; %compiling between rip times
            end
        end

        if (size(rem_riptimes,1) > 10) && (size(nrem_riptimes,1) > 10)

            %combine riptimes
            riptimes = sortrows([rem_riptimes; nrem_riptimes],1);

            ripnum = size(riptimes,1);
            tmp_nrem = [];
            tmp_rem = [];
            %%
            if ripnum > 1
                celldata = [];
                spikecounts = [];
                for cellcount = 1:ctxnum %get spikes for each cell
                    index = [day,ep,ctxidx(cellcount,:)] ;
                    if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
                        spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
                    else
                        spiketimes = [];
                    end
                    spiketimes = spiketimes(logical(~isExcluded(spiketimes,nonriptimes)));
                    spikebins = periodAssign(spiketimes, riptimes(:,[1 2])); %Assign spikes to align with each ripple event (same number = same rip event, number indicates ripple event)
                    if ~isempty(spiketimes)
                        validspikes = find(spikebins);
                        spiketimes = spiketimes(validspikes); %get spike times that happen during ripples
                        spikebins = spikebins(validspikes);
                        tmpcelldata = [spiketimes spikebins];
                    end
                    if ~isempty(spiketimes)
                        tmpcelldata(:,3) = cellcount; %keep count of how many cells active during rip event
                    else
                        tmpcelldata = [0 0 cellcount];
                    end
                    celldata = [celldata; tmpcelldata];
                    spikecount = zeros(1,size(riptimes,1));
                    for i = 1:length(spikebins)
                        spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
                    end
                    spikecounts = [spikecounts; spikecount]; %concatenating num spikes per cell, per event
                end
                cellcounts = sum((spikecounts > 0));
                eventindex = find(cellcounts >= cellcountthresh); %Find event indices that have more than cellthresh # of cells
                if ~isempty(eventindex)
                    for event = 1:length(eventindex)
                        rippletype = riptimes(eventindex(event),3);
                        cellsi = celldata(find(celldata(:,2)==eventindex(event)),3);
                        spktimes = celldata(find(celldata(:,2)==eventindex(event)),1);

                        [cellsi,ia] = unique(cellsi,'first');
                        spktimes = spktimes(ia);
                        [~, I] = sort(spktimes);
                        %                     [~,sortorder] = sort(ia);
                        I(:,2) = 1:length(I(:,1));
                        rankorder = [cellsi(I(:,1)) normalize(I(:,2),'range')];
                        templateranks = [];
                        for evs = 1:length(eventindex)
                            nanVec = nan(ctxnum,1);
                            rippletype2 = riptimes(eventindex(evs),3);

                            if (evs ~= event) && (rippletype == rippletype2)
                                cellsi2 = celldata(find(celldata(:,2)==eventindex(evs)),3);
                                spktimes2 = celldata(find(celldata(:,2)==eventindex(evs)),1);
                                [cellsi2,ia2] = unique(cellsi2,'first');
                                spktimes2 = spktimes2(ia2);
                                [~, I2] = sort(spktimes2);
                                %                             [~,sortorder2] = sort(ia2);
                                I2(:,2) = 1:length(I2(:,1));
                                rankorder2 = [cellsi2(I2(:,1)) normalize(I2(:,2),'range')];
                                nanVec(rankorder2(:,1)) = rankorder2(:,2);
                                tmpRank = nanVec;
                                templateranks = [templateranks tmpRank];
                            end
                        end
                        templatemeanranks = nanmean(templateranks,2);
                        ranksems = (nanstd(templateranks')./sqrt(length(templateranks(1,:))))';
                        currcells = templatemeanranks(rankorder(:,1));

                        if plotranks == 1
                            [~, sortidx] = sort(templatemeanranks);
                            sortmean = templatemeanranks(sortidx);
                            sortsem = ranksems(sortidx);
                            errorbar([1:length(sortmean)],...
                                sortmean',sortsem','r','LineStyle','none','Marker','o','MarkerFaceColor','auto');
                            title('rank order template')
                            ylabel('Normalized rank')
                            xlabel('Cell #')
                            set(gcf, 'renderer', 'painters')
                            keyboard;
                        end

%                         rho = corr(currcells,rankorder(:,2),'Type','Spearman','rows','complete');
                        rho = corr(currcells,rankorder(:,2),'rows','complete');
                        
                        if rippletype == 1
                            remRho = [remRho; rho];
                            tmp_rem = [tmp_rem; rho];
                        elseif rippletype == 2
                            nremRho = [nremRho; rho];
                            tmp_nrem = [tmp_nrem; rho];
                        end
                    end
                end
            end
            remRhoEp = [remRhoEp; mean(tmp_rem)];
            nremRhoEp = [nremRhoEp; mean(tmp_nrem)];
        end
    end
end
keyboard
[p h] = ranksum(remRho,nremRho);
[p2 h2] = signrank(remRhoEp,nremRhoEp);
figure
bar([mean(remRho) mean(nremRho)],'k')
hold on
errorbar([1:2],[mean(remRho) mean(nremRho)],...
    [(std(remRho)./sqrt(length(remRho))) ...
    (std(nremRho)./sqrt(length(nremRho)))],'k','LineStyle','none')
% ylim([0.28 0.315])
% yticks([0.28:0.01:0.31])
xticklabels({'REM','NREM'})
ylabel('Rank order correlation (r)')
title(['Rank order correlation p=' num2str(p)])
set(gcf, 'renderer', 'painters')

keyboard