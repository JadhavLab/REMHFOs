%jds_chainIsolatedCofiring
%Pairwise cofiring (z-scored coactivity) of PFC (or CA1) cells within chain
%vs isolated cortical ripples during REM. Compares the two distributions
%(rank-sum, dip tests) and plots boxplot and histograms.
%Set 'area' below to 'PFC' or 'CA1'.

clear all
close all
animalprefixlist = {'ZT2','JS17','JS15','JS14','JS12','JS13','JS34','BG1','JS21','KL8'};
day = 1;
area = 'PFC';
mean_isoremCorr = [];
mean_chainremCorr = [];
mean_diff = [];
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    load(sprintf('%s/%sctxrippletime_chainREM%02d.mat', dir, animalprefix, day));
    load(sprintf('%s/%sspikes%02d.mat', dir, animalprefix, day));
    load(sprintf('%s/%sremeps%02d.mat', dir, animalprefix, day));
    epochs = remeps;
    dat = [];
    for e = 1:length(epochs)
        epoch = epochs(e);

        chainrem_riptimes = [ctxripple{day}{epoch}.starttimeC ctxripple{day}{epoch}.endtimeC];
        isorem_riptimes = [ctxripple{day}{epoch}.starttimeNC ctxripple{day}{epoch}.endtimeNC];

        [ctxidx, hpidx] = jds_getallepcells(dir, animalprefix, day, epoch, []); %(tet, cell)
        ctxnum = length(ctxidx(:,1));
        hpnum = length(hpidx(:,1));

        if area == 'PFC'
            cellidx = ctxidx;
            numcells = ctxnum;
        elseif area == 'CA1'
            cellidx = hpidx;
            numcells = hpnum;
        end
        ncCorrMat = [];
        cCorrMat = [];
        celldata = [];
        spikecounts = [];
        if ((length(chainrem_riptimes) > 10) && (length(isorem_riptimes) > 10))
            for cellcount = 1:numcells %get spikes for each cell
                index = [day,epoch,cellidx(cellcount,:)] ;
                if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
                    spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
                else
                    spiketimes = [];
                end
                spikebins = periodAssign(spiketimes, chainrem_riptimes(:,[1 2])); %Assign spikes to align with each ripple event (same number = same rip event, number indicates ripple event)
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
                spikecount = zeros(1,size(chainrem_riptimes,1));
                for i = 1:length(spikebins)
                    spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
                end
                spikecounts = [spikecounts spikecount']; %concatenating num spikes per cell, per event
            end
            for i = 1:numcells
                tmpCo = [];
                for ii = 1:numcells
                    if (i ~= ii) && (ii > i)
                        n1 = spikecounts(:,i);
                        n2 = spikecounts(:,ii);
                        coactiveZ = coactivezscore(n1, n2);
                        tmpCo = [tmpCo; coactiveZ];
                    end
                end
                mean_chainremCorr = [mean_chainremCorr; (nanmean(tmpCo))];
            end
            %isolated REM ripples
            celldata = [];
            spikecounts = [];
            for cellcount = 1:numcells %get spikes for each cell
                index = [day,epoch,cellidx(cellcount,:)] ;
                if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
                    spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
                else
                    spiketimes = [];
                end
                spikebins = periodAssign(spiketimes, isorem_riptimes(:,[1 2])); %Assign spikes to align with each ripple event (same number = same rip event, number indicates ripple event)
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
                spikecount = zeros(1,size(isorem_riptimes,1));
                for i = 1:length(spikebins)
                    spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
                end
                spikecounts = [spikecounts spikecount']; %concatenating num spikes per cell, per event
            end
            for i = 1:numcells
                tmpCo = [];
                for ii = 1:numcells
                    if (i ~= ii) && (ii > i)
                        n1 = spikecounts(:,i);
                        n2 = spikecounts(:,ii);
                        coactiveZ = coactivezscore(n1, n2);
                        tmpCo = [tmpCo; coactiveZ];
                    end
                end
                mean_isoremCorr = [mean_isoremCorr; (nanmean(tmpCo))];
            end
        end
    end
end

[p h] = ranksum(mean_chainremCorr,mean_isoremCorr)
datacombinedCofiring = [mean_chainremCorr; mean_isoremCorr];
g1 = repmat({'Chain'},length(mean_chainremCorr),1);
g2 = repmat({'Isolated'},length(mean_isoremCorr),1);
g = [g1;g2];

figure;
h = boxplot(datacombinedCofiring,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
ylim([-0.4 1.8])
title(['RippleCofiring-p = ' num2str(p)])
ylabel('Abs. cofiring (z)')
set(gcf, 'renderer', 'painters')

[p1,dip,xl,xu]=dipTest(mean_chainremCorr(~isnan(mean_chainremCorr)))
[p2,dip,xl,xu]=dipTest(mean_isoremCorr(~isnan(mean_isoremCorr)))

figure
histogram(mean_chainremCorr,50,'DisplayStyle','stairs');
hold on
histogram(mean_isoremCorr,50,'DisplayStyle','stairs');

