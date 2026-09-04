clear all
close all
animalprefixlist = {'ZT2','JS17','JS15','JS14','JS12','JS13','JS34','BG1','JS21','KL8'};
day = 1;
area = 'PFC';
mean_nremCorr = [];
mean_remCorr = [];
mean_diff = [];
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    load(sprintf('%s/%sctxrippletime_REM%02d.mat', dir, animalprefix, day));
    rem_rip = ctxripple; clear ctxripple
    load(sprintf('%s/%sctxrippletime_SWS%02d.mat', dir, animalprefix, day));
    load(sprintf('%s/%sspikes%02d.mat', dir, animalprefix, day));
    load(sprintf('%s/%sremeps%02d.mat', dir, animalprefix, day));
    nrem_rip = ctxripple; clear ctxripple
    epochs = remeps;
    dat = [];
    for e = 1:length(epochs)
        epoch = epochs(e);
% 
        rem_riptimes = [rem_rip{day}{epoch}.starttime rem_rip{day}{epoch}.endtime];
        nrem_riptimes = [nrem_rip{day}{epoch}.starttime nrem_rip{day}{epoch}.endtime];

%         rem_riptimes = [rem_rip{day}{epoch}.starttimeC rem_rip{day}{epoch}.endtimeC];
%         nrem_riptimes = [rem_rip{day}{epoch}.starttimeNC rem_rip{day}{epoch}.endtimeNC];

        %might change to do all epochs combined (need tracked cells)
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
        if ((length(rem_riptimes) > 20) && (length(nrem_riptimes) > 20))
            for cellcount = 1:numcells %get spikes for each cell
                index = [day,epoch,cellidx(cellcount,:)] ;
                if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
                    spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
                else
                    spiketimes = [];
                end
                spikebins = periodAssign(spiketimes, rem_riptimes(:,[1 2])); %Assign spikes to align with each ripple event (same number = same rip event, number indicates ripple event)
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
                spikecount = zeros(1,size(rem_riptimes,1));
                for i = 1:length(spikebins)
                    spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
                end
                spikecounts = [spikecounts spikecount']; %concatenating num spikes per cell, per event
            end
            for i = 1:numcells
                tmpcoact = [];
                for ii = 1:numcells
                    if (i ~= ii) && (ii > i)
                        n1 = spikecounts(:,i);
                        n2 = spikecounts(:,ii);
                        coactiveZ = coactivezscore(n1, n2);
                        tmpcoact = [tmpcoact; coactiveZ];
                    end
                end
                mean_remCorr = [mean_remCorr; mean(tmpcoact)];
            end
            %do coord rips
            celldata = [];
            spikecounts = [];
            for cellcount = 1:numcells %get spikes for each cell
                index = [day,epoch,cellidx(cellcount,:)] ;
                if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
                    spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
                else
                    spiketimes = [];
                end
                spikebins = periodAssign(spiketimes, nrem_riptimes(:,[1 2])); %Assign spikes to align with each ripple event (same number = same rip event, number indicates ripple event)
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
                spikecount = zeros(1,size(nrem_riptimes,1));
                for i = 1:length(spikebins)
                    spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
                end
                spikecounts = [spikecounts spikecount']; %concatenating num spikes per cell, per event
            end
            for i = 1:numcells
                tmpcoact = [];
                for ii = 1:numcells
                    if (i ~= ii) && (ii > i)
                        n1 = spikecounts(:,i);
                        n2 = spikecounts(:,ii);
                        coactiveZ = coactivezscore(n1, n2);
                        tmpcoact = [tmpcoact; coactiveZ];
                    end
                end
                mean_nremCorr = [mean_nremCorr; mean(tmpcoact)];
            end
        end
    end
end

[p h] = ranksum(mean_remCorr,mean_nremCorr)
datacombinedCofiring = [mean_remCorr; mean_nremCorr];
g1 = repmat({'REM'},length(mean_remCorr),1);
g2 = repmat({'NREM'},length(mean_nremCorr),1);
g = [g1;g2];

figure;
h = boxplot(datacombinedCofiring,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-0.02 0.2])
title(['RippleCofiring-p = ' num2str(p)])
ylabel('Cofiring (z)')
set(gcf, 'renderer', 'painters')

[p1,dip,xl,xu]=dipTest(mean_remCorr(~isnan(mean_remCorr)))
[p2,dip,xl,xu]=dipTest(mean_nremCorr(~isnan(mean_nremCorr)))

figure
histogram(mean_remCorr,50,'DisplayStyle','stairs');
hold on
histogram(mean_nremCorr,50,'DisplayStyle','stairs');
keyboard

