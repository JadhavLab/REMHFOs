function out = jds_remRippleCA1CoactivitySpatialCorr(animalprefixlist)

day = 1;
%%
lowCorrSp = [];
highCorrSp = [];
spCorr = [];
for a = 1:length(animalprefixlist)
    
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
    
    spikes = loaddatastruct(dir, animalprefix, 'spikes', day); % get spikes
    % get ripple time
    load(sprintf('%s%sctxrippletime_chainREM0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%smapfields0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%sremeps0%d.mat',dir,animalprefix,day));
    epochs = remeps;
    
    for e = 1:length(epochs)
        ep = epochs(e);
        if ep == 1
            continue
        end
        rip = ctxripple{day}{ep};
        riptimes = [rip.starttimeC rip.endtimeC];
        if isempty(riptimes)
            continue
        end

        [hpidx, ctxidx] = jds_getallepcells(dir, animalprefix, day, ep, []);
        [hpidx2, ctxidx2] = jds_getallepcells(dir, animalprefix, day, ep-1, []); 

        hpnum = length(hpidx(:,1));
        
        if length(riptimes(:,1)) > 10
            celldata = [];
            CA1matrix = [];
            for cellcount = 1:hpnum %get spikes for each cell
                index = [day,ep,hpidx(cellcount,:)] ;
                if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
                    spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
                else
                    spiketimes = [];
                end
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
                CA1matrix = [CA1matrix spikecount']; %concatenating num spikes per cell, per event
            end
        else
            clear riptimes
            continue
        end
        for n=1:hpnum
            CA1_cell1 = CA1matrix(:,n);
            for nn=1:hpnum
                if n ~= nn
                    CA1_cell2 = CA1matrix(:,nn);
                    [b1, cellLoc] = ismember(hpidx(n,:),hpidx2,'rows','legacy');
                    [b2, cellLoc] = ismember(hpidx(nn,:),hpidx2,'rows','legacy');
                    coactiveZ = coactivezscore(CA1_cell1, CA1_cell2);
                    if (b1) && (b2)
                        if n < nn
                            mf1 = mapfields{day}{ep-1}{hpidx(n,1)}{hpidx(n,2)}.smoothedspikerate;
                            mf1(mf1 < 0) = 0;
                            mf2 = mapfields{day}{ep-1}{hpidx(nn,1)}{hpidx(nn,2)}.smoothedspikerate;
                            mf2(mf2 < 0) = 0;
                            if (~isempty(mf1)) && (~isempty(mf2))
                                R = corr2(mf1,mf2);

                                spCorr = [spCorr; [coactiveZ R]];

                                if coactiveZ > 0
                                    highCorrSp = [highCorrSp; R];
                                elseif coactiveZ < 0
                                    lowCorrSp = [lowCorrSp; R];
                                end
%                                 if (coactiveZ > 0.2) && (R > 0.5) %for plotting high or low pairs
%                                     keyboard
%                                 end
                            end
                        end
                    end
                end
            end
        end
    end
end
[p1 h1] = ranksum(lowCorrSp,highCorrSp)

datacombinedSpCorr = [lowCorrSp; highCorrSp];
g1 = repmat({'Low cofiring'},length(lowCorrSp),1);
g2 = repmat({'High cofiring'},length(highCorrSp),1);
g = [g1;g2];
keyboard