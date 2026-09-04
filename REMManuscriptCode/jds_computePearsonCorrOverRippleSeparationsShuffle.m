clear; close all;
animalprefixlist= {'ZT2','JS17','JS15','JS14','JS12','JS13','JS34','BG1','JS21','KL8'};

day=1;
compiledData = [];
rvalDataRem = [];
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    savedir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
    spikes = loaddatastruct(savedir, animalprefix, 'spikes', day); % get spikes
    load(sprintf('%s%sctxrippletime_REM0%d.mat',savedir,animalprefix,day));
    load(sprintf('%s%sremeps0%d.mat',savedir,animalprefix,day));
    rem = load(sprintf('%s%srem0%d.mat',savedir,animalprefix,day));
    rem = rem.rem;
    epochs = remeps;
    slp = rem;

    for e = 1:length(epochs)
        ep = epochs(e);
        sleeplist = [slp{day}{ep}.starttime slp{day}{ep}.endtime];
        [ctxidx, hpidx] = jds_getallepcells(savedir, animalprefix, day, ep, []); 
        ctxnum = length(ctxidx(:,1));
        hpnum = length(hpidx(:,1));
        cellidx = ctxidx;

        remrips = [ctxripple{day}{ep}.starttime ctxripple{day}{ep}.endtime];

        if size(remrips,1) > 1
            riptimes = remrips;
            ripnum = size(riptimes,1);
            celldata = [];
            pfc_matrix_corr = [];
            if size(riptimes,1) > 1
                for cellcount = 1:ctxnum 
                    index = [day,ep,ctxidx(cellcount,:)] ;
                    if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
                        spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
                    else
                        spiketimes = [];
                    end
                    spikebins = periodAssign(spiketimes, riptimes(:,[1 2])); 
                    if ~isempty(spiketimes)
                        validspikes = find(spikebins);
                        spiketimes = spiketimes(validspikes); 
                        spikebins = spikebins(validspikes);
                    end
                    spikecount = zeros(1,size(riptimes,1));
                    for i = 1:length(spikebins)
                        spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
                    end
                    spikecount(find(spikecount > 0)) = 1; %convert to 0 and 1
                    pfc_matrix_corr = [pfc_matrix_corr spikecount'];
                end
            end
            compiledData{a}{ep}.rvalMat = pfc_matrix_corr;
            for s = 1:size(sleeplist,1)
                tmpsleep = sleeplist(s,:);
                ripInBout = logical(isExcluded(riptimes(:,1),tmpsleep));
                subsetRips = riptimes(ripInBout,:);
                subsetRips = subsetRips(:,1) + (subsetRips(:,2)-subsetRips(:,1))/2; %middle
                ripDiffs = diff(subsetRips(:,1));
                subset_corr = pfc_matrix_corr(ripInBout,:);
                for n = 1:size(subset_corr,1)-1
                    if (sum(subset_corr(n,:)) > 0) && (sum(subset_corr(n+1,:)) > 0)
                        rval = corrcoef(subset_corr(n,:),subset_corr(n+1,:));
                        rval = rval(1,2);
                        if isnan(rval)
                            continue
                        end
                        if ripDiffs(n) <= 1
                            rvalDataRem = [rvalDataRem; [rval ripDiffs(n)]];
                        end
                    end
                end
            end
        end
    end
end
[rCorr pJI] = corrcoef(rvalDataRem)
rValDataRem = rCorr(1,2);

rValDataRem_s = [];
for shuf = 1:1000
    rvalDataRem_s = [];
    for a = 1:length(animalprefixlist)
        animalprefix = animalprefixlist{a};
        savedir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);

        %%
        %-----create the event matrix during SWRs-----%
        spikes = loaddatastruct(savedir, animalprefix, 'spikes', day); % get spikes
        load(sprintf('%s%sctxrippletime_REM0%d.mat',savedir,animalprefix,day));
        load(sprintf('%s%sremeps0%d.mat',savedir,animalprefix,day));
        rem = load(sprintf('%s%srem0%d.mat',savedir,animalprefix,day));
        rem = rem.rem;
        epochs = remeps;
        slp = rem;

        for e = 1:length(epochs)
            ep = epochs(e);
            sleeplist = [slp{day}{ep}.starttime slp{day}{ep}.endtime];

            remrips = [ctxripple{day}{ep}.starttime ctxripple{day}{ep}.endtime];

            if size(remrips,1) > 1

                riptimes = remrips;

                ripnum = size(riptimes,1);
                pfc_matrix_corr = compiledData{a}{ep}.rvalMat;
                pfc_matrix_corr = pfc_matrix_corr(randperm(size(pfc_matrix_corr,1)),:);
                for s = 1:size(sleeplist,1)
                    tmpsleep = sleeplist(s,:);
                    ripInBout = logical(isExcluded(riptimes(:,1),tmpsleep));
                    subsetRips = riptimes(ripInBout,:);
                    subsetRips = subsetRips(:,1) + (subsetRips(:,2)-subsetRips(:,1))/2; %middle
                    ripDiffs = diff(subsetRips(:,1));
                    subset_corr = pfc_matrix_corr(ripInBout,:);
                    for n = 1:size(subset_corr,1)-1
                        if (sum(subset_corr(n,:)) > 0) && (sum(subset_corr(n+1,:)) > 0)
                            rval = corrcoef(subset_corr(n,:),subset_corr(n+1,:));
                            rval = rval(1,2);
                            if isnan(rval)
                                continue
                            end
                            if ripDiffs(n) <= 1
                                rvalDataRem_s = [rvalDataRem_s; [rval ripDiffs(n)]];
                            end
                        end
                    end
                end
            end
        end
    end
    [rCorr_s pCorr_s] = corrcoef(rvalDataRem_s);
    rValDataRem_s(shuf) = rCorr_s(1,2);
    disp(['Shuffle number ' num2str(shuf) ' rval = ' num2str(rCorr_s(1,2))])
end
figure
histogram(rValDataRem_s,50,'DisplayStyle','stairs')
hold on
x = [rValDataRem rValDataRem];
y = [0 100];
plot(x,y,'-r')
keyboard;


