%jds_assemblyFieldMetricsShuffle
%Computes 2D spatial rate maps for PFC assembly activation events (run
%epochs, reactivation strength > 5) and quantifies their spatial
%information and sparsity against a null distribution from circularly
%shifted activation times (n=1000 shuffles, minimum 20 s shift). Plots
%z-scored spatial information and sparsity histograms (1.65 = one-tailed
%95th percentile) and data vs shuffle-median comparisons (sign-rank).

animalprefixlist = {'JS17','JS15','JS14','JS12','JS13','JS34','BG1','JS21','KL8','ZT2'};
day = 1;
epochs = 2:2:16;
area = 'PFC';
occThresh = 0.02;
minShift = 20;
nShuf = 1000;
spatialInfo = [];
spatialInfo_s = [];
spInfoZ = [];
sprsty = [];
sprsty_s = [];
sparsityZ = [];
fieldData = [];
aNum = 1;
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
    load(sprintf('%s%s%s_RTimeStrengthSleepNewSpk_20_%02d.mat',dir,animalprefix,area,day));
    pos = loaddatastruct(dir, animalprefix, 'pos', day); % get position
    for ep = 1:length(epochs)
        epoch = epochs(ep);
        posdata = pos{day}{epoch}.data;
        react = RtimeStrength{epoch+1}.reactivationStrengthRun;
        for assem = 1:length(react) %get spikes for each cell
            strength = react{assem}(:,2);
            assemtimes = react{assem}(:,1);
            assemtimes = assemtimes(strength > 5);
            out2d = jds_2DFieldsOptimized(assemtimes, posdata);
            assemblyfields{day}{epoch}{assem}.smoothedspikerate = ...
                out2d.smoothedspikerate;
            rateMap = out2d.smoothedspikerate;
            occMap = out2d.smoothedoccupancy;
            [spInfo, sparsity] = assemblyFieldStats(rateMap, occMap, occThresh);
            spatialInfo = [spatialInfo; spInfo];
            sprsty = [sprsty; sparsity];
            T = posdata(end,1) - posdata(1,1);
            shifts = floor(minShift + (T - 2*minShift) * rand(1, nShuf));
            shufVals = [];
            for s = 1:nShuf
                assemtimes_s = sort(mod(assemtimes - posdata(1,1) + shifts(s), T) + posdata(1,1));
                out2d_s = jds_2DFieldsOptimized(assemtimes_s, posdata);
                rateMap_s = out2d_s.smoothedspikerate;
                [spInfo_s, sparsity_s] = assemblyFieldStats(rateMap_s, occMap, occThresh);
                shufVals = [shufVals; [spInfo_s sparsity_s]];
            end
            spInfoZ = [spInfoZ; (spInfo - mean(shufVals(:,1)))/std(shufVals(:,1))];
            sparsityZ = [sparsityZ; (sparsity - mean(shufVals(:,2)))/std(shufVals(:,2))];
            spatialInfo_s = [spatialInfo_s; median(shufVals(:,1))];
            sprsty_s = [sprsty_s; median(shufVals(:,2))];
            fieldData{aNum}.field = rateMap;
            fieldData{aNum}.occ = occMap;
            fieldData{aNum}.spatialInfo = spInfo;
            fieldData{aNum}.sparsity = sparsity;
            fieldData{aNum}.spatialInfoZ = (spInfo - mean(shufVals(:,1)))/std(shufVals(:,1));
            fieldData{aNum}.sparsityZ = (sparsity - mean(shufVals(:,2)))/std(shufVals(:,2));
            fieldData{aNum}.spatialInfoShuf = shufVals(:,1);
            fieldData{aNum}.sparsityShuf = shufVals(:,2);
            fieldData{aNum}.anim = animalprefix;
            fieldData{aNum}.epoch = epoch;

            aNum = aNum + 1;
        end
    end
end

figure; hold on;
histogram(spInfoZ, 'BinWidth', 0.5, 'FaceColor', [.4 .4 .8]);
xline(1.65, '--k', '95th pct');
xlabel('Spatial info z-score'); ylabel('# assemblies');

figure
[p1 h1] = signrank(spatialInfo,spatialInfo_s)

datacombinedInfo = [spatialInfo; spatialInfo_s];
g1 = repmat({'Data'},length(spatialInfo),1);
g2 = repmat({'Shuffle'},length(spatialInfo_s),1);
g = [g1;g2];

h = boxplot(datacombinedInfo,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
title(['2D Assembly Map Spatial Info-p = ' num2str(p1)])
ylabel('Spatial Information (bits/event)')

figure; hold on;
histogram(sparsityZ, 'BinWidth', 0.5, 'FaceColor', [.4 .4 .8]);
xline(-1.65, '--k', '5th pct');
xlabel('Sparsity z-score'); ylabel('# assemblies');

figure
[p2 h2] = signrank(sprsty,sprsty_s)

datacombinedSparsity = [sprsty; sprsty_s];
g1 = repmat({'Data'},length(sprsty),1);
g2 = repmat({'Shuffle'},length(sprsty_s),1);
g = [g1;g2];

h = boxplot(datacombinedSparsity,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
title(['2D Assembly Map Sparsity-p = ' num2str(p2)])
ylabel('Sparsity')
