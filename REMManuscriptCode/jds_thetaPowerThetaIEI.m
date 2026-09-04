function jds_thetaPowerThetaIEI(animalprefixlist)


day = 1;
% epochs = [1:2:17];
daystring = sprintf('%02d',day);
savedata = 1;
smoothing_width = 0.002; %50 ms
samprate = 1500;
g1 = gaussian(smoothing_width*samprate, ceil(8*smoothing_width*samprate));
phasicIEI = [];
tonicIEI = [];
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);

    load(sprintf('%s%stetinfo.mat',dir,animalprefix));
    load(sprintf('%s%spos0%d.mat',dir,animalprefix,day));% get sws time
    rem = load(sprintf('%s%srem0%d.mat',dir,animalprefix,day));% get sws time
    rem = rem.rem;
    load(sprintf('%s%sphasicrembouts0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%sremeps0%d.mat',dir,animalprefix,day));% get sws time
    epochs = remeps;
    tets = tetinfo{1}{epochs(1)};

    ctxtets = []; %get all ctxriptet tetrodes
    for t = 1:length(tets)
        tmp = tets{t};
        if isfield(tmp, 'descrip')
            if isequal(tmp.descrip, 'riptet')
                %             if isequal(tmp.descrip, 'CA1Ref')
                ctxtets = [ctxtets; t];
            end
        end
    end

    for e = 1:length(epochs)
        ampdataall = [];
        epoch = epochs(e);

        remtime = rem{day}{epoch};
        phasiclist = [phasicrem{day}{epoch}.starttime phasicrem{day}{epoch}.endtime];
        remlist = [remtime.starttime remtime.endtime];

        if epoch <10
            epochstring = ['0',num2str(epoch)];
        else
            epochstring = num2str(epoch);
        end

        load(sprintf('%s%sripples0%d.mat',dir,animalprefix,day));
        rTets = find(~cellfun(@isempty,ripples{day}{epoch}));

        tetsNumRips = [];
        tetsAmpRips = [];
        for scan = 1:length(rTets)
            t = rTets(scan);
            numR = length(ripples{day}{epoch}{t}.startind);
            tetsNumRips = [tetsNumRips; numR];
            tetsAmpRips = [tetsAmpRips; mean(ripples{day}{epoch}{t}.maxthresh)];
        end
        %         [ripcnt idx] = max(tetsNumRips);
        [ripcnt idx] = max(tetsAmpRips);
% 
%         ctxtets = rTets(idx);

        for i = 1:length(ctxtets)
            ctxtet = ctxtets(i);

            if (ctxtet<10)
                ctxtetstring = ['0',num2str(ctxtet)];
            else
                ctxtetstring = num2str(ctxtet);
            end

            curreegfile = [dir,'/EEG/',animalprefix,'eeg', daystring,'-',epochstring,'-',ctxtetstring];
            load(curreegfile);
            eegdat = eeg{day}{epoch}{ctxtet}.data;
            [b,a] = butter(4,[5/(samprate/2) 12/(samprate/2)]);
%             eegdat =
            v = filtfilt(b,a,eegdat);

%             ampdatatmp = theta{day}{epoch}{ctxtet}.data(:,1);
            ampdatatmp = v;
            ampdataall = [ampdataall; ampdatatmp'];
        end
%         times = geteegtimes(theta{day}{epoch}{ctxtet}) ; % construct time array
        times = geteegtimes(eeg{day}{epoch}{ctxtet}) ; % construct time array

        phasicvec = list2vec(phasiclist,times);
        remvec = list2vec(remlist,times);

        tonicvec = remvec & ~phasicvec;
        toniclist = vec2list(tonicvec,times);

        ampdata = mean(ampdataall,1); %mean amplitude across all tetrodes
%         ampdata = conv(ampdata,g1,'same');
        ampdataDiff = diff(ampdata);

        %Zscore amplitude data and find epochs where peak is >=-2 (LFP is
        %flipped)
        z_amp = zscore(ampdataDiff);
        zci = find(diff(sign(z_amp))); %find the zero crossing indices
        startidx = zci(1:end-1); %look at pretty much every interval
        endidx = zci(2:end);
        indices = [startidx' endidx'];

        peak_idx = [];
        for l = 1:length(indices(:,1))
            idx = indices(l,:);
            ampvector = z_amp(idx(1):idx(2));
            if sum(ampvector) > 0
                [max_z maxidx] = max(ampvector);
                peak_idx = [peak_idx; idx(2)];
            end
        end

        thetatimes = times(peak_idx);

        inTonicPeaks = logical(isExcluded(thetatimes,toniclist));
        inPhasicPeaks = logical(isExcluded(thetatimes,phasiclist));

        thetaIEI_tonic = diff(thetatimes(inTonicPeaks));
        thetaIEI_tonic = thetaIEI_tonic(find(thetaIEI_tonic<1));
        thetaIEI_phasic = diff(thetatimes(inPhasicPeaks));
        thetaIEI_phasic = thetaIEI_phasic(find(thetaIEI_phasic<1));
        phasicIEI = [phasicIEI; nanmean(thetaIEI_phasic)];
        tonicIEI = [tonicIEI; nanmean(thetaIEI_tonic)];
    end
end

[p h] = ranksum(tonicIEI,phasicIEI)
datacombinedThetaPow = [tonicIEI; phasicIEI];
g1 = repmat({'Low Theta'},length(tonicIEI),1);
g2 = repmat({'High Theta'},length(phasicIEI),1);
g = [g1;g2];

figure; hold on
h = boxplot(datacombinedThetaPow,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-0.02 0.2])
% yticks([0.11:0.01:0.14])
title(['Theta Inter-peak-intervals-p = ' num2str(p)])
ylabel('IPI')
set(gcf, 'renderer', 'painters')
for i = 1:size(tonicIEI,1)
    if (~isnan(tonicIEI(i))) && (~isnan(phasicIEI(i)))
    x = [1 2];
    y = [tonicIEI(i) phasicIEI(i)];
    plot(x,y,'k')
    end
end
figure; hold on
bar([nanmean(tonicIEI) nanmean(phasicIEI)],'k')
errorbar(1:2,[nanmean(tonicIEI) nanmean(phasicIEI)],...
    [nanstd(tonicIEI)/sqrt(length(tonicIEI)) ...
    nanstd(phasicIEI)/sqrt(length(~isnan(phasicIEI)))],'k','LineStyle','none')
ylim([0.12 0.13])
yticks([0.12:0.002:0.13])
xticks([1:2])
xticklabels({'Low theta power','High theta power'})
title(['Theta Inter-peak-intervals-p = ' num2str(p)])

keyboard