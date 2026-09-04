%getripples_amp_freq_sleep_rem
close all;
clear all;
%%
animalprefixlist = {'ZT2','JS17','JS15','JS14','JS12','JS13','JS34','BG1','JS21','KL8'}; % animal prefix

animaltestday = 1;% animal experimental day
% eps = 1:2:17;% epochs
% eps = 11;
%%
%---- set parameters ----%
minstd = 3; % min std above mean for ripple detection
minrip = 1; % min number of tetrode detected ripples
% minenergy = 0; % min energy threshold
velfilter = 4; %velocity <= 4cm for ripple detection
matcheegtime = 0; % match EEG time?
mismatch = 0;
day = 1;

analyzeData = 1;
remAmp = [];
nremAmp = [];
remFreq = [];
nremFreq = [];
remLength = [];
nremLength = [];
%%
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};

    % set animal directory
    animaldir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    if analyzeData == 1
        load(sprintf('%s/%sremctxripple_amp_freq2_%02d.mat', animaldir, animalprefix, day));
        rippleampfreqrem = rippleampfreq; clear rippleampfreq
        load(sprintf('%s/%sctxripple_amp_freq_%02d.mat', animaldir, animalprefix, day));
        rippleampfreqnrem = rippleampfreq; clear rippleampfreq
        load(sprintf('%s/%sremeps%02d.mat', animaldir, animalprefix, day));
        for ep = 1:length(remeps)
            e = remeps(ep);
            tmpdata1 = rippleampfreqrem{day}{e};
            tmpdata2 = rippleampfreqnrem{day}{e};
            if (~isempty(tmpdata1)) && (~isempty(tmpdata2))
                if (length(tmpdata1(:,1))>10) && (length(tmpdata2(:,1))>10)
                    remAmp = [remAmp; mean(tmpdata1(:,4))];
                    nremAmp = [nremAmp; mean(tmpdata2(:,4))];
                    remFreq = [remFreq; mean(tmpdata1(:,5))];
                    nremFreq = [nremFreq; mean(tmpdata2(:,5))];
                    remLength = [remLength; mean(tmpdata1(:,2)-tmpdata1(:,1))];
                    nremLength = [nremLength; mean(tmpdata2(:,2)-tmpdata2(:,1))];
%                     if mean(tmpdata2(:,2)-tmpdata2(:,1)) >0.3
%                         keyboard
%                     end
                end
            end
        end
    else
        eegdir = [animaldir,'EEG/'];% EEG directory
        %%
        tetinfo = loaddatastruct(animaldir, animalprefix, 'tetinfo'); % get tetrode info
        %%
        % day loop
        for day = animaltestday
            load(sprintf('%s/%sctxrippletime_SWS%02d.mat', animaldir, animalprefix, day));
            rem_rip = ctxripple; 
            load(sprintf('%s/%sspikes%02d.mat', animaldir, animalprefix, day));
            load(sprintf('%s/%sremeps%02d.mat', animaldir, animalprefix, day));

            % loop for each epoch per day
            for id = 1:length(remeps)
                d = day;%day
                e = remeps(id);%epoch
                rem_riptimes = [rem_rip{day}{e}.starttime rem_rip{day}{e}.endtime];

                %combine riptimes
                riptimes = rem_riptimes;
                disp(['Animal: ',animalprefix,' Epoch:',num2str(e)])% display current animal, day and epoch

                % calculate using riptet only
                tetfilter = 'isequal($descrip, ''ctxriptet'')';% use riptet only
                tetlist =  evaluatefilter(tetinfo{d}{e}, tetfilter);
                tetlist = unique(tetlist(:,1))';
                if isequal(animalprefix,'JS14')
                    tetlist(find(tetlist == 17)) = [];
                end

                % get the mean and std for each frequency channel for z-scoring
                % note that only run it for the first time
%                 baselinespecgram_forref(animalprefix, d, e, tetlist, 'fpass',[0 400])% high frequency range, 0-400 Hz, for all tets

                ripples = loaddatastruct(animaldir, animalprefix, 'ctxripples01', d);% load ripple info

                r = ripples{d}{e}{tetlist(1)};
                % time range, 10ms bin
                times = r.timerange(1):0.001:r.timerange(end);
                %reset
                nrip = zeros(size(times));
                nstd=[];
                ripplestd = zeros(size(times));

                % tetrode loop
                if ~isempty(riptimes)
                    for t = 1:length(tetlist)
                        tmprip = ripples{d}{e}{tetlist(t)};
                        % get the indeces for the ripples with energy above minenergy
                        % and maxthresh above minstd
%                         rvalid = find((tmprip.energy >= minenergy) & (tmprip.maxthresh >= minstd));
                        rvalid = find(tmprip.maxthresh >= minstd);
                        rtimes = [tmprip.starttime(rvalid) tmprip.endtime(rvalid)];
                        tmpripplestd = [tmprip.maxthresh(rvalid) tmprip.maxthresh(rvalid)];
                        % create another parallel vector with bordering times for zeros
                        nrtimes = [(rtimes(:,1) - 0.00001) (rtimes(:,2) + 0.00001)];
                        rtimes = reshape(rtimes', length(rtimes(:)), 1);
                        rtimes(:,2) = 1;
                        tmpriplestd = [rtimes(:,1) tmpripplestd(:)];
                        nrtimes = [r.timerange(1) ; reshape(nrtimes', ...
                            length(nrtimes(:)), 1) ; r.timerange(2)];
                        nrtimes(:,2) = 0;
                        % create a new list with all of the times in it
                        tlist = sortrows([rtimes ; nrtimes]);
                        [junk, ind] = unique(tlist(:,1));
                        tlist = tlist(ind,:);

                        stdlist = sortrows([tmpriplestd ; nrtimes]);
                        stdlist =stdlist(ind,:);
                        nrip = nrip + interp1(tlist(:,1), tlist(:,2), times, 'nearest');
                        nstd(t,:) = interp1(stdlist(:,1), stdlist(:,2), times, 'nearest');  % carry forward amplitude of ripple
                    end

                    %find the start and end borders of each ripple
                    inripple = (nrip >= minrip);
                    startrippleind = find(diff(inripple) == 1)+1;
                    endrippleind = find(diff(inripple) == -1)+1;
                    ripplestdout = [];

                    if (endrippleind(1) < startrippleind(1))
                        endrippleind = endrippleind(2:end);
                    end
                    if (endrippleind(end) < startrippleind(end))
                        startrippleind = startrippleind(1:end-1);
                    end
                    startripple = times(startrippleind);
                    endripple = times(endrippleind);
                    %----- measure amplitude of each ripple-----%
                    % Get amplitude of "global" ripple: maximum across tetrodes
                    [max_nstd,tetid] = max(nstd,[],1);
                    ampripple = max_nstd(startrippleind);
                    riptet = tetid(startrippleind);

                    out = [startripple(:) endripple(:) ampripple(:)]; % amplitude of ripple
                    riptimes(:,4) = 0;
                    for r = 1:length(riptimes(:,1))
                        ripmidtmp = riptimes(r,1) + ((riptimes(r,2) - riptimes(r,1))/2);
                        idx = find((ripmidtmp > out(:,1)) & (ripmidtmp < out(:,2)));
                        if ~isempty(idx)
                            amptmp = out(idx,3);
                            riptimes(r,4) = amptmp;
                            riptetnew(r) = riptet(idx);
                        else
                            amptmp = NaN;
                            riptimes(r,4) = amptmp;
                            riptetnew(r) = NaN;
                            mismatch = mismatch+1;
                        end
                    end
                    riptimes(:,5) = 0;

                    %----- measure frequncy of each ripple-----%
                    for r = 1:length(riptimes(:,1))
                        riptime = riptimes(r,1:2);
                        if ~isnan(riptetnew(r))
                            riptimes(r,5) = ripple_frequency_fun(animalprefix, d, e, tetlist(riptetnew(r)), riptime);
                        else
                            riptimes(r,5) = NaN;
                        end
                    end
                else
                    riptimes = [];
                end
                rippleampfreq{d}{e} = riptimes;% save result
                clear out; clear freqripple; clear riptimes
            end
            save(sprintf('%s/%sctxripple_amp_freq2_%02d.mat', animaldir, animalprefix, d), 'rippleampfreq');% save files
            clear rippleampfreq;
        end
    end
end

if analyzeData == 1
    
    meanAmpInd = nanmean(remAmp)
    semAmpInd = nanstd(remAmp)./sqrt(length(find(~isnan(remAmp))))
    medianAmpInd = nanmedian(remAmp)

    meanAmpCoord = nanmean(nremAmp)
    semAmpCoord = nanstd(nremAmp)./sqrt(length(find(~isnan(nremAmp))))
    medianAmpCoord = nanmedian(nremAmp)

    meanFreqInd = nanmean(remFreq)
    semFreqInd = nanstd(remFreq)./sqrt(length(find(~isnan(remFreq))))
    medianFreqInd = nanmedian(remFreq)

    meanFreqCoord = nanmean(nremFreq)
    semFreqcoord = nanstd(nremFreq)./sqrt(length(find(~isnan(nremFreq))))
    medianFreqCoord = nanmedian(nremFreq)

    meanLengthInd = nanmean(remLength)
    semLengthInd = nanstd(remLength)./sqrt(length(find(~isnan(remLength))))
    medianLengthInd = nanmedian(remLength)

    meanLengthCoord = nanmean(nremLength)
    semLengthcoord = nanstd(nremLength)./sqrt(length(find(~isnan(nremLength))))
    medianLengthCoord = nanmedian(nremLength)

    [p h] = ranksum(remAmp,nremAmp)
    datacombinedRipAmp = [remAmp; nremAmp];
    g1 = repmat({'REM'},length(remAmp),1);
    g2 = repmat({'NREM'},length(nremAmp),1);
    g = [g1;g2];

    figure;
    h = boxplot(datacombinedRipAmp,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
    ylim([-0.02 0.2])
    title(['PFC - Ripple amplitude-p = ' num2str(p)])
    ylim([3 6])
    set(gcf, 'renderer', 'painters')
    ylabel('Amplitude (STD)')

    [p2 h2] = ranksum(remFreq,nremFreq)
    datacombinedRipFreq = [remFreq; nremFreq];
    g1 = repmat({'REM'},length(remFreq),1);
    g2 = repmat({'NREM'},length(nremFreq),1);
    g = [g1;g2];

    figure;
    h = boxplot(datacombinedRipFreq,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
    ylim([150 210])
    title(['PFC - Ripple frequency-p = ' num2str(p2)])
    set(gcf, 'renderer', 'painters')
    ylabel('Frequency (Hz)')

    [p3 h3] = ranksum(remLength,nremLength)
    datacombinedRipLength = [remLength; nremLength];
    g1 = repmat({'REM'},length(remLength),1);
    g2 = repmat({'NREM'},length(nremLength),1);
    g = [g1;g2];

    figure;
    h = boxplot(datacombinedRipLength,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
    ylim([0.03 0.15])
    title(['PFC - Ripple lengths-p = ' num2str(p3)])
    set(gcf, 'renderer', 'painters')
    ylabel('Lengths (s)')
end
keyboard

