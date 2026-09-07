%jds_chainIsolatedRipplesAmpFreqDurREM
%With analyzeData = 0: detects cortical REM ripples across riptet tetrodes,
%measures amplitude and frequency of each event, and saves
%<animal>ctxripple_REM_amp_freq_<day>.mat.
%With analyzeData = 1: loads that file, splits events into chain vs
%isolated (non-chain) REM ripples, and compares ripple amplitude
%(rank-sum, boxplot).

close all;
clear all;
animalprefixlist = {'JS21','JS14','KL8'}; % animal prefix
animaltestday = 1;% animal experimental day
eps = 1:2:17;% epochs
%---- set parameters ----%
minstd = 3; % min std above mean for ripple detection
minrip = 1; % min number of tetrode detected ripples
mismatch = 0;
day = 1;

analyzeData = 1;
nonchain = [];
chain = [];
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};

    % set animal directory
    animaldir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    if analyzeData == 1
        load(sprintf('%s/%sctxripple_REM_amp_freq_%02d.mat', animaldir, animalprefix, day));
        load(sprintf('%s/%sctxrippletime_chainREM%02d.mat', animaldir, animalprefix, day));
        for ep = 1:length(eps)
            ripAssign = [];
            e = eps(ep);
            tmpdata = rippleampfreq{day}{e};
            tmpdata2 = ctxripple{day}{e}.C_sep;
            chainAssign = [];
            for r = 1:length(tmpdata2)
                tmp = tmpdata2{r};
                tmp(:,2) = 1;
                chainAssign = [chainAssign; tmp];
            end
            tmpdata3 = ctxripple{day}{e}.starttimeNC;
            tmpdata3(:,2) = 0;
            ripAssign = sortrows([tmpdata3; chainAssign],1);
            if ~isempty(tmpdata)
                for rr = 1:length(ripAssign(:,1))
                    chainidx = ripAssign(find(ripAssign(rr,1) == tmpdata(:,1)),2);
                    ampTmp = tmpdata(find(ripAssign(rr,1) == tmpdata(:,1)),4);
                    if chainidx == 0
                        nonchain = [nonchain; ampTmp];
                    elseif chainidx == 1
                        chain = [chain; ampTmp];
                    end
                end
            end
        end
    else
        eegdir = [animaldir,'EEG/'];% EEG directory
        tetinfo = loaddatastruct(animaldir, animalprefix, 'tetinfo'); % get tetrode info
        % day loop
        for day = animaltestday
            load(sprintf('%s/%sctxrippletime_REM%02d.mat', animaldir, animalprefix, day));
            rip = ctxripple; clear ripple


            % loop for each epoch per day
            for id = 1:length(eps)
                d = day;%day
                e = eps(id);%epoch
                riptimes = [rip{day}{e}.starttime rip{day}{e}.endtime];

                %combine riptimes
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
                baselinespecgram_forref(animalprefix, d, e, tetlist, 'fpass',[0 400])% high frequency range, 0-400 Hz, for all tets

                ripples = loaddatastruct(animaldir, animalprefix, 'ctxripples', d);% load ripple info

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
                        % get the indeces for the ripples with maxthresh above minstd
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
            save(sprintf('%s/%sctxripple_REM_amp_freq_%02d.mat', animaldir, animalprefix, d), 'rippleampfreq');% save files
            clear rippleampfreq;
        end
    end
end

if analyzeData == 1
    meanAmpIso = nanmean(nonchain)
    semAmpIso = nanstd(nonchain)./sqrt(length(find(~isnan(nonchain))))
    medianAmpIso = nanmedian(nonchain)

    meanAmpChain = nanmean(chain)
    semAmpChain = nanstd(chain)./sqrt(length(find(~isnan(chain))))
    medianAmpChain = nanmedian(chain)

    [p h] = ranksum(nonchain,chain)
    datacombinedRipAmp = [nonchain; chain];
    g1 = repmat({'Isolated'},length(nonchain),1);
    g2 = repmat({'Chain'},length(chain),1);
    g = [g1;g2];

    figure;
    h = boxplot(datacombinedRipAmp,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
    title(['Ctx REM ripple amplitude-p = ' num2str(p)])
    ylim([2 12])
    set(gcf, 'renderer', 'painters')
end

