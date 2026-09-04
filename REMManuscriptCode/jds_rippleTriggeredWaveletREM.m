function out = jds_rippleTriggeredWaveletREM(animalprefixlist, days)
% Plot ripple triggered spectrogram
% -------------------------------------------------------------------------

%parse the options
Fs = 1500;
win = [0.5 0.5];
S0 = 2*(1/Fs);
DJ = 0.15;
DT = 1/Fs;
N = sum(win)*Fs;
J1 = fix(log(N*DT/S0)/log(2))/DJ;
allAnimSpecs = [];
allAnimSpecsEps = [];
numrips = 0;
for d = 1:length(days)
    day = days(d);
    daystring = sprintf('%02d',day);
    for a = 1:length(animalprefixlist)
        allEpSpecs = [];
        animalprefix = animalprefixlist{a};
        dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
        load(sprintf('%s%sctxrippletime_REM0%d.mat',dir,animalprefix,day));
        load(sprintf('%s%sremeps0%d.mat',dir,animalprefix,day));
        epochs = remeps;
        for i = 1:length(epochs)
            epoch = epochs(i);

            if epoch <10
                epochstring = ['0',num2str(epoch)];
            else
                epochstring = num2str(epoch);
            end

            load(sprintf('%s%stetinfo.mat',dir,animalprefix));

            rip = ctxripple{day}{epoch};
            
            if isempty(rip)
                continue
            end

            riptimes = [rip.starttime rip.endtime];
            if isempty(riptimes)
                continue
            end
            numrips = numrips+size(riptimes,1);
            
            load(sprintf('%s%sctxripples0%d.mat',dir,animalprefix,day));
            rTets = find(~cellfun(@isempty,ctxripples{day}{epoch}));

            tetsNumRips = [];
            for scan = 1:length(rTets)
                t = rTets(scan);
                numR = length(ctxripples{day}{epoch}{t}.startind);
                tetsNumRips = [tetsNumRips; numR];
            end
            [ripcnt idx] = max(tetsNumRips);

            tetlist = rTets(idx); %take the tetrode with the most ripples detected
            if ~isempty(riptimes)
                tets = tetinfo{day}{epoch};

                allTetSpecs = [];

                for ii = 1:length(tetlist)

                    tet = tetlist(ii);

                    if (tet<10)
                        tetstring = ['0',num2str(tet)];
                    else
                        tetstring = num2str(tet);
                    end

                    curreegfile = [dir,'/EEG/',animalprefix,'eegref', daystring,'-',epochstring,'-',tetstring];
                    load(curreegfile);

                    % Define EEG
                    e = eegref{day}{epoch}{tet}.data';
                    starttime = eegref{day}{epoch}{tet}.starttime;
                    endtime = (length(e)-1) * (1 / Fs);

                    % Define triggering events as the start of each ripple
                    triggers = riptimes(:,1)-starttime;

                    %Remove triggering events that are too close to the beginning or end
                    while triggers(1)<win(1)
                        triggers(1) = [];
                    end
                    while triggers(end)> endtime-win(2)
                        triggers(end) = [];
                    end

                    % Compute a z-scored spectrogram using the mean and std for the entire session
                    meanbase = [];
                    idxcnt = 0;
                    FsBase = Fs*(sum(win));
                    for s = 1:length(e)/(FsBase)
                        basetmp = e((idxcnt+1):(FsBase+idxcnt));
                        ptmp = wavelet(basetmp,1/Fs,0,DJ,S0,J1, 'MORLET');
                        ptmp = abs(ptmp);
                        idxcnt = idxcnt + FsBase;
                        meanbase = cat(3, meanbase, ptmp);
                    end

                    allmeanbase = mean(meanbase,3);

                    meanP = allmeanbase;
                    stdP = std(meanbase,[],3);

                    % Calculate the event triggered spectrogram
                    windoweddata = jds_createdatamatc_wavelet(e,triggers,Fs,win);
                    allwavsp = [];
                    for s = 1:length(windoweddata(:,1))
                        [temp,period] = wavelet(windoweddata(s, 1:end),1/Fs,0,DJ,S0,J1,'MORLET');
                        pow=abs(temp);
                        allwavsp(:,:,s)=(pow-meanP)./stdP;
                    end

                    allTetSpecs = cat(3, allTetSpecs, mean(allwavsp,3));
                    clear allwavsp
                end
            end
            allEpSpecs = cat(3, allEpSpecs, mean(allTetSpecs,3));
            allAnimSpecsEps = [allAnimSpecsEps; mean(mean(allTetSpecs,3),2)'];
        end
        clear triggers riptimes
        allAnimSpecs = cat(3, allAnimSpecs, mean(allEpSpecs,3));
    end
end

%plot the spectrogram
meanspect = mean(allAnimSpecs,3);
figure
imagesc(meanspect)
frequency=round(arrayfun(@(x) 1/x, period));
yLabels = arrayfun(@num2str,frequency,'UniformOutput',false);
yticklabels({yLabels{10},yLabels{20},yLabels{30},yLabels{40},yLabels{50},yLabels{60}})
ylabel('Frequency')
xlim([0 1500])
xticks([0:750:1500])
xticklabels({'-500','0','500'})
xlabel('Time From Ripple Onset (ms)')
colorbar
set(gcf, 'renderer', 'painters')
end

