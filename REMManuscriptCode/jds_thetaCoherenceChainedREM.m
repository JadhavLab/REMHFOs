function out = jds_thetaCoherenceChainedREM(animalprefixlist,Hz)

day = 1;
savedata = 1;
window = [4 4];

%theta freq window used = [6 12];

g1 = gaussian(3, 10);
movingwin = [1000 20]/1000;
params.Fs = 1500;
params.err = [2 0.05];
params.fpass = [0 40];
params.tapers = [2 3]; %DO NOT GO LOWER THAN THIS-MESSES UP
meanThetaCohEpsNC = [];
meanThetaCohEpsC = [];
for r = 1:2
    zCohEps = [];
    for a = 1:length(animalprefixlist)
        cohData = [];
        animalprefix = animalprefixlist{a};
        dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
        load(sprintf('%s%stetinfo.mat',dir,animalprefix));
        load(sprintf('%s%sctxrippletime_chainREM0%d.mat',dir,animalprefix,day));
        load(sprintf('%s%sremeps0%d.mat',dir,animalprefix,day));
        rem = load(sprintf('%s%srem0%d.mat',dir,animalprefix,day));
        rem = rem.rem;
        tets = tetinfo{1}{1};
        epochs = remeps;
        pfctets = [];
        ca1tets = [];
        for t = 1:length(tets)
            tmp = tets{t};
            if isfield(tmp, 'descrip')
                if isequal(tmp.descrip, 'ctxriptet')
                    pfctets = [pfctets; t];
                elseif isequal(tmp.descrip, 'riptet')
                    ca1tets = [ca1tets; t];
                end
            end
        end

        for e = 1:length(epochs)
            ep = epochs(e);
            if r == 1
                triggers = ctxripple{day}{ep}.starttimeNC;
            else
                triggers = ctxripple{day}{ep}.starttimeC;
            end


            if length(triggers) > 20
                if ep <10
                    epochstring = ['0',num2str(ep)];
                else
                    epochstring = num2str(ep);
                end
                meanThetaCohTet = [];
                tPair = 0;
                zCohTets = [];
                for t = 1:length(pfctets)
                    ptet = pfctets(t);
                    if (ptet<10)
                        ptetstring = ['0',num2str(ptet)];
                    else
                        ptetstring = num2str(ptet);
                    end

                    for tt = 1:length(ca1tets)
                        hptet = ca1tets(tt);
                        tPair = tPair + 1;
                        if (hptet<10)
                            hptetstring = ['0',num2str(hptet)];
                        else
                            hptetstring = num2str(hptet);
                        end

                        %----- Get the eeg data and time -----%
                        currpfcfile = [dir,'/EEG/',animalprefix,'eeg', '01' ,'-',epochstring,'-',ptetstring];
                        load(currpfcfile);
                        p_eeg = eeg;
                        currca1file = [dir,'/EEG/',animalprefix,'eeg', '01' ,'-',epochstring,'-',hptetstring];
                        load(currca1file);
                        hp_eeg = eeg;

                        eegpfc = p_eeg{day}{ep}{ptet}.data;
                        tvec = geteegtimes(p_eeg{day}{ep}{ptet});
                        eegca1 = hp_eeg{day}{ep}{hptet}.data;

                        %----- Do the cohereogram calc -----%
                        disp(['Doing day', num2str(day), ' epoch', num2str(ep), ' tets ',num2str(ptet), '/', num2str(hptet)]);
                        [Coh,Phi,~,~,~,t,freq] = cohgramc(eegpfc,eegca1,movingwin,params);
                        Coh = Coh';
                        t = t + tvec(1);
                        fs_c = round(1/(t(2)-t(1))); %length of time bins
                        Nfreq = length(freq); %number of different frequency bands
                        win = [-window(1):(1/fs_c):window(2)];
                        winidx = round(win.*fs_c);

                        %----- Find the mean and std for whole epoch for zscore later -----%
                        meanCohEpoch = mean(Coh,2);
                        stdCohEpoch = std(Coh,0,2);
                        %
                        freqbandindx = (freq >= Hz(1) & freq <= Hz(2))';
                        meanCohEpoch_freq = mean(meanCohEpoch(freqbandindx),1); %epoch mean and std within a certain freq range
                        stdCohEpoch_freq = std(mean(Coh(freqbandindx,:)));
                        tmpTrig = [];
                        allCoh = [];
                        for i = 1:length(triggers)
                            if triggers(i) > (t(1) + 2) && triggers(i) < (t(end) - 2) % throw away the first  and last 2s
                                [junk, trigidx] = min(abs(t - triggers(i)));
                                trialindx = [trigidx+winidx];

                                %%% COHEROGRAM %%%
                                allCoh(:,:,i) = Coh(:,trialindx); %all trials cat

                                %%% COHERENCE LEVEL WITHIN FREQ BAND (Hz) %%%
                                %----- take only the coherence within the frequency band -----%
                                cohfreqband = Coh(freqbandindx,trialindx);
                                meanwithinfreqband = mean(cohfreqband,1);
                                zscorefreqband = (meanwithinfreqband - meanCohEpoch_freq)./stdCohEpoch_freq;

                                tmpTrig = [tmpTrig; zscorefreqband];
                            end
                        end
                        meanThetaCohTet = [meanThetaCohTet; smoothvect(mean(tmpTrig),g1)]; %mean for pair all trigs
                        meanAllCohTet = mean(allCoh,3); %mean of trials
                        zscoreCohTet = bsxfun(@rdivide,(meanAllCohTet - meanCohEpoch),stdCohEpoch(:)); %z of all trials for pair
                        zCohTets = cat(3, zCohTets, zscoreCohTet); %cat all tet pairs
                        cohData{day}{ep}{tPair}.descrip = 'tets - PFC/CA1';
                        cohData{day}{ep}{tPair}.tetpair = [ptet hptet];
                        cohData{day}{ep}{tPair}.coh = meanAllCohTet;
                        cohData{day}{ep}{tPair}.cohBaseline = meanCohEpoch;
                        cohData{day}{ep}{tPair}.cohBaselineStd = stdCohEpoch;
                        cohData{day}{ep}{tPair}.window = win;
                        cohData{day}{ep}{tPair}.CohBand = Hz;
                        cohData{day}{ep}{tPair}.MeanCohBand = meanThetaCohTet;

                        clear Coh allCoh eegca1 eegpfc phi t tvec
                    end
                end
                if r == 1
                    meanThetaCohEpsNC = [meanThetaCohEpsNC; mean(meanThetaCohTet)];
                elseif r == 2
                    meanThetaCohEpsC = [meanThetaCohEpsC; mean(meanThetaCohTet)];
                end
                zCohEps = cat(3, zCohEps, mean(zCohTets,3)); %cat cohereograms per epoch (avg across all pairs)
            end
        end
        if (savedata == 1) && (r == 1)
            save(sprintf('%s%sCA1PFC_coherence_REMnonchained%02d.mat', dir,animalprefix,day), 'cohData');
        elseif (savedata == 1) && (r == 2)
            save(sprintf('%s%sCA1PFC_coherence_REMchained%02d.mat', dir,animalprefix,day), 'cohData');
        end
    end
    meanAnimCoh = mean(zCohEps,3);
    sd = 4;
    s = gaussian2(sd,(2*sd));
    cohgram = filter2(s,meanAnimCoh,'valid'); %'Valid' tag excludes edges that would be messed up by filter
    
    figure
    imagesc(win, 1:40, cohgram)
    set(gca,'YDir','normal')
    colorbar
    hold on
    if r == 1
        title('non-chained')
    else
        title('chained')
    end

    figure(2); hold on
    if r == 1
        ax1 = gca;
        ax1.FontSize = 14;
        pl1 = plot(win,mean(meanThetaCohEpsNC,1),'-k','LineWidth',1)
        boundedline(win,mean(meanThetaCohEpsNC),std(meanThetaCohEpsNC)./sqrt(size(meanThetaCohEpsNC,2)),'-k');
    elseif r == 2
        pl2 = plot(win,mean(meanThetaCohEpsC,1),'-r','LineWidth',1)
        boundedline(win,mean(meanThetaCohEpsC),std(meanThetaCohEpsC)./sqrt(size(meanThetaCohEpsC,2)),'-r');
        legend([pl1 pl2],{'NonChained','Chained'})
    end
end
keyboard