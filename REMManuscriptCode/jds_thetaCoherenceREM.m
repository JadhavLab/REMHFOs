function out = jds_thetaCoherenceREM(animalprefixlist,Hz)
% Calculates and plots theta coherence during REM PFC ripples compared to
% baseline.
% -------------------------------------------------------------------------
day = 1;
savedata = 1;
window = [4 4];

g1 = gaussian(3, 10);
movingwin = [1000 20]/1000;
params.Fs = 1500;
params.err = [2 0.05];
params.fpass = [0 40];
params.tapers = [2 3]; 
zCohEps = [];
meanThetaCohEps = [];

zCohEps_s = [];
meanThetaCohEps_s = [];
for a = 1:length(animalprefixlist)
    cohData = [];
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
    load(sprintf('%s%stetinfo.mat',dir,animalprefix));
    load(sprintf('%s%sctxrippletime_REM0%d.mat',dir,animalprefix,day));
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
        triggers = ctxripple{day}{ep}.starttime;
        
        if length(triggers) > 20
            if ep <10
                epochstring = ['0',num2str(ep)];
            else
                epochstring = num2str(ep);
            end
            meanThetaCohTet = [];
            meanThetaCohTet_s = [];
            tPair = 0;
            zCohTets = [];
            zCohTets_s = [];
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
                    shift = randi(length(eegca1),1);
                    eegca1_s = circshift(eegca1,shift);

                    %----- Do the cohereogram calc -----%
                    disp(['Doing day', num2str(day), ' epoch', num2str(ep), ' tets ',num2str(ptet), '/', num2str(hptet)]);
                    [Coh,Phi,~,~,~,t,freq] = cohgramc(eegpfc,eegca1,movingwin,params);
                    Coh = Coh';
                    [Coh_s,Phi,~,~,~,t,freq] = cohgramc(eegpfc,eegca1_s,movingwin,params); %shuffle
                    Coh_s = Coh_s';
                    t = t + tvec(1);
                    fs_c = round(1/(t(2)-t(1))); %length of time bins
                    Nfreq = length(freq); %number of different frequency bands
                    win = [-window(1):(1/fs_c):window(2)];
                    winidx = round(win.*fs_c);

                    %----- Find the mean and std for whole epoch for zscore later -----%
                    meanCohEpoch = mean(Coh,2);
                    stdCohEpoch = std(Coh,0,2);

                    meanCohEpoch_s = mean(Coh_s,2);
                    stdCohEpoch_s = std(Coh_s,0,2);
                    %
                    freqbandindx = (freq >= Hz(1) & freq <= Hz(2))';
                    meanCohEpoch_freq = mean(meanCohEpoch(freqbandindx),1); %epoch mean and std within a certain freq range
                    stdCohEpoch_freq = std(mean(Coh(freqbandindx,:)));

                    meanCohEpoch_freq_s = mean(meanCohEpoch_s(freqbandindx),1); %epoch mean and std within a certain freq range
                    stdCohEpoch_freq_s = std(mean(Coh_s(freqbandindx,:)));

                    tmpTrig = [];
                    allCoh = [];

                    tmpTrig_s = [];
                    allCoh_s = [];
                    for i = 1:length(triggers)
                        if triggers(i) > (t(1) + 2) && triggers(i) < (t(end) - 2) % throw away the first  and last 2s
                            [junk, trigidx] = min(abs(t - triggers(i)));
                            trialindx = [trigidx+winidx];

                            %%% COHEROGRAM %%%
                            allCoh(:,:,i) = Coh(:,trialindx); %all trials cat

                            allCoh_s(:,:,i) = Coh_s(:,trialindx); %all trials cat

                            %%% COHERENCE LEVEL WITHIN FREQ BAND (Hz) %%%
                            %----- take only the coherence within the frequency band -----%
                            cohfreqband = Coh(freqbandindx,trialindx);
                            meanwithinfreqband = mean(cohfreqband,1);
                            zscorefreqband = (meanwithinfreqband - meanCohEpoch_freq)./stdCohEpoch_freq;

                            cohfreqband_s = Coh_s(freqbandindx,trialindx);
                            meanwithinfreqband_s = mean(cohfreqband_s,1);
                            zscorefreqband_s = (meanwithinfreqband_s - meanCohEpoch_freq_s)./stdCohEpoch_freq_s;

                            tmpTrig = [tmpTrig; zscorefreqband];
                            tmpTrig_s = [tmpTrig_s; zscorefreqband_s];
                        end
                    end
                    meanThetaCohTet = [meanThetaCohTet; smoothvect(mean(tmpTrig),g1)]; %mean for pair all trigs
                    meanAllCohTet = mean(allCoh,3); %mean of trials
                    zscoreCohTet = bsxfun(@rdivide,(meanAllCohTet - meanCohEpoch),stdCohEpoch(:)); %z of all trials for pair
                    zCohTets = cat(3, zCohTets, zscoreCohTet); %cat all tet pairs

                    meanThetaCohTet_s = [meanThetaCohTet_s; smoothvect(mean(tmpTrig_s),g1)]; %mean for pair all trigs
                    meanAllCohTet_s = mean(allCoh_s,3); %mean of trials
                    zscoreCohTet_s = bsxfun(@rdivide,(meanAllCohTet_s - meanCohEpoch_s),stdCohEpoch_s(:)); %z of all trials for pair
                    zCohTets_s = cat(3, zCohTets_s, zscoreCohTet_s); %cat all tet pairs

                    cohData{day}{ep}{tPair}.descrip = 'tets - PFC/CA1';
                    cohData{day}{ep}{tPair}.tetpair = [ptet hptet];
                    cohData{day}{ep}{tPair}.coh = meanAllCohTet;
                    cohData{day}{ep}{tPair}.coh_s = meanAllCohTet_s;
                    cohData{day}{ep}{tPair}.window = win;
                    cohData{day}{ep}{tPair}.CohBand = Hz;
                    cohData{day}{ep}{tPair}.MeanCohBand = meanThetaCohTet;
                    cohData{day}{ep}{tPair}.MeanCohBand_s = meanThetaCohTet_s;
                end
            end
            meanThetaCohEps = [meanThetaCohEps; mean(meanThetaCohTet)];
            zCohEps = cat(3, zCohEps, mean(zCohTets,3)); %cat cohereograms per epoch (avg across all pairs)
            meanThetaCohEps_s = [meanThetaCohEps_s; mean(meanThetaCohTet_s)];
            zCohEps_s = cat(3, zCohEps_s, mean(zCohTets_s,3));
        end
    end
    if savedata == 1
        save(sprintf('%s%sCA1PFC_coherence_REMwithshuf%02d.mat', dir,animalprefix,day), 'cohData');
    end
end
meanAnimCoh = mean(zCohEps,3);
sd = 4;
s = gaussian2(sd,(2*sd));
cohgram = filter2(s,meanAnimCoh,'valid');

meanAnimCoh_s = mean(zCohEps_s,3);
cohgram_s = filter2(s,meanAnimCoh_s,'valid'); 

figure
imagesc(win, 1:40, cohgram)
set(gca,'YDir','normal')
colorbar
hold on

figure
imagesc(win, 1:40, cohgram_s)
set(gca,'YDir','normal')
colorbar
hold on

figure; hold on
pl1 = plot(win,mean(meanThetaCohEps,1),'-k','LineWidth',1)
pl2 = plot(win,mean(meanThetaCohEps_s,1),'-r','LineWidth',1)
boundedline(win,mean(meanThetaCohEps),std(meanThetaCohEps)./sqrt(size(meanThetaCohEps,2)),'-k');
boundedline(win,mean(meanThetaCohEps_s),std(meanThetaCohEps_s)./sqrt(size(meanThetaCohEps_s,2)),'-r');
keyboard