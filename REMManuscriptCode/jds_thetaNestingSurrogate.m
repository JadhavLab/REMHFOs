% Calculate theta nesting of gamma and ripple oscillations. Compared nested
% oscillation frequencies and amplitude of nested oscillations vs surrogate
% waveforms
% -------------------------------------------------------------------------
clear all

%% Working with actual LFP signals (example .mat file available at github)
animalprefixlist = {'KL8','JS17','JS15','JS14','JS12','JS13','JS34','BG1','JS21','ZT2'};
days = 1;
modFreqRipAll = [];
modFreqRipAllSurr = [];
modFreqGamAll = [];
modFreqGamAllSurr = [];
sigTetRip = [];
sigTetRipWavs = [];
sigTetRipWavsSurr = [];
sigTetGam = [];
sigTetGamWavs = [];
sigTetGamWavsSurr = [];
gamPref = [];
ripPref = [];
gamRipMI = [];

ripDataSurr = [];
gamDataSurr = [];
for d = 1:length(days)
    day = days(d);
    daystring = sprintf('%02d',day);
    for a = 1:length(animalprefixlist)
        allEpSpecs = [];
        animalprefix = animalprefixlist{a};
        dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
        rem = load(sprintf('%s%srem0%d.mat',dir,animalprefix,day));
        load(sprintf('%s%sphasicrembouts0%d.mat',dir,animalprefix,day));

        rem = rem.rem;
        remPh = phasicrem;
        load(sprintf('%s%sremeps0%d.mat',dir,animalprefix,day));
        epochs = remeps;
        for ep = 1:length(epochs)
            epoch = epochs(ep);

            if epoch <10
                epochstring = ['0',num2str(epoch)];
            else
                epochstring = num2str(epoch);
            end

            if remPh{day}{epoch}.total_duration == 0
                continue
            end

            remList = [rem{day}{epoch}.starttime rem{day}{epoch}.endtime];
            phasicList = [remPh{day}{epoch}.starttime remPh{day}{epoch}.endtime];

            load(sprintf('%s%stetinfo.mat',dir,animalprefix));

            tets = tetinfo{day}{epoch};

            ctxtets = []; %get all ctxriptet tetrodes
            for t = 1:length(tets)
                tmp = tets{t};
                if isfield(tmp, 'area')
                    if isequal(tmp.area, 'PFC')
                        ctxtets = [ctxtets; t];
                    end
                end
            end

            tetlist = ctxtets;

            allTetSpecs = [];

            for t = 1:length(tetlist)
                tetPhaseModRip = [];
                tetPhaseModGam = [];
                tetRipMeanWav = [];
                tetGamMeanWav = [];
                tetRipMeanWavSurr = [];
                tetGamMeanWavSurr = [];
                tet = tetlist(t);

                if (tet<10)
                    tetstring = ['0',num2str(tet)];
                else
                    tetstring = num2str(tet);
                end

                curreegfile = [dir,'/EEG/',animalprefix,'eegref', daystring,'-',epochstring,'-',tetstring];
                load(curreegfile);

                % Define EEG
                e = eegref{day}{epoch}{tet}.data';
                times = geteegtimes(eegref{day}{epoch}{tet});
                tmpMI = [];
                for r = 1:size(remList,1)
                    remEp = remList(r,:);
                    stIdx = lookup(remEp(1),times);
                    endIdx = lookup(remEp(2),times);
                    lfp = e(stIdx:endIdx);

                    data_length = length(lfp);
                    srate = 1500;
                    dt = 1/srate;
                    tt = (1:data_length)*dt;

                    %% Define the amplitude- and phase-frequencies

                    PhaseFreqVector=2:1:12; %theta
                    AmpFreqVector=40:1:250; %gamma to ripple

                    PhaseFreq_BandWidth=1;
                    AmpFreq_BandWidth=2;

                    %% Define phase bins

                    nbin = 18; % number of phase bins
                    position=zeros(1,nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
                    winsize = 2*pi/nbin;
                    for j=1:nbin
                        position(j) = -pi+(j-1)*winsize;
                    end

                    %% Compute MI and comodulogram
                    Pf1 = 2;
                    Pf2 = 12;

                    Af1 = 150;
                    Af2 = 250;

                    [MIrip,MeanAmp] = ModIndex_v1(lfp,srate,Pf1,Pf2,Af1,Af2,position);
                    tetPhaseModRip = [tetPhaseModRip; MeanAmp/sum(MeanAmp)];
                    Af1 = 40;
                    Af2 = 100;

                    [MIgam,MeanAmp] = ModIndex_v1(lfp,srate,Pf1,Pf2,Af1,Af2,position);
                    tetPhaseModGam = [tetPhaseModGam; MeanAmp/sum(MeanAmp)];

                    tmpMI = [tmpMI; [MIgam MIrip]];
                end
                gamRipMI = [gamRipMI; mean(tmpMI,1)];
                
                [maxVal maxRipIdx] = max(mean(tetPhaseModRip,1));
                [maxVal maxGamIdx] = max(mean(tetPhaseModGam,1));
                currripfile = [dir,'/EEG/',animalprefix,'ripple', daystring,'-',epochstring,'-',tetstring];
                load(currripfile);
                ripamp = double(ripple{day}{epoch}{tet}.data(:,1));
                currgammafile = [dir,'/EEG/',animalprefix,'gamma', daystring,'-',epochstring,'-',tetstring];
                load(currgammafile);
                gammaamp = double(gamma{day}{epoch}{tet}.data(:,1));
                currthetafile = [dir,'/EEG/',animalprefix,'theta', daystring,'-',epochstring,'-',tetstring];
                load(currthetafile);
                thetaPh = double(theta{day}{epoch}{tet}.data(:,2));
                thetaPh = thetaPh/10000; %to radians
                if maxGamIdx < length(position)
                    prefPhGam = position(maxGamIdx:maxGamIdx+1);
                else
                    prefPhGam = position(maxGamIdx-1:maxGamIdx);
                end
               
                if maxRipIdx < length(position)
                    prefPhRip = position(maxRipIdx:maxRipIdx+1);
                else
                    prefPhRip = position(maxRipIdx-1:maxRipIdx);
                end
                
                goodphases = logical(thetaPh > prefPhGam(1) & thetaPh < prefPhGam(2)); %logical vector of phases that fall in bin
                prefPhTimes = vec2list(goodphases,times); %list of phase bin times
                inPhasic = logical(isExcluded(prefPhTimes(:,1),phasicList)); %times in phasic/high theta times
                inPhasicTimes = prefPhTimes(inPhasic,:);
                inPhasicTimes_surr = jitter(inPhasicTimes(:,1));
                [sorted sortI] = sort(inPhasicTimes_surr);
                tDiff = inPhasicTimes(:,2)-inPhasicTimes(:,1);
                tDiff = tDiff(sortI);
                inPhasicTimes_surr = [inPhasicTimes_surr (inPhasicTimes_surr + tDiff)];
                phaseIdxG = [];
                phaseIdxSurr = [];
                for t = 1:size(inPhasicTimes,1)
                    idxst = lookup(inPhasicTimes(t,1), times);
                    idxend = lookup(inPhasicTimes(t,2), times);
                    phaseIdxG = [phaseIdxG; [idxst idxend]];

                    idxstSurr = lookup(inPhasicTimes_surr(t,1), times);
                    idxendSurr = lookup(inPhasicTimes_surr(t,2), times);
                    phaseIdxSurr = [phaseIdxSurr; [idxstSurr idxendSurr]];
                end
                for r = 1:size(phaseIdxG,1)
                    tmpidx = phaseIdxG(r,:);
                    [M Igam] = max(gammaamp(tmpidx(1):tmpidx(2)));
                    Igam = Igam + tmpidx(1);
                    gamLfp = e(Igam-300:Igam+300);
                    tetGamMeanWav = [tetGamMeanWav; gamLfp];
                    
                    tmpidxSurr = phaseIdxSurr(r,:);
                    [M IgamSurr] = max(gammaamp(tmpidxSurr(1):tmpidxSurr(2)));
                    IgamSurr = IgamSurr + tmpidxSurr(1);
                    gamLfpSurr = e(IgamSurr-300:IgamSurr+300);
                    tetGamMeanWavSurr = [tetGamMeanWavSurr; gamLfpSurr];
                end
                if size(phaseIdxG,1) > 20
                    tetGamMeanWavAll = mean(tetGamMeanWav,1);
                    tetGamMeanWavAllSurr = mean(tetGamMeanWavSurr,1);
                    gamWin = tetGamMeanWavAll(301-75:301+75);
                    gamWinSurr = tetGamMeanWavAllSurr(301-75:301+75);
                    numCyclesGam = sum(islocalmax(gamWin));
                    if numCyclesGam >= 5
                        sigTetGam = [sigTetGam; 1];
                        sigTetGamWavs = [sigTetGamWavs; tetGamMeanWavAll];
                        modFreqGam = numCyclesGam*10;
                        modFreqGamAll = [modFreqGamAll; modFreqGam];
                        gamPref = [gamPref; position(maxGamIdx)];
                    else
                        sigTetGam = [sigTetGam; 0];
                    end

                    
                    [pks, pklocs] = findpeaks(gamWin);                               
                    [troughs,trlocs] = findpeaks(-gamWin);    
                    tmpdata = [];
                    for p = 1:length(pklocs)
                        tmp = lookup(pklocs(p),trlocs);
                        tmpdata = [tmpdata; abs(gamWin(pklocs(p))*0.195-gamWinSurr(trlocs(tmp))*0.195)];
                    end
                    
                    [pks, pklocs] = findpeaks(gamWinSurr);                               
                    [troughs,trlocs] = findpeaks(-gamWinSurr);  
                    tmpdata2 = [];
                    for p = 1:length(pklocs)
                        tmp = lookup(pklocs(p),trlocs);
                        tmpdata2 = [tmpdata2; abs(gamWinSurr(pklocs(p))*0.195-gamWinSurr(trlocs(tmp))*0.195)];
                    end
                    gamDataSurr = [gamDataSurr; [mean(tmpdata) mean(tmpdata2)]];
                end

                goodphases = logical(thetaPh > prefPhRip(1) & thetaPh < prefPhRip(2)); %logical vector of phases that fall in bin
                prefPhTimes = vec2list(goodphases,times); %list of phase bin times
                inPhasic = logical(isExcluded(prefPhTimes(:,1),phasicList)); %times in phasic/high theta times
                inPhasicTimes = prefPhTimes(inPhasic,:);
                inPhasicTimes_surr = jitter(inPhasicTimes(:,1));
                [sorted sortI] = sort(inPhasicTimes_surr);
                tDiff = inPhasicTimes(:,2)-inPhasicTimes(:,1);
                tDiff = tDiff(sortI);
                inPhasicTimes_surr = [inPhasicTimes_surr (inPhasicTimes_surr + tDiff)];
                phaseIdx = [];
                phaseIdxSurr = [];
                for t = 1:size(inPhasicTimes,1)
                    idxst = lookup(inPhasicTimes(t,1), times);
                    idxend = lookup(inPhasicTimes(t,2), times);
                    phaseIdx = [phaseIdx; [idxst idxend]];

                    idxstSurr = lookup(inPhasicTimes_surr(t,1), times);
                    idxendSurr = lookup(inPhasicTimes_surr(t,2), times);
                    phaseIdxSurr = [phaseIdxSurr; [idxstSurr idxendSurr]];
                end
                for r = 1:size(phaseIdx,1)
                    tmpidx = phaseIdx(r,:);
                    [M Irip] = max(ripamp(tmpidx(1):tmpidx(2)));
                    Irip = Irip + tmpidx(1);
                    ripLfp = e(Irip-300:Irip+300);
                    tetRipMeanWav = [tetRipMeanWav; ripLfp];

                    tmpidxSurr = phaseIdxSurr(r,:);
                    [M IripSurr] = max(ripamp(tmpidxSurr(1):tmpidxSurr(2)));
                    IripSurr = IripSurr + tmpidxSurr(1);
                    ripLfpSurr = e(IripSurr-300:IripSurr+300);
                    tetRipMeanWavSurr = [tetRipMeanWavSurr; ripLfpSurr];
                end

                if size(phaseIdx,1) > 20
                    tetRipMeanWavAll = mean(tetRipMeanWav,1);
                    tetRipMeanWavAllSurr = mean(tetRipMeanWavSurr,1);
                    ripWin = tetRipMeanWavAll(301-75:301+75);
                    ripWinSurr = tetRipMeanWavAllSurr(301-75:301+75);
                    numCyclesRip = sum(islocalmax(ripWin));
                    if numCyclesRip >= 5
                        sigTetRip = [sigTetRip; 1];
                        sigTetRipWavs = [sigTetRipWavs; tetRipMeanWavAll];
                        modFreqRip = numCyclesRip*10;
                        modFreqRipAll = [modFreqRipAll; modFreqRip];
                        ripPref = [ripPref; position(maxRipIdx)];
                    else
                        sigTetRip = [sigTetRip; 0];
                    end
                   
                    [pks, pklocs] = findpeaks(ripWin);                               
                    [troughs,trlocs] = findpeaks(-ripWin);    
                    tmpdata = [];
                    for p = 1:length(pklocs)
                        tmp = lookup(pklocs(p),trlocs);
                        tmpdata = [tmpdata; abs(ripWin(pklocs(p))*0.195-ripWin(trlocs(tmp))*0.195)];
                    end
                    
                    [pks, pklocs] = findpeaks(ripWinSurr);                               
                    [troughs,trlocs] = findpeaks(-ripWinSurr);  
                    tmpdata2 = [];
                    for p = 1:length(pklocs)
                        tmp = lookup(pklocs(p),trlocs);
                        tmpdata2 = [tmpdata2; abs(ripWinSurr(pklocs(p))*0.195-ripWinSurr(trlocs(tmp))*0.195)];
                    end

                    ripDataSurr = [ripDataSurr; [mean(tmpdata) mean(tmpdata2)]];
                end
            end
        end
    end
end
%% Phase pref difference
[p,U2] = watsons_U2_approx_p(ripPref, gamPref)

datacombinedPref = [ripPref; gamPref];
g1 = repmat({'Theta pref - Ripple'},length(ripPref),1);
g2 = repmat({'Theta pref - Gamma'},length(gamPref),1);
g = [g1;g2];

figure
h = boxplot(datacombinedPref,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
title(['Theta phase preference-p = ' num2str(p)])
set(gcf, 'renderer', 'painters')

%% Frequency of nested oscillation
p2 = ranksum(modFreqRipAll,modFreqGamAll);

datacombinedWavFreq = [modFreqRipAll; modFreqGamAll];
g1 = repmat({'Wav freq - Ripple'},length(modFreqRipAll),1);
g2 = repmat({'Wav freq - Gamma'},length(modFreqGamAll),1);
g = [g1;g2];

figure
h = boxplot(datacombinedWavFreq,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-0.2 1.2])
title(['Mean wav freq-p = ' num2str(p2)])
set(gcf, 'renderer', 'painters')

%% Nested ripple oscillation vs surrogate
p3 = signrank(ripDataSurr(:,1),ripDataSurr(:,2));

datacombinedRipplePeaks = [ripDataSurr(:,1); ripDataSurr(:,2)];
g1 = repmat({'Ripple peaks'},length(ripDataSurr(:,1)),1);
g2 = repmat({'Surrogate peaks'},length(ripDataSurr(:,2)),1);
g = [g1;g2];

figure
h = boxplot(datacombinedRipplePeaks,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-0.2 1.2])
title(['Mean wav peak-p = ' num2str(p3)])
ylim([-2 6])
set(gcf, 'renderer', 'painters')

%% Nested gamma oscillation vs surrogate

p4 = signrank(gamDataSurr(:,1),gamDataSurr(:,2));

datacombinedGammaPeaks = [gamDataSurr(:,1); gamDataSurr(:,2)];
g1 = repmat({'Gamma peaks'},length(gamDataSurr(:,1)),1);
g2 = repmat({'Surrogate peaks'},length(gamDataSurr(:,2)),1);
g = [g1;g2];

figure
h = boxplot(datacombinedGammaPeaks,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-0.2 1.2])
title(['Mean wav peak-p = ' num2str(p4)])
ylim([-5 35])
set(gcf, 'renderer', 'painters')

keyboard
