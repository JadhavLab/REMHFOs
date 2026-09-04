% Calculate and plot phase amplitude coupling in PFC - theta, gamma, ripple
% coupling. Calculate PAC modulation index (MI) for discreet gamma and
% ripples bands
% -------------------------------------------------------------------------

%% Working with actual LFP signals (example .mat file available at github)
animalprefixlist = {'KL8','JS17','JS15','JS14','JS12','JS13','JS34','BG1','JS21','ZT2'};
days = 1;
allComod = [];
allComodMean = [];
allPhaseModRip = [];
allPhaseModGam = [];
phaseModRipMI = [];
PhaseModGamMI = [];
numChans = 0;
for d = 1:length(days)
    day = days(d);
    daystring = sprintf('%02d',day);
    for a = 1:length(animalprefixlist)
        allEpSpecs = [];
        animalprefix = animalprefixlist{a};
        dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
        rem = load(sprintf('%s%srem0%d.mat',dir,animalprefix,day));
        load(sprintf('%s%sphasicrembouts0%d.mat',dir,animalprefix,day));
        rem = phasicrem; %constrain analysis to high theta power bouts in REM (putative phasic rem) 
        load(sprintf('%s%sremeps0%d.mat',dir,animalprefix,day));
        epochs = remeps;
        for ep = 1:length(epochs)
            epoch = epochs(ep);

            if epoch <10
                epochstring = ['0',num2str(epoch)];
            else
                epochstring = num2str(epoch);
            end

            remList = [rem{day}{epoch}.starttime rem{day}{epoch}.endtime];

            load(sprintf('%s%stetinfo.mat',dir,animalprefix));

            load(sprintf('%s%sctxripples0%d.mat',dir,animalprefix,day));
            rTets = find(~cellfun(@isempty,ctxripples{day}{epoch}));

            tetsNumRips = [];
            for scan = 1:length(rTets)
                t = rTets(scan);
                numR = length(ctxripples{day}{epoch}{t}.startind);
                tetsNumRips = [tetsNumRips; numR];
            end
            [ripcnt idx] = max(tetsNumRips);

            tetlist = rTets(idx);
%             tetlist = rTets;
            if ~isempty(rem{day}{epoch}.starttime)
                tets = tetinfo{day}{epoch};

                allTetSpecs = [];

                for t = 1:length(tetlist)
                    numChans = numChans + 1;
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
                        AmpFreqVector=150:2:250; %gamma to ripple

                        PhaseFreq_BandWidth=1;
                        AmpFreq_BandWidth=2;

                        %% Define phase bins

                        nbin = 18; % number of phase bins
                        position=zeros(1,nbin); % this variable will get the beginning (not the center) of each phase bin (in rads)
                        winsize = 2*pi/nbin;
                        for j=1:nbin
                            position(j) = -pi+(j-1)*winsize;
                        end

                        %% Filtering and Hilbert transform

                        'CPU filtering'
                        tic
                        Comodulogram=single(zeros(length(PhaseFreqVector),length(AmpFreqVector)));
                        AmpFreqTransformed = zeros(length(AmpFreqVector), data_length);
                        PhaseFreqTransformed = zeros(length(PhaseFreqVector), data_length);

                        for ii=1:length(AmpFreqVector)
                            Af1 = AmpFreqVector(ii);
                            Af2=Af1+AmpFreq_BandWidth;
                            AmpFreq=eegfilt(lfp,srate,Af1,Af2); % filtering
                            AmpFreqTransformed(ii, :) = abs(hilbert(AmpFreq)); % getting the amplitude envelope
                        end

                        for jj=1:length(PhaseFreqVector)
                            Pf1 = PhaseFreqVector(jj);
                            Pf2 = Pf1 + PhaseFreq_BandWidth;
                            PhaseFreq=eegfilt(lfp,srate,Pf1,Pf2); % filtering
                            PhaseFreqTransformed(jj, :) = angle(hilbert(PhaseFreq)); % getting the phase time series
                        end
                        toc

                        %% Compute MI and comodulogram

                        'Comodulation loop'

                        counter1=0;
                        for ii=1:length(PhaseFreqVector)
                            counter1=counter1+1;

                            Pf1 = PhaseFreqVector(ii);
                            Pf2 = Pf1+PhaseFreq_BandWidth;

                            counter2=0;
                            for jj=1:length(AmpFreqVector)
                                counter2=counter2+1;

                                Af1 = AmpFreqVector(jj);
                                Af2 = Af1+AmpFreq_BandWidth;
                                [MI,MeanAmp]=ModIndex_v2(PhaseFreqTransformed(ii, :), AmpFreqTransformed(jj, :), position);
                                Comodulogram(counter1,counter2)=MI;
                            end
                        end
                        allComodMean = [allComodMean; mean(Comodulogram)];
                        allComod = cat(3, allComod, Comodulogram);

                        Pf1 = 2;
                        Pf2 = 12;


                        Af1 = 150;
                        Af2 = 250;

                        [MI,MeanAmp] = ModIndex_v1(lfp,srate,Pf1,Pf2,Af1,Af2,position);
                        allPhaseModRip = [allPhaseModRip; MeanAmp/sum(MeanAmp)];
                        phaseModRipMI = [phaseModRipMI; MI];

                        Af1 = 40;
                        Af2 = 100;

                        [MI,MeanAmp] = ModIndex_v1(lfp,srate,Pf1,Pf2,Af1,Af2,position);
                        allPhaseModGam = [allPhaseModGam; MeanAmp/sum(MeanAmp)];
                        PhaseModGamMI = [PhaseModGamMI; MI];
                    end
                end
            end
        end
    end
end
keyboard
figure
contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,mean(allComod,3)',30,'lines','none')
set(gca,'fontsize',14)
ylabel('Amplitude Frequency (Hz)')
xlabel('Phase Frequency (Hz)')
colormap(inferno)
colorbar

figure
allPhaseModRipMean = mean(allPhaseModRip);
allPhaseModRipSem = std(allPhaseModRip)./sqrt(size(allPhaseModRip,1));
xaxis = [10:20:720];
figure; hold on
pl1 = plot(xaxis,[allPhaseModRipMean allPhaseModRipMean],'-r')
boundedline(xaxis,[allPhaseModRipMean allPhaseModRipMean],...
    [allPhaseModRipSem allPhaseModRipSem],'-r');

allPhaseModGamMean = mean(allPhaseModGam);
allPhaseModGamSem = std(allPhaseModGam)./sqrt(size(allPhaseModGam,1));
pl2 = plot(xaxis,[allPhaseModGamMean allPhaseModGamMean],'-b')
boundedline(xaxis,[allPhaseModGamMean allPhaseModGamMean],...
    [allPhaseModGamSem allPhaseModGamSem],'-b');
legend([pl1 pl2],{'Ripple','Gamma'})

xlim([90 630])
set(gca,'xtick',90:90:720)
xlabel('Phase (Deg)')
ylabel('Amplitude')
title(['MI = ' num2str(MI)])

[p h] = ranksum(phaseModRipMI,PhaseModGamMI);
datacombinedPAC = [phaseModRipMI; PhaseModGamMI];
g1 = repmat({'Ripple'},length(phaseModRipMI),1);
g2 = repmat({'Gamma'},length(PhaseModGamMI),1);
g = [g1;g2];

figure;
h = boxplot(datacombinedPAC,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
title(['Theta phase amplitude coupling-p = ' num2str(p)])
ylabel('Modulation index (MI)')
set(gcf, 'renderer', 'painters')

keyboard
