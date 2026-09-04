%jds_REM_phaseAmpCoupling_shuf

%% Working with actual LFP signals (example .mat file available at github)
animalprefixlist = {'KL8','JS17','JS15','JS14','JS12','JS13','JS34','BG1','JS21','ZT2'};
days = 1;
allComod = [];
allPhaseModRip = [];
allPhaseModGam = [];
phaseModRipMI = [];
PhaseModGamMI = [];
allComod_s = [];
allPhaseModRip_s = [];
allPhaseModGam_s = [];
phaseModRipMI_s = [];
PhaseModGamMI_s = [];
for d = 1:length(days)
    day = days(d);
    daystring = sprintf('%02d',day);
    for a = 1:length(animalprefixlist)
        allEpSpecs = [];
        animalprefix = animalprefixlist{a};
        dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
%         rem = load(sprintf('%s%srem0%d.mat',dir,animalprefix,day));
        load(sprintf('%s%sphasicrembouts0%d.mat',dir,animalprefix,day));
%         rem = rem.rem;
        rem = phasicrem;
        load(sprintf('%s%sremeps0%d.mat',dir,animalprefix,day));
        epochs = remeps;
        %         epochs = 1:length(ctxripple{day});
        for ep = 1:length(epochs)
            epoch = epochs(ep);

            if epoch <10
                epochstring = ['0',num2str(epoch)];
            else
                epochstring = num2str(epoch);
            end

%             remList = [rem{day}{epoch}.starttime rem{day}{epoch}.endtime];
            remList = [rem{day}{epoch}.starttime-5 rem{day}{epoch}.endtime+5];

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

                        PhaseFreqVector=2:2:12; %theta
                        AmpFreqVector=40:2:250; %gamma to ripple

                        PhaseFreq_BandWidth=2;
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
                        
                        ampRand = randi(length(AmpFreqTransformed),[1 1]);
                        phaseRand = randi(length(PhaseFreqTransformed),[1 1]);
                        AmpFreqTransformed_shuf = circshift(AmpFreqTransformed, [0 ampRand]);
                        PhaseFreqTransformed_shuf = circshift(PhaseFreqTransformed, [0 phaseRand]);

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

                        
                        %% Compute MI and comodulogram for shuf

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
                                [MI,MeanAmp]=ModIndex_v2(PhaseFreqTransformed_shuf(ii, :), AmpFreqTransformed_shuf(jj, :), position);
                                Comodulogram_s(counter1,counter2)=MI;
                            end
                        end
                        allComod_s = cat(3, allComod_s, Comodulogram_s);

                        Pf1 = 2;
                        Pf2 = 12;


                        Af1 = 150;
                        Af2 = 250;

                        [MI,MeanAmp] = ModIndex_v1_shuf(lfp,srate,Pf1,Pf2,Af1,Af2,position);
                        allPhaseModRip_s = [allPhaseModRip_s; MeanAmp/sum(MeanAmp)];
                        phaseModRipMI_s = [phaseModRipMI_s; MI];

                        Af1 = 40;
                        Af2 = 100;

                        [MI,MeanAmp] = ModIndex_v1_shuf(lfp,srate,Pf1,Pf2,Af1,Af2,position);
                        allPhaseModGam_s = [allPhaseModGam_s; MeanAmp/sum(MeanAmp)];
                        PhaseModGamMI_s = [PhaseModGamMI_s; MI];
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
title('Data')
colormap(inferno)
colorbar

figure
contourf(PhaseFreqVector+PhaseFreq_BandWidth/2,AmpFreqVector+AmpFreq_BandWidth/2,mean(allComod_s,3)',30,'lines','none')
set(gca,'fontsize',14)
ylabel('Amplitude Frequency (Hz)')
xlabel('Phase Frequency (Hz)')
title('Shuffle')
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

xlim([0 720])
set(gca,'xtick',0:360:720)
xlabel('Phase (Deg)')
ylabel('Amplitude')
title(['MI = ' num2str(MI)])

[p h] = ranksum(PhaseModGamMI,PhaseModGamMI_s);
datacombinedPAC_gamma = [PhaseModGamMI; PhaseModGamMI_s];
g1 = repmat({'Gamma'},length(PhaseModGamMI),1);
g2 = repmat({'Gamma-shuffle'},length(PhaseModGamMI_s),1);
g = [g1;g2];

figure;
h = boxplot(datacombinedPAC_gamma,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-0.02 0.2])
title(['Theta phase amplitude coupling-p = ' num2str(p)])
ylabel('Modulation index (MI) - Gamma')
set(gcf, 'renderer', 'painters')

%%

[p2 h2] = ranksum(phaseModRipMI,phaseModRipMI_s);
datacombinedPAC_ripple = [phaseModRipMI; phaseModRipMI_s];
g1 = repmat({'Ripple'},length(phaseModRipMI),1);
g2 = repmat({'Ripple-shuffle'},length(phaseModRipMI_s),1);
g = [g1;g2];

figure;
h = boxplot(datacombinedPAC_ripple,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-0.02 0.2])
title(['Theta phase amplitude coupling-p = ' num2str(p2)])
ylabel('Modulation index (MI) - Ripple')
set(gcf, 'renderer', 'painters')

keyboard
