function out = jds_gammaTriggeredWaveletREM(animalprefixlist, days)
% out = calcripspectrum(index, excludeperiods, eeg,ripples,cellinfo, options)
%  Computes the spectrogram around the middle of each decoded event.
%  Options:
%       appendindex-- Determines whether index is included in output vector
%           Default: 1
%       fpass-- Determines the frequency range for computing spectrum.
%           Default: [2 350]
%       average_trials-- Determines if events are averaged or not.
%           Default: 0
%       spectrum_window-- Determines the sliding window used to compute
%           the event triggered spectrogram. Default: [0.1 0.01]
%       event_window--Determines the size of the window around each
%           triggering event. Default: [0.2 0.2]

%  out is a structure with the following fields
%       S-- This is a MxNxT matrix of the spectrogram for each tetrode
%           M is time relative to triggering event, N is frequency, T is event
%       F-- Frequency vector
%       T-- time relative to triggering event
%       fit-- This is the fit based on the spectrum computed for the entire
%           epoch to normalize S. To reconstruct S without normalization,
%           add log10(frequency)*fit(2)+fit(1) to S.
%       index-- Only if appendindex is set to 1 (default)




%parse the options
Fs = 1500;
win = [0.5 0.5];
S0 = 2*(1/Fs);
DJ = 0.15;
DT = 1/Fs;
N = sum(win)*Fs;
J1 = fix(log(N*DT/S0)/log(2))/DJ;
allAnimSpecsC = [];
gamModC = [];
allAnimSpecsNC = [];
gamModNC = [];
for r = 1:2
    for d = 1:length(days)
        day = days(d);
        daystring = sprintf('%02d',day);
        for a = 1:length(animalprefixlist)
            allEpSpecs = [];
            animalprefix = animalprefixlist{a};
            dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
            %         load(sprintf('%s%sctxrippletime0%d.mat',dir,animalprefix,day));
            load(sprintf('%s%sctxgammatime2_nRip0%d.mat',dir,animalprefix,day));
            load(sprintf('%s%sremeps0%d.mat',dir,animalprefix,day));
%             epochs = 1:2:17;
            epochs = remeps;
            for i = 1:length(epochs)
                epoch = epochs(i);

                if epoch <10
                    epochstring = ['0',num2str(epoch)];
                else
                    epochstring = num2str(epoch);
                end

                load(sprintf('%s%stetinfo.mat',dir,animalprefix));

                gam = ctxgamma{day}{epoch};

                if isempty(gam)
                    continue
                end

                if r == 1
                    gamtimes = [gam.starttimeC gam.endtimeC];
                elseif r == 2
                    gamtimes = [gam.starttimeNC gam.endtimeNC];
                end
                if isempty(gamtimes)
                    continue
                end
                %
                %             iri = diff(gamtimes(:,1));
                %             keepidx = [1;find(iri>=0.5)+1];
                %
                %             if isempty(keepidx)
                %                 continue
                %             end
                %             gamtimes = gamtimes(keepidx,:);

                ctxgammas = load(sprintf('%s%sctxripples0%d.mat',dir,animalprefix,day));
                ctxgammas = ctxgammas.ctxripples;
                rTets = find(~cellfun(@isempty,ctxgammas{day}{epoch}));

                tetsNumGam = [];
                for scan = 1:length(rTets)
                    t = rTets(scan);
                    numR = length(ctxgammas{day}{epoch}{t}.startind);
                    tetsNumGam = [tetsNumGam; numR];
                end
                [ripcnt idx] = max(tetsNumGam);

                tetlist = rTets(idx);
                %             tetlist = rTets;
                if ~isempty(gamtimes)
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
                        triggers = gamtimes(:,1)-starttime;

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
                        for s = 1:length(e)/Fs
                            basetmp = e((idxcnt+1):(Fs+idxcnt));
                            ptmp = wavelet(basetmp,1/Fs,0,DJ,S0,J1, 'MORLET');
                            ptmp = abs(ptmp);
                            idxcnt = idxcnt + Fs;
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
                frequency=round(arrayfun(@(x) 1/x, period));
                gamIdx = find(frequency > 40 & frequency < 100);
                tmpEp = mean(allTetSpecs,3);
                if r == 1
                    gamModC = [gamModC; mean(tmpEp(gamIdx,:))];
                elseif r == 2
                    gamModNC = [gamModNC; mean(tmpEp(gamIdx,:))];
                end
                allEpSpecs = cat(3, allEpSpecs, mean(allTetSpecs,3));
            end
            clear triggers gamtimes
            meanspectAnim = mean(allEpSpecs,3);

            %         imagesc(meanspectAnim)
            %         ax = gca;
            %         ax.FontSize = 14
            %         frequency=round(arrayfun(@(x) 1/x, period));
            %         yLabels = arrayfun(@num2str,frequency,'UniformOutput',false);
            %         yticklabels({yLabels{10},yLabels{20},yLabels{30},yLabels{40},yLabels{50},yLabels{60}})
            %         ylabel('Frequency')
            %         xlim([0 1500])
            %         xticks([0:750:1500])
            %         xticklabels({'-500','0','500'})
            %         xlabel('Time From Ripple Onset (ms)')
            %         colorbar
            %         set(gcf, 'renderer', 'painters')
            %         saveg=0;
            %         if saveg==1
            %             figfile = ['/Volumes/JUSTIN/FigureWorking/REM_Prelim/wavelet/','REMctxrip',animalprefix,'Day',num2str(day),'Tetrode',num2str(tet)];
            %             print('-djpeg', figfile);
            %             print('-depsc', figfile);
            %         end
            %         close
            if r == 1
                allAnimSpecsC = cat(3, allAnimSpecsC, mean(allEpSpecs,3));
            elseif r == 2
                allAnimSpecsNC = cat(3, allAnimSpecsNC, mean(allEpSpecs,3));
            end
        end
    end
end
meanspectC = mean(allAnimSpecsC,3);
meanspectNC = mean(allAnimSpecsNC,3);

figure
imagesc(meanspectC)
ax = gca;
ax.FontSize = 14
frequency=round(arrayfun(@(x) 1/x, period));
yLabels = arrayfun(@num2str,frequency,'UniformOutput',false);
yticklabels({yLabels{10},yLabels{20},yLabels{30},yLabels{40},yLabels{50},yLabels{60}})
ylabel('Frequency')
xlim([0 1500])
xticks([0:750:1500])
xticklabels({'-500','0','500'})
xlabel('Time From Gamma Onset (ms)')
title('Ripple coupled')
colorbar
set(gcf, 'renderer', 'painters')

figure
imagesc(meanspectNC)
ax = gca;
ax.FontSize = 14
frequency=round(arrayfun(@(x) 1/x, period));
yLabels = arrayfun(@num2str,frequency,'UniformOutput',false);
yticklabels({yLabels{10},yLabels{20},yLabels{30},yLabels{40},yLabels{50},yLabels{60}})
ylabel('Frequency')
xlim([0 1500])
xticks([0:750:1500])
xticklabels({'-500','0','500'})
xlabel('Time From Gamma Onset (ms)')
title('Ripple uncoupled')

colorbar
set(gcf, 'renderer', 'painters')

saveg = 0;
if saveg==1
    figfile = ['/Volumes/JUSTIN/FigureWorking/CA1ReactivationSuppression/EPS/Wavelet/PFC/','AllAnimWaveletCA1tetIndPFCripTrig'];
    print('-djpeg', figfile);
    print('-depsc', figfile);
end

sigbins = [];
for r = 1:length(gamModC(1,:))
    [p h] = signrank(gamModC(:,r),gamModNC(:,r));
    if p < 0.05
        sig = 1;
    else
        sig = 0;
    end
    sigbins = [sigbins; [r sig p]];
end

comps = bwconncomp(sigbins(:,2));
list = comps.PixelIdxList;
tbins = 1:1500;

figure; hold on
pl1 = plot(mean(gamModC),'-r');
boundedline([1:1500],nanmean(gamModC),nanstd(gamModC)./sqrt(size(gamModC,1)),'-r');
pl2 = plot(mean(gamModNC),'-k');
boundedline([1:1500],nanmean(gamModNC),nanstd(gamModNC)./sqrt(size(gamModNC,1)),'-k');

for l = 1:length(list)
    startidx = tbins(list{l}(1));
    endidx = tbins(list{l}(end));
    x = [startidx endidx];
    y = ones(1,size(x,2))*0.25;
    plot(x,y,'*','MarkerSize',8,'MarkerEdgeColor','k')
end

legend([pl1 pl2],{'Ripple coupled','Ripple uncoupled'})
xlim([0 1500])
xticks([0:750:1500])
xticklabels({'-500','0','500'})
ylabel('Gamma power (z)')
xlabel('Time from gamma (ms)')
set(gcf, 'renderer', 'painters')

meanGamPowC = mean(gamModC')';
meanGamPowNC = mean(gamModNC')';

[p h] = ranksum(meanGamPowC,meanGamPowNC);
datacombinedGamPower = [meanGamPowC; meanGamPowNC];
g1 = repmat({'Coupled'},length(meanGamPowC),1);
g2 = repmat({'Uncoupled'},length(meanGamPowNC),1);
g = [g1;g2];

figure;
h = boxplot(datacombinedGamPower,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-0.02 0.2])
title(['Gamma power-p = ' num2str(p)])
ylabel('Power (z)')
set(gcf, 'renderer', 'painters')
keyboard;
end

