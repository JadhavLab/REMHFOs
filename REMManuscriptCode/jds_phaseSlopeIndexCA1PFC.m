function jds_phaseSlopeIndexCA1PFC(animalprefixlist)
prePost = 2;
day = 1;
winSz = prePost*1500;
allPSI = [];
allPSI_s = [];
allPSImean = [];
allPSImean_s = [];
allPSIcurve = [];
allPSIcurve_s = [];
allPSIcurvemean = [];
allPSIcurvemean_s = [];

ca1pfclfp = [];
addpath('/Users/justinshin/Desktop/Code/FieldTrip/')
ft_defaults
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
%         if r == 1
%             triggers = ctxripple{day}{ep}.starttimeNC;
%         else
            ctxrip = ctxripple{day}{ep}.C_sep;
            tmp = [];
            for c = 1:length(ctxrip)
                mid = ctxrip{c}(1,1) + ((ctxrip{c}(1,end) - ctxrip{c}(1,1))/2); %middle of chain
                tmp = [tmp; mid]; 
                %                     tmp = [tmp; ctxrip{c}]; %all
            end
            triggers = tmp;
%         end

        if length(triggers) >= 5

            if ep <10
                epochstring = ['0',num2str(ep)];
            else
                epochstring = num2str(ep);
            end
            tetPSI = [];
            tetPSI_s = [];
            tetPSIcurve = [];
            tetPSIcurve_s = [];
            tPair = 0;
            for tet = 1:length(pfctets)
                ptet = pfctets(tet);
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
                    eegca1Shuf = circshift(eegca1,shift);

                    %lookup to get indices
                    react_idx = [];
                    for t = 1:length(triggers)
                        idxtmp = lookup(triggers(t), tvec);
                        react_idx = [react_idx; idxtmp];
                    end
                    react_idxShuf = react_idx(randperm(length(react_idx)));

                    hpEeg = [];
                    ctxEeg = [];
                    data = [];
                    dataS = [];
                    plfpTmp = [];
                    hlfpTmp = [];
                    for t = 1:length(react_idx)
                        %concatenate eeg snippets
                        tmpP = eegpfc(react_idx(t)-winSz:react_idx(t)+(winSz-1));
                        tmpH = eegca1(react_idx(t)-winSz:react_idx(t)+(winSz-1));
                        tmpH_s = eegca1(react_idxShuf(t)-winSz:react_idxShuf(t)+(winSz-1));
                        evTime = tvec(react_idx(t)-winSz:react_idx(t)+(winSz-1));
                        data.trial{t}(1,:) = tmpP;
                        data.trial{t}(2,:) = tmpH;
                        data.time{t} = evTime;

                        dataS.trial{t}(1,:) = tmpP;
                        dataS.trial{t}(2,:) = tmpH_s;
                        dataS.time{t} = evTime;

                        plfpTmp = [plfpTmp; tmpP];
                        hlfpTmp = [hlfpTmp; tmpH];
                    end

                    ca1pfclfp.pfc = mean(plfpTmp,1);
                    ca1pfclfp.ca1 = mean(hlfpTmp,1);
                    ca1pfclfp.pfcTet = ptet;
                    ca1pfclfp.ca1Tet = hptet;
%%
                    data.label = {'a';'b'};
                    cfg = [];
                    cfg.length = 2;
                    data = ft_redefinetrial(cfg, data);

                    % spectral decomposition
                    cfg = [];
                    cfg.method = 'mtmfft';
                    cfg.output = 'fourier';
                    cfg.tapsmofrq = 2;
                    cfg.foilim = [0 100];
                    freq = ft_freqanalysis(cfg, data);

                    %perform PSI calculation per tet pair
                    cfg = [];
                    cfg.method = 'psi';
%                     cfg.method = 'wpli_debiased';
                    cfg.bandwidth = 5;
                    psi = ft_connectivityanalysis(cfg, freq);
                    spctrm = psi.psispctrm;
                    thetIdx = find(psi.freq > 6 & psi.freq < 12); %get only theta freq bins
                    meanPsiPfcLead = mean(spctrm(1,2,thetIdx));
%                     meanPsiCa1Lead = mean(spctrm(2,1,thetIdx));
%                     if tt == 1
%                         cfg = [];
%                         cfg.parameter = 'psispctrm';
%                         ft_connectivityplot(cfg, psi);
%                         keyboard
%                         close all
%                     end
                    
                    tetPSI = [tetPSI; meanPsiPfcLead]; %if CA1 leading, positive
                    allPSI = [allPSI; meanPsiPfcLead];
                    tetPSIcurve = [tetPSIcurve; squeeze(spctrm(1,2,:))'];
                    allPSIcurve = [allPSIcurve; squeeze(spctrm(1,2,:))'];
%%
                    dataS.label = {'a';'b'};
                    cfg = [];
                    cfg.length = 2;
                    dataS = ft_redefinetrial(cfg, dataS);

                    % spectral decomposition
                    cfg = [];
                    cfg.method = 'mtmfft';
                    cfg.output = 'fourier';
                    cfg.tapsmofrq = 2;
                    cfg.foilim = [0 100];
                    freq = ft_freqanalysis(cfg, dataS);

                    %perform PSI calculation per tet pair
                    cfg = [];
                    cfg.method = 'psi';
%                     cfg.method = 'wpli_debiased';
                    cfg.bandwidth = 5;
                    psi_s = ft_connectivityanalysis(cfg, freq);
                    spctrm_s = psi_s.psispctrm;
                    thetIdx_s = find(psi_s.freq > 6 & psi_s.freq < 12); %get only theta freq bins
                    meanPsiPfcLead_s = mean(spctrm_s(1,2,thetIdx_s));
%                     meanPsiCa1Lead = mean(spctrm_s(2,1,thetIdx_s));

                    tetPSI_s = [tetPSI_s; meanPsiPfcLead_s];
                    allPSI_s = [allPSI_s; meanPsiPfcLead_s];
                    tetPSIcurve_s = [tetPSIcurve_s; squeeze(spctrm_s(1,2,:))'];
                    allPSIcurve_s = [allPSIcurve_s; squeeze(spctrm_s(1,2,:))'];
                end
            end
            allPSImean = [allPSImean; mean(tetPSI)];
            allPSImean_s = [allPSImean_s; mean(tetPSI_s)];
            allPSIcurvemean = [allPSIcurvemean; mean(tetPSIcurve)];
            allPSIcurvemean_s = [allPSIcurvemean_s; mean(tetPSIcurve_s)];
        end
    end
end
g1 = gaussian(3,3);
allPSIcurvemeanSm = [];
for i = 1:size(allPSIcurvemean,1)
    tmp = conv(allPSIcurvemean(i,:),g1,'same');
    allPSIcurvemeanSm = [allPSIcurvemeanSm; tmp];
end

figure; boundedline(freq.freq,mean(allPSIcurvemeanSm,1),std(allPSIcurvemeanSm)./sqrt(size(allPSIcurvemeanSm,1)),'-k');
hold on
% boundedline(freq.freq,mean(allPSIcurvemean_s,1),std(allPSIcurvemean_s)./sqrt(size(allPSIcurvemean_s,1)),'-r');
ylim([-0.2 0.2])
x = [0 50];
y = [0 0];
plot(x,y,'--k')
xlim([0 50])
xlabel('Frequency')
ylabel('Phase slope index')
ttest(allPSImean)
[h p] = ttest(allPSImean)
title(['PSI - ttest p=' num2str(p)])
set(gcf, 'renderer', 'painters')

figure
bar(mean(allPSImean))
hold on
% plot(ones(1,length(allPSImean)),allPSImean,'bo')
errorbar(1,mean(allPSImean), std(allPSImean)/sqrt(length(allPSImean)),'-k','LineStyle','none')
set(gcf, 'renderer', 'painters')

x = repmat(1,length(allPSImean),1);  % create the x data needed to overlay the swarmchart on the boxchart
scatter(x(:),allPSImean(:),'filled','MarkerFaceAlpha',0.6','jitter','on','jitterAmount',0.05);
set(gcf, 'renderer', 'painters')

[h p] = ttest2(allPSImean,allPSImean_s)

datacombinedPSI = [allPSImean; allPSImean_s;];
g1 = repmat({'Data'},length(allPSImean),1);
g2 = repmat({'Shuf'},length(allPSImean_s),1);

g = [g1;g2];

figure;
h = boxplot(datacombinedPSI,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-0.02 0.2])
ylabel('Phase slope index')
set(gcf, 'renderer', 'painters')
title(['PSI - ttest p=' num2str(p)])
keyboard