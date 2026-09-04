function jds_thetaPhaseLockingPFCRipplesShifters(animalprefixlist)

day = 1;

daystring = '01';

%%

PPC_shift = [];
PPC_nonshift = [];
PPC_shiftKappa = [];
PPC_nonshiftKappa = [];
pdf_shift = [];
pdf_nonshift = [];
phaseShift = [];
phaseNonshift = [];
allpdf = [];

for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};

    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
    %%
    %-----create the event matrix during SWRs-----%
    load(sprintf('%s%sspikes%02d.mat',dir,animalprefix,day)); % get spikes
    load(sprintf('%s%sctxrippletime_chainREM%02d.mat',dir,animalprefix,day));
    load(sprintf('%s%sctxripples%02d.mat',dir,animalprefix,day));
    rem = load(sprintf('%s%srem%02d.mat',dir,animalprefix,day));
    rem = rem.rem;
    load(sprintf('%s%sremeps%02d.mat',dir,animalprefix,day));
    load(sprintf('%s%stetinfo.mat',dir,animalprefix));
    load(sprintf('%s%sCA1remshiftershighthresh%02d.mat',dir,animalprefix,day));
    epochs = remeps;
    tets = tetinfo{1}{epochs(1)};
    %
    %     reftet = [];
    %     for t = 1:length(tets)
    %         tmp = tets{t};
    %         if isfield(tmp, 'descrip')
    %             if isequal(tmp.descrip, 'CA1Ref')
    %                 reftet = [reftet; t];
    %             end
    %         end
    %     end
    %
    for ep = 1:length(epochs)

        eps = [(epochs(ep) - 1) epochs(ep)];

        epsleep = eps(2);
        if epsleep == 1
            continue
        end
        eprun = eps(1);

        rTets = find(~cellfun(@isempty,ctxripples{day}{epsleep})); %for phase locking of CA1 cells to PFC

        tetsNumRips = [];
        for scan = 1:length(rTets)
            t = rTets(scan);
            numR = length(ctxripples{day}{epsleep}{t}.startind);
            tetsNumRips = [tetsNumRips; numR];
        end
        [ripcnt idx] = max(tetsNumRips);
        reftet = rTets(idx);

        if (epsleep <10) && (isequal(animalprefix, 'ZT2'))
            epochstring = ['0',num2str(epsleep)];
        else
            epochstring = num2str(epsleep);
        end

        [ctxidx, hpidx] = matchidx_acrossep_singleday(dir, animalprefix, day, eps, []); %(tet, cell)
        cellnum = size(hpidx,1);

        shifters = shiftList{day}{epsleep}.cellidx;

        ctxrip = ctxripple{day}{epsleep}.C_sep;
        tmp = [];
        for c = 1:length(ctxrip)
            midtmp = ctxrip{c}(1,1) + ((ctxrip{c}(end,2) - ctxrip{c}(1,1))/2);
            stend = [(midtmp-2) (midtmp+2)];
            tmp = [tmp; stend]; %middle of chain extended +-2 sec
        end

        thetalist = tmp;
        if isempty(shifters)
            continue
        end

        for cellcount = 1:cellnum %get spikes for each cell
            index = [day,epsleep,hpidx(cellcount,[1 2])];
            cidx1 = find(hpidx(cellcount,1) == shifters(:,1));
            cidx2 = find(hpidx(cellcount,2) == shifters(:,2));
            cidx3 = intersect(cidx1,cidx2);
            if ~isempty(cidx3)
                isShift = shifters(cidx3,3);
                if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)})
                    spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
                else
                    spiketimes = [];
                end

                goodspikes = isExcluded(spiketimes, thetalist);

                if (reftet<10)
                    reftetstring = ['0',num2str(reftet)];
                else
                    reftetstring = num2str(reftet);
                end

                curreegfile = [dir,'/EEG/',animalprefix,'thetagnd', daystring,'-',epochstring,'-',reftetstring];
                load(curreegfile);
                thetagnd_sleep = thetagnd; clear thetagnd;

                phasedata = thetagnd_sleep{day}{epsleep}{reftet}.data(:,2);

                t = geteegtimes(thetagnd_sleep{day}{epsleep}{reftet});

                if ~isempty(goodspikes)
                    sph = phasedata(lookup(spiketimes, t));
                    sph = double(sph(logical(goodspikes))) / 10000;  % If no spikes, this will be empty
                else
                    sph = [];
                end

                if length(sph)>10

                    % Rayleigh and Modulation: Originally in lorenlab Functions folder
                    stats = rayleigh_test(sph); % stats.p and stats.Z, and stats.n
                    [m, ph] = modulation(sph);
                    phdeg = ph*(180/pi);
                    % Von Mises Distribution - From Circular Stats toolbox
                    [thetahat, kappa] = circ_vmpar(sph); % Better to give raw data. Can also give binned data.
                    thetahat_deg = thetahat*(180/pi);
                    
                    alpha = linspace(-pi, pi, 50)';
                    allphases = sph*(180/pi);

                    [pdf] = circ_vmpdf(alpha,thetahat,kappa);

                    [prayl, zrayl] = circ_rtest(sph); % Rayleigh test for non-uniformity of circular data
                    cellPPC = [];
%                     if prayl < 0.05
                        for p = 1:length(allphases)
                            for p2 = 1:length(allphases)
                                if (p ~= p2) && (p2 > p)
                                    phaseDiff = allphases(p2) - allphases(p);
                                    %calculates the pairwise phase consistency
                                    %(PPC)
                                    PPCtmp = cosd(phaseDiff);
                                    cellPPC = [cellPPC; PPCtmp];
                                end
                            end
                        end
                        if isShift == 1
                            PPC_shift = [PPC_shift; mean(cellPPC)];
                            PPC_shiftKappa = [PPC_shiftKappa; kappa];
                            pdf_shift = [pdf_shift; pdf'];
                            phaseShift = [phaseShift; phdeg];
                        elseif isShift == 0
                            PPC_nonshift = [PPC_nonshift; mean(cellPPC)];
                            PPC_nonshiftKappa = [PPC_nonshiftKappa; kappa];
                            pdf_nonshift = [pdf_nonshift; pdf'];
                            phaseNonshift = [phaseNonshift; phdeg];
                        end
                        allpdf = [allpdf; pdf'];
%                     end
                end
            end
        end
    end
end

[p,U2] = watsons_U2_approx_p(phaseShift, phaseNonshift)

meanPCCexc = mean(PPC_shift)
semPCCexc = std(PPC_shift)/sqrt(length(PPC_shift))
medianPCCexc = median(PPC_shift)

meanPCCinh = mean(PPC_nonshift)
semPCCinh = std(PPC_nonshift)/sqrt(length(PPC_nonshift))
medianPCCinh = median(PPC_nonshift)

% keyboard;

datameansPPC = [mean(PPC_shift) mean(PPC_nonshift)];
datasemsPPC = [(std(PPC_shift)/sqrt(length(PPC_shift)))...
    (std(PPC_nonshift)/sqrt(length(PPC_nonshift)))];

[p h] = ranksum(PPC_shift,PPC_nonshift)

datacombinedPPC = [PPC_shift; PPC_nonshift];
g1 = repmat({'Shift'},length(PPC_shift),1);
g2 = repmat({'Non'},length(PPC_nonshift),1);
g = [g1;g2];

% boxplot(datacombinedCoact,g);
figure
h = boxplot(datacombinedPPC,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
ylim([-0.15 0.3])
title(['Pairwise phase consistency-p = ' num2str(p)])
set(gcf, 'renderer', 'painters')

[p2 h2] = ranksum(PPC_shiftKappa,PPC_nonshiftKappa)

datacombinedKappa = [PPC_shiftKappa; PPC_nonshiftKappa];
g1 = repmat({'Shift'},length(PPC_shiftKappa),1);
g2 = repmat({'Non'},length(PPC_nonshiftKappa),1);
g = [g1;g2];

% boxplot(datacombinedCoact,g);
figure
h = boxplot(datacombinedKappa,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
ylim([-0.2 1.6])
title(['Kappa -p = ' num2str(p2)])
set(gcf, 'renderer', 'painters')

kappaMn = [mean(PPC_shiftKappa) mean(PPC_nonshiftKappa)];
kappaSems =  [(std(PPC_shiftKappa)/sqrt(length(PPC_shiftKappa)))...
    (std(PPC_nonshiftKappa)/sqrt(length(PPC_nonshiftKappa)))];

figure
bar([1:2],datameansPPC,'k')
hold on
er = errorbar([1:2],datameansPPC,datasemsPPC);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('Pairwise Phase Consistency (PPC)')
title('CA1 Cell Phase Locking Strength (PPC)')
xticklabels({'Shift','Non'}); xtickangle(45)

figure
bar([1:2],kappaMn,'k')
hold on
er = errorbar([1:2],kappaMn,kappaSems);
er.Color = [0 0 0]; er.LineWidth = 2; er.LineStyle = 'none';
ylabel('Kappa')
title('CA1 Cell Phase Locking Strength (Kappa)')
xticklabels({'Shift','Non'}); xtickangle(45)

figure;
datacombined = [PPC_shift; PPC_nonshift];
g1 = repmat({'Shift'},length(PPC_shift),1);
g2 = repmat({'Non'},length(PPC_nonshift),1);
g = [g1;g2];

boxplot(datacombined,g,'PlotStyle','compact');

pdfShift = [pdf_shift pdf_shift(:,(2:end))];
pdfNonshift = [pdf_nonshift pdf_nonshift(:,(2:end))];
allpdf = [allpdf allpdf(:,(2:end))];

[p,U2] = watsons_U2_approx_p(phaseShift, phaseNonshift)
figure; 
shadedErrorBar([1:99],mean(allpdf),nanstd(allpdf)./sqrt(size(allpdf,1)),'-k',0); 
title(['shift-non hp REMPFCriptheta p=' num2str(p)])

figure;
shadedErrorBar([1:99],mean(pdfShift),nanstd(pdfShift)./sqrt(size(pdfShift,1)),'-b',0); hold on
shadedErrorBar([1:99],mean(pdfNonshift),nanstd(pdfNonshift)./sqrt(size(pdfNonshift,1)),'-r',0);
xlim([1 99]); xticks([1 25 50 75 99]); xticklabels({'-180','0','180','360','540'})
xlabel('Degrees'); ylabel('Probability')
title('CA1 REM peri-ripple PFC-Theta Phase Locking PDF')

keyboard

