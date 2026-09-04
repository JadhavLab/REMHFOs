function jds_rippleThetaPhaseLockingREM(animalprefixlist)

savedata = 0;
day = 1;

daystring = '01';

isoRippleDeg = [];
chainRippleDeg = [];
isoRippleRad = [];
chainRippleRad = [];
pdfIso = [];
pdfChain = [];
isoPPC = [];
chainPPC = [];
for a = 1:length(animalprefixlist)
    shiftList = [];
    animalprefix = animalprefixlist{a};

    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
    %%
    %-----create the event matrix during SWRs-----%
    load(sprintf('%s%sremeps%02d.mat',dir,animalprefix,day));
    rem = load(sprintf('%s%srem%02d.mat',dir,animalprefix,day));
    rem = rem.rem;
    epochs = remeps;
    load(sprintf('%s%stetinfo.mat',dir,animalprefix));
    load(sprintf('%s%scellinfo.mat',dir,animalprefix));
%     load(sprintf('%s%sctxrippletime_chainSWS0%d.mat',dir,animalprefix,day));
%     swsrips = ctxripple; clear ctxripple
    load(sprintf('%s%sctxrippletime_chainREM0%d.mat',dir,animalprefix,day));

    tets = tetinfo{1}{epochs(1)};

    load(sprintf('%s%sripples0%d.mat',dir,animalprefix,day));

    for ep = 1:length(epochs)

        eps = [(epochs(ep) - 1) epochs(ep)];

        epsleep = eps(2);
        if epsleep == 1
            continue
        end

        rTets = find(~cellfun(@isempty,ripples{day}{epsleep}));

        tetsNumRips = [];
        for scan = 1:length(rTets)
            t = rTets(scan);
            numR = length(ripples{day}{epsleep}{t}.startind);
            tetsNumRips = [tetsNumRips; numR];
        end
        [ripcnt idx] = max(tetsNumRips);
        thetatet = rTets(idx);

        if (eps(2) <10) && (isequal(animalprefix, 'ZT2'))
            epochstringrem = ['0',num2str(epsleep)];
        else
            epochstringrem = num2str(epsleep);
        end

        if (thetatet<10)
            tetstring = ['0',num2str(thetatet)];
        else
            tetstring = num2str(thetatet);
        end

%         isotimes = swsrips{day}{epsleep}.starttimeC +...
%             ((swsrips{day}{epsleep}.endtimeC - swsrips{day}{epsleep}.starttimeC)/2);
        isotimes = ctxripple{day}{epsleep}.starttimeNC +...
            ((ctxripple{day}{epsleep}.endtimeNC - ctxripple{day}{epsleep}.starttimeNC)/2);
        chaintimes = ctxripple{day}{epsleep}.starttimeC +...
            ((ctxripple{day}{epsleep}.endtimeC - ctxripple{day}{epsleep}.starttimeC)/2);

        curreegfileRem = [dir,'/EEG/',animalprefix,'thetagnd', daystring,'-',epochstringrem,'-',tetstring];
        load(curreegfileRem);
        thetagnd_rem = thetagnd; clear thetagnd;

        phasedata = thetagnd_rem{day}{epsleep}{thetatet}.data(:,2)*-1;

        t = geteegtimes(thetagnd_rem{day}{epsleep}{thetatet});

        if length(isotimes)~=0
            sphIso = double(phasedata(lookup(isotimes, t)))/10000;
%             sphIso = double(sphIso(logical(goodspikes))) / 10000;  % If no spikes, this will be empty
        else
            sphIso = [];
        end

        if length(chaintimes)~=0
            sphChain = double(phasedata(lookup(chaintimes, t)))/10000;
%             sphChain = double(sphChain(logical(goodspikes2))) / 10000;  % If no spikes, this will be empty
        else
            sphChain = [];
        end

        if (length(sphIso)>20) && (length(sphChain)>20)

            % Rayleigh and Modulation: Originally in lorenlab Functions folder
            stats = rayleigh_test(sphIso); % stats.p and stats.Z, and stats.n
            [mIso, phIso] = modulation(sphIso);
            phdegIso = sphIso*(180/pi);
            % Von Mises Distribution - From Circular Stats toolbox
            [thetahatIso, kappaIso] = circ_vmpar(sphIso); % Better to give raw data. Can also give binned data.
            thetahat_degIso = thetahatIso*(180/pi);

            [praylIso, zraylIso] = circ_rtest(sphIso); % Rayleigh test for non-uniformity of circular data

            stats2 = rayleigh_test(sphChain); % stats.p and stats.Z, and stats.n
            [mChain, phChain] = modulation(sphChain);
            phdegChain = sphChain*(180/pi);
            % Von Mises Distribution - From Circular Stats toolbox
            [thetahatChain, kappaChain] = circ_vmpar(sphChain); % Better to give raw data. Can also give binned data.
            thetahat_degChain = thetahatChain*(180/pi);

            [praylChain, zraylChain] = circ_rtest(sphChain);

            % Make finer polar plot and overlay Von Mises Distribution Fit.
            % Circ Stats Box Von Mises pdf uses a default of 100 angles/nbin
            % -------------------------------------------------------------
            nbins = 50;
            bins = -pi:(2*pi/nbins):pi;
            countIso = histc(sphIso, bins);
            countChain = histc(sphChain, bins);
            isoPPCtmp = [];
            for p = 1:length(phdegIso)
                for p2 = 1:length(phdegIso)
                    if (p ~= p2) && (p2 > p)
                        phaseDiff = phdegIso(p2) - phdegIso(p);
                        %calculates the pairwise phase consistency
                        %(PPC)
                        PPCtmp = cosd(phaseDiff);
                        isoPPCtmp = [isoPPCtmp; PPCtmp];
                    end
                end
            end
            isoPPC = [isoPPC; mean(isoPPCtmp)];

            chainPPCtmp = [];
            for p = 1:length(phdegChain)
                for p2 = 1:length(phdegChain)
                    if (p ~= p2) && (p2 > p)
                        phaseDiff = phdegChain(p2) - phdegChain(p);
                        %calculates the pairwise phase consistency
                        %(PPC)
                        PPCtmp = cosd(phaseDiff);
                        chainPPCtmp = [chainPPCtmp; PPCtmp];
                    end
                end
            end
            chainPPC = [chainPPC; mean(chainPPCtmp)];

            % Make Von Mises Fit
            alpha = linspace(-pi, pi, 50)';
            [pdfIsoTmp] = circ_vmpdf(alpha,thetahatIso,kappaIso);
            [pdfChainTmp] = circ_vmpdf(alpha,thetahatChain,kappaChain);
            isoRippleDeg = [isoRippleDeg; sphIso*(180/pi)];
            chainRippleDeg = [chainRippleDeg; sphChain*(180/pi)];
            isoRippleRad = [isoRippleRad; sphIso];
            chainRippleRad = [chainRippleRad; sphChain];
            pdfIso = [pdfIso; pdfIsoTmp'];
            pdfChain = [pdfChain; pdfChainTmp'];
        end
    end
end

% [pval_remALL, k_remALL, K_remALL] = circ_kuipertest(phaseExc_rem, phaseInh_rem,200,1)
[p1,U2_1] = watsons_U2_approx_p(isoRippleRad*(pi/180), chainRippleRad*(pi/180))

p1 = signrank(isoPPC,chainPPC)

datacombinedPPC = [chainPPC; isoPPC];
g1 = repmat({'Chain PPC'},length(chainPPC),1);
g2 = repmat({'Isolated PPC'},length(isoPPC),1);
g = [g1;g2];

figure; hold on
h = boxplot(datacombinedPPC,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
title(['Theta PPC-p = ' num2str(p1)])
xlim([0.5 2.5])
ylim([-0.1 0.8])
yticks(0:0.2:0.8)
ylabel('Ripple to theta PPC')
set(gcf, 'renderer', 'painters')

for i = 1:length(chainPPC)
    x = [1 2];
    y = [chainPPC(i) isoPPC(i)];
    if chainPPC(i) > isoPPC(i)
        plot(x,y,'-r')
    else
        plot(x,y,'-k')
    end
end

pdfIso = [pdfIso pdfIso(:,(2:end))];
pdfChain = [pdfChain pdfChain(:,(2:end))];

figure;
shadedErrorBar([1:99],mean(pdfIso),nanstd(pdfIso)./sqrt(size(pdfIso,1)),'-b',0); hold on
shadedErrorBar([1:99],mean(pdfChain),nanstd(pdfChain)./sqrt(size(pdfChain,1)),'-r',0);
xlim([1 99]); xticks([1 25 50 75 99]); xticklabels({'-180','0','180','360','540'})
xlabel('Degrees'); ylabel('Probability')
title('Iso-Chain theta phase locking')
set(gcf, 'renderer', 'painters')

wakeCellPh2 = [wakeCellPh; wakeCellPh+360; wakeCellPh+720];

remCellPh2 = [remCellPh; remCellPh+360; remCellPh+720];

bins1=0:45:1080;
[CA1wake t]=histcounts(wakeCellPh2,bins1);
[CA1rem t]=histcounts(remCellPh2,bins1);

figure; hold on
plot(bins1(1:end-1),CA1wake/sum(CA1wake),'r','linewidth',2)
plot(bins1(1:end-1),CA1rem/sum(CA1rem),'b','linewidth',2)
legend({'wake','rem'})
xlim([0 720])
xticks([0:180:720])
title(['Watsons U = ' num2str(U2_1) ' p = ' num2str(p1)])
set(gcf, 'renderer', 'painters')

phaseShiftCell2 = [phaseShiftCell; phaseShiftCell+360;...
    [phaseShiftCell(:,1)+360 phaseShiftCell(:,2)];...
    [phaseShiftCell(:,1) phaseShiftCell(:,2)+360]];

phaseNonshiftCell2 = [phaseNonshiftCell; phaseNonshiftCell+360;...
    [phaseNonshiftCell(:,1)+360 phaseNonshiftCell(:,2)];...
    [phaseNonshiftCell(:,1) phaseNonshiftCell(:,2)+360]];

figure
scatter(phaseShiftCell2(:,2),phaseShiftCell2(:,1),'r.')
hold on
scatter(phaseNonshiftCell2(:,2),phaseNonshiftCell2(:,1),'k.')
xlim([0 720])
xticks([0:180:720])
ylim([0 720])
yticks([0:180:720])
set(gcf, 'renderer', 'painters')
legend({'shifting','nonshifting'})

figure
histogram((shifts),40)
xlabel('Shift (Deg)')
ylabel('Count')
set(gcf, 'renderer', 'painters')
xlim([-180 180])
xticks([-180:90:180])
keyboard

