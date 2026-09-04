function jds_thetaPhaseShiftersREM(animalprefixlist)

savedata = 0;
day = 1;

daystring = '01';

%%

phaseNonshiftCell = [];
phaseShiftCell = [];
wakeCell = [];
remCell = [];
wakeCellPh = [];
remCellPh = [];
shifts = [];
for a = 1:length(animalprefixlist)
    shiftList = [];
    animalprefix = animalprefixlist{a};
    
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
    %%
    %-----create the event matrix during SWRs-----%
    load(sprintf('%s%sspikes%02d.mat',dir,animalprefix,day)); % get spikes
    load(sprintf('%s%sthetatime%02d.mat',dir,animalprefix,day));
    load(sprintf('%s%sremeps%02d.mat',dir,animalprefix,day));
    rem = load(sprintf('%s%srem%02d.mat',dir,animalprefix,day));
    rem = rem.rem;
    epochs = remeps;
    load(sprintf('%s%stetinfo.mat',dir,animalprefix));
    load(sprintf('%s%scellinfo.mat',dir,animalprefix));
    
    tets = tetinfo{1}{epochs(1)}; 

    load(sprintf('%s%sripples0%d.mat',dir,animalprefix,day));
    
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
    
    for ep = 1:length(epochs)

        eps = [(epochs(ep) - 1) epochs(ep)];

        epsleep = eps(2);
        if epsleep == 1
            continue
        end
        eprun = eps(1);

        rTets = find(~cellfun(@isempty,ripples{day}{epsleep}));

        tetsNumRips = [];
        for scan = 1:length(rTets)
            t = rTets(scan);
            numR = length(ripples{day}{epsleep}{t}.startind);
            tetsNumRips = [tetsNumRips; numR];
        end
        [ripcnt idx] = max(tetsNumRips);
        reftet = rTets(idx);

        if (eprun <10) && (isequal(animalprefix, 'ZT2'))
            epochstring = ['0',num2str(eprun)];
        else
            epochstring = num2str(eprun);
        end
        
        if (eps(2) <10) && (isequal(animalprefix, 'ZT2'))
            epochstringrem = ['0',num2str(eps(2))];
        else
            epochstringrem = num2str(eps(2));
        end
        
        thetalist = thetatime{day}{eprun};
        thetalist = [thetalist.starttime thetalist.endtime];
        remlist = rem{day}{eps(2)};
        remlist = [remlist.starttime remlist.endtime];

        [ctxidx, hpidx] = matchidx_acrossep_singleday(dir, animalprefix, day, eps, []); %(tet, cell)
        cellnum = size(hpidx,1);
        idx = [];
        for cellcount = 1:cellnum %get spikes for each cell
            index = [day,eprun,hpidx(cellcount,[1 2])];
            index2 = [day,epsleep,hpidx(cellcount,[1 2])];
            if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)})
                spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
            else
                spiketimes = [];
            end
            if ~isempty(spikes{index2(1)}{index2(2)}{index2(3)}{index2(4)})
                spiketimes2 = spikes{index2(1)}{index2(2)}{index2(3)}{index2(4)}.data(:,1);
            else
                spiketimes2 = [];
            end
            
            goodspikes = isExcluded(spiketimes, thetalist);
            goodspikes2 = isExcluded(spiketimes2, remlist);
                        
            if (reftet<10)
                reftetstring = ['0',num2str(reftet)];
            else
                reftetstring = num2str(reftet);
            end
            
            curreegfile = [dir,'/EEG/',animalprefix,'thetagnd', daystring,'-',epochstring,'-',reftetstring];
            load(curreegfile);
            thetagnd_run = thetagnd; clear thetagnd;
            
            curreegfileRem = [dir,'/EEG/',animalprefix,'thetagnd', daystring,'-',epochstringrem,'-',reftetstring];
            load(curreegfileRem);
            thetagnd_rem = thetagnd; clear thetagnd;
            
            phasedata = thetagnd_run{day}{eprun}{reftet}.data(:,2)*-1;
            
            phasedata2 = thetagnd_rem{day}{epsleep}{reftet}.data(:,2)*-1;
            
            t = geteegtimes(thetagnd_run{day}{eprun}{reftet});
            t2 = geteegtimes(thetagnd_rem{day}{epsleep}{reftet});
            
            if length(goodspikes)~=0
                sph = phasedata(lookup(spiketimes, t));
                sph = double(sph(logical(goodspikes))) / 10000;  % If no spikes, this will be empty
            else
                sph = [];
            end
            
            if length(goodspikes2)~=0
                sph2 = phasedata2(lookup(spiketimes2, t2));
                sph2 = double(sph2(logical(goodspikes2))) / 10000;  % If no spikes, this will be empty
            else
                sph2 = [];
            end
            
            if (length(sph)>10) && (length(sph2)>10)
                
                % Rayleigh and Modulation: Originally in lorenlab Functions folder
                stats = rayleigh_test(sph); % stats.p and stats.Z, and stats.n
                [m, ph] = modulation(sph);
                phdeg = ph*(180/pi);
                % Von Mises Distribution - From Circular Stats toolbox
                [thetahat, kappa] = circ_vmpar(sph); % Better to give raw data. Can also give binned data.
                thetahat_deg = thetahat*(180/pi);
                
                
                [prayl, zrayl] = circ_rtest(sph); % Rayleigh test for non-uniformity of circular data
                
                stats2 = rayleigh_test(sph2); % stats.p and stats.Z, and stats.n
                [m2, ph2] = modulation(sph2);
                phdeg2 = ph2*(180/pi);
                % Von Mises Distribution - From Circular Stats toolbox
                [thetahat2, kappa2] = circ_vmpar(sph2); % Better to give raw data. Can also give binned data.
                thetahat_deg2 = thetahat2*(180/pi);
                
                [prayl2, zrayl2] = circ_rtest(sph2); 
                
                % Make finer polar plot and overlay Von Mises Distribution Fit.
                % Circ Stats Box Von Mises pdf uses a default of 100 angles/nbin
                % -------------------------------------------------------------
                nbins = 50;
                bins = -pi:(2*pi/nbins):pi;
                count = histc(sph, bins);
                count2 = histc(sph2, bins);
                
                % Make Von Mises Fit
                alpha = linspace(-pi, pi, 50)';
%                 alpha = linspace(-pi, pi, 100)';
                [pdf] = circ_vmpdf(alpha,thetahat,kappa);
                [pdf2] = circ_vmpdf(alpha,thetahat2,kappa2);

                if (prayl < 0.05) && (prayl2 < 0.05) %get proportion shifters here
%                     if (phdeg2 < 140) || (phdeg2 > 320) %shifters
                    phDiff= phdeg2-phdeg;
                    normDiff = mod(phDiff,360);
                    if normDiff > 180
                        normDiff = normDiff - 360;
                    elseif normDiff < -180
                        normDiff = normDiff + 360;
                    end
                    if abs(normDiff) > 90
                        phaseShiftCell = [phaseShiftCell; [phdeg phdeg2]];
                        shift = 1;
                    else
                        phaseNonshiftCell = [phaseNonshiftCell; [phdeg phdeg2]];
                        shift = 0;
                    end
                    idx = [idx; [hpidx(cellcount,[1 2]) shift phdeg phdeg2 abs(normDiff)]];
                    shifts = [shifts; normDiff];
                    wakeCell = [wakeCell; pdf'];
                    wakeCellPh = [wakeCellPh; phdeg];
                    remCell = [remCell; pdf2'];
                    remCellPh = [remCellPh; phdeg2];
                end
            end
        end
        shiftList{day}{epsleep}.cellidx = idx;
    end
    if savedata == 1
        save(sprintf('%s%sCA1remshiftershighthreshPFCtheta%02d.mat',dir,animalprefix,day), 'shiftList');
    end
end

% [pval_remALL, k_remALL, K_remALL] = circ_kuipertest(phaseExc_rem, phaseInh_rem,200,1)
[p1,U2_1] = watsons_U2_approx_p(wakeCellPh*(pi/180), remCellPh*(pi/180))

pdfwake = [wakeCell wakeCell(:,(2:end))];
pdfrem = [remCell remCell(:,(2:end))];

figure;
shadedErrorBar([1:99],mean(pdfwake),nanstd(pdfwake)./sqrt(size(pdfwake,1)),'-b',0); hold on
shadedErrorBar([1:99],mean(pdfrem),nanstd(pdfrem)./sqrt(size(pdfrem,1)),'-r',0);
xlim([1 99]); xticks([1 25 50 75 99]); xticklabels({'-180','0','180','360','540'})
xlabel('Degrees'); ylabel('Probability')
title('REM-WAKE theta phase locking')
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

