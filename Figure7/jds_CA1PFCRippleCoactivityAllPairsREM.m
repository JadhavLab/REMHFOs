function jds_CA1PFCRippleCoactivityAllPairsREM(animalprefixlist)
%JDS_CA1PFCRIPPLECOACTIVITYALLPAIRSREM CA1-PFC pairwise ripple cofiring, REM vs NREM.
%   jds_CA1PFCRippleCoactivityAllPairsREM(animalprefixlist) computes z-scored
%   ripple cofiring for all CA1-PFC cell pairs within REM chain cortical
%   ripples and within NREM cortical ripples, compares the two distributions
%   (rank-sum, dip test), and plots histograms, boxplots, and the
%   NREM-vs-REM cofiring correlation.
%
%   animalprefixlist - cell array of animal prefix strings

day = 1;

ncRippleCo = [];
cRippleCo = [];

for a = 1:length(animalprefixlist)

    animalprefix = char(animalprefixlist(a));
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    %Load reactivation strength file for all assemblies and epochs
    load(sprintf('%s%sCA1_RTimeStrengthSleepNewSpk_20_%02d.mat',dir,animalprefix,day));
    CA1_R = RtimeStrength; clear RtimeStrength
    load(sprintf('%s%sPFC_RTimeStrengthSleepNewSpk_20_%02d.mat',dir,animalprefix,day));
    PFC_R = RtimeStrength; clear RtimeStrength
    load(sprintf('%s%sremeps%02d.mat',dir,animalprefix,day));


    %Load ripples
    coordrips = load(sprintf('%s%sctxrippletime_SWS%02d.mat',dir,animalprefix,day));
    indrips = load(sprintf('%s%sctxrippletime_chainREM%02d.mat',dir,animalprefix,day));

    epochs = remeps;

    for e = 1:length(epochs)

        ep = epochs(e);
        if ep == 1
            continue
        end
        CA1assemblytmp = CA1_R{ep}.reactivationStrength;
        CA1num = length(CA1assemblytmp);
        PFCassemblytmp = PFC_R{ep}.reactivationStrength;
        PFCnum = length(PFCassemblytmp);
        ncrips = [indrips.ctxripple{day}{ep}.starttimeC indrips.ctxripple{day}{ep}.endtimeC];
        crips = [coordrips.ctxripple{day}{ep}.starttime coordrips.ctxripple{day}{ep}.endtime];
        if (isempty(crips)) || (isempty(ncrips))
            continue
        end
        if (length(crips(:,1)) > 10) && (length(ncrips(:,1)) > 10)
            load(sprintf('%s%sspikes%02d.mat',dir,animalprefix,day));
            numcrips = length(crips(:,1));
            numncrips = length(ncrips(:,1));

            CA1idx = CA1_R{ep}.cellidx;
            PFCidx = PFC_R{ep}.cellidx;
            Chigh = CA1idx;
            Phigh = PFCidx;
            for c = 1:length(Chigh(:,1))
                for pp = 1:length(Phigh(:,1))

                    if (~isempty(spikes{day}{ep}{Chigh(c,1)}{Chigh(c,2)})) &&...
                            (~isempty(spikes{day}{ep}{Phigh(pp,1)}{Phigh(pp,2)}))
                        cell1spks = spikes{day}{ep}{Chigh(c,1)}{Chigh(c,2)}.data(:,1);
                        cell2spks = spikes{day}{ep}{Phigh(pp,1)}{Phigh(pp,2)}.data(:,1);

                        spkbins1_c = periodAssign(cell1spks, ncrips);
                        spkbins2_p = periodAssign(cell2spks, ncrips);

                        ripnum1_c = unique(spkbins1_c);
                        activeinrip1_c = ripnum1_c(find(ripnum1_c ~= 0));

                        ripnum2_p = unique(spkbins2_p);
                        activeinrip2_p = ripnum2_p(find(ripnum2_p ~= 0));

                        common_hprip = length(find(ismember(activeinrip1_c, activeinrip2_p)));

                        %calculate zscored ripple coactivity

                        nAB_hp = common_hprip;
                        nA_hp = length(activeinrip1_c);
                        nB_hp = length(activeinrip2_p);

                        coact_hprip = (nAB_hp - (nA_hp*nB_hp/numncrips))/...
                            sqrt(nA_hp*nB_hp*(numncrips - nA_hp)*(numncrips - nB_hp)/...
                            (numncrips^2*(numncrips-1)));

                        ncRippleCo = [ncRippleCo; coact_hprip];
                    end
                end
            end
            for c = 1:length(Chigh(:,1))
                for pp = 1:length(Phigh(:,1))
                    if (~isempty(spikes{day}{ep}{Chigh(c,1)}{Chigh(c,2)})) &&...
                            (~isempty(spikes{day}{ep}{Phigh(pp,1)}{Phigh(pp,2)}))
                        cell1spks = spikes{day}{ep}{Chigh(c,1)}{Chigh(c,2)}.data(:,1);
                        cell2spks = spikes{day}{ep}{Phigh(pp,1)}{Phigh(pp,2)}.data(:,1);

                        spkbins1_c = periodAssign(cell1spks, crips);
                        spkbins2_p = periodAssign(cell2spks, crips);

                        ripnum1_c = unique(spkbins1_c);
                        activeinrip1_c = ripnum1_c(find(ripnum1_c ~= 0));

                        ripnum2_p = unique(spkbins2_p);
                        activeinrip2_p = ripnum2_p(find(ripnum2_p ~= 0));

                        common_hprip = length(find(ismember(activeinrip1_c, activeinrip2_p)));

                        %calculate zscored ripple coactivity

                        nAB_hp = common_hprip;
                        nA_hp = length(activeinrip1_c);
                        nB_hp = length(activeinrip2_p);

                        coact_hprip = (nAB_hp - (nA_hp*nB_hp/numcrips))/...
                            sqrt(nA_hp*nB_hp*(numcrips - nA_hp)*(numcrips - nB_hp)/...
                            (numcrips^2*(numcrips-1)));

                        cRippleCo = [cRippleCo; coact_hprip];
                    end
                end
            end
        end
    end
end

[p3 h3] = ranksum(cRippleCo(~isnan(cRippleCo)),ncRippleCo(~isnan(ncRippleCo)))
datacombinedSpatialCorr = [cRippleCo; ncRippleCo];
g1 = repmat({'NREM'},length(cRippleCo),1);
g2 = repmat({'REM'},length(ncRippleCo),1);
g = [g1;g2];

figure
histogram(ncRippleCo,50,'DisplayStyle','stairs');
hold on
x = [nanmedian(ncRippleCo) nanmedian(ncRippleCo)];
y = [0 500];
histogram(cRippleCo,50,'DisplayStyle','stairs');

[p1,dip,xl,xu]=dipTest(ncRippleCo(~isnan(ncRippleCo)))
[p2,dip,xl,xu]=dipTest(cRippleCo(~isnan(cRippleCo)))

title(['PFC-CA1 cofiring - dipTest p'])
legend({['REM chain rips p = ' num2str(p1)],['NREM rips p = ' num2str(p2)]})
xlabel('Co-firing (z)')
ylabel('Count (pairs)')
set(gcf, 'renderer', 'painters')

figure;
h = boxplot(datacombinedSpatialCorr,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
title(['CA1PFC - Ripple Cofiring-p = ' num2str(p3)])
ylabel('Ripple Cofiring (z)')
set(gcf, 'renderer', 'painters')

[r p] = corrcoef([cRippleCo ncRippleCo],'rows','complete')
scatter(cRippleCo,ncRippleCo,'.k')
hold on
lsline
set(gcf, 'renderer', 'painters')
ylabel('REM cofiring')
xlabel('NREM cofiring')
title(['Co-firing correlation r=' num2str(r(1,2)) ' p=' num2str(p(1,2))])
