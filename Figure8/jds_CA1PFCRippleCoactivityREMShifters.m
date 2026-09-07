function jds_CA1PFCRippleCoactivityREMShifters(animalprefixlist)
%JDS_CA1PFCRIPPLECOACTIVITYREMSHIFTERS Ripple cofiring of theta-phase shifters vs nonshifters.
%   jds_CA1PFCRippleCoactivityREMShifters(animalprefixlist) computes z-scored
%   CA1-PFC pairwise cofiring within REM chain cortical ripples separately
%   for CA1 cells classified as REM theta-phase shifters vs nonshifters, and
%   compares the two distributions (rank-sum) and the proportion of pairs
%   with positive cofiring.
%
%   animalprefixlist - cell array of animal prefix strings

day = 1;

shifterRippleCo = [];
nonshifterRippleCo = [];
shiftCofiringCorr = [];
nonshiftCofiringCorr = [];

for a = 1:length(animalprefixlist)

    animalprefix = char(animalprefixlist(a));
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    %Load REM epochs and chain ripples
    load(sprintf('%s%sremeps%02d.mat',dir,animalprefix,day));
    load(sprintf('%s%sctxrippletime_chainREM%02d.mat',dir,animalprefix,day));

    load(sprintf('%s%sCA1remshiftershighthresh%02d.mat',dir,animalprefix,day));

    epochs = remeps;

    for e = 1:length(epochs)

        ep = epochs(e);
        if ep == 1
            continue
        end
        [ctxidx, hpidx] = jds_getallepcells(dir, animalprefix, day, ep, []);

        shifters = shiftList{day}{ep}.cellidx;

        ctxrips = [ctxripple{day}{ep}.starttimeC ctxripple{day}{ep}.endtimeC];
        if (isempty(ctxrips))
            continue
        end
        if (length(ctxrips(:,1)) > 20)
            load(sprintf('%s%sspikes%02d.mat',dir,animalprefix,day));
            numrips = length(ctxrips(:,1));

            CA1idx = hpidx;
            PFCidx = ctxidx;

            Chigh = CA1idx;
            Phigh = PFCidx;
            for c = 1:length(Chigh(:,1))
                cidx1 = find(Chigh(c,1) == shifters(:,1));
                cidx2 = find(Chigh(c,2) == shifters(:,2));
                cidx3 = intersect(cidx1,cidx2);
                if ~isempty(cidx3)
                    tmpcell = [];
                    isShift = shifters(cidx3,3);
                    shiftMag = abs(shifters(cidx3,4) - shifters(cidx3,5));
                    for pp = 1:length(Phigh(:,1))
                        if (~isempty(spikes{day}{ep}{Chigh(c,1)}{Chigh(c,2)})) &&...
                                (~isempty(spikes{day}{ep}{Phigh(pp,1)}{Phigh(pp,2)}))
                            cell1spks = spikes{day}{ep}{Chigh(c,1)}{Chigh(c,2)}.data(:,1);
                            cell2spks = spikes{day}{ep}{Phigh(pp,1)}{Phigh(pp,2)}.data(:,1);

                            spkbins1_c = periodAssign(cell1spks, ctxrips);
                            spkbins2_p = periodAssign(cell2spks, ctxrips);

                            ripnum1_c = unique(spkbins1_c);
                            activeinrip1_c = ripnum1_c(find(ripnum1_c ~= 0));

                            ripnum2_p = unique(spkbins2_p);
                            activeinrip2_p = ripnum2_p(find(ripnum2_p ~= 0));

                            common_hprip = length(find(ismember(activeinrip1_c, activeinrip2_p)));

                            %calculate zscored ripple coactivity

                            nAB_hp = common_hprip;
                            nA_hp = length(activeinrip1_c);
                            nB_hp = length(activeinrip2_p);

                            coact_rip = (nAB_hp - (nA_hp*nB_hp/numrips))/...
                                sqrt(nA_hp*nB_hp*(numrips - nA_hp)*(numrips - nB_hp)/...
                                (numrips^2*(numrips-1)));

                            if isShift == 1
                                shifterRippleCo = [shifterRippleCo; coact_rip];
                            elseif isShift == 0
                                nonshifterRippleCo = [nonshifterRippleCo; coact_rip];
                            end
                            tmpcell = [tmpcell; coact_rip];
                            
                        end
                    end
                    if isShift == 1
                        shiftCofiringCorr = [shiftCofiringCorr; [nanmean(tmpcell) shiftMag]];
                    elseif isShift == 0
                        nonshiftCofiringCorr = [nonshiftCofiringCorr; [nanmean(tmpcell) shiftMag]];
                    end
                end
            end
        end
    end
end

nonshifterRippleCo = nonshifterRippleCo(find(~isnan(nonshifterRippleCo)));
shifterRippleCo = shifterRippleCo(find(~isnan(shifterRippleCo)));

nonshifterRippleCoHigh = nonshifterRippleCo(find(nonshifterRippleCo>0));
shifterRippleCoHigh = shifterRippleCo(find(shifterRippleCo>0));

[p1 h1] = ranksum(nonshifterRippleCo,shifterRippleCo)
datacombinedRemShifter = [nonshifterRippleCo; shifterRippleCo];
g1 = repmat({'Nonshifter'},length(nonshifterRippleCo),1);
g2 = repmat({'Shifter'},length(shifterRippleCo),1);
g = [g1;g2];

figure
histogram(shifterRippleCo,50,'DisplayStyle','stairs');
hold on
x = [nanmedian(shifterRippleCo) nanmedian(shifterRippleCo)];
y = [0 500];
histogram(nonshifterRippleCo,50,'DisplayStyle','stairs');

title(['PFC-CA1 shifters'])
legend({'Shifters','Non-shifters'})
xlabel('Co-firing (z)')
ylabel('Count')
set(gcf, 'renderer', 'painters')

figure;
h = boxplot(datacombinedRemShifter,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
title(['CA1PFC - Ripple Cofiring-p = ' num2str(p1)])
ylabel('Ripple Cofiring (z)')
set(gcf, 'renderer', 'painters')
ylim([-4 5])

shiftProp = length(find(shifterRippleCo > 0))/length(shifterRippleCo);
nonshiftProp = length(find(nonshifterRippleCo > 0))/length(nonshifterRippleCo);

figure
bar([nonshiftProp shiftProp])
yticks([0.4 0.45])
xticklabels({'NonShifters','Shifters'})
ylabel('Proportion CA1PFC cofiring > 0')
set(gcf, 'renderer', 'painters')
ylim([0.35 0.5])
yticks([0.35:0.05:0.5])
title(['ztest for proportions-p = ' num2str(0.00964)])

