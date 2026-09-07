function jds_CA1PFCRippleCoactivityREMCA1RipAlignedMUA(animalprefixlist)
%JDS_CA1PFCRIPPLECOACTIVITYREMCA1RIPALIGNEDMUA Ripple-aligned CA1 activity by cofiring level.
%   jds_CA1PFCRippleCoactivityREMCA1RipAlignedMUA(animalprefixlist) computes,
%   for each CA1 cell, its mean z-scored cofiring with PFC cells within REM
%   chain cortical ripples, splits cells into high (>0), low (<0), and
%   non-cofiring groups, and plots ripple-aligned smoothed spiking for each
%   group with a rank-sum comparison of peri-ripple (+/-100 ms) activity.
%
%   animalprefixlist - cell array of animal prefix strings

day = 1;

pret = 1;
binsize = 0.03;
postt = 1;

peakbins = find(abs(-pret:binsize:postt)<=0.1);

nstd = round(binsize*2/binsize); 
g1 = gaussian(nstd, 3*nstd+1);

ripTrigMuaLow = [];
ripTrigMuaHigh = [];
ripTrigMuaNon = [];
for a = 1:length(animalprefixlist)

    animalprefix = char(animalprefixlist(a));
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    %Load reactivation strength file for all assemblies and epochs
    load(sprintf('%s%sCA1_RTimeStrengthSleepNewSpk_20_%02d.mat',dir,animalprefix,day));
    CA1_R = RtimeStrength; clear RtimeStrength
    load(sprintf('%s%sPFC_RTimeStrengthSleepNewSpk_20_%02d.mat',dir,animalprefix,day));
    PFC_R = RtimeStrength; clear RtimeStrength
    load(sprintf('%s%sremeps%02d.mat',dir,animalprefix,day));
    load(sprintf('%s%sswsALL%02d.mat',dir,animalprefix,day));
    rem = load(sprintf('%s%srem%02d.mat',dir,animalprefix,day));
    rem = rem.rem;


    %Load ripples
    indrips = load(sprintf('%s%sctxrippletime_chainREM%02d.mat',dir,animalprefix,day));

    epochs = remeps;

    for e = 1:length(epochs)

        ep = epochs(e);
        if ep == 1
            continue
        end
        firstsws = [sws{day}{ep}.starttime(1) sws{day}{ep}.endtime(1)];
        lastsws = [sws{day}{ep}.starttime(end) sws{day}{ep}.endtime(end)];
        allrem = [rem{day}{ep}.starttime rem{day}{ep}.endtime];
        remdur = rem{day}{ep}.total_duration;
        if ((firstsws(2)-firstsws(1) < 30)) || ((lastsws(2)-lastsws(1) < 30))
            continue
        end
        CA1assemblytmp = CA1_R{ep}.reactivationStrength;
        CA1num = length(CA1assemblytmp);
        PFCassemblytmp = PFC_R{ep}.reactivationStrength;
        PFCnum = length(PFCassemblytmp);
        ncrips = [indrips.ctxripple{day}{ep}.starttimeC indrips.ctxripple{day}{ep}.endtimeC];

        %use the second ripple of each chain as the cofiring event set;
        %alignment is still to all chain ripple start times
        ctxrip = indrips.ctxripple{day}{ep}.C_sep;
        tmp = [];
        for c = 1:length(ctxrip)
            tmp = [tmp; [ctxrip{c}(2,1) ctxrip{c}(2,2)]];
        end
        alignRips = ncrips;
        ncrips = tmp;

        if isempty(ncrips)
            continue
        end
        if length(ncrips(:,1)) > 10
            load(sprintf('%s%sspikes%02d.mat',dir,animalprefix,day));
            numncrips = length(ncrips(:,1));

            CA1idx = CA1_R{ep}.cellidx;
            PFCidx = PFC_R{ep}.cellidx;

            Chigh = CA1idx;
            Phigh = PFCidx;
            for c = 1:length(Chigh(:,1))
                tmp = [];
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

                        tmp = [tmp; coact_hprip];
                    end
                end
                tmpCell = [];
                for i=1:length(alignRips(:,1))
                    
                    currrip = alignRips(i,1);

                    % PFC
                    currspks =  cell1spks(find( (cell1spks>=(currrip-pret)) & (cell1spks<=(currrip+postt)) ));
                    currspks = currspks-(currrip);
                    histspks = histc(currspks,[-pret:binsize:postt]);
                   
                    if isempty(histspks) % Only happens if spikeu is empty
                        histspks = zeros(size([-pret:binsize:postt]));
                    end
                    try
                        histspks = smoothvect(histspks, g1);
                    catch
                        disp('Stopped in DFAsj_getripalignspiking');
                    end
                    histspks = histspks(:).';
                    tmpCell = [tmpCell; histspks];
                end
                zTmpCell = zscore(mean(tmpCell));
                if nanmean(tmp) > 0 
                    ripTrigMuaHigh = [ripTrigMuaHigh; zTmpCell];
                elseif nanmean(tmp) < 0
                    ripTrigMuaLow = [ripTrigMuaLow; zTmpCell];
                elseif isnan(nanmean(tmp))
                    ripTrigMuaNon = [ripTrigMuaNon; zTmpCell];
                end             
            end
        end
    end
end

figure; hold on
ax1 = gca;
ax1.FontSize = 14;
pl1 = plot([-pret:binsize:postt],mean(ripTrigMuaLow),'-b','LineWidth',1)
boundedline([-pret:binsize:postt],mean(ripTrigMuaLow),std(ripTrigMuaLow)./sqrt(size(ripTrigMuaLow,1)),'-b');
pl2 = plot([-pret:binsize:postt],mean(ripTrigMuaHigh),'-r','LineWidth',1)
boundedline([-pret:binsize:postt],mean(ripTrigMuaHigh),std(ripTrigMuaHigh)./sqrt(size(ripTrigMuaHigh,1)),'-r');
pl3 = plot([-pret:binsize:postt],mean(ripTrigMuaNon),'-k','LineWidth',1)
boundedline([-pret:binsize:postt],mean(ripTrigMuaNon),std(ripTrigMuaNon)./sqrt(size(ripTrigMuaNon,1)),'-k');
xlim([-0.5 0.5])

cofiring = [ripTrigMuaLow; ripTrigMuaNon];

figure; hold on
ax1 = gca;
ax1.FontSize = 14;
pl2 = plot([-pret:binsize:postt],mean(cofiring),'-r','LineWidth',1)
boundedline([-pret:binsize:postt],mean(cofiring),std(cofiring)./sqrt(size(cofiring,1)),'-r');
pl3 = plot([-pret:binsize:postt],mean(ripTrigMuaHigh),'-k','LineWidth',1)
boundedline([-pret:binsize:postt],mean(ripTrigMuaHigh),std(ripTrigMuaHigh)./sqrt(size(ripTrigMuaHigh,1)),'-k');
xlim([-0.5 0.5])

meanCofiring = mean(cofiring(:,peakbins),2);
meanHigh = mean(ripTrigMuaHigh(:,peakbins),2);
[p1 h] = ranksum(meanCofiring, meanHigh)

title(['CA1 High cofiring vs other-chainrips 100ms surr p=' num2str(p1)])


