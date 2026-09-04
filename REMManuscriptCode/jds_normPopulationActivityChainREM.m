function jds_normPopulationActivityChainREM(animalprefixlist,area)

% --------------- Parameters ---------------
pret=200; postt=200;

%% --------------------------------------------------
%  %Align Single Cell Firing Rate to event
% -------------------------------------------------
day = 1;
binsize = 5;
g1 = gaussian(5, 5);

crip_compiled_data = [];
crip_compiled_data_ep = [];

ncrip_compiled_data = [];
ncrip_compiled_data_ep = [];
allPSDc = [];
allPSDnc = [];
anim_data = [];
for r = 1:2
    for a = 1:length(animalprefixlist)
        anim_data_tmp = [];
        animalprefix = animalprefixlist{a};
        dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
        %     dir = sprintf('/Volumes/JUSTIN/NovelFamiliarNovel/%s_direct/',animalprefix);

        load(sprintf('%s%sctxrippletime_chainctxREM0%d.mat',dir,animalprefix,day));
        load(sprintf('%s%sremeps0%d.mat',dir,animalprefix,day));

        load(sprintf('%s%s_spikematrix_ev_allepochallcell5_%02d.mat',dir,animalprefix,day));
        rem = load(sprintf('%s%sremctx0%d.mat',dir,animalprefix,day));
%         rem = load(sprintf('%s%sswsALL0%d.mat',dir,animalprefix,day));
        rem = rem.rem;
%         rem = rem.sws;
        %     epochs = find(~cellfun('isempty',observation_matrix));
        epochs = remeps;
        for e = 1:length(epochs)
            ep_tmp = [];
            ep = epochs(e);

            if (mod(ep,2) == 0 || ep == 1)
                eps = [ep ep+1];
            else
                eps = [ep ep-1];
            end

            % chain
            

            %% nonchain
            if r == 1
%                 ctxrip.starttime = [ctxripple{day}{ep}.starttimeC ctxripple{day}{ep}.endtimeC];
                ctxrip = ctxripple{day}{ep}.C_sep;
                tmp = [];
                for c = 1:length(ctxrip)
                    tmp = [tmp; [ctxrip{c}(1,1) ctxrip{c}(1,2)]]; %first
                    %                 tmp = [tmp; ctxrip{c}]; %all
                end
                ctxrip = [];
                ctxrip.starttime = tmp;
            else
                ctxrip.starttime = [ctxripple{day}{ep}.starttimeNC ctxripple{day}{ep}.endtimeNC];
            end
            %         ctxrip.starttime = [ctxgamma{day}{ep}.starttimeNC ctxgamma{day}{ep}.endtimeNC];

            if strcmp(area,'CA1')
                datamat = observation_matrix{ep}.hpdata;
            else
                datamat = observation_matrix{ep}.ctxdata;
            end

            timevect = observation_matrix{ep}.timeeeg(1)*1000:5:observation_matrix{ep}.timeeeg(end)*1000;
            if ~isempty(rem{day}{ep}.starttime)
                remlist = [rem{day}{ep}.starttime rem{day}{ep}.endtime];
                remlength = rem{day}{ep}.total_duration; %in seconds
                [~,swsvec] = wb_list2vec(remlist,timevect/1000);
                swsidx = find(swsvec == 1);

                meanpopFR = sum(sum(datamat(:,swsidx)))/remlength;

                swsconfine = datamat(:,swsidx);
                cellid = [];
                for i = 1:length(swsconfine(:,1))
                    numspks = sum(swsconfine(i,:));
                    meanfr = numspks/remlength;
                    if meanfr ~= 0
                        cellid(i) = i;
                    end
                end
                cellid = find(cellid ~= 0)';

                if length(ctxrip.starttime) > 10
                    criptimes = ctxrip.starttime*1000;
                    criptimes = criptimes(:,1) + (criptimes(:,2)-criptimes(:,1))/2;
                    activitysummed = sum(datamat(cellid,:));
                    activitysummed = smoothvect(activitysummed, g1);
                    for t = 1:length(criptimes)
                        currtrig = criptimes(t); %maybe get rid of first trigger?
                        peakidx = lookup(currtrig,timevect);
                        if ((peakidx + postt)<length(activitysummed)) && ((pret-100) > 0)
                            currspks = activitysummed(peakidx-pret:peakidx+postt);
                            binnedFR = currspks./0.005;
                            normFR = binnedFR./meanpopFR;
                            if r == 1
                                crip_compiled_data = [crip_compiled_data; normFR];
                            else
                                ncrip_compiled_data = [ncrip_compiled_data; normFR];
                            end
                            anim_data_tmp = [anim_data_tmp; normFR];
                            ep_tmp = [ep_tmp; normFR];
                        end
                    end
                    Params = struct();
                    Params.WindowSec = size(ep_tmp,2)/200;
                    Params.SamplFreq = 200;
                    [PSDcurve,FreqList] = REM_PFCrippleMUA_PSD(mean(ep_tmp,1),Params);
                    if r == 1
                        allPSDc = [allPSDc; zscore(PSDcurve)];
                    else
                        allPSDnc = [allPSDnc; zscore(PSDcurve)];
                    end
                end
                if r == 1
                    crip_compiled_data_ep = [crip_compiled_data_ep; zscore(mean(ep_tmp,1))];
                else
                    ncrip_compiled_data_ep = [ncrip_compiled_data_ep; zscore(mean(ep_tmp,1))];
                end
            end
        end
        anim_data{a} = anim_data_tmp;
    end
end

figure; hold on %plot bounded line
ax1 = gca;
ax1.FontSize = 14;
pl1 = plot(-pret*binsize:binsize:postt*binsize,nanmean(crip_compiled_data),'-k','LineWidth',1)
boundedline(-pret*binsize:binsize:postt*binsize,nanmean(crip_compiled_data),...
    nanstd(crip_compiled_data)./sqrt(length(crip_compiled_data(:,1))),'-k');
pl2 = plot(-pret*binsize:binsize:postt*binsize,nanmean(ncrip_compiled_data),'-r','LineWidth',1)
boundedline(-pret*binsize:binsize:postt*binsize,nanmean(ncrip_compiled_data),...
    nanstd(ncrip_compiled_data)./sqrt(length(ncrip_compiled_data(:,1))),'-r');

% xlim([-250 250])
set(gcf, 'renderer', 'painters')
ylabel('Normalized M.U.A.')
xlabel('Time from event (ms)')

figure; hold on
for i = 1:length(anim_data)
    plot(zscore(mean(anim_data{i})))
end

figure; hold on
ax1 = gca;
ax1.FontSize = 14;
pl1 = plot(FreqList,mean(allPSDc),'-b','LineWidth',1)
boundedline(FreqList,mean(allPSDc),std(allPSDc)./sqrt(length(allPSDc(:,1))),'-b');
pl1 = plot(FreqList,mean(allPSDnc),'-k','LineWidth',1)
boundedline(FreqList,mean(allPSDnc),std(allPSDnc)./sqrt(length(allPSDnc(:,1))),'-k');

xlim([4 12])
% thetIdx = find(FreqList == 8);
thetIdx = find(FreqList >=6 & FreqList <= 10);
[p h] = ranksum(mean(allPSDc(:,thetIdx)'),mean(allPSDnc(:,thetIdx)'))
title(['MUA PSD Chained/Nonchained 6-10Hz mean p = ' num2str(p)])
set(gcf, 'renderer', 'painters')

keyboard