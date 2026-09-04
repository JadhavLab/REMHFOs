function jds_normPopulationActivityREM(animalprefixlist,area)
% Plots ripple aligned, normalized population activity. Ripple aligned
% activity is normalized to baseline activity during the specific sleep
% state. ALso calculated PSD of multiunit activity aligned to ripples.
% -------------------------------------------------------------------------

pret=400; postt=400;

day = 1;
binsize = 5;
g1 = gaussian(3, 3);

crip_compiled_data = []; 
crip_compiled_data_ep = []; 
PSDall = [];
anim_data = [];
animPsd = [];
for a = 1:length(animalprefixlist)
    anim_data_tmp = [];
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    load(sprintf('%s%sctxrippletime_REM0%d.mat',dir,animalprefix,day));

    load(sprintf('%s%s_spikematrix_ev_allepochallcell5_%02d.mat',dir,animalprefix,day));
    rem = load(sprintf('%s%srem0%d.mat',dir,animalprefix,day));
    rem = rem.rem;
    load(sprintf('%s%sremeps0%d.mat',dir,animalprefix,day));
    epochs = remeps;
    for e = 1:length(epochs)
        ep_tmp = [];
        ep = epochs(e);

        ctxrip = ctxripple{day}{ep};
        
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

            if length(ctxrip.starttime) > 20
                criptimes = [ctxrip.starttime*1000 ctxrip.endtime*1000];
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
                        crip_compiled_data = [crip_compiled_data; normFR];
                        anim_data_tmp = [anim_data_tmp; normFR];
                        ep_tmp = [ep_tmp; normFR];
                    end
                end
                Params = struct();
                Params.WindowSec = length(ep_tmp(1,:))/200;
                Params.SamplFreq = 200;
                [PSDcurve,FreqList] = REM_PFCrippleMUA_PSD(mean(ep_tmp),Params);
                PSDall = [PSDall; zscore(PSDcurve)];
                animPsd = [animPsd; a];
            end
            crip_compiled_data_ep = [crip_compiled_data_ep; zscore(mean(ep_tmp,1))];
        end
    end
    anim_data{a} = anim_data_tmp;
end

figure; hold on %plot bounded line
pl1 = plot(-pret*binsize:binsize:postt*binsize,nanmean(crip_compiled_data),'-k','LineWidth',1)
boundedline(-pret*binsize:binsize:postt*binsize,nanmean(crip_compiled_data),...
    nanstd(crip_compiled_data)./sqrt(length(crip_compiled_data(:,1))),'-k');

set(gcf, 'renderer', 'painters')
ylabel('Normalized M.U.A.')
xlabel('Time from event (ms)')
xlim([-1000 1000])

figure; hold on
pl2 = plot(FreqList,mean(PSDall,1),'-k','LineWidth',1)
boundedline(FreqList,mean(PSDall,1),std(PSDall,1)./sqrt(length(PSDall(:,1))),'-k');
xlim([4 12])
ylabel('zscore power')
xlabel('Frequency')
title('REM PFC All ripple triggered MUA PSD')
set(gcf, 'renderer', 'painters')

keyboard