function out = jds_pEMG(animalprefixlist)

day = 1;

fs = 1500;  % Sampling frequency (Hz)
win_size = 0.5;  % Window size in seconds
samples_per_win = round(win_size * fs);
smoothing_width = 1;
% Frequency band of interest
f_low = 300;
f_high = 600;
shoulder_low = 275;
shoulder_high = 625;
kernel = gaussian(2, 10);

% Create Butterworth bandpass filter (with shoulders)
[b_filt, a_filt] = butter(4, [shoulder_low, shoulder_high] / (fs/2), 'bandpass');

allWakeEmg = [];
allRemEmg = [];

allWakeEmgMn = [];
allRemEmgMn = [];

rem2wakeEMG = [];
rem2wakenremEMG = [];

for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
    load(sprintf('%s%stetinfo.mat',dir,animalprefix));
    load(sprintf('%s%swaking0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%sswsALL0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%sremeps0%d.mat',dir,animalprefix,day));
    rem = load(sprintf('%s%srem0%d.mat',dir,animalprefix,day));
    rem = rem.rem;
    tets = tetinfo{1}{1};
    epochs = remeps;
    pfctets = [];
    reftet = [];
    for t = 1:length(tets)
        tmp = tets{t};
        if isfield(tmp, 'descrip')
            if isequal(tmp.descrip, 'riptet')
                pfctets = [pfctets; t];
            elseif isequal(tmp.descrip, 'CA1Ref')
                reftet = [reftet; t];
            end
        end
    end

    for e = 1:length(epochs)
        pemgtmp = [];
        ep = epochs(e);

        wakelist = [waking{day}{ep}.starttime waking{day}{ep}.endtime];
        remlist = [rem{day}{ep}.starttime rem{day}{ep}.endtime];
        swslist = [sws{day}{ep}.starttime sws{day}{ep}.endtime];
        allbouts = [[wakelist zeros(size(wakelist,1),1)];...
            [swslist ones(size(swslist,1),1)];... 
            [remlist ones(size(remlist,1),1)*2]];

        allbouts = sortrows(allbouts,1);


        if ep <10
            epochstring = ['0',num2str(ep)];
        else
            epochstring = num2str(ep);
        end
        for t = 1:length(pfctets)
            ptet = pfctets(t);
            if (ptet<10)
                ptetstring = ['0',num2str(ptet)];
            else
                ptetstring = num2str(ptet);
            end
            tmptets = pfctets;
            tmptets(find(pfctets == ptet)) = [];

            rand = randi(length(tmptets));
            %                 ref = tmptets(rand);
            ref = reftet;

            if (ref<10)
                reftetstring = ['0',num2str(ref)];
            else
                reftetstring = num2str(ref);
            end

            %----- Get the eeg data and time -----%
            currpfcfile = [dir,'/EEG/',animalprefix,'eeg', '01' ,'-',epochstring,'-',ptetstring];
            load(currpfcfile);
            p_eeg = eeg;
            currreffile = [dir,'/EEG/',animalprefix,'eeg', '01' ,'-',epochstring,'-',reftetstring];
            load(currreffile);
            ref_eeg = eeg;

            tvec = geteegtimes(p_eeg{day}{ep}{ptet});
            eegpfc = p_eeg{day}{ep}{ptet}.data;
            eegref = ref_eeg{day}{ep}{ref}.data;

            filtered_pfc = filtfilt(b_filt, a_filt, eegpfc)';
            filtered_ref = filtfilt(b_filt, a_filt, eegref)';

            num_bins = floor(size(filtered_pfc,2) / samples_per_win);
            emg_scores = zeros(1, num_bins);
            for b = 1:num_bins
                idx_start = (b-1)*samples_per_win + 1;
                idx_end = b*samples_per_win;
                bin_pfc = filtered_pfc(idx_start:idx_end);
                bin_ref = filtered_ref(idx_start:idx_end);

                r = corr(bin_pfc', bin_ref');

                % Mean of pairwise correlations = EMG score
                emg_scores(b) = r;
            end
            pemgtmp = [pemgtmp; emg_scores];
        end
        time_axis = tvec(1)+ (0:num_bins-1) * win_size;
        epEmg = normalize(mean(pemgtmp,1),'range');
        epEmgSm = smoothvect(epEmg,kernel);

        for w = 1:size(wakelist)
            tmpwk = wakelist(w,:);
            thiswk = logical(isExcluded(time_axis,tmpwk));
            emgtmpMn = mean(epEmgSm(thiswk));
            emgtmp = epEmgSm(thiswk)';
            allWakeEmg = [allWakeEmg; emgtmp];
            allWakeEmgMn = [allWakeEmgMn; emgtmpMn];
        end

        for w = 1:size(remlist)
            tmprem = remlist(w,:);
            thisrem = logical(isExcluded(time_axis,tmprem));
            emgtmpMn = mean(epEmgSm(thisrem));
            emgtmp = epEmgSm(thisrem)';
            allRemEmg = [allRemEmg; emgtmp];
            allRemEmgMn = [allRemEmgMn; emgtmpMn];
        end

        if ~isempty(remlist)
            for rr = 1:length(remlist(:,1))
                idx = find(remlist(rr,2) == allbouts(:,2));
                endidx = lookup(remlist(rr,2),time_axis);
                endEMGvec = epEmgSm(endidx-20:endidx+10);
                if allbouts(idx+1,3) == 0
                    rem2wakeEMG = [rem2wakeEMG; endEMGvec];
                end
                rem2wakenremEMG = [rem2wakenremEMG; endEMGvec];
            end
        end
    end
end
rem2wakeEMGmean = mean(rem2wakeEMG);
rem2wakeEMGsem = std(rem2wakeEMG)./sqrt(size(rem2wakeEMG,1));
figure; hold on
ax1 = gca;
ax1.FontSize = 14;
pl1 = plot(-20:10,mean(rem2wakeEMG),'-k','LineWidth',1)
boundedline(-20:10,rem2wakeEMGmean,rem2wakeEMGsem,'-k');
xticks([-20 0 10])
xticklabels({'-10', '0', '5'})
title('REM to Wake iEMG')

rem2wakenremEMGmean = mean(rem2wakenremEMG);
rem2wakenremEMGsem = std(rem2wakenremEMG)./sqrt(size(rem2wakenremEMG,1));
figure; hold on
ax1 = gca;
ax1.FontSize = 14;
pl1 = plot(-20:10,mean(rem2wakenremEMG),'-k','LineWidth',1)
boundedline(-20:10,rem2wakenremEMGmean,rem2wakenremEMGsem,'-k');
xticks([-20 0 10])
xticklabels({'-10', '0', '5'})
title('REM to Wake/nrem iEMG')

[p h] = ranksum(allRemEmgMn,allWakeEmgMn);
datacombinedemg = [allRemEmgMn; allWakeEmgMn];
g1 = repmat({'REM'},length(allRemEmgMn),1);
g2 = repmat({'Wake'},length(allWakeEmgMn),1);
g = [g1;g2];

figure;
h = boxplot(datacombinedemg,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
% ylim([-0.02 0.2])
title(['EMG-p = ' num2str(p)])
set(gcf, 'renderer', 'painters')

figure
histogram(allRemEmg,50)
hold on
histogram(allWakeEmg,50)
xlim([0 1])

figure
[xpdf fpdf] = ksdensity(allRemEmg);
xpdf = xpdf/max(xpdf);
plot(fpdf,xpdf)
hold on
[xpdf fpdf] = ksdensity(allWakeEmg);
xpdf = xpdf/max(xpdf);
plot(fpdf,xpdf)
xlim([0.1 1])
ylabel('Normalized probability')
xlabel('EMG')
keyboard