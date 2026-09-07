%jds_compareNREMCA1SuppressionREMCofiringShifters
%Same suppression-vs-REM-cofiring analysis as
%jds_compareNREMCA1SuppressionREMCofiring, run separately for CA1 cells
%classified as REM theta-phase shifters vs nonshifters. Plots suppression
%by cofiring quartile for each group and tests whether the two
%correlation coefficients differ (Fisher r-to-z).

clear all
close all
savedir = '/Volumes/JUSTIN/SingleDay/ProcessedDataREM/';
load([savedir 'CA1nremca1allripmodsigdata.mat'])

day = 1;

animalprefixlist = {'KL8','JS17','JS15','JS14','JS12','JS13','JS34','BG1','JS21','ZT2'};

cofiringModShift = [];
cofiringModNonshift = [];
cofiringFRchangeShift = [];
cofiringFRchangeNonshift = [];

for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
    idxs = find(allripplemod_idx(:,1) == a);
    animIdx = allripplemod_idx(idxs,2:end);
    animhists = allripmodhists(idxs,:);
    animmod = allripmodMI(idxs);

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
    load(sprintf('%s/%sCA1remshiftershighthresh%02d.mat',dir,animalprefix,day));

    epochs = remeps;

    for e = 1:length(epochs)

        ep = epochs(e);
        if ep == 1
            continue
        end

        shifters = shiftList{day}{ep}.cellidx;

        epidx = find(animIdx(:,2) == ep);
        animIdxEp = animIdx(epidx,3:4);
        animhistsEp = animhists(epidx,:);
        animmodEp = animmod(epidx);
        
        rips = [indrips.ctxripple{day}{ep}.starttimeC indrips.ctxripple{day}{ep}.endtimeC];
        if isempty(rips)
            continue
        end
        firstsws = [sws{day}{ep}.starttime(1) sws{day}{ep}.endtime(1)];
        lastsws = [sws{day}{ep}.starttime(end) sws{day}{ep}.endtime(end)];
        allrem = [rem{day}{ep}.starttime rem{day}{ep}.endtime];
        remdur = rem{day}{ep}.total_duration;
        if length(rips(:,1)) > 10
            load(sprintf('%s%sspikes%02d.mat',dir,animalprefix,day));
            numncrips = length(rips(:,1));

            CA1idx2 = CA1_R{ep}.cellidx;
            CA1idx = animIdxEp;
            PFCidx = PFC_R{ep}.cellidx;

            for c = 1:length(CA1idx(:,1))
                tmp = [];
                idx1 = find(CA1idx(c,1) == CA1idx2(:,1));
                idx2 = find(CA1idx(c,2) == CA1idx2(:,2));
                idx3 = intersect(idx1,idx2);
                if isempty(idx3)
                    continue
                end
                idx1 = find(CA1idx(c,1) == shifters(:,1));
                idx2 = find(CA1idx(c,2) == shifters(:,2));
                idx3 = intersect(idx1,idx2);
                if isempty(idx3)
                    continue
                end
                isShift = shifters(idx3,3);
                cell1spks = spikes{day}{ep}{CA1idx(c,1)}{CA1idx(c,2)}.data(:,1);
                firstnremFR = (sum(isExcluded(cell1spks, firstsws)))/(firstsws(2)-firstsws(1));
                lastnremFR = (sum(isExcluded(cell1spks, lastsws)))/(lastsws(2)-lastsws(1));
                remFR = (sum(isExcluded(cell1spks, allrem)))/remdur;
                FRchange = lastnremFR - firstnremFR;
                if FRchange < -5
                    continue
                end
                for pp = 1:length(PFCidx(:,1))
                    if (~isempty(spikes{day}{ep}{CA1idx(c,1)}{CA1idx(c,2)})) &&...
                            (~isempty(spikes{day}{ep}{PFCidx(pp,1)}{PFCidx(pp,2)}))
                        cell1spks = spikes{day}{ep}{CA1idx(c,1)}{CA1idx(c,2)}.data(:,1);
                        cell2spks = spikes{day}{ep}{PFCidx(pp,1)}{PFCidx(pp,2)}.data(:,1);

                        spkbins1_c = periodAssign(cell1spks, rips);
                        spkbins2_p = periodAssign(cell2spks, rips);

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
                if isShift == 1
                    cofiringModShift = [cofiringModShift; [nanmean(tmp) animmodEp(c)]];
                    cofiringFRchangeShift = [cofiringFRchangeShift; [FRchange animmodEp(c)]];
                else
                    cofiringModNonshift = [cofiringModNonshift; [nanmean(tmp) animmodEp(c)]];
                    cofiringFRchangeNonshift = [cofiringFRchangeNonshift; [FRchange animmodEp(c)]];
                end
            end
        end
    end
end

cofiringModShift = cofiringModShift(~isnan(cofiringModShift(:,1)),:);
cofiringModNonshift = cofiringModNonshift(~isnan(cofiringModNonshift(:,1)),:);

[r1 p1] = corrcoef(cofiringModShift)
[r2 p2] = corrcoef(cofiringModNonshift)

quartile_sep = floor(length(cofiringModShift(:,1))/4);
spklatquar = sortrows(cofiringModShift,1);
vals = [];
cnt = 1;
for s = 1:4
    if s < 4
        tmp = spklatquar(cnt:quartile_sep*s,1:2);
        tmp(:,3) = s;
    else
        tmp = spklatquar(cnt:end,1:2);
        tmp(:,3) = s;
    end
    vals{s} = tmp;
    cnt = cnt + quartile_sep;
    clear tmp
end

pWRS1 = ranksum(vals{1, 1}(:,2),vals{1, 4}(:,2));

v = cellfun(@mean,vals,'UniformOutput',false);

v2 = cellfun((@(x) std(x)./sqrt(length(x(:,1)))),vals,'UniformOutput',false);

data_sems = vertcat(v2{:});

data_means = vertcat(v{:});

X = [1:4];

figure; hold on
errorbar(X, data_means(:,2), data_sems(:,2),'r','LineWidth',3,'Capsize',0);
hold on
xlim([0.5 4.5])
xticks([1:4])

quartile_sep = floor(length(cofiringModNonshift(:,1))/4);
spklatquar = sortrows(cofiringModNonshift,1);
vals = [];
cnt = 1;
for s = 1:4
    if s < 4
        tmp = spklatquar(cnt:quartile_sep*s,1:2);
        tmp(:,3) = s;
    else
        tmp = spklatquar(cnt:end,1:2);
        tmp(:,3) = s;
    end
    vals{s} = tmp;
    cnt = cnt + quartile_sep;
    clear tmp
end

pWRS2 = ranksum(vals{1, 1}(:,2),vals{1, 4}(:,2));

v = cellfun(@mean,vals,'UniformOutput',false);

v2 = cellfun((@(x) std(x)./sqrt(length(x(:,1)))),vals,'UniformOutput',false);

data_sems = vertcat(v2{:});

data_means = vertcat(v{:});

errorbar(X, data_means(:,2), data_sems(:,2),'b','LineWidth',3,'Capsize',0);
title(['shift - p=' num2str(p1(1,2)) ' r=' num2str(r1(1,2)) ' non - p' num2str(p2(1,2)) ' r=' num2str(r2(1,2))])
ylabel('CA1 suppression')
xlabel(['shiftQ1vQ4 - p=' num2str(pWRS1) ' nonQ1vQ4 - p' num2str(pWRS2)])
set(gcf, 'renderer', 'painters')

%test for significant difference of correlation coefficients

r1 = r1(1,2);   % Correlation coefficient for group 1
r2 = r2(1,2);   % Correlation coefficient for group 2
n1 = size(cofiringModShift,1);     % Sample size for group 1
n2 = size(cofiringModNonshift,1);     % Sample size for group 2

% Step 1: Apply Fisher's r-to-z transformation
z1 = 0.5 * log((1 + r1) / (1 - r1));
z2 = 0.5 * log((1 + r2) / (1 - r2));

% Step 2: Compute the standard error
SE = sqrt(1 / (n1 - 3) + 1 / (n2 - 3));

% Step 3: Compute the Z statistic
Z = (z1 - z2) / SE;

% Step 4: Determine the p-value (two-tailed test)
p_value = 2 * (1 - normcdf(abs(Z)));  % normcdf calculates the CDF of the standard normal distribution

% Output results
fprintf('Z-statistic: %.4f\n', Z);
fprintf('P-value: %.4f\n', p_value);

% Step 5: Decision based on significance level (alpha = 0.05)
alpha = 0.05;
if p_value < alpha
    fprintf('The difference between the correlation coefficients is statistically significant (reject H0).\n');
else
    fprintf('The difference between the correlation coefficients is not statistically significant (fail to reject H0).\n');
end


