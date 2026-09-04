function jds_rippleTriggeredAssemblyStrengthREMxVal(animalprefixlist,area,state)

%To do: plot the peak strength within specified ripple events. compare
%ripple types for different assemblies - might need to change binsize for
%raster generation
day = 1;

% bins = 400; %for 5 ms
% bins = 20; %for 100 ms
bins = 50; %for 20ms
% peakbins = find((-bins:bins)<=20 & (-bins:bins)>=0);
peakbins = find(abs(-bins:bins)<=25); %for +-500ms
% peakbins = find((-bins:bins)>=0 & (-bins:bins)<=25); %for +-500ms


rvalsAll = [];
rvalsAll_s = [];
rvalsSpearAll = [];
rvalsSpearAll_s = [];
seqFitFirstHalf = [];
seqFitSecondHalf = [];
slopediff = [];
slopediff_s = [];
shufnum = 1000;
g1 = gaussian(3, 10);
firstripinchain = 1;
dataComb = [];
ripthresh = [10];% 20 30 40 50 60 70 80 90 100];
for rips = 1:length(ripthresh)
    rvals = [];
    rvals_s = [];
    rvalsSpear = [];
    rvalsSpear_s = [];
    allData = [];
    rnum = ripthresh(rips);
    for s = 1:shufnum
        allevents_ctxriptrigRem1 = [];
        allevents_ctxriptrigRem2 = [];
        allevents_ctxriptrigRem_s1 = [];
        allevents_ctxriptrigRem_s2 = [];
        for a = 1:length(animalprefixlist)

            animalprefix = char(animalprefixlist(a));
            dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

            %Load reactivation strength file for all assemblies and epochs
            load(sprintf('%s%s%s_RTimeStrength%sNewSpk_20_%02d.mat',dir,animalprefix,area,state,day));
            %Load ripples
%             load(sprintf('%s%sctxgammatime2_chainREM%02d.mat',dir,animalprefix,day));
            load(sprintf('%s%sctxrippletime_chainREM%02d.mat',dir,animalprefix,day));
%% For matching
            ripcnts = load(sprintf('%s%sctxrippletime_chainREM%02d.mat',dir,animalprefix,day));
%%
            remripple = ctxripple;
%             remripple = ctxgamma;
            %     ctxripple = ripple;
            %     load(sprintf('%s%sctxrippletime_SWS%02d.mat',dir,animalprefix,day));
            load(sprintf('%s%sremeps%02d.mat',dir,animalprefix,day));

            %Which epochs to analyze
            %     epochs = find(~cellfun(@isempty,RtimeStrength));
            epochs = remeps;
            epochs = epochs(find(epochs>1));

            for e = 1:length(epochs)
                ep = epochs(e);
                assemblytmp = RtimeStrength{ep}.reactivationStrength;
                %         nremripstarts = ctxripple{day}{ep}.starttime;
                if firstripinchain == 1
                    tmp = remripple{day}{ep}.C_sep;
                    remripstarts = [];
                    for c = 1:length(tmp)
                        remripstarts = [remripstarts; tmp{c}(1,1)];
%                         remripstarts = [remripstarts; tmp{c}(1,1)+(tmp{c}(end,2)-tmp{c}(1,1))/2];
                    end
                else
%                     remripstarts = remripple{day}{ep}.starttimeC;
%                     rip_starttime = remripstarts*1000;
% 
%                     iri = diff(rip_starttime);
%                     keepidx = [1;find(iri>=1000)+1]; %events separated by 1 s
%                     
%                     remripstarts = remripstarts(keepidx);
                end

                remchaincnt = length(ripcnts.ctxripple{day}{ep}.C_sep);

                randvals = randperm(length(remripstarts));
                halfEvs = floor(length(randvals)/2);
                rips1 = remripstarts(randvals(1:halfEvs));
                rips2 = remripstarts(randvals(halfEvs+1:end));
                rips1 = sort(rips1);
                rips2 = sort(rips2);
                shiftvals = [];
                %Do PFC ripples
                if remchaincnt > rnum
                    if ~isempty(assemblytmp)
                        for ii = 1:length(assemblytmp)
                            react_idx = [];
                            for t = 1:length(rips1)
                                idxtmp = lookup(rips1(t), assemblytmp{ii}(:,1));
                                react_idx = [react_idx; idxtmp];
                            end
                            atmp = [];
                            atmp_s = [];
                            strengthstmp = assemblytmp{ii}(:,2);
                            strengthstmp(strengthstmp<0) = 0;
                            shift = randi(length(strengthstmp),1);
                            shiftvals = [shiftvals; shift];
                            strengthstmpshuf = circshift(strengthstmp,shift);
                            for r = 1:length(react_idx)
                                if ((react_idx(r) + bins) < length(strengthstmp)) && ((react_idx(r) - bins) > 1)
                                    tmp = strengthstmp((react_idx(r) - bins):(react_idx(r) + bins)); %get vector of reactivation strenths for specified time period
                                    atmp = [atmp; tmp'];
                                end
                            end
                            allevents_ctxriptrigRem1 = [allevents_ctxriptrigRem1; zscore(smoothvect(mean(atmp,1),g1))];
                            for r = 1:length(react_idx)
                                if ((react_idx(r) + bins) < length(strengthstmpshuf)) && ((react_idx(r) - bins) > 1)
                                    tmp = strengthstmpshuf((react_idx(r) - bins):(react_idx(r) + bins)); %get vector of reactivation strenths for specified time period
                                    atmp_s = [atmp_s; tmp'];
                                end
                            end
                            allevents_ctxriptrigRem_s1 = [allevents_ctxriptrigRem_s1; zscore(smoothvect(mean(atmp_s,1),g1))];
                        end
                        %Do CA1 ripples
                        for ii = 1:length(assemblytmp)
                            react_idx = [];
                            for t = 1:length(rips2)
                                idxtmp = lookup(rips2(t), assemblytmp{ii}(:,1));
                                react_idx = [react_idx; idxtmp];
                            end
                            atmp = [];
                            atmp_s = [];
                            strengthstmp = assemblytmp{ii}(:,2);
                            strengthstmp(strengthstmp<0) = 0;
                            strengthstmpshuf = circshift(strengthstmp,shiftvals(ii));
                            for r = 1:length(react_idx)
                                if ((react_idx(r) + bins) < length(strengthstmp)) && ((react_idx(r) - bins) > 1)
                                    tmp = strengthstmp((react_idx(r) - bins):(react_idx(r) + bins));
                                    atmp = [atmp; tmp'];
                                end
                            end
                            allevents_ctxriptrigRem2 = [allevents_ctxriptrigRem2; zscore(smoothvect(mean(atmp,1),g1))];
                            for r = 1:length(react_idx)
                                if ((react_idx(r) + bins) < length(strengthstmpshuf)) && ((react_idx(r) - bins) > 1)
                                    tmp = strengthstmpshuf((react_idx(r) - bins):(react_idx(r) + bins));
                                    atmp_s = [atmp_s; tmp'];
                                end
                            end
                            allevents_ctxriptrigRem_s2 = [allevents_ctxriptrigRem_s2; zscore(smoothvect(mean(atmp_s,1),g1))];
                        end
                    end
                end
            end
        end
        allData{s}.firstHalfRipples = allevents_ctxriptrigRem1;
        allData{s}.secondHalfRipples = allevents_ctxriptrigRem2;
        allData{s}.firstHalfRipplesShuf = allevents_ctxriptrigRem_s1;
        allData{s}.secondHalfRipplesShuf = allevents_ctxriptrigRem_s2;
        [M1 I1] = max(allevents_ctxriptrigRem1(:,peakbins)');
        [M2 I2] = max(allevents_ctxriptrigRem2(:,peakbins)');
        [B1,II1] = sort(I1);
        p1 = polyfit(B1*20,1:length(B1),1);
        seqFitFirstHalf = [seqFitFirstHalf; B1];
        B2 = I2(II1);
        p2 = polyfit(B2*20,1:length(B2),1);
        slDiff = p1(1)-p2(1);
        slopediff = [slopediff; slDiff];
        seqFitSecondHalf = [seqFitSecondHalf; B2];
        rvalsSpear = [rvalsSpear; corr(B1',B2','Type','Spearman')];
        R = corrcoef(I1,I2);
        rvals = [rvals; R(1,2)];
        disp(num2str(R(1,2)))
        allData{s}.rval = R(1,2);
        allData{s}.rvalSpear = corr(B1',B2','Type','Spearman');
        [M1s I1s] = max(allevents_ctxriptrigRem_s1(:,peakbins)');
        [M2s I2s] = max(allevents_ctxriptrigRem_s2(:,peakbins)');
        [B1s,II1s] = sort(I1s);
        p1s = polyfit(B1s*20,1:length(B1s),1);
        B2s = I2s(II1s);
        p2s = polyfit(B2s*20,1:length(B2s),1);
        slDiff_s = p1s(1)-p2s(1);
        slopediff_s = [slopediff_s; slDiff_s];
        rvalsSpear_s = [rvalsSpear_s; corr(B1s',B2s','Type','Spearman')];
        R_s = corrcoef(I1s,I2s);
        rvals_s = [rvals_s; R_s(1,2)];
        allData{s}.rval_s = R_s(1,2);
        allData{s}.rvalSpear_s = corr(B1s',B2s','Type','Spearman');
    end
    rvalsAll{rips} = rvals;
    rvalsAll_s{rips} = rvals_s;
    rvalsSpearAll{rips} = rvalsSpear;
    rvalsSpearAll_s{rips} = rvalsSpear_s;
    dataComb.data = allData;
    dataComb.descrip = 'split first chain ripple cross validation, pearson correlation';
    dataComb.rvals = rvals;
    dataComb.rvals_s = rvals_s;
    dataComb.rvalsSpear = rvalsSpear;
    dataComb.rvalsSpear_s = rvalsSpear_s;
end
figure
bar([mean(dataComb.rvals) mean(dataComb.rvals_s)],'k');
hold on; errorbar(1:2,[mean(dataComb.rvals) mean(dataComb.rvals_s)],[std(dataComb.rvals)./sqrt(1000) std(dataComb.rvals_s)./sqrt(1000)],'k','LineStyle','none');
set(gcf, 'renderer', 'painters')
ylabel('r')
p = ranksum(dataComb.rvals  ,dataComb.rvals_s)
title(['First NREM ripple chain assembly seq xval p=' num2str(p)])
set(gcf, 'renderer', 'painters')
figure
histogram(dataComb.rvals,50)
hold on
histogram(dataComb.rvals_s,50)
ylabel('Count')
title(['First NREM ripple chain assembly seq xval p=' num2str(p)])
set(gcf, 'renderer', 'painters')

figure
bar([mean(slopediff) mean(slopediff_s)],'k');
hold on; errorbar(1:2,[mean(slopediff) mean(slopediff_s)],[std(slopediff)./sqrt(1000) std(slopediff_s)./sqrt(1000)],'k','LineStyle','none');
set(gcf, 'renderer', 'painters')
keyboard
