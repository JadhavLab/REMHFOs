function [autocorrSWS, autocorrREM] = jds_xCorrREM(animalprefixlist)
% Calculate and plot cross correlation or autocorrelation (depending on
% input)
% -------------------------------------------------------------------------

day =1;
bin = 0.01;
tmax = 1;
tbins = -(tmax/bin):(tmax/bin)-1;
autobins = 92:109;
sw1 = bin*3; 
autocorrREM = [];
autocorrSWS = [];
for r = 1:2
    for a = 1:length(animalprefixlist)

        animalprefix = char(animalprefixlist(a));

        dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
        if r == 1
            load(sprintf('%s%sctxrippletime_SWS0%d.mat',dir,animalprefix,day));
        elseif r == 2
            load(sprintf('%s%sctxrippletime_REM0%d.mat',dir,animalprefix,day));
        end
        rem = load(sprintf('%s%srem0%d.mat',dir,animalprefix,day));
        load(sprintf('%s%sremeps0%d.mat',dir,animalprefix,day));
        rem = rem.rem;
        epochs = remeps;
        for ep=1:length(epochs)

            epoch = epochs(ep);

            ctx = ctxripple{day}{epoch};

            if ~isempty(rem{day}{epoch}.starttime)

                if (~isempty(ctx))
                    if length(ctx.starttime) > 10

                        ctxmidtimes = ctx.starttime;

                        ctxmidtimes_shuf = ctx.starttime;

                        xc = spikexcorr(ctxmidtimes, ctxmidtimes, bin, tmax);

                        if ~isempty(xc.c1vsc2)
                            Zcrosscov = zscore(xc.c1vsc2);
                            prob = xc.c1vsc2./sum(xc.c1vsc2);
                            nstd=round(sw1/(xc.time(2) - xc.time(1))); 
                            g1 = gaussian(nstd, nstd);
                            timebase = xc.time;
                            bins_run = find(abs(timebase) <= tmax); 

                            Zcrosscov_sm = smoothvect(Zcrosscov, g1);
                            Zcrosscov_sm(autobins) = NaN;
                            if r == 1
                                autocorrSWS = [autocorrSWS; Zcrosscov_sm];
                            elseif r == 2
                                autocorrREM = [autocorrREM; Zcrosscov_sm];
                            end
                        end
                    end
                end
            end
        end
    end
end

%%
figure; hold on
ylim([-0.5 0.5])

sigbins = [];
for r = 1:length(autocorrSWS(1,:))
    if ~isempty(find(autobins == r))
        sig = 0;
        p = NaN;
    else
%         [h p] = ttest2(autocorrSWS(:,r),autocorrREM(:,r));
        [p h] = signrank(autocorrSWS(:,r),autocorrREM(:,r));
        if p < 0.05
            sig = 1;
        else
            sig = 0;
        end
    end
    sigbins = [sigbins; [r sig p]]
end
comps = bwconncomp(sigbins(:,2));
list = comps.PixelIdxList;

figure; hold on
boundedline(tbins,nanmean(autocorrSWS),nanstd(autocorrSWS)./sqrt(size(autocorrSWS,1)),'-k');
boundedline(tbins,nanmean(autocorrREM),nanstd(autocorrREM)./sqrt(size(autocorrREM,1)),'-r');

for l = 1:length(list)
    startidx = tbins(list{l}(1));
    endidx = tbins(list{l}(end));
    x = [startidx endidx];
    y = ones(1,size(x,2))*0.25;
    plot(x,y,'*','MarkerSize',8,'MarkerEdgeColor','k')
end
xlim([-50 50])
xticks([-50:25:50])
set(gcf, 'renderer', 'painters')
keyboard
