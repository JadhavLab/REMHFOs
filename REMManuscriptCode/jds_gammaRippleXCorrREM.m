function [c1vsc2, xcorr_sm_all, Zcrosscov, Zcrosscov_sm_all] = jds_gammaRippleXCorrREM(animalprefixlist)

day =1;
bin = 0.01;
tmax = 1;
tbins = -(tmax/bin):(tmax/bin)-1;
sw1 = bin*3; % for smoothing corrln. Necessary?
shufnum = 10;
shuf = 1;
ripgammacorr = [];
ripgammacorr_shuf = [];
for a = 1:length(animalprefixlist)

    animalprefix = char(animalprefixlist(a));

    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
    load(sprintf('%s%sctxrippletime_REM0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%sctxgammatime2_REM0%d.mat',dir,animalprefix,day));
    %             load(sprintf('%s%sctxrippletime_REM0%d.mat',dir,animalprefix,day));
    %             load(sprintf('%s%sctxgammatime_REM0%d.mat',dir,animalprefix,day));
    rem = load(sprintf('%s%srem0%d.mat',dir,animalprefix,day));
    load(sprintf('%s%sremeps0%d.mat',dir,animalprefix,day));
    %     load(sprintf('%s%sswsALL0%d.mat',dir,animalprefix,day));
    rem = rem.rem;
    epochs = remeps;
    for ep=1:length(epochs)

        epoch = epochs(ep);

        ctx = ctxripple{day}{epoch};
        gam = ctxgamma{day}{epoch};

        if ~isempty(rem{day}{epoch}.starttime)

            if (~isempty(ctx))
                if length(ctx.starttime) > 10

                    ctxriptimes = ctx.starttime;
                    ctxgamtimes = gam.starttime;
                    ctxriptimes_shuf = ctx.starttime;
                    ctxgamtimes_shuf = gam.starttime;
                    iri = diff(ctxriptimes_shuf); %use distribution of IRIs to jitter ctxriptimes
                    iri2 = diff(ctxgamtimes_shuf);
                    iri(find(iri > 5)) = [];
                    iri2(find(iri2 > 5)) = [];
                    if shuf
                        for s = 1:shufnum
                            %%
                            %jitter ctx
                            sizeR = [1 length(ctxriptimes_shuf)] ;
                            num1 = floor(length(ctxriptimes_shuf)/2);
                            R = zeros(sizeR);  % set all to zero
                            ix = randperm(numel(R)); % randomize the linear indices
                            ix = ix(1:num1); % select the first
                            R(ix) = 1; % set the corresponding positions to 1
                            R(find(R == 0)) = -1;
                            if ~isempty(iri)
                                jittimes = (datasample(iri,length(ctxriptimes_shuf)).*R')+(randi(5,1,length(R))'.*R');
                            else
                                jittimes = (randi(5,1,length(R))'.*R');
                            end
                            %                             jittimes = (rand(1,length(R))')+(randi(10,1,length(R))'.*R');

                            ctxriptimes_shuf = sort(ctxriptimes + jittimes);
                            %%
                            %                             jitter hp
                            sizeR = [1 length(ctxgamtimes_shuf)] ;
                            num1 = floor(length(ctxgamtimes_shuf)/2);
                            R = zeros(sizeR);  % set all to zero
                            ix = randperm(numel(R)); % randomize the linear indices
                            ix = ix(1:num1); % select the first
                            R(ix) = 1; % set the corresponding positions to 1
                            R(find(R == 0)) = -1;
                            if ~isempty(iri2)
                                jittimes = (datasample(iri2,length(ctxgamtimes_shuf)).*R')+(randi(5,1,length(R))'.*R');
                            else
                                jittimes = (randi(5,1,length(R))'.*R');
                            end
                            ctxgamtimes_shuf = sort(ctxgamtimes + jittimes);
                            %%
                            xc_shuf = spikexcorr(ctxgamtimes_shuf, ctxriptimes_shuf, bin, tmax);

                            if ~isempty(xc_shuf.c1vsc2)
                                Zcrosscov_shuf = zscore(xc_shuf.c1vsc2);

                                nstd=round(sw1/(xc_shuf.time(2) - xc_shuf.time(1))); % will be 3 std
%                                 g1 = gaussian(nstd, nstd);
                                g1 = gaussian(0.1*10, ceil(8*0.1*10))
                                timebase = xc_shuf.time;
                                bins_run = find(abs(timebase) <= tmax); % +/- Corrln window

                                %                     xcorr_sm = smoothvect(xc.c1vsc2, g1);
                                Zcrosscov_sm_shuf = smoothvect(Zcrosscov_shuf, g1);% smoothed
                                %                     Zcrosscov_sm = Zcrosscov;

                                ripgammacorr_shuf = [ripgammacorr_shuf; Zcrosscov_sm_shuf];
                            end
                        end
                    end
                    xc = spikexcorr(ctxgamtimes, ctxriptimes, bin, tmax);

                    if ~isempty(xc.c1vsc2)
                        Zcrosscov = zscore(xc.c1vsc2);
                        prob = xc.c1vsc2./sum(xc.c1vsc2);
                        nstd=round(sw1/(xc.time(2) - xc.time(1))); % will be 3 std
%                         g1 = gaussian(nstd, nstd);
                        g1 = gaussian(0.1*10, ceil(8*0.1*10));
                        timebase = xc.time;
                        bins_run = find(abs(timebase) <= tmax); % +/- Corrln window

                        Zcrosscov_sm = smoothvect(Zcrosscov, g1);% smoothed
                        ripgammacorr = [ripgammacorr; Zcrosscov_sm];
                    end
                end
            end
        end
    end
end

%%
figure; hold on
boundedline(tbins,nanmean(ripgammacorr),nanstd(ripgammacorr)./sqrt(size(ripgammacorr,1)),'-k');
upperBnd = prctile(ripgammacorr_shuf,95);
lowerBnd = prctile(ripgammacorr_shuf,5);
xticks([-100:50:100])
xticklabels({'-1','-0.5','0','0.5','1'})
xlabel('Lag(s)')
ylabel('Cross correlation (z)')
set(gcf, 'renderer', 'painters')
x = [0 0];
y = [-0.6 1.2];
hold on
plot(x,y,'--r')
title('PFC REM Ripple-Gamma xcorr')
%Compare to shuffled times
if shufnum > 0
    figure(1); hold on;
    plot([-(tmax/bin):(tmax/bin)-1],upperBnd,'--r')
    plot([-(tmax/bin):(tmax/bin)-1],lowerBnd,'--r')
end

set(gcf, 'renderer', 'painters')
keyboard
