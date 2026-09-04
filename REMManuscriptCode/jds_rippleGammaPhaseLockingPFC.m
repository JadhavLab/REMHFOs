function jds_rippleGammaPhaseLockingPFC(animalprefixlist)
% Phase locking of PFC units to extracted gamma and ripple oscillations
% -------------------------------------------------------------------------
days = 1;
allcellphaseripple_rem = [];
allcellphasegamma_rem = [];
for r = 1:2
    for a = 1:length(animalprefixlist)
        animalprefix = animalprefixlist{a};

        dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

        for d = 1:length(days)
            day = days(d);

            spikes = loaddatastruct(dir, animalprefix, 'spikes', day); % get spikes
            if r == 1
                load(sprintf('%s%sctxrippletime_REM0%d.mat',dir,animalprefix,day));
                evs = ctxripple;
            else
                load(sprintf('%s%sctxgammatime_REM0%d.mat',dir,animalprefix,day));
                evs = ctxgamma;
            end
            load(sprintf('%s%sremeps0%d.mat',dir,animalprefix,day));
            
            epochs = remeps;

            for ep = 1:length(epochs)
                epoch = epochs(ep);
                [ctxidx, hpidx] = jds_getallepcells(dir, animalprefix, day, ep, []); %(tet, cell)
                ctxnum = length(ctxidx(:,1));
                hpnum = length(hpidx(:,1));


                if epoch <10
                    epochstring = ['0',num2str(epoch)];
                else
                    epochstring = num2str(epoch);
                end

                eventtimes = [evs{day}{epoch}.starttime evs{day}{epoch}.endtime];
                if ~isempty(eventtimes)
                    if length(eventtimes(:,1)) > 10
                        for cellcount = 1:ctxnum %get spikes for each cell
                            index = [day,epoch,ctxidx(cellcount,:)] ;
                            if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)})
                                if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
                                    spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
                                else
                                    spiketimes = [];
                                end
                            else
                                spiketimes = [];
                            end
                            spikebins = periodAssign(spiketimes, eventtimes(:,[1 2])); %Assign spikes to align with each ripple event (same number = same rip event, number indicates ripple event)
                            if ~isempty(spiketimes)
                                validspikes = find(spikebins);
                                spiketimes = spiketimes(validspikes); %get spike times that happen during ripples
                            end
                            if ~isempty(spiketimes)
                                hptet = ctxidx(cellcount,1);

                                if (hptet<10)
                                    hptetstring = ['0',num2str(hptet)];
                                else
                                    hptetstring = num2str(hptet);
                                end
                                if r == 1
                                    curreegfile = [dir,'/EEG/',animalprefix,'ripple', ['0' num2str(day)],'-',epochstring,'-',hptetstring];
                                    g = load(curreegfile);
                                    starttime = g.ripple{day}{epoch}{hptet}.starttime;
                                    spikeidxs = floor((spiketimes(:,1)-starttime)*1500)+1;
                                    phasedata = g.ripple{day}{epoch}{hptet}.data(:,2)*-1;
                                else
                                    curreegfile = [dir,'/EEG/',animalprefix,'gamma', ['0' num2str(day)],'-',epochstring,'-',hptetstring];
                                    g = load(curreegfile);
                                    starttime = g.gamma{day}{epoch}{hptet}.starttime;
                                    spikeidxs = floor((spiketimes(:,1)-starttime)*1500)+1;
                                    phasedata = g.gamma{day}{epoch}{hptet}.data(:,2)*-1;
                                end

                                spikephasetmp = double(phasedata(spikeidxs))/10000;

                                [ripplehat, kappa] = circ_vmpar(spikephasetmp); % Better to give raw data. Can also give binned data.
                                [m, ph] = modulation(spikephasetmp);
                                % Make finer polar plot and overlay Von Mises Distribution Fit.
                                % Circ Stats Box Von Mises pdf uses a default of 100 angles/nbin
                                % -------------------------------------------------------------
                                nbins = 50;
                                bins = -pi:(2*pi/nbins):pi;
                                count = histc(spikephasetmp, bins);

                                % Make Von Mises Fit
                                alpha = linspace(-pi, pi, 50)';
                                [pdf] = circ_vmpdf(alpha,ripplehat,kappa);
                                stats = rayleigh_test(spikephasetmp);
                                if r == 1
                                    allcellphaseripple_rem = [allcellphaseripple_rem; ph];
                                else
                                    allcellphasegamma_rem = [allcellphasegamma_rem; ph];
                                end
                            end
                        end
                    end
                end
                clear eventtimes
            end
        end
    end
end
allcellphaseripple_remDup = [allcellphaseripple_rem*(180/pi); (allcellphaseripple_rem*(180/pi))+360];
allcellphasegamma_remDup = [allcellphasegamma_rem*(180/pi); (allcellphasegamma_rem*(180/pi))+360];
figure; hold on
histogram(allcellphaseripple_remDup,50,'DisplayStyle','bar')
histogram(allcellphasegamma_remDup,50,'DisplayStyle','bar')
ylabel('Count')
xticks([0:180:720])
xlim([0 720])
ylim([0 110])
xlabel('Phase')
legend({'Ripple','Gamma'})
x = [360 360];
y = [0 110];
plot(x,y)
set(gcf, 'renderer', 'painters')

[p1,U2_1] = watsons_U2_approx_p(allcellphaseripple_rem*(180/pi), allcellphasegamma_rem*(180/pi))
title(['PFC phase locking to REM ripples and gamma - p=' num2str(p1) ' U2=' num2str(U2_1)])
keyboard;
