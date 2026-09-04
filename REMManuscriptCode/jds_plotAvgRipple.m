function jds_plotAvgRipple(animalprefixlist)
% Plot grand ripple average waveform
% -------------------------------------------------------------------------

pre = 150;
post = 150;
day = 1;
scaling = 0.195;
alleventseeg = [];
alleventsamp = [];
for a = 1:length(animalprefixlist)
    animAmp = [];
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    load(sprintf('%s%sctxripples0%d.mat',dir,animalprefix,day));% get ripple time
    rem = load(sprintf('%s%srem0%d.mat',dir,animalprefix,day));% get immmobility time
    rem = rem.rem;
    load(sprintf('%s%stetinfo.mat',dir,animalprefix));
    load(sprintf('%s%sremeps0%d.mat',dir,animalprefix,day));
    epochs = remeps;

    daystring = ['0',num2str(day)];

    samprate = 1500;
    for ep=1:length(epochs)
        epoch = epochs(ep);
        if epoch <10
            epochstring = ['0',num2str(epoch)];
            epochstring2 = ['0',num2str(epoch)];
        else
            epochstring = num2str(epoch);
            epochstring2 = num2str(epoch);
        end

        if strcmp(animalprefix,'ER1')
            epochstring = num2str(epoch);
            if epoch <10
                epochstring2 = ['0',num2str(epoch)];
            else
                epochstring2 = num2str(epoch);
            end
        end

        tets = tetinfo{1}{1};

        %to use all tetrodes
        ctxtets = [];
        for t = 1:length(tets)
            if (t == 17) && (strcmp(animalprefix,'JS14'))
                continue
            end
            tmp = tets{t};
            if isfield(tmp, 'descrip')
                if isequal(tmp.descrip, 'riptet')
                    ctxtets = [ctxtets; t];
                end
            end
        end
        ctxtimetet = ctxtets(1);
        remEp = rem{day}{epoch};
        remlist = [remEp.starttime remEp.endtime];

        %cortical ripples
        if (ctxtimetet<10)
            ctxtimetetstring = ['0',num2str(ctxtimetet)];
        else
            ctxtimetetstring = num2str(ctxtimetet);
        end

        curreegfile = [dir,'/EEG/',animalprefix,'ripple', daystring,'-',epochstring,'-',ctxtimetetstring];
        load(curreegfile);

        time1 = geteegtimes(ripple{day}{epoch}{ctxtimetet}) ; % construct time array

        for t = 1:length(ctxtets)
            tet = ctxtets(t);
            alltetrips = [];
            if (tet<10)
                ctxtetstring = ['0',num2str(tet)];
            else
                ctxtetstring = num2str(tet);
            end

            ctxrips = ripples{day}{epoch}{tet}.starttime;
            ctxinds = [ripples{day}{epoch}{tet}.startind ripples{day}{epoch}{tet}.endind];

            inRemRips = logical(isExcluded(ctxrips(:,1),remlist));
            ctxrips = ctxrips(inRemRips,:);
            ctxinds = ctxinds(inRemRips,:);

            ripeegfile = [dir,'/EEG/',animalprefix,'ripple', daystring,'-',epochstring,'-',ctxtetstring];
            eegfile = [dir,'/EEG/',animalprefix,'eegref', daystring,'-',epochstring2,'-',ctxtetstring];
            load(ripeegfile);
            load(eegfile);
            eegData = eegref{day}{epoch}{tet}.data;
            ramp = ripple{day}{epoch}{tet}.data(:,1);
            
            for r = 1:size(ctxinds,1)
                tmpInd = ctxinds(r,:);
                [M maxInd] = max(ramp(tmpInd(1):tmpInd(2)));
                peakInd = tmpInd(1) + maxInd;
                if tmpInd(1)-pre > 0 && tmpInd(2)+post < length(ramp)
                    ripAmp = double(ramp(peakInd-pre:peakInd+post))*double(scaling);
                    eegtmp = double(eegData(peakInd-pre:peakInd+post))*double(scaling);
                    alleventseeg = [alleventseeg; eegtmp'];
                    alleventsamp = [alleventsamp; ripAmp'];
                    animAmp = [animAmp;ripAmp'];
                    alltetrips = [alltetrips; eegtmp'];
                end
            end
        end
        clear remlist
    end
    subplot(2,5,a)
    boundedline(-pre:post,nanmean(animAmp),nanstd(animAmp)./sqrt(size(animAmp,1)),'-k');
    title(animalprefix)
end
figure; hold on;
pl1 = plot(-pre:post, mean(alleventsamp)); 
boundedline(-pre:post,nanmean(alleventsamp),nanstd(alleventsamp)./sqrt(size(alleventsamp,1)),'-k');
xticks([-150 0 150])
xlim([-150 150])
ylim([-12 12])
xticklabels({'-100','0','100'})
xlabel('Time from ripple (ms)')
ylabel('Ripple amplitude')

figure; hold on;
pl1 = plot(-pre:post, mean(alleventseeg)); 
boundedline(-pre:post,nanmean(alleventseeg),nanstd(alleventseeg)./sqrt(size(alleventseeg,1)),'-k');
xticks([-150 0 150])
xlim([-150 150])
xticklabels({'-100','0','100'})
xlabel('Time from ripple (ms)')
ylabel('Voltage (uV)')
title('PFC ripple grand average')

keyboard;

