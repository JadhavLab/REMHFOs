function jds_plotAvgSpindle(animalprefixlist)

% SET DATA
%--------------------add your animals' directories--------------------%
pre = 3000;
post = 3000;
day = 1;
scaling = 0.195;

for rips = 1:2
    alleventseeg = [];
    for a = 1:length(animalprefixlist)
        animAmp = [];
        animalprefix = animalprefixlist{a};
        dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
        %     dir = sprintf('/Volumes/JUSTIN/Inference/%s_direct/',animalprefix);

        load(sprintf('%s%sctxspindletime_chainSWS0%d.mat',dir,animalprefix,day));% get ripple time
        ctxripples = ctxspindle;
        load(sprintf('%s%sswsALL0%d.mat',dir,animalprefix,day));% get immmobility time
        rem = sws;
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
                    if isequal(tmp.descrip, 'ctxriptet')
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

            curreegfile = [dir,'/EEG/',animalprefix,'eeg', daystring,'-',epochstring,'-',ctxtimetetstring];
            load(curreegfile);

            time1 = geteegtimes(eeg{day}{epoch}{ctxtimetet}) ; % construct time array
            allTetEeg = [];
            for t = 1:length(ctxtets)
                tet = ctxtets(t);
                alltetrips = [];
                if (tet<10)
                    ctxtetstring = ['0',num2str(tet)];
                else
                    ctxtetstring = num2str(tet);
                end

                eegfile = [dir,'/EEG/',animalprefix,'eeg', daystring,'-',epochstring2,'-',ctxtetstring];
                load(eegfile);
                eegData = eeg{day}{epoch}{tet}.data;
                allTetEeg = [allTetEeg; eegData'];
            end
            meanEeg = mean(allTetEeg);
            [b,a] = butter(3,[10/(samprate/2) 16/(samprate/2)]);
            spinAmp = filtfilt(b,a,meanEeg);

            if rips == 1
                rtimes = [ctxspindle{day}{epoch}.starttimeNC ctxspindle{day}{epoch}.endtimeNC];
            elseif rips == 2
                rtimes = [ctxspindle{day}{epoch}.starttimeC ctxspindle{day}{epoch}.endtimeC];
            end
            ctxinds = [lookup(rtimes(:,1),time1) lookup(rtimes(:,2),time1)];
            for r = 1:size(ctxinds,1)
                tmpInd = ctxinds(r,:);
                [M maxInd] = max(spinAmp(tmpInd(1):tmpInd(2)));
                peakInd = tmpInd(1) + maxInd;
                if tmpInd(1)-pre > 0 && tmpInd(2)+post < length(meanEeg)
                    eegtmp = double(meanEeg(peakInd-pre:peakInd+post))*double(scaling);
                    alleventseeg = [alleventseeg; eegtmp];
                end
            end
            clear remlist
        end
    end
    figure(rips); hold on
    pl1 = plot(-pre:post, mean(alleventseeg));
    boundedline(-pre:post,nanmean(alleventseeg),nanstd(alleventseeg)./sqrt(size(alleventseeg,1)),'-k');
    if rips == 1
        title('Avg Isolated spindle')
    elseif rips == 2
        title('Avg train spindle')
    end
end
keyboard
figure; hold on;
subplot(2,1,1)
pl1 = plot(-pre:post, mean(alleventsamp));
boundedline(-pre:post,nanmean(alleventsamp),nanstd(alleventsamp)./sqrt(size(alleventsamp,1)),'-k');
xticks([-150 0 150])
xlim([-150 150])
ylim([-12 12])
xticklabels({'-100','0','100'})
xlabel('Time from ripple (ms)')
ylabel('Ripple amplitude')

subplot(2,1,2)
pl2 = plot(-pre:post, mean(alleventsenv));
boundedline(-pre:post, nanmean(alleventsenv),nanstd(alleventsenv)./sqrt(size(alleventsenv,1)),'-k');
xticks([-150 0 150])
xlim([-150 150])
xticklabels({'-100','0','100'})
xlabel('Time from ripple (ms)')
ylabel('Ripple envelope')
set(gcf, 'renderer', 'painters')

figure; hold on;
pl1 = plot(-pre:post, mean(alleventsamp));
boundedline(-pre:post,nanmean(alleventsamp),nanstd(alleventsamp)./sqrt(size(alleventsamp,1)),'-k');
xticks([-150 0 150])
xlim([-150 150])
ylim([-12 12])
xticklabels({'-100','0','100'})
xlabel('Time from ripple (ms)')
ylabel('Ripple amplitude (uV)')


figure; hold on;
% subplot(2,1,1)
pl1 = plot(-pre:post, mean(alleventseeg));
boundedline(-pre:post,nanmean(alleventseeg),nanstd(alleventseeg)./sqrt(size(alleventseeg,1)),'-k');
xticks([-150 0 150])
xlim([-150 150])
% ylim([-12 12])
xticklabels({'-100','0','100'})
xlabel('Time from ripple (ms)')
ylabel('Voltage (uV)')
title('PFC ripple grand average')
%
% subplot(2,1,2)
% pl1 = plot(-pre:post, mean(alleventsamp));
% boundedline(-pre:post,nanmean(alleventsamp),nanstd(alleventsamp)./sqrt(size(alleventsamp,1)),'-k');
% xticks([-150 0 150])
% xlim([-150 150])
% ylim([-12 12])
% xticklabels({'-100','0','100'})
% xlabel('Time from ripple (ms)')
% ylabel('Ripple amplitude')


keyboard;

