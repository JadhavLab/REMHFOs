function jds_extractHighThetaPhasicREM(animalprefixlist)


day = 2;
epochs = [1:2:17];
daystring = sprintf('%02d',day);
savedata = 1;
smoothing_width = 0.005; %50 ms
W = rectwin(11);
samprate = 1500;
g1 = gaussian(smoothing_width*samprate, ceil(8*smoothing_width*samprate));
totalPhasic = [];
totalRem = [];
for a = 1:length(animalprefixlist)
    animalprefix = animalprefixlist{a};
    dir = sprintf('/Volumes/JUSTIN/Inference/%s_direct/', animalprefix);
    
    load(sprintf('%s%stetinfo.mat',dir,animalprefix));
    rem = load(sprintf('%s%srem0%d.mat',dir,animalprefix,day));% get sws time
    rem = rem.rem;
    

    tets = tetinfo{1}{epochs(1)};
   
    ctxtets = []; %get all ctxriptet tetrodes
    for t = 1:length(tets)
        tmp = tets{t};
        if isfield(tmp, 'descrip')
            if isequal(tmp.descrip, 'ctxriptet')
                ctxtets = [ctxtets; t];
            end
        end
    end
    
    for e = 1:length(epochs)
        ampdataall = [];
        envdataall = [];
        epoch = epochs(e);

        remtime = rem{day}{epoch};
        remlist = [remtime.starttime remtime.endtime];
        
        if epoch <10
            epochstring = ['0',num2str(epoch)];
        else
            epochstring = num2str(epoch);
        end
% 
%         load(sprintf('%s%sctxripples0%d.mat',dir,animalprefix,day));
%         rTets = find(~cellfun(@isempty,ctxripples{day}{epoch}));
% 
%         tetsNumRips = [];
%         for scan = 1:length(rTets)
%             t = rTets(scan);
%             numR = length(ctxripples{day}{epoch}{t}.startind);
%             tetsNumRips = [tetsNumRips; numR];
%         end
%         [ripcnt idx] = max(tetsNumRips);
% 
%         ctxtets = rTets(idx);

        for i = 1:length(ctxtets)
            ctxtet = ctxtets(i);
            
            if (ctxtet<10)
                ctxtetstring = ['0',num2str(ctxtet)];
            else
                ctxtetstring = num2str(ctxtet);
            end
            
            curreegfile = [dir,'/EEG/',animalprefix,'theta', daystring,'-',epochstring,'-',ctxtetstring];
            load(curreegfile);
            
            ampdatatmp = theta{day}{epoch}{ctxtet}.data(:,1);
            envdatatmp = theta{day}{epoch}{ctxtet}.data(:,3);
            ampdataall = [ampdataall; ampdatatmp'];
            envdataall = [envdataall; envdatatmp'];
        end
        times = geteegtimes(theta{day}{epoch}{ctxtet}) ; % construct time array
        
        ampdata = mean(ampdataall,1); %mean amplitude across all tetrodes
        ampdata = conv(ampdata,g1,'same');
        envdata = mean(envdataall,1);
        envdata = conv(envdata,W,'same');
        ampdataDiff = diff(ampdata);
        %Get the mean env of theta during REM here for later
        %thresholding
        
        inRemAmp = logical(isExcluded(times,remlist));
        
        meanRemEnv = mean(envdata(inRemAmp));

        highThet = envdata > (meanRemEnv + 4*std(meanRemEnv));

        evList = vec2list(highThet,times);

        longIdx = find(evList(:,2)-evList(:,1) > 0.9);
        finalEpochs = evList(longIdx,:); %greater than 900 ms
% 
%         %Zscore amplitude data and find epochs where peak is >=-2 (LFP is
%         %flipped)
%         z_amp = zscore(ampdataDiff);
%         zci = find(diff(sign(z_amp))); %find the zero crossing indices
%         startidx = zci(1:end-1); %look at pretty much every interval
%         endidx = zci(2:end);
%         indices = [startidx' endidx'];
%         
%         peak_idx = [];
%         for l = 1:length(indices(:,1))
%             idx = indices(l,:);
%             ampvector = z_amp(idx(1):idx(2));
%             [max_z minidx] = max(ampvector);
%             
%             if max_z > 0
%                 peak_idx = [peak_idx; idx(1)];
%             end
%         end
%         
%         thetatimes = times(peak_idx);
% 
%         thetaIEI = diff(thetatimes);
%         
%         thetaIEIsm = conv(thetaIEI,W,'same'); %smoothing with this amplifies values but shouldnt matter downstream
%         lowIEI = prctile(thetaIEIsm,10);
% 
%         candEvs = thetaIEIsm < lowIEI;
%         evList = vec2list(candEvs,thetatimes);
%         longIdx = find(evList(:,2)-evList(:,1) > 0.9);
%         candEvs = evList(longIdx,:); %greater than 900 ms
% 
%         finalEpochs = [];
%         for i = 1:size(candEvs,1)
%             tmpInterval = candEvs(i,:);
%             tmpIdx = logical(isExcluded(times,tmpInterval));
%             ampTmp = mean(ampdata(tmpIdx));
% %             ampTmp = mean(envdata(tmpIdx));
%             if ampTmp > meanRemEnv
%                 finalEpochs = [finalEpochs; candEvs(i,:)];
%             end
%         end

        %Constrain by REM
        if (~isempty(remtime.starttime)) && (~isempty(finalEpochs))
            totalRem = [totalRem; remtime.total_duration];
            remlist = [remtime.starttime remtime.endtime];
            [~,remvec] = wb_list2vec(remlist,times);
            [~,phasicvec] = wb_list2vec(finalEpochs,times);
            
            phasic_new = remvec & phasicvec;
            
            phasictimes = vec2list(phasic_new,times);

            phasictimes = phasictimes(find(phasictimes(:,2)-phasictimes(:,1) > 0.9),:);
            
            if ~isempty(phasictimes)
                totalPhasic = [totalPhasic; sum(phasictimes(:,2)-phasictimes(:,1))];
                phasicrem{day}{epoch}.starttime = phasictimes(:,1);
                phasicrem{day}{epoch}.endtime = phasictimes(:,2);
                phasicrem{day}{epoch}.total_duration = sum(phasictimes(:,2)-phasictimes(:,1));
                phasicrem{day}{epoch}.descrip = 'phasic rem extracted from mean of ctxriptets, Mizugeki method';
            else
                totalPhasic = [totalPhasic; 0];
                phasicrem{day}{epoch}.starttime = [];
                phasicrem{day}{epoch}.endtime = [];
                phasicrem{day}{epoch}.total_duration = 0;
                phasicrem{day}{epoch}.descrip = 'phasic rem extracted from mean of ctxriptets, Mizugeki method';
            end
        else
            if ~isempty(remtime.starttime)
                totalRem = [totalRem; remtime.total_duration];
            else
                totalRem = [totalRem; 0];
            end
            totalPhasic = [totalPhasic; 0];
            phasicrem{day}{epoch}.starttime = [];
            phasicrem{day}{epoch}.endtime = [];
            phasicrem{day}{epoch}.total_duration = 0;
            phasicrem{day}{epoch}.descrip = 'phasic rem extracted from mean of ctxriptets, Mizugeki method';
        end
        clear remlist
    end
    
    if savedata == 1
        save(sprintf('%s%sphasicrembouts%02d.mat', dir,animalprefix,day), 'phasicrem');
    end
    clear phasicrem
end
keyboard