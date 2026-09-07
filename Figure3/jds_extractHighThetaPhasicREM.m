function jds_extractHighThetaPhasicREM(animalprefixlist)
%JDS_EXTRACTHIGHTHETAPHASICREM Extract high-theta (phasic) REM bouts.
%   jds_extractHighThetaPhasicREM(animalprefixlist) averages the
%   theta-filtered LFP envelope across all ctxriptet tetrodes, keeps
%   periods more than 4 SD above the mean REM theta envelope that last
%   longer than 900 ms and fall within REM, and saves them to
%   <animal>phasicrembouts<day>.mat.
%
%   animalprefixlist - cell array of animal prefix strings
%

day = 2;
epochs = [2];
daystring = sprintf('%02d',day);
savedata = 1;
W = rectwin(11);
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
        envdataall = [];
        epoch = epochs(e);

        remtime = rem{day}{epoch};
        remlist = [remtime.starttime remtime.endtime];
        
        if epoch <10
            epochstring = ['0',num2str(epoch)];
        else
            epochstring = num2str(epoch);
        end

        for i = 1:length(ctxtets)
            ctxtet = ctxtets(i);
            
            if (ctxtet<10)
                ctxtetstring = ['0',num2str(ctxtet)];
            else
                ctxtetstring = num2str(ctxtet);
            end
            
            curreegfile = [dir,'/EEG/',animalprefix,'theta', daystring,'-',epochstring,'-',ctxtetstring];
            load(curreegfile);
            
            envdatatmp = theta{day}{epoch}{ctxtet}.data(:,3);
            envdataall = [envdataall; envdatatmp'];
        end
        times = geteegtimes(theta{day}{epoch}{ctxtet}) ; % construct time array

        envdata = mean(envdataall,1); %mean envelope across all tetrodes
        envdata = conv(envdata,W,'same');
        %Get the mean env of theta during REM here for later
        %thresholding

        inRemAmp = logical(isExcluded(times,remlist));
        
        meanRemEnv = mean(envdata(inRemAmp));

        highThet = envdata > (meanRemEnv + 4*std(envdata(inRemAmp)));

        evList = vec2list(highThet,times);

        longIdx = find(evList(:,2)-evList(:,1) > 0.9);
        finalEpochs = evList(longIdx,:); %greater than 900 ms

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
