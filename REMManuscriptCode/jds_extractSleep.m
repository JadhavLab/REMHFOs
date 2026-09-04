function wb_extractsleep2_withspike(animaldir, animalprefix, varargin)
%-------Reads in these files:
%  delta (gnd)
%  theta (gnd)
%  position
%  spindle (gnd)

%  animaldir - example '/data99/user/animaldatafolder/', a folder
%                containing processed matlab data for the animal
%
%  fileprefix	- folder name where the day's data is stored
%

eegdir = [animaldir,'EEG/'];

%-----------------------------------------------%
%----             initial setup         --------%
%-----------------------------------------------%
detect_tet = [];  % to process the data from a certain tetrode, default empty: process all
manual_day = [];  %  to process the data from a certain day, default empty: process all
sleep_waking_only = 0; % separate sleep and waking only

ctx_flag = 0;  % cortical data = 1; hippocampal = 0;
rem_flag = 1;  % detect REM
spindle_flag = 0;  % detect Spindle
sws_flag = 1;  % detect SWS

plot_TDratio = 1; % Plot result
smoothing_width = 1;    % SD of gaussian used to smooth filtered eeg magnitude, in s

velocity_thresh = 4;   % velocity threshhold for moving and immobility,  in cm/s
velocity_thresh_waking = 4; % velocity threshhold for moving and sleep,  in cm/s
time_immobile_waking = 7;  % minimal time for immobility, in s
time_immobile = 60;   % minimal time for sleep, in s


rem_thresh = [];     % Theta:Delta threshold for separating NREM with REM
spindle_thresh = 1;  %  Spindle:Delta threshold for separating Spindle with SWS

mindur_spindle = 0.5;  % minimal period for spindle, in s
mindur_rem = 10; % minimal period for REM, in s


mindur_sleep_extractsleep = [];
ignorelist = {};

%-----------------------------------------------%
%---        load customized parameters      ----%
%-----------------------------------------------%

for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'manual_day'
                manual_day = varargin{option+1};
            case 'detect_tet'
                detect_tet = varargin{option+1};
            case 'sleep_waking_only'
                sleep_waking_only = varargin{option+1};
            case 'rem'
                rem_flag = varargin{option+1};
            case 'spindle'
                spindle_flag = varargin{option+1};
            case 'sws'
                sws_flag = varargin{option+1};
            case 'cortex'
                ctx_flag = varargin{option+1};
            case 'immobility'
                time_immobile = varargin{option+1}(1);
                velocity_thresh = varargin{option+1}(2);
            case 'immobility_waking'
                time_immobile_waking = varargin{option+1}(1);
                velocity_thresh_waking = varargin{option+1}(2);
            case 'smoothing_width'
                smoothing_width = varargin{option+1};
            case 'plot_TDratio'
                plot_TDratio = varargin{option+1};
            case 'rem_thresh'
                rem_thresh = varargin{option+1};
            case 'spindle_max'
                spindle_max = varargin{option+1};
            case 'spindle_thresh'
                spindle_thresh = varargin{option+1};
            case 'mindur_rem'
                mindur_rem = varargin{option+1};
            case 'mindur_spindle'
                mindur_spindle = varargin{option+1};
            case 'ignorelist'
                mindur_sleep_extractsleep = varargin{option+1}{1}; % select to ignore some dataset
                ignorelist = varargin{option+1}{2};
                animnum = varargin{option+1}{3};
            otherwise
                error(['Option ',varargin{option},' unknown.']);  % errors for undefined parameters
        end
    else
        error('Options must be strings, followed by the variable');
    end
end
%%
%-----------------------------------------------%
%---          Sleep cycle detection         ----%
%-----------------------------------------------%
% load these data structures
tetinfo = loaddatastruct(animaldir, animalprefix, 'tetinfo');
days = find(~cellfun('isempty',tetinfo));
% if epoch manually specified, pick it
if ~isempty(manual_day)
    days = manual_day;
end

% iteration struction: day -- epoch -- tetrode
for day = 1%days % day iterations
    task = loaddatastruct(animaldir,animalprefix,'task',day); % get task type, only run for sleep task
    pos = loaddatastruct(animaldir, animalprefix, 'pos', day); % get animal [time, positions, velocity]

    if ctx_flag
        multi = loaddatastruct(animaldir, animalprefix, 'multi', day);
    end
    for ep = 1:length(pos{day}) %CHANGE
       
        %%
        %-----------------------------------------------%
        %---        Sleep, immobility and waking    ----%
        %-----------------------------------------------%
        sleepperiods = [];
        if ~strcmp(task{day}{ep}.type,'sleep') % only run for sleep task
            continue
        end
        if ~isempty(pos{day}{ep})
            if ep > length(pos{day})
                disp(sprintf('task on day %d has more epochs than pos, skipping',day))
                continue
            end
            if size(pos{day}{ep}.data,2) > 5  % get velocity
                velocity = pos{day}{ep}.data(:,9);% smoothed velocity
            else
                velocity = pos{day}{ep}.data(:,5);
            end
            postimes = pos{day}{ep}.data(:,1);  % get time stamps

            %**********************sleep period*******************%
            immobile = vec2list(velocity < velocity_thresh,postimes); % generate [start end] list of immobile epochs
            immobile = [immobile(:,1)+time_immobile, immobile(:,2)]; % at least immobile for 60 seconds
            sleepperiods = immobile(find( (immobile(:,2)-immobile(:,1)) > 0),:); % v < velocity_thresh and last for at least time_immobile s is sleep
            immobile_spindle = [immobile(:,1)+5, immobile(:,2)]; % v < velocity_thresh and last for at least time_immobile+5 s is sleep
            sleepperiods_spindle = immobile_spindle(find( (immobile_spindle(:,2)-immobile_spindle(:,1)) > 0),:);
            % choose to ignore some customized periods
            if ~isempty(ignorelist)
                % identify valid periods
                sp_valid = [];
                % 1. exceed mindur_sleep_extractsleep
                spnums = find( (sleepperiods(:,2) - sleepperiods(:,1))  > mindur_sleep_extractsleep)';
                for sp = spnums
                    if rowfind([animnum day ep],ignorelist{animnum}) == 0
                        % 2. are not on the ignorelist
                        sp_valid = [sp_valid  sp];  % add
                    end
                end
                disp(sprintf('ignorelist:  %d of %d sleep periods pass test',length(sp_valid),size(sleepperiods,1)))
                sleepperiods = sleepperiods(sp_valid,:);
            end

            % display the amount of sleep periods
            sleepduration = round(sum(sleepperiods(:,2)-sleepperiods(:,1)));
            epochduration = round(pos{day}{ep}.data(end,1) - pos{day}{ep}.data(1,1));
            disp(sprintf('d %d e %d: %d s of sleep (%d s total)',day,ep,sleepduration,epochduration));

            % fit into save file
            if ~isempty(sleepperiods)
                slp.starttime = sleepperiods(:,1);
                slp.endtime = sleepperiods(:,2);
            else
                slp.starttime = [];
                slp.endtime = [];
            end

            if ~isempty(pos{day}{ep})
                slp.timerange = [pos{day}{ep}.data(1,1) pos{day}{ep}.data(end,1)];
            else
                disp(sprintf('no position data for this epoch! assigning [0 0] for sleep.timerange d%de%d',day,ep))
                slp.timerange = [ 0 0 ];
            end
            slp.time_immobile = time_immobile;    % in seconds
            slp.velocity_thresh = velocity_thresh;
            % install output in proper place
            sleep{day}{ep} = slp;
            clear slp;


            %*********************waking period******************%
            %----immobility----%
            % v < velocity_thresh_waking for at least time_immobile_waking
            immobile_w = vec2list(velocity < velocity_thresh_waking,postimes);
            longinds = (immobile_w(:,2) - immobile_w(:,1)) > time_immobile_waking ;
            nonwakingperiods = [immobile_w(longinds,1) + time_immobile_waking    immobile_w(longinds,2)]; % non-waking periods = immobile period
            %----waking----%
            wakingvec = ~list2vec(nonwakingperiods,postimes);% waking = non-immobile
            wakingperiods = vec2list(wakingvec,postimes);
            nwakingvec = list2vec(nonwakingperiods,postimes);
            imobw_vec = (nwakingvec) &...
                ~list2vec(sleepperiods,postimes);
            imobwperiod = vec2list(imobw_vec,postimes);

            % display the amount of waking periods
            wakingduration = round(sum(wakingperiods(:,2)-wakingperiods(:,1)));
            epochduration = round(pos{day}{ep}.data(end,1) - pos{day}{ep}.data(1,1));
            disp(sprintf('d %d e %d: %d s of waking (%d s total)',day,ep,wakingduration,epochduration));

            % fit into save file
            if ~isempty(wakingperiods)
                wak.starttime = wakingperiods(:,1);
                wak.endtime = wakingperiods(:,2);
            else
                wak.starttime = [];
                wak.endtime = [];
            end

            if ~isempty(pos{day}{ep})
                wak.timerange = [pos{day}{ep}.data(1,1) pos{day}{ep}.data(end,1)];
            else
                disp(sprintf('no position data for this epoch! assigning [0 0] for sleep.timerange d%de%d',day,ep))
                wak.timerange = [ 0 0 ];
            end
            wak.time_immobile_waking = time_immobile_waking;    % in seconds
            wak.velocity_thresh_waking = velocity_thresh_waking;
            % install output in proper place
            waking{day}{ep} = wak;
            clear wak;

            if ~isempty(imobwperiod)
                nwak.starttime = imobwperiod(:,1);
                nwak.endtime = imobwperiod(:,2);
            else
                nwak.starttime = [];
                nwak.endtime = [];
            end

            if ~isempty(pos{day}{ep})
                nwak.timerange = [pos{day}{ep}.data(1,1) pos{day}{ep}.data(end,1)];
            else
                disp(sprintf('no position data for this epoch! assigning [0 0] for waking immobility.timerange d%de%d',day,ep))
                nwak.timerange = [ 0 0 ];
            end
            nwak.time_immobile_waking = time_immobile_waking;    % in seconds
            nwak.velocity_thresh_waking = velocity_thresh_waking;
            nwaking{day}{ep} = nwak;
            clear nwak;

            if sleep_waking_only  % only need sleep and waking periods? okay, the end.
                continue
            end

            %%
            %-----------------------------------------------%
            %----    REM versus NREM (Spindle +SWS) --------%
            %-----------------------------------------------%
            % if no preselected tetrode, determine the ripple tets
            if isempty(detect_tet)
                detect_tet = [];
                for ttt = 1:length(tetinfo{day}{ep})
                    if ~isempty(tetinfo{day}{ep}{ttt})
                        if isfield(tetinfo{day}{ep}{ttt},'descrip')
                            if strcmp(tetinfo{day}{ep}{ttt}.descrip,'riptet')
                                detect_tet = [detect_tet ttt];
                            end
                        end
                    end
                end
            end

            TD = [];
            SD = [];
            times_filteeg = [];
            tmpflist1 = sprintf('%s%stheta%02d-%02d-%02d.mat', eegdir,animalprefix, day, ep, detect_tet(1));
            load(tmpflist1);
            times_filteeg = geteegtimes(theta{day}{ep}{detect_tet(1)}) ;
            times_filteeg = times_filteeg(:)';
            total_neurons = 0;
            for tet = detect_tet  

                % load EEG-gnd data
                tmpflist1 = sprintf('%s%stheta%02d-%02d-%02d.mat', eegdir,animalprefix, day, ep, tet);
                tmpflist2 = sprintf('%s%sdelta%02d-%02d-%02d.mat', eegdir,animalprefix, day, ep, tet);
                tmpflist3 = sprintf('%s%sspindlegnd%02d-%02d-%02d.mat', eegdir,animalprefix, day, ep, tet);
                load(tmpflist1);
                load(tmpflist2);
                load(tmpflist3);
                
                % convert the envelope fields to double
                tenv = double(theta{day}{ep}{tet}.data(:,3));
                denv = double(delta{day}{ep}{tet}.data(:,3));
                spenv = double(spindlegnd{day}{ep}{tet}.data(:,3));

                % verify match the Fs between filtered eeg
                fs1 = theta{day}{ep}{tet}.samprate; % most files downsample to 150 Hz
                fs2 = delta{day}{ep}{tet}.samprate;
                fs3 = spindlegnd{day}{ep}{tet}.samprate;
                if (fs1 ~= fs2) || (fs1 ~= fs3) || (fs2 ~= fs3)
                    error('sampling rates for filtered eegs need to match')
                end

                % smooth the magnitude envelope with a 1s Gaussian kernal
                samprate = theta{day}{ep}{tet}.samprate;
                kernel = gaussian(smoothing_width*samprate, ceil(8*smoothing_width*samprate));
                tenv = smoothvect(tenv, kernel);
                tenv = tenv(:)';
                denv = smoothvect(denv, kernel);
                denv = denv(:)';
                spenv = smoothvect(spenv, kernel);
                spenv = spenv(:)';

                % calculate Theta:Delta ratio and Spindle:Delta ratio
                tdratio = tenv ./ denv ;
                sdratio = spenv ./ denv ;

                % check that eeg times vectors (if multiple tetrodes) are
                % equivalent; use the same time vector for all tetrode
                if isempty(times_filteeg)
                    times_filteeg = geteegtimes(theta{day}{ep}{tet});
                    TD = [TD ; tdratio];
                    SD = [SD ; sdratio];
                else
                    times_filteeg2 = geteegtimes(theta{day}{ep}{tet});
                    times_filteeg2 = times_filteeg2(:)';
                    % now check the times vectors are matched
                    if all(times_filteeg2 == times_filteeg)
                        if length(tdratio) - size(TD,2) == 1   % correction in error for Bond
                            tdratio(end) = [];
                            sdratio(end) = [];
                            sleepvec(end) = [];
                        end
                        TD = [TD ; tdratio];
                        SD = [SD ; sdratio];
                    else
                        %if not, interpolate
                        TD = [TD ; interp1(times_filteeg2,tdratio,times_filteeg,'nearest')];
                        SD = [SD ; interp1(times_filteeg2,sdratio,times_filteeg,'nearest')];
                    end
                end
            end
            % Take the mean over the detection tetrodes
            TD_mean = mean(TD,1);
            if isempty(rem_thresh)
                %                 TD_sort = sort(TD_mean);
                %                 TD_percent = round(length(TD_mean).*0.90);
                %                 rem_thresh = TD_sort(TD_percent);
                rem_thresh = mean(TD_mean)+ 1*std(TD_mean);
            end
            SD_mean = mean(SD,1);
            indTD = find(isnan(TD_mean) == 1);
            TD_mean(indTD) = TD_mean(indTD-1);
            indSD = find(isnan(SD_mean) == 1);
            SD_mean(indSD) = SD_mean(indSD-1);
            % transform sleepperiods list into 0-1 vector, 1 sleep, 0 not
            % sleep
            [stime,sleepvec] = wb_list2vec(sleepperiods,times_filteeg);
            [spindletime,sleepvec_spindle] = wb_list2vec(sleepperiods_spindle,times_filteeg);
            [imtime,imvec] = wb_list2vec(imobwperiod,times_filteeg);

            if length(sleepvec) - size(TD,2) == 1   % correction for error in Bond
                sleepvec(end) = [];
                sleepvec_spindle(end) = [];
            end

            % rem_vec = (sleepvec) & (TD_mean > rem_thresh) & (SD_mean < spindle_max);
            % if you know how to define spindle_max, use the above one
            rem_vec = (sleepvec) & (TD_mean > rem_thresh);
            spindle_vec = (sleepvec_spindle) &...
                (SD_mean > spindle_thresh) &...
                ~rem_vec ;   % notice that spindles must fall outside of REM times
            remlist = vec2list(rem_vec,times_filteeg);
            spindlelist = vec2list(spindle_vec,times_filteeg);

            % eliminate periods that are shorter than minimum durations
            remlist(( (remlist(:,2)-remlist(:,1)) < mindur_rem),:) = []; % the REM period should be at least 10-s long
            [remtime,rem_vec] = wb_list2vec(remlist,times_filteeg);
            nrem_vec = (sleepvec) & ~rem_vec; % sleep = REM + NREM
            nremlist = vec2list(nrem_vec,times_filteeg);

            spindlelist(( (spindlelist(:,2)-spindlelist(:,1)) < mindur_spindle),:) = [];% the Spindle period should be at least 0.5-s long
            spindle_vec = list2vec(spindlelist,times_filteeg);
            if spindle_flag
                sws_vec = nrem_vec & spindle_vec';   %NREM = Spindle+ SWS
            else
                sws_vec = nrem_vec;
            end
            swslist = vec2list(sws_vec,times_filteeg);

            % fit REM into save file
            if rem_flag
                if ~isempty(remlist)
                    r.starttime = remlist(:,1);
                    r.endtime = remlist(:,2);
                    r.TD_mean = TD_mean;
                    r.timevec = times_filteeg;
                    remduration = round(sum(remlist(:,2)-remlist(:,1)));
                    r.total_duration = remduration;
                else
                    r.starttime = [];
                    r.endtime = [];
                    r.td_mean = [];
                    remduration = 0;
                end
                r.detect_tet = detect_tet;
                r.timerange = [0 length(tenv)/fs1] + theta{day}{ep}{tet}.starttime;
                r.samprate = theta{day}{ep}{tet}.samprate;
                r.rem_thresh = rem_thresh;
                r.baseline = 'relying on static threshold a la Buzsaki';
                r.std = 'relying on static threshold a la Buzsaki';
                r.mindur_rem = mindur_rem;
                % install output in proper place
                rem{day}{ep} = r;
                % print total duration
                disp(sprintf('d %d e %d t %d: %d s of REM',day,ep,tet,remduration))
                clear r;
            end
            % fit Spindle into save file
            if spindle_flag
                if ~isempty(spindlelist)
                    sp.starttime = spindlelist(:,1);
                    sp.endtime = spindlelist(:,2);
                    sp.TD_mean = TD_mean;
                    sp.timevec = times_filteeg;
                    spindleduration = round(sum(spindlelist(:,2)-spindlelist(:,1)));
                    sp.total_duration = spindleduration;
                else
                    sp.starttime = [];
                    sp.endtime = [];
                    sp.td_mean = [];
                    spindleduration = 0;
                end
                sp.detect_tet = detect_tet;
                sp.timerange = [0 length(tenv)/fs1] + theta{day}{ep}{tet}.starttime;
                sp.samprate = theta{day}{ep}{tet}.samprate;
                sp.spindle_thresh = spindle_thresh;
                sp.baseline = 'relying on static threshold a la Buzsaki';
                sp.std = 'relying on static threshold a la Buzsaki';
                sp.mindur_spindle = mindur_spindle;
                % install output in proper place
                spindles{day}{ep} = sp;
                % print total duration
                disp(sprintf('d %d e %d t %d: %d s of spindle',day,ep,tet,spindleduration))
                clear sp;
            end
            if sws_flag
                if ~isempty(swslist)
                    sw.starttime = swslist(:,1);
                    sw.endtime = swslist(:,2);
                    sw.TD_mean = TD_mean;
                    sw.timevec = times_filteeg;
                    swsduration = round(sum(swslist(:,2)-swslist(:,1)));
                    sw.total_duration = swsduration;
                else
                    sw.starttime = [];
                    sw.endtime = [];
                    sw.td_mean = [];
                    swsduration = 0;
                end
                sw.detect_tet = detect_tet;
                sw.timerange = [0 length(tenv)/fs1] + theta{day}{ep}{tet}.starttime;
                sw.samprate = theta{day}{ep}{tet}.samprate;
                sw.spindle_thresh = spindle_thresh;
                sw.baseline = 'relying on static threshold a la Buzsaki';
                sw.std = 'relying on static threshold a la Buzsaki';
                sw.mindur_sws = 0;
                % install output in proper place
                sws{day}{ep} = sw;
                % print total duration
                disp(sprintf('d %d e %d t %d: %d s of sws',day,ep,tet,swsduration))
                clear sw;
            end

            %%
            %-----------------------------------------------%
            %---------------   Plot Results ----------------%
            %-----------------------------------------------%
            if plot_TDratio

                % time plot of ratio trace
                xvec = times_filteeg-times_filteeg(1);
                stime = stime-times_filteeg(1);
                imtime = imtime-times_filteeg(1);
                remtime = remtime -times_filteeg(1);

                if length(xvec) - size(TD,2) == 1   % correction for error in Bond
                    xvec(end) = [];
                end

                whitebg([0 0 0])
                % fill the non-period bin with NaN
                svec = sleepvec;             % sleep
                svec(svec == 0) = nan;
                rvec = rem_vec;              % rem
                rvec(rvec==0) = nan;
                spvec = spindle_vec;         % spindle
                spvec(spvec==0) = nan;

                if length(rvec) - size(TD,2) == 1   % correction for error in Bond
                    rvec(end) = [];
                    spvec(end) = [];
                end

                figure;
                %***********velocity, TD ratio and periods*********%
                subplot(2,1,1);
                minTD = min(TD_mean);
                maxTD = max(TD_mean);
                minRect = minTD-3*std(TD_mean);
                maxRect = maxTD + 3*std(TD_mean);

                for i =1:size(sleepperiods,1),
                    rectangle('Position',[stime(i,1),minRect,...
                        (stime(i,2)-stime(i,1)),...
                        maxRect-minRect],'FaceColor',[.7 .7 .7],...
                        'EdgeColor','none')
                    hold on
                end

                for i =1:size(remlist,1),
                    rectangle('Position',[remtime(i,1),minRect,...
                        remtime(i,2)-remtime(i,1),maxRect-minRect],...
                        'FaceColor',[.8 .8 1],'EdgeColor','none')
                    hold on
                end

                %                 for i =1:size(imobwperiod,1),
                %                     rectangle('Position',[imtime(i,1),minRect,...
                %                         imtime(i,2)-imtime(i,1),maxRect-minRect],...
                %                         'FaceColor',[1,.5,.5],'EdgeColor','none')
                %                     hold on
                %                 end
                %                 ntet = 0;
                %                 nscale_n1 = 0;
                %                 for tet = detect_tet  % tetrode iterations
                %                     if max(length(spks{day}{ep})) >= tet
                %                         spksvec = spikvec{tet};
                %                         if  strcmp(tetinfo{day}{ep}{tet}.area,'PFC')
                %                             if ~isempty(spksvec)
                %                                 ntet = ntet+1;
                %                                 plot(spksvec(:,1)-times_filteeg(1),(spksvec(:,2)+nscale_n1)./total_neurons*maxRect,'r.')
                %                                 nscale_n1 = nscale_n1 + max(spksvec(:,2));
                %                             end
                %                         else
                %                             if ~isempty(spksvec)
                %                                 ntet = ntet+1;
                %                                 plot(spksvec(:,1)-times_filteeg(1),(spksvec(:,2)+nscale_n1)./total_neurons*maxRect,'w.')
                %                                 nscale_n1 = nscale_n1 + max(spksvec(:,2));
                %                             end
                %                         end
                %                     end
                %                 end
                % velocity is scaled for visualization purpose
                plot(postimes-times_filteeg(1),velocity./5+6,'Color','w','linewidth',1); hold on
                plot(postimes-times_filteeg(1),velocity_thresh./5+6.*ones(1,length(postimes)),'w--')
                plot(xvec,TD_mean,'Color',[.5 .5 1],'linewidth',2); hold on
                plot(xvec,rem_thresh.*ones(1,length(xvec)),'b--','linewidth',1); hold on
                axis tight;  ylim([0 10]);
                title(sprintf('%s  day: %d epoch: %d',animalprefix,day,ep),'fontsize',20,'fontweight','bold')
                xlabel('Time (s)','fontsize',15,'fontweight','bold')
                %***********histogram of TD ratio (on left side)*********%
                subplot(2,2,3)
                [N X] = hist(TD_mean,0:.1:5); h = bar(X,N); xlim([0 5])
                set(h(1),'facecolor',[.5 .5 1]); set(h(1),'edgecolor','none');
                patch=findobj(h,'Type','patch');
                hold on
                plot(rem_thresh.*ones(1,max(N)),1:max(N),'b--','linewidth',2)
                title('Theta:Delta Threshold','fontsize',20,'fontweight','bold')
                %***********histogram of SD ratio (on right side)*********%
                subplot(2,2,4)
                [N X] = hist(SD_mean,0:.1:5); h = bar(X,N); xlim([0 5])
                set(h(1),'facecolor',[.5 1 .5 ]); set(h(1),'edgecolor','none');
                patch=findobj(h,'Type','patch');
                hold on
                plot(spindle_thresh.*ones(1,max(N)),1:max(N),'g--','linewidth',2)
                title('Spindle:Delta Threshold','fontsize',20,'fontweight','bold')
                figfile = (sprintf('%s/%s_d%d_e%d.fig',animaldir,animalprefix,day,ep));
                %                 saveas(gcf,figfile,'fig');
            end
        end
    end
end
%     %%
%     %-----------------------------------------------%
%     %-----        Save results for each day     ----%
%     %-----------------------------------------------%
% save sleep
save(sprintf('%s/%ssleep%02d.mat', animaldir, animalprefix, day), 'sleep');
clear sleep

% save waking
save(sprintf('%s/%swaking%02d.mat', animaldir, animalprefix, day), 'waking');
clear waking

% save waking and immobility
save(sprintf('%s/%snonwaking%02d.mat', animaldir, animalprefix, day), 'nwaking');
clear nwaking

if ~sleep_waking_only
    % save REM
    if rem_flag && exist('rem','var')
        disp(sprintf('saving rem.. (rem_thresh = %d, mindur_rem = %d, smoothing_width = %d',...
            rem_thresh,mindur_rem,smoothing_width))
        save(sprintf('%s/%srem%02d.mat', animaldir, animalprefix, day), 'rem');
        clear rem
    end
    % save Spindle
    if spindle_flag
        disp(sprintf('saving spindles.. (spindle_thresh = %d, mindur_spindle = %d, smoothing_width = %d',...
            spindle_thresh,mindur_spindle,smoothing_width))
        save(sprintf('%s/%sspindles%02d.mat', animaldir, animalprefix, day), 'spindles');
        clear spindles
    end
    % save SWS
    if sws_flag
        disp(sprintf('saving sws.. (sws_thresh = %d, mindur_sws = %d, smoothing_width = %d',...
            rem_thresh,0,smoothing_width))
        save(sprintf('%s/%ssws%02d.mat', animaldir, animalprefix, day), 'sws');
        clear sws
    end
end
cd(animaldir)

end




