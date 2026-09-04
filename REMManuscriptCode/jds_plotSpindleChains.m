function jds_plotSpindleChains(animalprefixlist)

%function extractspindles(directoryname, fileprefix, day, tetrode_list, %%min_suprathresh_duration, nstd, maxpeakval, options)
%
%	Reads in the spindle files from the specified day and tetrodes and
%	extracts all of the spindle from that tetrodes.
%
%	assumes position data stores in pos file in animdirectory
%
%directoryname - example '/data99/user/animaldatafolder/', a folder
%                containing processed matlab data for the animal
%
%fileprefix	- folder name where the day's data is stored
%
%day		- the day to process
%
%tetrode_list	- the tetrode(s) to process.
%			-1 indicates all tetrodes should be processed
%			 0 indicates that the tetrode with the most cells will
%			 	be processed
%			 #(s) indicate the set of tetrodes to proces
%min_suprathresh_duration
%		- the time (in seconds) which the signal
%       must remain above threshold to be counted as as spindle; this guards
%       against short spikes in signal (e.g. noise) as being counted as
%       spindles. Set min_suprathreshold_duration to some small value, like
%       0.015 s.
%
%nstd		- the number of standard dev that spindle must be from mean to
%			be detected. Start with 2.
%
%
%
% Outputs:
%spindles 	- structue with various fields, including the following which
%			describe each ripple.
%	starttime - time of beginning of spindle
%	endtime	  - time of end of spindle
%	midtime   - time of midpoint of energy of event
%	peak	  - peak height of waveform)
%	maxthresh - the largest threshold in stdev units at which this spindle
%			would still be detected.
%	energy	  - total sum squared energy of waveform
%	startind  - index of start time in spindle structure
%	endind    - index of end time in spindle structure
%	midind    - index of middle time in spindle structure
%	posind    - index into pos structure for the midpoint of this ripple
% 	posinterp - interpolation factor value between adjacent position
% 		    elements
epochs = 1:2:17;
day = 1;

for a = 1:length(animalprefixlist)
    fileprefix = animalprefixlist{a};
    directoryname = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', fileprefix);

    load(sprintf('%s%stetinfo.mat',directoryname,fileprefix));
    load(sprintf('%s%sremeps0%d.mat',directoryname,fileprefix,day));
    tetrode = evaluatefilter(tetinfo{1}{1},'(isequal($descrip, ''ctxriptet''))');
    load(sprintf('%s%sctxspindletime_chainSWS0%d.mat',directoryname,fileprefix,day));
    % define the standard deviation for the Gaussian smoother which we
    % apply before thresholding (this reduces sensitivity to spurious
    % flucutations in the spindle envelope)
    smoothing_width = 0.004; % 4 ms

    d = day;

    if tetrode == 0
        %select eeg with most cells
        numcells = findnumcell(directoryname, fileprefix, d, 2);
        [i tet] = max(numcells) ;
    else
        tet = tetrode;
    end

    dsz = '';
    if (d < 10)
        dsz = '0';
    end

    % load the positionf file
    eval(['load ', directoryname, '/', fileprefix, 'pos', dsz, num2str(d), '.mat']); %load pos file

    % move to the EEG directory
    pushd([directoryname,'/EEG/']);
    if (tet == -1)
        % set up to handle all tetrodes
        tet = 1:1000;
    end
    sigma = 2.5;

    % Create a Gaussian smoothing window
    sz = 48 * sigma; % Size of the window (adjustable)
    if mod(sz, 2) == 0
        sz = sz + 1; % Ensure odd size for proper centering
    end

%     kernel = fspecial('gaussian', [sz, sz], sigma);
    kernel = gausswin(48, 2.5);
    for eps = 1:length(remeps)
        allrenv = [];
        e = remeps(eps);
        % go through each tetrode
        for t = 1:length(tet)

            trode = tet(t);
            load(sprintf('%seeg%02d-%02d-%02d.mat',fileprefix, day, e, trode));
            % convert the spindle envelope field to double
%             temprenv = double(spindlegnd{d}{e}{trode}.data(:,3));
            temprenv = double(eeg{d}{e}{trode}.data);
            samprate = 1500;
%             kernel = gaussian(smoothing_width*samprate, ceil(8*smoothing_width*samprate));
%             temprenv = smoothvect(temprenv, kernel);
            allrenv = [allrenv; temprenv'];
        end
        times_filteeg = geteegtimes(eeg{d}{e}{trode});
        [b,a] = butter(3,[10/(samprate/2) 16/(samprate/2)]);
        v = filtfilt(b,a,mean(allrenv));
        v = zscore(v);
%         for s = 1:length(ctxspindle{d}{e}.C_sep)
%             tmp = ctxspindle{d}{e}.C_sep{s};
%             if size(tmp,1) == 3
%                 first_st = lookup(tmp(1,1), times_filteeg);
%                 last_end = lookup(tmp(end,2), times_filteeg);
%                 figure; hold on
%                 plot(first_st-750:last_end+750, (v(first_st-750:last_end+750)));
%                 for ss = 1:size(tmp,1)
%                     idx = lookup(tmp(ss,:),times_filteeg);
%                     plot(idx(1):idx(2), (v(idx(1):idx(2))),'-r');
%                 end
%                 x = [first_st first_st+1500];
%                 y = [3 3];
%                 plot(x,y,'k','LineWidth',3)
%                 xlim([first_st-750 last_end+750 ]);
%                 ylim([-6 6]);
%                 keyboard
%             end
%         end

        for s = 1:length(ctxspindle{d}{e}.starttimeNC)
            tmp = [ctxspindle{d}{e}.starttimeNC(s) ctxspindle{d}{e}.endtimeNC(s)];
            first_st = lookup(tmp(1,1), times_filteeg);
            last_end = lookup(tmp(1,2), times_filteeg);
            figure; hold on
            plot(first_st-2500:last_end+2500, (v(first_st-2500:last_end+2500)));
            for ss = 1:size(tmp,1)
                idx = lookup(tmp(ss,:),times_filteeg);
                plot(idx(1):idx(2), (v(idx(1):idx(2))),'-r');
            end
            x = [first_st first_st+1500];
            y = [3 3];
            plot(x,y,'k','LineWidth',3)
            xlim([first_st-2500 last_end+2500 ]);
            ylim([-6 6]);
            keyboard
            close
        end
    end
end


