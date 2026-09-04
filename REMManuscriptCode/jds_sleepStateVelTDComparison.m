function jds_sleepStateVelTDComparison(animalprefixlist, varargin)
% Compare velocity and theta/delta ratio (TD) across NREM and REM sleep.
%
% INPUT
%   animalprefixlist : cell array of animal ID strings
%   varargin         : optional name-value pairs
%
% The function produces:
%   • Box-plots of velocity and TD (z-scored) in NREM vs. REM
%   • Mean ± SEM traces of TD around NREM→REM and REM→NREM transitions
%
% -------------------------------------------------------------------------

% Accumulators across animals
nremVel      = [];  remVel      = [];
nremTD       = [];  remTD       = [];
nrem2remTD   = [];  rem2nremTD  = [];

for a = 1:numel(animalprefixlist)
    %% ---------- Defaults (override with varargin) ----------
    animalprefix           = animalprefixlist{a};
    animaldir              = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/', animalprefix);
    eegdir                = fullfile(animaldir, 'EEG');

    detect_tet            = [];   % tetrodes to use for TD (empty → auto-detect rip-tets)
    manual_day            = [];   % restrict to specific day(s) if desired
    smoothing_width       = 1;    % Gaussian σ (s) for envelope smoothing
    velocity_thresh       = 4;    % immobility threshold (cm/s)
    velocity_thresh_wake  = 4;    % immobility threshold within wake (cm/s)
    time_immobile_wake    = 7;    % min immobility to count as non-wake (s)
    time_immobile_sleep   = 60;   % min immobility to count as sleep (s)
    rem_thresh            = [];   % TD threshold for REM (auto if empty)
    mindur_rem            = 10;   % minimum REM bout (s)

    %% ---------- Parse optional parameters ----------
    for k = 1:2:numel(varargin)
        opt = lower(varargin{k});
        val = varargin{k+1};

        switch opt
            case 'manual_day',           manual_day           = val;
            case 'detect_tet',           detect_tet           = val;
            case 'immobility'
                time_immobile_sleep     = val(1);
                velocity_thresh         = val(2);
            case 'immobility_waking'
                time_immobile_wake      = val(1);
                velocity_thresh_wake    = val(2);
            case 'smoothing_width',      smoothing_width      = val;
            case 'rem_thresh',           rem_thresh           = val;
            case 'mindur_rem',           mindur_rem           = val;
            case 'ignorelist'
               
            otherwise
                error('Unknown option "%s".', varargin{k});
        end
    end

    %% ---------- Load tetrode info & loop over days ----------
    tetinfo = loaddatastruct(animaldir, animalprefix, 'tetinfo');
    days    = find(~cellfun('isempty', tetinfo));
    if ~isempty(manual_day), days = manual_day; end

    for day = days
        load(fullfile(animaldir, sprintf('%sremeps%02d.mat', animalprefix, day)), 'remeps');
        task = loaddatastruct(animaldir, animalprefix, 'task', day);
        pos  = loaddatastruct(animaldir, animalprefix, 'pos',  day);

        for ep = remeps
            if ~strcmp(task{day}{ep}.type, 'sleep') || isempty(pos{day}{ep}), continue; end

            %% ---------- Basic position/velocity vectors ----------
            posdata   = pos{day}{ep}.data;
            velocity  = posdata(:, size(posdata,2) > 5 + 4);  % smoothed velocity column (either 9 or 5)
            postimes  = posdata(:,1);

            %% ---------- Identify sleep and wake periods ----------
            sleepperiods = vec2list(velocity < velocity_thresh, postimes);
            sleepperiods = [sleepperiods(:,1) + time_immobile_sleep, sleepperiods(:,2)];
            sleepperiods = sleepperiods((sleepperiods(:,2)-sleepperiods(:,1)) > 0, :);

            % Wake immobility → exclude from wake
            imm_wake      = vec2list(velocity < velocity_thresh_wake, postimes);
            long_imm      = (imm_wake(:,2) - imm_wake(:,1)) > time_immobile_wake;
            nonwake       = [imm_wake(long_imm,1) + time_immobile_wake, imm_wake(long_imm,2)];
            waking_vec    = ~list2vec(nonwake, postimes);
            imobw_periods = vec2list(list2vec(nonwake, postimes) & ~list2vec(sleepperiods, postimes), postimes);

            %% ---------- Auto-detect rip-tetrodes if none specified ----------
            if isempty(detect_tet)
                for t = 1:numel(tetinfo{day}{ep})
                    info = tetinfo{day}{ep}{t};
                    if ~isempty(info) && isfield(info, 'descrip') && strcmp(info.descrip, 'riptet')
                        detect_tet(end+1) = t;
                    end
                end
            end
            if isempty(detect_tet), warning('No riptet found for day %d ep %d.', day, ep); continue; end

            %% ---------- Compute TD ratio across chosen tetrodes ----------
            TD_all     = [];
            eeg_times  = [];

            for tet = detect_tet
                load(fullfile(eegdir, sprintf('%stheta%02d-%02d-%02d.mat',  animalprefix, day, ep, tet)),  'theta');
                load(fullfile(eegdir, sprintf('%sdelta%02d-%02d-%02d.mat',  animalprefix, day, ep, tet)),  'delta');

                tenv = smoothvect(double(theta{day}{ep}{tet}.data(:,3)), ...
                                   gaussian(smoothing_width*theta{day}{ep}{tet}.samprate, ...
                                   ceil(8*smoothing_width*theta{day}{ep}{tet}.samprate)));
                denv = smoothvect(double(delta{day}{ep}{tet}.data(:,3)), ...
                                   gaussian(smoothing_width*delta{day}{ep}{tet}.samprate, ...
                                   ceil(8*smoothing_width*delta{day}{ep}{tet}.samprate)));

                tdratio = tenv ./ denv;

                this_times = geteegtimes(theta{day}{ep}{tet})';
                if isempty(eeg_times)
                    eeg_times = this_times;
                    TD_all    = tdratio;
                else
                    % interpolate if slightly mis-aligned
                    TD_all = [TD_all; interp1(this_times, tdratio, eeg_times, 'nearest', 'extrap')];
                end
            end

            TD_mean = mean(TD_all, 1, 'omitnan');
            TD_z    = zscore(TD_mean);

            if isempty(rem_thresh), rem_thresh = mean(TD_mean) + std(TD_mean); end

            %% ---------- Split REM vs. NREM ----------
            sleep_vec  = list2vec(sleepperiods,   eeg_times);
            rem_vec    = sleep_vec & (TD_mean > rem_thresh);
            rem_list   = vec2list(rem_vec, eeg_times);
            rem_list   = rem_list((rem_list(:,2)-rem_list(:,1)) >= mindur_rem, :);

            nrem_vec   = sleep_vec & ~list2vec(rem_list, eeg_times);
            nrem_list  = vec2list(nrem_vec, eeg_times);

            %% ---------- Gather metrics ----------
            for r = 1:size(rem_list, 1)
                st = lookup(rem_list(r,1), eeg_times);
                en = lookup(rem_list(r,2), eeg_times);

                % TD around transition: NREM→REM
                nrem2remTD = [nrem2remTD; TD_z(st-3000 : st+7500)]; 
                % TD values within REM
                remTD      = [remTD; mean(TD_z(st:en))];          

                % Velocity
                sv = lookup(rem_list(r,1), postimes);
                ev = lookup(rem_list(r,2), postimes);
                remVel     = [remVel;  mean(velocity(sv:ev), 'omitnan')]; 
            end

            for n = 1:size(nrem_list, 1)
                st = lookup(nrem_list(n,1), eeg_times);
                en = lookup(nrem_list(n,2), eeg_times);

                % TD around transition: REM→NREM
                rem2nremTD = [rem2nremTD; TD_z(en-7500 : en+3000)];
                % TD values within NREM
                nremTD     = [nremTD; mean(TD_z(st:en))];         

                % Velocity
                sv = lookup(nrem_list(n,1), postimes);
                ev = lookup(nrem_list(n,2), postimes);
                nremVel    = [nremVel; mean(velocity(sv:ev), 'omitnan')]; 
            end
        end
    end
end

%% ---------- Plot TD around transitions ----------
figure; hold on; set(gca,'FontSize',14);
t1 = -3000:7500;
boundedline(t1, mean(nrem2remTD), std(nrem2remTD)./sqrt(size(nrem2remTD,1)), '-k');
plot(t1, mean(nrem2remTD), 'k', 'LineWidth',1);
xticks([-3000 0 7500]); xticklabels({'-2','0','5'});
title('NREM \rightarrow REM transition (TD z)');

figure; hold on; set(gca,'FontSize',14);
t2 = -7500:3000;
boundedline(t2, mean(rem2nremTD), std(rem2nremTD)./sqrt(size(rem2nremTD,1)), '-k');
plot(t2, mean(rem2nremTD), 'k', 'LineWidth',1);
xticks([-7500 0 3000]); xticklabels({'-5','0','2'});
title('REM \rightarrow NREM transition (TD z)');

%% ---------- Box-plots ----------
[pVel,~] = signrank(nremVel, remVel);
figure; boxplot([nremVel; remVel], [repmat({'NREM'},numel(nremVel),1); repmat({'REM'},numel(remVel),1)], ...
                'OutlierSize',7,'Symbol','k+'); title(sprintf('Velocity (p = %.3f)', pVel));
set(gca,'FontSize',12); set(gcf,'Renderer','painters');

[pTD,~]  = signrank(nremTD, remTD);
figure; boxplot([nremTD; remTD],  [repmat({'NREM'},numel(nremTD),1);  repmat({'REM'},numel(remTD),1)], ...
                'OutlierSize',7,'Symbol','k+'); title(sprintf('TD ratio z (p = %.3f)', pTD));
set(gca,'FontSize',12); set(gcf,'Renderer','painters');
end
