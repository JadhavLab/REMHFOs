
clear all
close all
runscript = 0;
savedata = 1; % save data option - only works if runscript is also on

savedir = '/Volumes/JUSTIN/Inference/';

[y, m, d] = datevec(date);

respRange = 200;

savefig1=0;
val=6;savefile = [savedir 'Allanim_250ctxcoordripplemod_by250mscrit_sleep_VTA_alldata_largewin']; area = 'VTA'; clr = 'b';
% val=7;savefile = [savedir 'Allanim_Coord250ctxripplemod400_by250mscrit_sleep_CA1_alldata_largewin']; area = 'CA1'; clr = 'r';
% val=8;savefile = [savedir 'Allanim_noncoordoord250ctxripplemod400_by250mscrit_sleep_CA1INT_alldata_largewin']; area = 'CA1'; clr = 'r';

% Plot options
plotanimidx =  []; % To pick animals for plotting
plotdays = []; % If you only load data when runscript=0 and savedata=0, then this field will supplant days


%If runscript, run Datafilter and save data
if runscript == 1
    
    %Animal selection
    %-----------------------------------------------------
%     animals = {'ZT2','JS34','JS17','JS21','JS14','JS15','ER1','KL8'};%
    animals = {'TH605'};%

    %-----------------------------------------------------
    
     runepochfilter = 'isequal($type, ''sleep'')';
%     runepochfilter = 'isequal($type, ''run'')';
    
    % Cell filter
    % -----------
    switch val
        
        case 6
            cellfilter = 'strcmp($area, ''VTA'') && ($numspikes > 100)'; % PFC cells with spiking criterion
        case 7
            cellfilter = '(strcmp($area, ''CA1'')|| strcmp($area, ''iCA1'')) && ($numspikes > 100) && ($meanrate < 7) ';
        case 8
            cellfilter = '(strcmp($area, ''CA1'')|| strcmp($area, ''iCA1'')) && ($numspikes > 100) && ($meanrate > 7) ';
    end
    
    % Time filter
    % -----------
    riptetfilter = '(isequal($descrip, ''ctxriptet''))';
    
    % Iterator
    % --------
    iterator = 'singlecellanal';
    
    % Filter creation
    % ----------------
    modf = createfilter('animal',animals,'epochs',runepochfilter, 'cells',...
        cellfilter, 'iterator', iterator);
    disp('Done Filter Creation');
    
    % Set analysis function
    % ----------------------
    
    switch val
        
        % NOTE CHANGE HERE TO MINRIP=1
%         case 6
%             modf = setfilterfunction(modf,'DFAsj_getripalignspiking5',{'spikes', 'rippletime_noncoordSWS', 'tetinfo', 'pos'},'dospeed',0,'lowsp_thrs',4,'minrip',1); % Default stdev is 3
%         case 7
%             modf = setfilterfunction(modf,'DFAsj_getripalignspiking5',{'spikes', 'rippletime_noncoordSWS', 'tetinfo', 'pos'},'dospeed',0,'lowsp_thrs',4,'minrip',1); % Default stdev is 3

        case 6
            modf = setfilterfunction(modf,'DFAsj_getripalignspiking5',{'spikes', 'ctxrippletime_coordSWS', 'tetinfo', 'pos'},'dospeed',0,'lowsp_thrs',4,'minrip',1); % Default stdev is 3
        case 7
            modf = setfilterfunction(modf,'DFAsj_getripalignspiking5',{'spikes', 'ctxrippletime_REM', 'tetinfo', 'pos'},'dospeed',0,'lowsp_thrs',4,'minrip',1); % Default stdev is 3
        case 8
            modf = setfilterfunction(modf,'DFAsj_getripalignspiking5',{'spikes', 'ctxrippletime_REM', 'tetinfo', 'pos'},'dospeed',0,'lowsp_thrs',4,'minrip',1); % Default stdev is 3
    end
    
    % Run analysis
    % ------------
    modf = runfilter(modf);
    disp('Finished running filter script');
    %--------------------- Finished Filter Function Run -------------------
    
    if savedata == 1
        clear figopt1 runscript plotdays plotanimidx savedata
        save(savefile, '-v7.3');
    end
    
else
    load(savefile);
end

if ~exist('savedata')
    return
end


% -------------------------  Filter Format Done -------------------------

% ----------------------------------
% Whether to gather data or to load previously gathered data
% --------------------------------------------------------------------
gatherdata = 1; savegatherdata = 1; consolidate_epochs = 0;

[y, m, d] = datevec(date);

switch val
    
%     case 6
%         
%         gatherdatafile = [savedir 'Allanim_noncoord250ca1ripplemod_by250mscrit_sleep_PFC_alldata_largewin_sepeps_gather_X6'];
%         
%     case 7
%         
%         gatherdatafile = [savedir 'Allanim_noncoord250ca1ripplemod_by250mscrit_sleep_CA1_alldata_largewin_sepeps_gather_X6'];

    case 6

        gatherdatafile = [savedir 'Allanim_250ctxcoordripplemod_by250mscrit_sleep_VTA_alldata_largewin_sepeps_gather_X6'];

    case 7
% 
        gatherdatafile = [savedir 'Allanim_Coord250ctxripplemod400_by250mscrit_sleep_CA1_alldata_largewin_sepeps_gather_X6'];

    case 8
% 
        gatherdatafile = [savedir 'Allanim_noncoord250ca1ripplemod400_by250mscrit_sleep_CA1INT_alldata_largewin_sepeps_gather_X6'];

end

if gatherdata
    
    cnt=0; % Count how many cells will be kept based on nspikes in output: >0
    allanimindex=[]; alldataraster=[]; alldatahist = []; all_Nspk=[];
    
    for an = 1:length(modf)
        for i=1:length(modf(an).output{1})
            % Check for empty output - If Cell defined in rpoch and Nspks in ripple response wndow > 0
            if (modf(an).output{1}(i).Nspikes > 0)
                cnt=cnt+1;
                anim_index{an}(cnt,:) = modf(an).output{1}(i).index;
                % Only indexes
                animindex=[an modf(an).output{1}(i).index]; % Put animal index in front
                allanimindex = [allanimindex; animindex]; % Collect all Anim Day Epoch Tet Cell Index
                % Data
                alldataraster{cnt} = modf(an).output{1}(i).rip_spks_cell; % Only get raster and histogram response
                alldatahist{cnt} = modf(an).output{1}(i).rip_spkshist_cell;
                all_Nspk(cnt) = modf(an).output{1}(i).Nspikes;
                alldataraster_rdm{cnt} = modf(an).output{1}(i).rdm_spks_cell; % Only get raster and histogram response
                alldatahist_rdm{cnt} = modf(an).output{1}(i).rdm_spkshist_cell;
                % trialResps: Summed Nspks/trial in respective window
                alldatatrialResps{cnt} = modf(an).output{1}(i).trialResps;
                alldatatrialResps_bck{cnt} = modf(an).output{1}(i).trialResps_bck;
                alldatatrialResps_rdm{cnt} = modf(an).output{1}(i).trialResps_rdm;
                % Nspikes summed across trials in response and bckgnd window
                all_Nspk_resp(cnt) = sum(modf(an).output{1}(i).trialResps);
                all_Nspk_bck(cnt) = sum(modf(an).output{1}(i).trialResps_bck);
                % Properties
                allcellfr(cnt) = modf(an).output{1}(i).cellfr;
                
                %end
                if cnt==1
                    pret =  modf(an).output{1}(i).pret;
                    postt = modf(an).output{1}(i).postt;
                    binsize = modf(an).output{1}(i).binsize;
                    rwin = modf(an).output{1}(i).rwin;
                    bckwin = modf(an).output{1}(i).bckwin;
                    bins_resp  = modf(an).output{1}(i).bins_resp;
                    bins_bck = modf(an).output{1}(i).bins_bck;
                    timeaxis = modf(an).output{1}(i).timeaxis;
                end
            end
        end
        
    end
    
    % Consolidate single cells across epochs. Multiple methods: see also DFSsj_getcellinfo and DFSsj_xcorrmeasures2
    % ----------------------------------------------------------------------------
    if consolidate_epochs == 1
        allripplemod = struct;
        
        % ---------
        dummyindex=allanimindex;  % all anim-day-epoch-tet-cell indices
        cntcells=0;
        for i=1:length(alldatahist)
            animdaytetcell=allanimindex(i,[1 2 4 5]);
            ind=[];
            while rowfind(animdaytetcell,dummyindex(:,[1 2 4 5]))~=0          % collect all rows (epochs)
                ind = [ind rowfind(animdaytetcell,dummyindex(:,[1 2 4 5]))];        % finds the first matching row
                dummyindex(rowfind(animdaytetcell,dummyindex(:,[1 2 4 5])),:)=[0 0 0 0 0]; % after adding index, remove the corresponding row
                % so you could find the next one
            end
            
            % Gather everything for the current cell across epochs
            currhist=[]; currraster=[]; currNspk=0; currNspk_resp=0; currNspk_bck=0; curr_cellfr=[];
            currhist_rdm=[]; currraster_rdm=[];
            currtrialResps=[]; currtrialResps_rdm=[]; currtrialResps_bck=[];
            for r=ind
                currNspk = currNspk + all_Nspk(r);
                currNspk_resp = currNspk_resp + all_Nspk_resp(r);
                currNspk_bck = currNspk_bck + all_Nspk_bck(r);
                currhist = [currhist; alldatahist{r}];
                currraster = [currraster, alldataraster{r}];
                currhist_rdm = [currhist_rdm; alldatahist_rdm{r}];
                currraster_rdm = [currraster_rdm, alldataraster_rdm{r}];
                currtrialResps = [currtrialResps, alldatatrialResps{r}];
                currtrialResps_rdm = [currtrialResps_rdm, alldatatrialResps_rdm{r}];
                currtrialResps_bck = [currtrialResps_bck, alldatatrialResps_bck{r}];
                curr_cellfr = [curr_cellfr; allcellfr(r)];
            end
            
            % Condition for Nspk. Version 1 had a min of 50 for entire window. Increase it to 100,
            % and can also add a condition for spikes in (resp+bck) window. Need to try a few values
            if (currNspk >= 50)
                
                cntcells = cntcells + 1;
                allripplemod_idx(cntcells,:)=animdaytetcell;
                allripplemod(cntcells).index=animdaytetcell;
                allripplemod(cntcells).hist=currhist*(1000/binsize); % Convert to firing rate in Hz
                allripplemod(cntcells).raster=currraster;
                allripplemod(cntcells).Nspk=currNspk;
                allripplemod(cntcells).hist_rdm=currhist_rdm*(1000/binsize); % Convert to firing rate in Hz
                allripplemod(cntcells).raster_rdm=currraster_rdm;
                % Trial Resps
                allripplemod(cntcells).trialResps = currtrialResps';
                allripplemod(cntcells).trialResps_rdm = currtrialResps_rdm';
                allripplemod(cntcells).trialResps_bck = currtrialResps_bck';
                % Properties
                allripplemod(cntcells).cellfr = nanmean(curr_cellfr);
            end
        end
        
        %Don't consolidate epochs here
        
    elseif consolidate_epochs == 0
        
        allripplemod = struct;
        
        % ---------
        dummyindex=allanimindex;  % all anim-day-epoch-tet-cell indices
        cntcells=0;
        
        
        % Gather everything for the current cell per epoch
        for r=1:length(allanimindex(:,1))
            currhist=[]; currraster=[]; curr_cellfr=[]; currhist_rdm=[]; currraster_rdm=[];
            currtrialResps=[]; currtrialResps_rdm=[]; currtrialResps_bck=[];
            
            animdayepochtetcell = allanimindex(r,:);
            currNspk = all_Nspk(r);
            currNspk_resp = all_Nspk_resp(r);
            currNspk_bck = all_Nspk_bck(r);
            currhist = alldatahist{r};
            currraster = alldataraster{r};
            currhist_rdm = alldatahist_rdm{r};
            currraster_rdm = alldataraster_rdm{r};
            currtrialResps = alldatatrialResps{r};
            currtrialResps_rdm = alldatatrialResps_rdm{r};
            currtrialResps_bck = alldatatrialResps_bck{r};
            curr_cellfr = allcellfr(r);
            
            % Condition for Nspk. Version 1 had a min of 50 for entire window. Increase it to 100,
            % and can also add a condition for spikes in (resp+bck) window. Need to try a few values
            if (currNspk >= 50) && (size(currhist,1) > 1)
                
                cntcells = cntcells + 1;
                
                allripplemod_idx(cntcells,:)=animdayepochtetcell;
                allripplemod(cntcells).index=animdayepochtetcell;
                allripplemod(cntcells).hist=currhist*(1000/binsize); % Convert to firing rate in Hz
                allripplemod(cntcells).raster=currraster;
                allripplemod(cntcells).Nspk=currNspk;
                allripplemod(cntcells).hist_rdm=currhist_rdm*(1000/binsize); % Convert to firing rate in Hz
                allripplemod(cntcells).raster_rdm=currraster_rdm;
                % Trial Resps
                allripplemod(cntcells).trialResps = currtrialResps';
                allripplemod(cntcells).trialResps_rdm = currtrialResps_rdm';
                allripplemod(cntcells).trialResps_bck = currtrialResps_bck';
                % Properties
                allripplemod(cntcells).cellfr = nanmean(curr_cellfr);
            end
        end
    end
   
    % Calculations/ Stats. Stats between response and bck
    % Similar to ...getrip4
    % -----------------------------------------------------------
    for i=1:cntcells
        
        curr_cellfr = allripplemod(i).cellfr;
        currhist = allripplemod(i).hist; %currraster = allripplemod(i).raster;
        currhist_rdm = allripplemod(i).hist_rdm; %currraster_rdm = allripplemod(i).raster_rdm;
        
        % Get the bckresp, rresp and rdmresp again - using firing rate matrices
        rresp = currhist(:,bins_resp);
        bckresp = currhist(:,bins_bck);
        rresp_rdm = currhist_rdm(:,bins_resp);
        
        % Bck
        avgbckresp_trial = mean(bckresp,2); avgrespbck_trial = avgbckresp_trial; % Mean in bck for each ripple
        %trialRespsh_bck = sum(bckresp,2);  % Nspikes/ trial in resp window
        avgbckhist = mean(bckresp); % Avg bck histogram
        mean_bckresp = mean(mean(bckresp)); % Single value
        distr_bckresp = bckresp(:); %All values taken by bins in background
        
        % Response
        avgresp_trial = mean(rresp,2); % Mean for each ripple
        %trialRespsh = sum(rresp,2); % Nspikes/ trial in resp window
        avgresphist = mean(rresp,1); % Avg resp histogram
        mean_rresp = mean(mean(rresp)); % Single value
        
        % RandomResponse
        avgresprdm_trial = mean(rresp_rdm,2); % Mean for each random ripple
        %trialRespsh_rdm = sum(rresp_rdm,2); % Nspikes/ trial in resp window
        avgrdmhist = mean(rresp_rdm,1); % Avg resp histogram
        mean_rdmresp = mean(mean(rresp_rdm)); % Single value
        
        sig_shuf = 0; sig_ttest = 0;
        
        % 0) Simple t-test
        % -----------------
        [sig_ttest, p] = ttest2(avgbckresp_trial, avgresp_trial); % WAS COMMENTED. NEED TO RE_RUN TO GET t-test!! RUN OUTSIDE THIS CODE
        
        
        % % 1) Significance test - USE SHUFFLE BETWEEN MEAN RESP AND MEAN BCK FOR EACH TRIAL
        % % ---------------------------------------------------------------------------
        % Get the actual mean difference between resp and back
        Dm = mean(avgresp_trial) - mean(avgbckresp_trial); % DONT want absolute value. want to shuffle.
        
        % Get significance by comparing Dm to Dshuf. One-tailed test
        % --------------------------------------------------------------
        if Dm>=0
            type = 'exc'; peakresp = max(avgresphist); % Peak in response histogram
        else
            type = 'inh'; peakresp = min(avgresphist); % Trough in response histogram
        end
        
        % Get %tage change over baseline and peak/trough: Save with sign- +ve or -ve modln. Not abs
        % ----------------------------------------------
        modln_raw = Dm;
        modln = 100*(Dm)/mean(avgbckresp_trial); % Mean %tage change above/below baseline
        peakchange = peakresp - mean(avgbckhist); % mean(avgbckhist) is same as mean(avgbckresp_trial)
        modln_peak = 100*(peakchange)/mean(avgbckresp_trial); % Peak/trough %tage change above/below baseline
        modln_div = 100*peakresp/mean(avgbckresp_trial);
        
        % Save
        % -----
        allripplemod(i).Dm = Dm;
        allripplemod(i).sig_ttest = sig_ttest; % Sig or not
        allripplemod(i).modln_peak = modln_peak; % %tage peak change over baseline
        allripplemod(i).modln_div = modln_div;
        allripplemod(i).modln = modln; % %tage mean change over baseline
        allripplemod(i).modln_raw = modln_raw; % Raw change in firing rate. The statistic
        allripplemod(i).type = type; % exc or inh
        allripplemod(i).anim = allripplemod(i).index(1); allanim(i) = allripplemod(i).index(1);
        allripplemod(i).days = allripplemod(i).index(2); alldays(i) = allripplemod(i).index(2);
        
        % Mean resp from histogram in respective window
        allripplemod(i).avghisttrialResps = avgresp_trial;
        %  allripplemod(i).avghisttrialResps_rdm = avgresprdm_trial;
        allripplemod(i).avghisttrialResps_bck = avgrespbck_trial;
        
        % Properties
        allripplemod(i).cellfr = curr_cellfr;
        allripplemod(i).Nrip = length(avgresp_trial);
        
    end
    
    b=gaussian(20,61);
    
    if respRange == 200
    swronset=1050;
        varRange=[1:200]+swronset;
    elseif respRange == 400
        swronset=1050;
        varRange=[-200:200]+swronset;
    else
        keyboard
    end
    
    for ii=1:cntcells
        ii
        curRast=allripplemod(ii).raster;
        numRips=length(curRast);
        % creating a raster matrix of 0/1 with 1ms resolution
        % it goes from -1000 to +1000ms relative to ripples
        
        curRastMat=zeros(numRips,2200);
        
        for i=1:numRips,curRastMat(i,round(curRast{i}+1101))=1;end
        
        respVarShufs=[];
        numRuns=1000;
        allShufMeans=[];
        % shuffling each trial in the raster separately in time,
        % cyclically. Doing that numRuns times (currently 1000).
        for runs=1:numRuns
            curRastMatShuf=zeros(numRips,2200);
            
            for i=1:numRips
                shiftVal=round(rand(1)*2200);
                curRastMatShuf(i,1+mod(round(curRast{i}+1101+shiftVal),2200))=1;
            end
            
            % shuffled "psth"
            meanRespShuf=filtfilt(b,1,mean(curRastMatShuf(:,50:end-50),1));
            
            allShufMeans=[allShufMeans; meanRespShuf];
            
            % for each shuffled psth, calculate the variance (mean squared
            % distance from mean)
            respVarShufs=[respVarShufs var(meanRespShuf(varRange))];
            
        end
        
        meanResp=filtfilt(b,1,mean(curRastMat(:,50:2150)));
        meanRespRange=meanResp(varRange);
        meanRespShuf=(mean(allShufMeans));
        meanRespShufRange= meanRespShuf(varRange);
       
        %new measure: instead of mean squared distance from its own mean,
        %mean squared distance from the mean shuffled psth's
        respVar2=mean((meanRespRange-meanRespShufRange).^2);
        % this is used to get a positive deviation from mean
        varPSTH=(meanResp-meanRespShuf).^2;
        rasterShufP2=1-sum(respVar2>respVarShufs)/numRuns;
        
        allripplemod(ii).rasterShufP2=rasterShufP2; %significant modulation?
        
        varRespAmp2=respVar2/mean(respVarShufs);
        
        allripplemod(ii).varRespAmp2=varRespAmp2;
        allripplemod(ii).varPSTH=varPSTH;
    end
    
    
    
    % Save
    % -----
    if savegatherdata == 1
        save(gatherdatafile, '-v7.3');
        return;
    end
    
else % gatherdata=0
    
    load(gatherdatafile);
    
end % end gather data





