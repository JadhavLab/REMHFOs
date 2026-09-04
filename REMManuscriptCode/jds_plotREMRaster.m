% Plot example raster - Plots CA1-PFC cells, CA1-PFC theta amplitude and
% PFC ripples
% -------------------------------------------------------------------------

clear; close all;
prefix='KL8';
day=1;
epochs=15;
win = 30;

% Get Data Location
% -------------------------------
switch prefix
    case 'KL8'
        directoryname = '/Volumes/JUSTIN/SingleDay/KL8_direct/';
        dire = '/Volumes/JUSTIN/SingleDay/KL8_direct/';
        animdirect = directoryname;
        riptetlist = 15;  % No Need if no ripples
        ctxriptetlist = 30;
        maineegtet = 54;  % CA1 tet % No need if no EEG
        peegtet = 30; % PFCtet % No need if no EEG
end

currdir = pwd;
if (directoryname(end) == '/')
    directoryname = directoryname(1:end-1);
end
if (dire(end) == '/')
    dire = dire(1:end-1);
end

if (day < 10)
    daystring = ['0',num2str(day)];
else
    daystring = num2str(day);
end

saveg = 0;

for e = 1:length(epochs)
    epoch = epochs(e);
    % Get data files
    spikefile = sprintf('%s/%sspikes%02d.mat', directoryname, prefix, day);
    load(spikefile);
    tetinfofile = sprintf('%s/%stetinfo.mat', directoryname, prefix);
    load(tetinfofile);
    cellinfofile = sprintf('%s/%scellinfo.mat', directoryname, prefix);
    load(cellinfofile);
  
    pripfile = sprintf('%s/%sctxrippletime_REM%02d.mat', directoryname, prefix, day);
    load(pripfile);

    thetaamp = [];
    ripamp = [];
    for c = 1:length(ctxriptetlist)
        tet = ctxriptetlist(c);
        thetafile = sprintf('%s/EEG/%sthetagnd%02d-%02d-%02d.mat', directoryname, prefix, day,epoch,tet);
        load(thetafile);
        thetaamptmp = thetagnd{day}{epoch}{tet}.data(:,1);
        thetaamp = [thetaamp; thetaamptmp'];

        teegTheta = geteegtimes(thetagnd{day}{epoch}{tet});

        ripfile = sprintf('%s/EEG/%sripple%02d-%02d-%02d.mat', directoryname, prefix, day,epoch,tet);
        load(ripfile);
        ripamptmp = ripple{day}{epoch}{tet}.data(:,1);
        ripamp = [ripamp; ripamptmp'];

        teegRip = geteegtimes(ripple{day}{epoch}{tet});
    end
    thetaamp = mean(thetaamp,1);
    ripamp = mean(ripamp,1);

    thetaampHp = [];
    ripampHp = [];
    for c = 1:length(riptetlist)
        tet = riptetlist(c);
        thetafile = sprintf('%s/EEG/%sthetagnd%02d-%02d-%02d.mat', directoryname, prefix, day,epoch,tet);
        load(thetafile);
        thetaamptmp = thetagnd{day}{epoch}{tet}.data(:,1);
        thetaampHp = [thetaampHp; thetaamptmp'];

        teegTheta = geteegtimes(thetagnd{day}{epoch}{tet});

        ripfile = sprintf('%s/EEG/%sripple%02d-%02d-%02d.mat', directoryname, prefix, day,epoch,tet);
        load(ripfile);
        ripamptmp = ripple{day}{epoch}{tet}.data(:,1);
        ripampHp = [ripampHp; ripamptmp'];

        teegRip = geteegtimes(ripple{day}{epoch}{tet});
    end
    thetaampHp = mean(thetaampHp,1);
    ripampHp = mean(ripampHp,1);

    rem = load(sprintf('%s/%srem%02d.mat', directoryname, prefix, day));
    rem = rem.rem;
    remtimes = rem{day}{epoch}.starttime;

    % Get cells
    % ---------
    % CA1 cells (black)
    filterString = 'strcmp($tag2, ''CA1Pyr'') && ($numspikes > 100)';

    cellindices = evaluatefilter(cellinfo{day}{epoch}, filterString);
    cellsi = [repmat([day epoch], size(cellindices,1),1 ), cellindices]; % day-epoch-tet-cell for CA1 cells
    usecellsi = 1:size(cellsi,1);

    % PFC cells
    filterString = 'strcmp($area, ''PFC'') && ($numspikes > 100)';
    pcellindices = evaluatefilter(cellinfo{day}{epoch}, filterString);
    cellsp = [repmat([day epoch], size(pcellindices,1),1 ), pcellindices]; % day-epoch-tet-cell for PFC cells
    usecellsp = 1:size(cellsp,1);
 
    for i=1:size(cellsi,1)
        eval(['spiketimei{',num2str(i),'}= spikes{cellsi(',num2str(i),',1)}{cellsi(',num2str(i),',2)}'...
            '{cellsi(',num2str(i),',3)}{cellsi(',num2str(i),',4)}.data(:,1);']);
        eval(['spikeposi{',num2str(i),'}= spikes{cellsi(',num2str(i),',1)}{cellsi(',num2str(i),',2)}'...
            '{cellsi(',num2str(i),',3)}{cellsi(',num2str(i),',4)}.data(:,2:3);']);
        eval(['spikeposidxi{',num2str(i),'}= spikes{cellsi(',num2str(i),',1)}{cellsi(',num2str(i),',2)}'...
            '{cellsi(',num2str(i),',3)}{cellsi(',num2str(i),',4)}.data(:,7);']);
    end

    for i=1:size(cellsp,1)
        eval(['spiketimep{',num2str(i),'}= spikes{cellsp(',num2str(i),',1)}{cellsp(',num2str(i),',2)}'...
            '{cellsp(',num2str(i),',3)}{cellsp(',num2str(i),',4)}.data(:,1);']);
        eval(['spikeposp{',num2str(i),'}= spikes{cellsp(',num2str(i),',1)}{cellsp(',num2str(i),',2)}'...
            '{cellsp(',num2str(i),',3)}{cellsp(',num2str(i),',4)}.data(:,2:3);']);
        eval(['spikeposidxp{',num2str(i),'}= spikes{cellsp(',num2str(i),',1)}{cellsp(',num2str(i),',2)}'...
            '{cellsp(',num2str(i),',3)}{cellsp(',num2str(i),',4)}.data(:,7);']);
    end

    p_riptimes = [ctxripple{day}{epoch}.starttime ctxripple{day}{epoch}.endtime]; %Noncoordinated ctx ripples
    p_rip_starttime = p_riptimes(:,1);
    p_rip_endtime = p_riptimes(:,2);

    % ------------------------------
    % Figure Parametersand Font Sizes
    % ------------------------------
    forppr = 0;

    set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);

    if forppr==1
        set(0,'defaultaxesfontsize',16);
        tfont = 18; % title font
        xfont = 16;
        yfont = 16;
    else
        set(0,'defaultaxesfontsize',24);
        tfont = 28;
        xfont = 20;
        yfont = 20;
    end
    clr = {'b',[0.8500 0.3250 0.0980],'g','y',[0.4940 0.1840 0.5560],'r','b','g','y','b',[0.8500 0.3250 0.0980],'g','y',[0.4940 0.1840 0.5560],'r','b','g','y','b',[0.8500 0.3250 0.0980],'g','y',[0.4940 0.1840 0.5560],'r','b','g','y'};

    clr1='k';
    clr2='r';

    figdir = '/Volumes/JUSTIN/RastersREM/';

%         winst = 16608.5; %EXAMPLE PLOT
%         winend = 16618.5; % secs
    %%
    % epochend = 6840;

    for rast = 1:length(remtimes(:,1))
%         winst = remtimes(rast,1); %- 2; %start 2 seconds before start of REM
%         winend = remtimes(rast,1) + 10; % full 10 seconds
%         winst = 20080; 
%         winend = 20085; 
        winst = 18925.9; 
        winend = 18935.9; 
        epochend = rem{day}{epoch}.endtime(end);
        epochend = winend;

        ii = 1;

        rastnum
        rastnum = rastnum + 1;

        figure(1); xlim([0 win]); hold on;
        redimscreen;

        taxis = winst:winend; 
        taxis = taxis - winst;

        winst_ms = winst*1000;
        winend_ms = winend*1000;

        eind1Rip = lookup(winst, teegRip);
        eind2Rip = lookup(winend, teegRip);
        taxisEEGRip = teegRip(eind1Rip:eind2Rip);
        taxisEEGRip = taxisEEGRip - winst;

        eind1Theta = lookup(winst, teegTheta);
        eind2Theta = lookup(winend, teegTheta);
        taxisEEGTheta = teegTheta(eind1Theta:eind2Theta);
        taxisEEGTheta = taxisEEGTheta - winst;

        baseline = 0;

        % First PFC Spikes on bottom

        cnt = 0;
        activepfc = 0;
        [B, I] = sort(cellfun(@length,spiketimep),'descend');
        for c=usecellsp
            cc = I(c);
            eval(['currspkt = spiketimep{',num2str(cc),'};']);
            currspkt = currspkt;
            currspkt = currspkt(find(currspkt>=winst & currspkt<=winend ));

            if ~isempty(currspkt)
                currspkt = currspkt - winst;
            end
            figure(1); xlim([0 win]); hold on;
            if ~isempty(currspkt)
                activepfc = activepfc+1;
                cnt=cnt+1;
                if size(currspkt,2)~=1, currspkt=currspkt'; end
                plotraster(currspkt,(baseline+2*(cnt-1))*ones(size(currspkt)),1.8,[],'Color','r','LineWidth',1);
            end
        end

        baseline = baseline + (activepfc)*2;
        baseline = baseline+1;

        % Now, CA1 spikes
        % ---------------
        cnt = 0;
        activeca1cnt = 0;
        [B, I] = sort(cellfun(@length,spiketimei),'descend');
        for c=usecellsi
            cc = I(c);
            eval(['currspkt = spiketimei{',num2str(cc),'};']);
            currspkt = currspkt;
            currspkt = currspkt(find(currspkt>=winst & currspkt<=winend ));

            % If spikes, subtract from subtract from start time and bin
            if ~isempty(currspkt)
                currspkt = currspkt - winst;
            end

            figure(1); hold on;
            if ~isempty(currspkt)
                activeca1cnt = activeca1cnt+1;
                cnt=cnt+1;
                if size(currspkt,2)~=1, currspkt=currspkt'; end
                plotraster(currspkt,(baseline+2*(cnt-1))*ones(size(currspkt)),1.8,[],'Color','k','LineWidth',1);
            end
        end
        cnt = 0;
        baseline = baseline + (activepfc)*2;
        baseline = baseline+4;

        cripsinwin = p_rip_starttime(find(p_rip_starttime>=winst & p_rip_starttime<=winend ));
        cripsinwinend = p_rip_endtime(find(p_rip_endtime>=winst & p_rip_endtime<=winend ));

        if ~isempty(cripsinwin)
            cripsinwin = cripsinwin - winst;
            cripsinwinend = cripsinwinend - winst;
        end

        plotraster(cripsinwin,(baseline+2*(cnt-1))*ones(size(cripsinwin)),1.8,[],'Color','m','LineWidth',1);
        plotraster(cripsinwinend,(baseline+2*(cnt-1))*ones(size(cripsinwinend)),1.8,[],'Color','k','LineWidth',1);
        
        baseline = baseline+4;

        curreeg = double(ripamp(eind1Rip:eind2Rip));

        %% Plot PFC AMP
        eegscale = max(curreeg)-min(curreeg);
        downeeg = baseline; upeeg = downeeg+8;
        plotscale = 8;
        curreeg = downeeg + (plotscale/2) + curreeg.*(plotscale/eegscale);
        plot(taxisEEGRip,curreeg,'k-','LineWidth',1);
        
        baseline = baseline+4;
        curreeg = double(thetaamp(eind1Theta:eind2Theta));

        eegscale = max(curreeg)-min(curreeg);
        downeeg = baseline; upeeg = downeeg+8;
        plotscale = 8;
        curreeg = downeeg + (plotscale/2) + curreeg.*(plotscale/eegscale);
        plot(taxisEEGTheta,curreeg,'b-','LineWidth',1);

        %% Plot CA1 theta AMP
        curreeg = double(ripampHp(eind1Rip:eind2Rip));

        baseline = baseline+4;
        curreeg = double(thetaampHp(eind1Theta:eind2Theta));

        eegscale = max(curreeg)-min(curreeg);
        downeeg = baseline; upeeg = downeeg+8;
        plotscale = 8;
        curreeg = downeeg + (plotscale/2) + curreeg.*(plotscale/eegscale);
        plot(taxisEEGTheta,curreeg,'r-','LineWidth',1);
        %%
        winsecs = [0:10:winend-winst]; 
        secs = [winst:10:winend];
        secs = round(secs); 
        msecs = [winst_ms:10000:winend_ms];
        xlabel('Time (secs)','FontSize',18,'Fontweight','normal');
        title([prefix ' - Day',num2str(day) 'Ep' num2str(epoch) ' Window Time: ' num2str(roundn(winst,-1))...
            ' - ' num2str(roundn(winend,-1)) '  hpctet - ' num2str(maineegtet)...
            ' pfctet - ' num2str(peegtet)],'FontSize',18,'Fontweight','normal');

        baseline = baseline+1;
        set(gca,'XLim',[0 winend-winst]);
        set(gca,'YLim',[0 baseline+10]);

        keyboard; % Pause after each plot

        if saveg==1
            figfile = [figdir,prefix,'Day',num2str(day),'Ep',num2str(epoch),'RasteregNo',num2str(ii),'Window',num2str(win),'Overlap',num2str(overlap)];
            print('-djpeg', figfile);
        end

        ii = ii+1;
        close all

    end
end

keyboard;

