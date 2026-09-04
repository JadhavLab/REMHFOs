% based on ripTrigFiringOnsetAnalysis3
clear all
close all

scrsz = get(0,'ScreenSize');
savedirX = '/Volumes/JUSTIN/SingleDay/ProcessedDataREM/';
% savedirX = '/Volumes/JUSTIN/Inference/';

% savedirX = '/Volumes/OneTouch/Taste/EJ7/20220519_LinearTrackRec_direct/';
plotIndividualCells=0;
areas1={'CA1','PFC'};

usedata = 1; % 1=alldata, 2 = only wtrack, 3 = only y-track, 4 = skip iCA1)

switch usedata
    case 1
        filestr = '';
    case 2
        filestr = '_Wtr_noiCA1';
    case 3
        filestr = '_Ytr';
    case 4
        filestr = '_NOiCA1';
end



figdir = '/data25/sjadhav/HPExpt/Figures/RevFig1/';

savefig1=0;
set(0,'defaultaxesfontsize',16);
% twin = 801:1300;
twin = 501:1500;

histOns=[];
histOnsW=[];
allxcorrLags={};
allallripmodhists={};
allallripmodhistsW={};
allallripmodhistssig={};
allallripmodhistsWsig={};
allanimdayvecRip={};
allanimdayvecRipW={};
allallripmodhistssigInd={};
allallripmodhiststype={};

allallsiginds={};

% how to smooth psths
b=gaussian(20,61);
animIdx = {};

for areaInd=1:2
    area1=areas1{areaInd};
    switch area1
        
        case 'CA1'
            
            load([savedirX 'Allanim_250SWRremripplemodwin_0mscrit_sleep_CA1_alldata_largewin_sepeps_gather_X6'])
%             load([savedirX 'Allanim_noncoord250ctxripplemod_by250mscrit_sleep_PFC_alldata_largewin_sepeps_gather_X6.mat'])
            allripplemod_idx=[];
            for w=1:length(allripplemod),allripplemod_idx=[allripplemod_idx;allripplemod(w).index];end
            
        case 'PFC'
%             load([savedirX 'Allanim_250ctxswsripplemod_0mscrit_sleep_PFC_alldata_largewin_sepeps_gather_X6.mat'])
            load([savedirX 'Allanim_250swrnremcoordripplemod400_0mscrit_sleep_CA1_alldata_largewin_sepeps_gather_X6']);
%             load([savedirX 'Allanim_250ctxremripplemod_by250mscrit_sleep_PFC_alldata_largewin_sepeps_gather_X6.mat']);
            allripplemod_idx=[];
            for w=1:length(allripplemod),allripplemod_idx=[allripplemod_idx;allripplemod(w).index];end

        case 'VTA'
            if areaInd == 1
                load([savedirX 'Allanim_250ctxcoordripplemod_by250mscrit_sleep_VTA_alldata_largewin_sepeps_gather_X6.mat']);
                allripplemod_idx=[];
                for w=1:length(allripplemod),allripplemod_idx=[allripplemod_idx;allripplemod(w).index];end
            else
                load([savedirX 'Allanim_250ctxnoncoordripplemod_by250mscrit_sleep_VTA_alldata_largewin_sepeps_gather_X6.mat']);
                allripplemod_idx=[];
                for w=1:length(allripplemod),allripplemod_idx=[allripplemod_idx;allripplemod(w).index];end
            end
    end
    allinds=unique([allripplemod_idx],'rows');
    allripmodhists=[];
    allripmodhistssig=[];
    allripmodhistssigInd=[];
    allripmodonset3=[];
    allripmodhiststype=[];
    allanimidx = [];
    
    allsiginds=[];
    
    for i=1:length(allripplemod)
        allripmodhists=[allripmodhists; zscore(filtfilt(b,1,mean(rast2mat_lrg(allripplemod(i).raster))))];
        
        
        %For CA2 tag file
        doflag=0;
        switch area1
            case 'CA1'                
                if allripplemod(i).rasterShufP2<0.05
                   doflag=1;
                end
                
            case 'PFC'
                if allripplemod(i).rasterShufP2<0.05
                    doflag=1;
                end
            case 'VTA'
                if allripplemod(i).rasterShufP2<0.05
                    doflag=1;
                end
        end
        if doflag==1
            curranim = allripplemod(i).index(1);
            allanimidx = [allanimidx; curranim];
                switch usedata
                    case 1 % all data
                        allripmodhistssig=[allripmodhistssig; zscore(filtfilt(b,1,mean(rast2mat_lrg(allripplemod(i).raster))))];
                        allripmodhistssigInd=[allripmodhistssigInd; allripplemod(i).index];
                        % 1 for swr-exc, 0 for swr-inh
                        allripmodhiststype=[allripmodhiststype strcmp(allripplemod(i).type, 'exc')];
                        allsiginds = [allsiginds; allripplemod(i).index];
                        
                    case 2 % only W track
                        %                         if curranim<=3
                        %                             allripmodhistssig=[allripmodhistssig; zscore(filtfilt(b,1,mean(rast2mat(allripplemod(i).raster))))];
                        %                             allripmodhistssigInd=[allripmodhistssigInd; allripplemod(i).index];
                        %                             % 1 for swr-exc, 0 for swr-inh
                        %                             allripmodhiststype=[allripmodhiststype strcmp(allripplemod(i).type, 'exc')];
                        %                             allsiginds = [allsiginds; allripplemod(i).index];
                        %                         end
                        
                        % For no iCA1
                        if curranim<=3
                            if areaInd==1 % if CA1 file
                                if strcmp(allripplemod(i).iCA1tag,'n')
                                    allripmodhistssig=[allripmodhistssig; zscore(filtfilt(b,1,mean(rast2mat_lrg(allripplemod(i).raster))))];
                                    allripmodhistssigInd=[allripmodhistssigInd; allripplemod(i).index];
                                    % 1 for swr-exc, 0 for swr-inh
                                    allripmodhiststype=[allripmodhiststype strcmp(allripplemod(i).type, 'exc')];
                                    allsiginds = [allsiginds; allripplemod(i).index];
                                end
                            else % Do usual if PFC file
                                allripmodhistssig=[allripmodhistssig; zscore(filtfilt(b,1,mean(rast2mat_lrg(allripplemod(i).raster))))];
                                allripmodhistssigInd=[allripmodhistssigInd; allripplemod(i).index];
                                % 1 for swr-exc, 0 for swr-inh
                                allripmodhiststype=[allripmodhiststype strcmp(allripplemod(i).type, 'exc')];
                                allsiginds = [allsiginds; allripplemod(i).index];
                            end
                        end
                        
                    case 3 % only Y track
                        if curranim>=4
                            allripmodhistssig=[allripmodhistssig; zscore(filtfilt(b,1,mean(rast2mat_lrg(allripplemod(i).raster))))];
                            allripmodhistssigInd=[allripmodhistssigInd; allripplemod(i).index];
                            % 1 for swr-exc, 0 for swr-inh
                            allripmodhiststype=[allripmodhiststype strcmp(allripplemod(i).type, 'exc')];
                            allsiginds = [allsiginds; allripplemod(i).index];
                        end
                        
                    case 4 % skip iCA1 if its the CA1 file,
                        if areaInd==1 % if CA1 file
                            if strcmp(allripplemod(i).iCA1tag,'n')
                                allripmodhistssig=[allripmodhistssig; zscore(filtfilt(b,1,mean(rast2mat_lrg(allripplemod(i).raster))))];
                                allripmodhistssigInd=[allripmodhistssigInd; allripplemod(i).index];
                                % 1 for swr-exc, 0 for swr-inh
                                allripmodhiststype=[allripmodhiststype strcmp(allripplemod(i).type, 'exc')];
                                allsiginds = [allsiginds; allripplemod(i).index];
                            end
                        else % Do usual if PFC file
                            allripmodhistssig=[allripmodhistssig; zscore(filtfilt(b,1,mean(rast2mat_lrg(allripplemod(i).raster))))];
                            allripmodhistssigInd=[allripmodhistssigInd; allripplemod(i).index];
                            % 1 for swr-exc, 0 for swr-inh
                            allripmodhiststype=[allripmodhiststype strcmp(allripplemod(i).type, 'exc')];
                            allsiginds = [allsiginds; allripplemod(i).index];
                        end
                end
          
            %else
            %   allripmodonset3=[allripmodonset3 nan];
            %allripmodhists=[allripmodhists; nan(1,size(allripmodhists,2))];
            
            %    end;
        end
    end
    
    allallripmodhists{end+1}=allripmodhists;
    allallripmodhistssig{end+1}=allripmodhistssig;
    allallripmodhistssigInd{end+1}=allripmodhistssigInd;
    allallripmodhiststype{end+1}=allripmodhiststype;
    animIdx{end+1} = allanimidx;
    %allallsiginds{end+1} = allsiginds;
    
    cntcellsRip=length(allripplemod);
    
    %     if areaInd==2
    %     if plotIndividualCells
    %         for kk=1:length(allripplemod)
    %         figure('Position',[300 200 scrsz(3)/3 600])
    %
    %
    %             imagesc(rast2mat(allripplemod(kk).raster));colormap(gray)
    %
    %
    %         keyboard
    %         end
    %     end
    %  end
end
clear allripplemod
CA1psths=allallripmodhistssig{1};
CA1inds=allallripmodhistssigInd{1};
PFCpsths=allallripmodhistssig{2};
PFCinds=allallripmodhistssigInd{2};

%% RATE OF SWR-MODULATED CELLS

CA1sleepSWRmodrate=size(CA1psths,1)./size(allallripmodhists{1},1)
PFCsleepSWRmodrate=size(PFCpsths,1)./size(allallripmodhists{2},1)


%% Try plotting non ripmod ca1 cells: have to set condition for rastershuf to >= 0.05 above
% xaxis=-499:500;
% CA1 non ripmod cells
% [tmp timingindCA1a]=max(CA1psths(:,251:750)');
% [A2 B2]=sort(timingindCA1a);
% figure;
% subplot(3,1,1)
% imagesc(xaxis,1:length(B2), CA1psths(B2,:));caxis([-3 3])
% xlabel('Time (ms)')
% ylabel('#Cell')
% title('CA1 SWR-nonmod cells')
% xlim([-400 400])
%
% PFC non ripmod cells
% [tmp timingindPFCa]=max(PFCpsths(:,251:750)');
% [A2 B2]=sort(timingindPFCa);
% figure;
% subplot(3,1,1)
% imagesc(xaxis,1:length(B2), PFCpsths(B2,:));caxis([-3 3])
% xlabel('Time (ms)')
% ylabel('#Cell')
% title('PFC SWR-nonmod cells')
% xlim([-400 400])

%% ripple-triggered spiking in all regions. Currently sorted by response amplitude, can
% also sort by days using:
% [sortedDays2 sortedDaysS2]=sort(allanimdayvecRip{2}(:,2),'descend');

% SJ - Comment Fig
% % sleep
% figure('Position',[900 100 scrsz(3)/5 scrsz(4)])
% xaxis=-500:10:500;
% subplot(2,1,1);
% [sorted1 sortedS1]=sort(mean(CA1psths(:,50:70),2),'descend');
% imagesc(xaxis,1:size(CA1psths,1),CA1psths(sortedS1,:));
% xlim(zoomWindow)
% title('Ripple-triggered PSTHs, all rip-mod cells, CA1');xlabel('Time (ms)');ylabel('Cells')
% subplot(2,1,2);
% [sorted2 sortedS2]=sort(mean(PFCpsths(:,50:70),2),'descend');
% imagesc(xaxis,1:size(PFCpsths,1),PFCpsths(sortedS2,:));
% xlim(zoomWindow)
% title('Ripple-triggered PSTHs, all rip-mod cells, PFC');xlabel('Time (ms)');ylabel('Cells')
% annotation('textbox', [0.5, 1, 0, 0], 'string', 'Sleep','fontsize',20)
%% unsorted
zoomWindow=[-1000 1000];
figure('Position',[900 100 scrsz(3)/5 scrsz(4)])
xaxis=-1049:1050;
subplot(2,1,1);
imagesc(xaxis,1:size(CA1psths,1),CA1psths);
xlim(zoomWindow)
title('Ripple-triggered PSTHs, all rip-mod cells, PFCexcREM');xlabel('Time (ms)');ylabel('Cells')
caxis([-4 4])
subplot(2,1,2);
imagesc(xaxis,1:size(PFCpsths,1),PFCpsths);
xlim(zoomWindow)
caxis([-4 4])
title('Ripple-triggered PSTHs, all rip-mod cells, PFCexcNREM');xlabel('Time (ms)');ylabel('Cells')

%annotation('textbox', [0.5, 1, 0, 0], 'string', 'Sleep','fontsize',20)
annotation('textbox', [0.5, 1, 0, 0], 'string', 'Run','fontsize',20)

% figfile = [figdir,'RipTrig_PoplnZscore_unsorted_CA1andPFC',filestr]
if savefig1==1
    print('-dpdf', figfile); print('-dpng', figfile, '-r300'); saveas(gcf,figfile,'fig'); print('-depsc2', figfile); print('-djpeg', figfile);
end


%%  exc/inh

%excinh=sign(mean(PFCpsths(:,500:650),2)-mean(PFCpsths(:,250:500),2));
excinhPFC=allallripmodhiststype{2};
excinhCA1=allallripmodhiststype{1};
%excinhPFCinds=allallsiginds{2};

excPFCpsths=PFCpsths(excinhPFC>0,:);
inhPFCpsths=PFCpsths(excinhPFC==0,:);

excCA1psths=CA1psths(excinhCA1>0,:);
inhCA1psths=CA1psths(excinhCA1==0,:);

%excinh1=sign(mean(CA1psths(:,50:65),2)-mean(CA1psths(:,25:50),2));
% FOR DEMETRIS: INDICES OF SIG SWR-EXC/INH
PFCindsExc=PFCinds(excinhPFC>0,:);
PFCindsInh=PFCinds(excinhPFC==0,:);

% SJ - Comment Fig
% figure('Position',[900 100 scrsz(3)/5 scrsz(4)]);
%
% subplot(3,1,1)
% imagesc(xaxis,1:size(CA1psths,1),CA1psths)
% caxis([-4 4]);xlim([-400 400])
% title('CA1');xlabel('Time (ms)')
%
%
%
% subplot(3,1,2)
% imagesc(xaxis,1:size(excPFCpsths,1),excPFCpsths)
% title('sig rip-exc');xlabel('Time (ms)')
%
% caxis([-4 4]);xlim([-400 400])
% subplot(3,1,3)
% imagesc(xaxis,1:size(inhPFCpsths,1),inhPFCpsths)
% title('sig rip-inh');xlabel('Time (ms)')
%
% caxis([-4 4]);xlim([-400 400])
%%
figure
hold on
% shadedErrorBar(xaxis,mean(inhCA1psths,1),std(inhCA1psths)./sqrt(size(inhCA1psths,1)),'-k',1);
% shadedErrorBar(xaxis,mean(inhPFCpsths,1),std(inhPFCpsths)./sqrt(size(inhPFCpsths,1)),'-r',1);
shadedErrorBar(xaxis,mean(CA1psths,1),std(CA1psths)./sqrt(size(CA1psths,1)),'-b',1);
for a = 1:8
    plot(xaxis, mean(CA1psths(find(animIdx{1}==a),:),1),'-k')
end
xlim([-1000 1000])
ylabel('Mean z-scored psth')
xlabel('Time (ms)')
set(gcf, 'renderer', 'painters')
figure
hold on
shadedErrorBar(xaxis,mean(PFCpsths,1),std(PFCpsths)./sqrt(size(PFCpsths,1)),'-r',1);
for a = 1:8
    plot(xaxis, mean(PFCpsths(find(animIdx{2}==a),:),1),'-k')
end
xlim([-1000 1000])
ylabel('Mean z-scored psth')
xlabel('Time (ms)')
set(gcf, 'renderer', 'painters')

figure
hold on
% shadedErrorBar(xaxis,mean(CA1psths,1),std(CA1psths)./sqrt(size(CA1psths,1)),'-k',1);
shadedErrorBar(xaxis,mean(excCA1psths,1),std(excCA1psths)./sqrt(size(excCA1psths,1)),'-r');
% shadedErrorBar(xaxis,mean(inhPFCpsths,1),std(inhPFCpsths)./sqrt(size(inhPFCpsths,1)),'-r',0);
shadedErrorBar(xaxis,mean(excPFCpsths,1),std(excPFCpsths)./sqrt(size(excPFCpsths,1)),'-b');
xlim([-1000 1000])
ylabel('Mean z-scored psth')
xlabel('Time (ms)')
legend('PFCexc REM','PFCexc NREM');
% legend('CA1','PFCexc','PFCinh');

figure
hold on
shadedErrorBar(xaxis,mean(CA1psths,1),std(CA1psths)./sqrt(size(CA1psths,1)),'-k',1);
% shadedErrorBar(xaxis,mean(inhCA1psths,1),std(inhCA1psths)./sqrt(size(inhCA1psths,1)),'-r');
shadedErrorBar(xaxis,mean(PFCpsths,1),std(PFCpsths)./sqrt(size(PFCpsths,1)),'-r',0);
% shadedErrorBar(xaxis,mean(inhPFCpsths,1),std(inhPFCpsths)./sqrt(size(inhPFCpsths,1)),'-b');
xlim([-1000 1000])
ylabel('Mean z-scored psth')
xlabel('Time (ms)')
legend('PFCexc REM','PFCexc NREM');
% legend('CA1','PFCexc','PFCinh');


% Use jbfill to avoid shaded area errors:
figure; hold on;
plot(xaxis,mean(CA1psths,1),'k-');
jbfill(xaxis,mean(CA1psths,1)+(std(CA1psths)./sqrt(size(CA1psths,1))),...
    mean(CA1psths,1)-(std(CA1psths)./sqrt(size(CA1psths,1))),'k','k',1,1);
plot(xaxis,mean(excPFCpsths,1),'r-');
jbfill(xaxis,mean(excPFCpsths,1)+(std(excPFCpsths)./sqrt(size(excPFCpsths,1))),...
    mean(excPFCpsths,1)-(std(excPFCpsths)./sqrt(size(excPFCpsths,1))),'r','r',1,1);
plot(xaxis,mean(inhPFCpsths,1),'b-');
jbfill(xaxis,mean(inhPFCpsths,1)+(std(inhPFCpsths)./sqrt(size(inhPFCpsths,1))),...
    mean(inhPFCpsths,1)-(std(inhPFCpsths)./sqrt(size(inhPFCpsths,1))),'b','b',1,1);
xlim([-1000 1000])
ylabel('Mean z-scored psth')
xlabel('Time (ms)')


% figfile = [figdir,'RipTrig_MeanZscorePSTHN2',filestr]
if savefig1==1,
    print('-dpdf', figfile); print('-dpng', figfile, '-r300'); saveas(gcf,figfile,'fig'); print('-depsc2', figfile); print('-djpeg', figfile);
end


% forFigIndsexc=[1 2 17 2; 2 2 12 1];
% forFigIndsinh=[2 1 9 1; 4 9 21 2; ];
% forFigIndsneu=[2 3 10 1];
% 
% if usedata ==1
%     exc1=find(ismember(PFCindsExc,forFigIndsexc(1,:),'rows')); % OR, exc1=find(ismember(excPFCinds,forFigInds(1,:),'rows'));
%     exc2=find(ismember(PFCindsExc,forFigIndsexc(2,:),'rows')); % OR, exc1=find(ismember(excPFCinds,forFigInds(1,:),'rows'));
%     inh1=find(ismember(PFCindsInh,forFigIndsinh(1,:),'rows'));
%     inh2=find(ismember(PFCindsInh,forFigIndsinh(2,:),'rows'));
% end

%% another look at the SWR-modulated PFC cells, and CA1 cells
% CA1 all ripmod cells
[tmp timingindCA1a]=min(inhCA1psths(:,twin)');
[tmp timingindCA1b]=max(excCA1psths(:,twin)');
[A1 B1]=sort(timingindCA1a);
[A2 B2]=sort(timingindCA1b);
figure;
subplot(2,1,1)
imagesc(xaxis,1:length([B1';B2']), [inhCA1psths(B1,:); excCA1psths(B2,:)]);caxis([-3 3])
xlabel('Time (ms)')
ylabel('#Cell')
title(sprintf('CA1 all CA1rip-mod cells - %dINH %dEXC',length(inhCA1psths(:,1)), length(excCA1psths(:,1))))
xlim([-500 500])
hold on
subplot(2,1,2); hold on
boundedline(xaxis,mean(excCA1psths,1),std(excCA1psths)./sqrt(size(excCA1psths,1)),'-r');
boundedline(xaxis,mean(inhCA1psths,1),std(inhCA1psths)./sqrt(size(inhCA1psths,1)),'-b');
xlim([-500 500])
ylabel('Mean z-scored psth')
xlabel('Time (ms)')
set(gcf, 'renderer', 'painters')

figure;
subplot(2,1,1)
[tmp timingindPFCa]=min(inhPFCpsths(:,twin)');
[tmp timingindPFCb]=max(excPFCpsths(:,twin)');
[A1 B1]=sort(timingindPFCa);
[A2 B2]=sort(timingindPFCb);
subplot(2,1,1)
imagesc(xaxis,1:length([B1';B2']), [inhPFCpsths(B1,:); excPFCpsths(B2,:)]);caxis([-3 3])
xlabel('Time (ms)')
ylabel('#Cell')
title(sprintf('PFC all CA1rip-mod cells - %dINH %dEXC',length(inhPFCpsths(:,1)), length(excPFCpsths(:,1))))
xlim([-500 500])
hold on
subplot(2,1,2); hold on
boundedline(xaxis,mean(excPFCpsths,1),std(excPFCpsths)./sqrt(size(excPFCpsths,1)),'-r');
boundedline(xaxis,mean(inhPFCpsths,1),std(inhPFCpsths)./sqrt(size(inhPFCpsths,1)),'-b');
xlim([-500 500])
ylabel('Mean z-scored psth')
xlabel('Time (ms)')
set(gcf, 'renderer', 'painters')

keyboard
if savefig1==1,
    print('-dpdf', figfile); print('-dpng', figfile, '-r300'); saveas(gcf,figfile,'fig'); print('-depsc2', figfile); print('-djpeg', figfile);
end

% With CA1 (also looking only at SWR-excited)
excinhCA1=allallripmodhiststype{1};
excCA1psths=CA1psths(excinhCA1>0,:);
inhCA1psths=CA1psths(excinhCA1==0,:);
[tmp timingindCA1]=max(excCA1psths(:,twin)');
[A2 B2]=sort(timingindCA1);
figure;
subplot(4,1,1)
imagesc(xaxis,1:length(B2),excCA1psths(B2,:));caxis([-3 3])
xlabel('Time (ms)')
ylabel('#Cell')
title('CA1 NC-ctxRip-excited cells')
xlim([-1000 1000])

%look at CA1 inh
[tmp timingindCA1inh]=min(inhCA1psths(:,twin)');
[A2 B2]=sort(timingindCA1inh);
subplot(4,1,2)
imagesc(xaxis,1:length(B2),inhCA1psths(B2,:));caxis([-3 3])
xlabel('Time (ms)')
ylabel('#Cell')
title('CA1 NC-ctxRip-inhibited cells')
xlim([-1000 1000])


% another look at the SWR-excited PFC cells
[tmp timingind]=max(excPFCpsths(:,twin)');
[A B]=sort(timingind);
subplot(4,1,3);imagesc(xaxis,1:length(B),excPFCpsths(B,:));caxis([-3 3])
xlabel('Time (ms)')
ylabel('#Cell')
title('PFC NC-ctxRip-excited cells')



% if usedata==1
%     % find indices after sorting
%     exc1s = find(B==exc1),
%     exc2s = find(B==exc2),
% end

% another look at the SWR-inhibited PFC cells
[tmp timingind2]=min(inhPFCpsths(:,twin)');
[A2 B2]=sort(timingind2);
xlim([-1000 1000])
subplot(4,1,4);imagesc(xaxis,1:length(B2),inhPFCpsths(B2,:));caxis([-3 3])
xlim([-1000 1000])
xlabel('Time (ms)')
ylabel('#Cell')
title('PFC NC-ctxRip-inhibited cells')

% if usedata==1
%     % find indices after sorting
%     inh1s = find(B2==inh1),
%     inh2s = find(B2==inh2),
% end

% figfile = [figdir,'RipTrig_PoplnZscore',filestr]
% if savefig1==1,
%     print('-dpdf', figfile); print('-dpng', figfile, '-r300'); saveas(gcf,figfile,'fig'); print('-depsc2', figfile); print('-djpeg', figfile);
% end


%STACKED PLOTS

[tmp timingindCA1inh]=min(inhCA1psths(:,twin)');
[A2 B2]=sort(timingindCA1inh);

% figure; yyaxis left; imagesc(xaxis,1:length(B2),inhCA1psths(B2,:));caxis([-3 3])
% yyaxis right
% hold on; shadedErrorBar(xaxis,mean(inhCA1psths,1),std(inhCA1psths)./sqrt(size(inhCA1psths,1)),'-b',1)
% x = [-1000 1000]; y = [0 0];
% plot(x,y,'--k','LineWidth',2)
% xlim([-1000 1000])
% yyaxis right; ylabel('Z-scored Response')
% yyaxis left
% ylabel('Cell #')
% xlabel('Time from ripple onset (ms)')

figure; yyaxis left; imagesc(xaxis,1:length(B2),inhCA1psths(B2,:));caxis([-3 3])
yyaxis right
hold on; boundedline(xaxis,mean(inhCA1psths,1),std(inhCA1psths)./sqrt(size(inhCA1psths,1)),'-b')
x = [-1000 1000]; y = [0 0];
plot(x,y,'--k','LineWidth',2)
xlim([-1000 1000])
yyaxis right; ylabel('Z-scored Response')
yyaxis left
ylabel('Cell #')
xlabel('Time from ripple onset (ms)')
%%

[tmp timingindCA1exc]=max(excCA1psths(:,twin)');
[A2 B2]=sort(timingindCA1exc);

% figure; yyaxis left; imagesc(xaxis,1:length(B2),excCA1psths(B2,:));caxis([-3 3])
% yyaxis right
% hold on; shadedErrorBar(xaxis,mean(excCA1psths,1),std(excCA1psths)./sqrt(size(excCA1psths,1)),'-r',1)
% x = [-1000 1000]; y = [0 0];
% plot(x,y,'--k','LineWidth',2)
% xlim([-1000 1000])
% yyaxis right; ylabel('Z-scored Response')
% yyaxis left
% ylabel('Cell #')
% xlabel('Time from ripple onset (ms)')

figure; yyaxis left; imagesc(xaxis,1:length(B2),excCA1psths(B2,:));caxis([-3 3])
yyaxis right
hold on; boundedline(xaxis,mean(excCA1psths,1),std(excCA1psths)./sqrt(size(excCA1psths,1)),'-r')
x = [-1000 1000]; y = [0 0];
plot(x,y,'--k','LineWidth',2)
xlim([-1000 1000])
yyaxis right; ylabel('Z-scored Response')
yyaxis left
ylabel('Cell #')
xlabel('Time from ripple onset (ms)')


[tmp timingind]=max(excPFCpsths(:,twin)');
[A B]=sort(timingind);

% figure; yyaxis left; imagesc(xaxis,1:length(B),excPFCpsths(B,:));caxis([-3 3])
% yyaxis right
% hold on; shadedErrorBar(xaxis,mean(excPFCpsths,1),std(excPFCpsths)./sqrt(size(excPFCpsths,1)),'-r',1)
% x = [-1000 1000]; y = [0 0];
% plot(x,y,'--k','LineWidth',2)
% xlim([-1000 1000])
% yyaxis right; ylabel('Z-scored Response')
% yyaxis left
% ylabel('Cell #')
% xlabel('Time from ripple onset (ms)')

figure; yyaxis left; imagesc(xaxis,1:length(B),excPFCpsths(B,:));caxis([-3 3])
yyaxis right
hold on; boundedline(xaxis,mean(excPFCpsths,1),std(excPFCpsths)./sqrt(size(excPFCpsths,1)),'-r')
x = [-1000 1000]; y = [0 0];
plot(x,y,'--k','LineWidth',2)
xlim([-1000 1000])
yyaxis right; ylabel('Z-scored Response')
yyaxis left
ylabel('Cell #')
xlabel('Time from ripple onset (ms)')

[tmp timingind]=min(inhPFCpsths(:,twin)');
[A B]=sort(timingind);

figure; yyaxis left; imagesc(xaxis,1:length(B),inhPFCpsths(B,:));caxis([-3 3])
yyaxis right
hold on; boundedline(xaxis,mean(inhPFCpsths,1),std(inhPFCpsths)./sqrt(size(inhPFCpsths,1)),'-b')
x = [-1000 1000]; y = [0 0];
plot(x,y,'--k','LineWidth',2)
xlim([-1000 1000])
yyaxis right; ylabel('Z-scored Response')
yyaxis left
ylabel('Cell #')
xlabel('Time from ripple onset (ms)')

figure
shadedErrorBar(xaxis,mean(inhCA1psths,1),std(inhCA1psths)./sqrt(size(inhCA1psths,1)),'-b',1);
hold on
shadedErrorBar(xaxis,mean(excCA1psths,1),std(excCA1psths)./sqrt(size(excCA1psths,1)),'-k',0);
xlim([-1000 1000])
title('CA1 mod');

figure
shadedErrorBar(xaxis,mean(inhPFCpsths,1),std(inhPFCpsths)./sqrt(size(inhPFCpsths,1)),'-k',1);
hold on
shadedErrorBar(xaxis,mean(excPFCpsths,1),std(excPFCpsths)./sqrt(size(excPFCpsths,1)),'-r',0);
xlim([-1000 1000])
title('PFC mod');

figure
% shadedErrorBar(xaxis,mean(PFCpsths,1),std(PFCpsths)./sqrt(size(PFCpsths,1)),'-r',1);
boundedline(xaxis,mean(PFCpsths,1),std(PFCpsths)./sqrt(size(PFCpsths,1)),'-r');
hold on
% shadedErrorBar(xaxis,mean(CA1psths,1),std(CA1psths)./sqrt(size(CA1psths,1)),'-b',1);
boundedline(xaxis,mean(CA1psths,1),std(CA1psths)./sqrt(size(CA1psths,1)),'-b');

xlim([-1000 1000])
title('Overall Mod');
%% unsorted 2
% SJ - Comment Fig

% figure;subplot(2,1,1);imagesc(xaxis,1:size(excPFCpsths,1),excPFCpsths);caxis([-3 3])
% xlabel('Time (ms)')
% ylabel('#Cell')
% title('PFC SWR-excited cells')
% subplot(2,1,2);imagesc(xaxis,1:size(inhPFCpsths,1),inhPFCpsths);caxis([-3 3])
% xlabel('Time (ms)')
% ylabel('#Cell')
% title('PFC SWR-inhibited cells')
% % With CA1 (also looking only at SWR-excited)
% excinhCA1=allallripmodhiststype{1};
%
% excCA1psths=CA1psths(excinhCA1>0,:);
% inhCA1psths=CA1psths(excinhCA1==0,:);
%
%
% [tmp timingindCA1]=max(excCA1psths(:,251:750)');
% [A2 B2]=sort(timingindCA1);
% figure;
% subplot(2,1,1)
% imagesc(xaxis,1:length(B2),excCA1psths(B2,:));caxis([-3 3])
% xlabel('Time (ms)')
% ylabel('#Cell')
% title('CA1 SWR-excited cells')
% xlim([-400 400])
%
% subplot(2,1,2)
% imagesc(xaxis,1:length(B),excPFCpsths(B,:));caxis([-3 3])
% xlabel('Time (ms)')
% ylabel('#Cell')
% xlim([-400 400])

keyboard

%% peak timing
bins1=-1000:2:1000;
[tmp timingindCA1exc]=max(excCA1psths(:,twin)');
[tmp timingindCA1inh]=min(inhCA1psths(:,twin)');

[tmp timingindPFCexc]=max(excPFCpsths(:,twin)');
[tmp timingindPFCinh]=min(inhPFCpsths(:,twin)');

% [pfc1 t]=hist(timingindPFCexc-250,bins1);
% [pfc2 t]=hist(timingindPFCinh-250,bins1);
% 
% [CA11 t]=hist(timingindCA1exc-250,bins1);
% [CA12 t]=hist(timingindCA1inh-250,bins1);


[pfc1 t]=hist(timingindPFCexc-500,bins1);
[pfc2 t]=hist(timingindPFCinh-500,bins1);

[CA11 t]=hist(timingindCA1exc-500,bins1);
[CA12 t]=hist(timingindCA1inh-500,bins1);
figure
plot(bins1,CA11/sum(CA11),'r','linewidth',2)
hold on;
plot(bins1,CA12/sum(CA12),'b','linewidth',2)
plot(bins1,pfc1/sum(pfc1),'k','linewidth',2);
plot(bins1,pfc2/sum(pfc2),'m','linewidth',2);
% plot(bins1,pfc2/sum(pfc2),'linewidth',2)
ylabel('fraction of units')
xlabel('Peak time relative to PFCripple onset')

%%
figure
plot(bins1,(CA11+CA12)./sum(CA11+CA12),'k','linewidth',2)
hold on;
plot(bins1,(pfc1+pfc2)./sum(pfc1+pfc2),'r','linewidth',2);
ylabel('fraction of units')
xlabel('Peak time relative to PFCripple onset')
p5peak = ranksum([timingindCA1exc timingindCA1inh],[timingindPFCinh timingindPFCexc])
title(sprintf(['CA1 pyr-int IndPFCrip timing - p=' num2str(p5peak)]))
set(gcf, 'renderer', 'painters')

%%

p1peak = ranksum(timingindCA1exc,timingindPFCexc) %p=1.0758e-30 CA1exc-PFCexc JDS NCPFCrip
p2peak = ranksum(timingindCA1inh,timingindPFCexc) %p=0.002 CA1inh-PFCexc JDS NCPFCrip
p3peak = ranksum(timingindCA1exc,timingindCA1inh) %p=2.6907e-21 CA1exc-CA1inh JDS NCPFCrip
p4peak = ranksum(timingindCA1exc,timingindPFCexc) %p= CA1exc-PFCexc JDS
p5peak = ranksum([timingindCA1exc timingindCA1inh],[timingindPFCinh timingindPFCexc]) %p= CA1exc-CA1inh JDS
p6peak = ranksum(timingindCA1inh,timingindPFCinh)
% 
% 
% 
%% timing of start of rise to peak
risestartsPFCexc=[];
risestartsPFCinh=[];
risestartsCA1=[];

for i=1:size(excPFCpsths,1)
    [a maxind]=max(excPFCpsths(i,251:750));
    maxind=maxind+250;
    aa=excPFCpsths(i,1:maxind);
    risestart=find(aa<mean(excPFCpsths(i,:))+std(excPFCpsths(i,:)),1,'last');
    risestartsPFCexc=[risestartsPFCexc risestart];
end
for i=1:size(inhPFCpsths,1)
    [a maxind]=min(inhPFCpsths(i,251:750));
    maxind=maxind+250;
    aa=inhPFCpsths(i,1:maxind);
    risestart=find(aa>mean(inhPFCpsths(i,:))-std(inhPFCpsths(i,:)),1,'last');
    risestartsPFCinh=[risestartsPFCinh risestart];
end
for i=1:size(excCA1psths,1)
    [a maxind]=max(excCA1psths(i,251:750));
    maxind=maxind+250;
    aa=excCA1psths(i,1:maxind);
    risestart=find(aa<mean(excCA1psths(i,:))+std(excCA1psths(i,:)),1,'last');
    risestartsCA1=[risestartsCA1 risestart];
end

%figure;plot(risestarts);hold on;plot(risestartsCA1,'r')
[risePFCexc t]=hist(risestartsPFCexc-500,bins1);
[risePFCinh t]=hist(risestartsPFCinh-500,bins1);
[riseCA1 t]=hist(risestartsCA1-500,bins1);

figure
plot(bins1,riseCA1/sum(riseCA1),'k','linewidth',2)
hold on;
plot(bins1,risePFCexc/sum(risePFCexc),'r','linewidth',2);
% plot(bins1,risePFCinh/sum(risePFCinh),'linewidth',2)
ylabel('fraction of units')
xlabel('Rise time relative to SWR onset')
legend('CA1','PFC swr-exc');%,'PFC swr-inh')
% 
% figfile = [figdir,'RipTrig_RiseTime',filestr]
% if savefig1==1,
%     print('-dpdf', figfile); print('-dpng', figfile, '-r300'); saveas(gcf,figfile,'fig'); print('-depsc2', figfile); print('-djpeg', figfile);
% end
% 
% mean(risestartsPFCexc-500), sem(risestartsPFCexc-500),
% mean(risestartsPFCinh-500), sem(risestartsPFCinh-500),
% mean(risestartsCA1-500), sem(risestartsCA1-500),
% 
p1rise = ranksum(risestartsCA1,risestartsPFCexc) %p=
p2rise = ranksum(risestartsCA1,risestartsPFCinh) %p=
p3rise = ranksum(risestartsPFCexc,risestartsPFCinh) %p=
% 

inhPFCsuppression = mean(inhPFCpsths(:,801:1300)');
inhCA1suppression = mean(inhCA1psths(:,801:1300)');
p1suppression = ranksum(inhPFCsuppression,inhCA1suppression)
% 
% %% modulation duration
% risedurPFCexc=[];
% risedurPFCinh=[];
% risedurCA1=[];
% 
% for i=1:size(excPFCpsths,1)
%     [a maxind]=max(excPFCpsths(i,251:750));
%     maxind=maxind+250;
%     aa=excPFCpsths(i,1:maxind);
%     risestart=find(aa<mean(excPFCpsths(i,:))+std(excPFCpsths(i,:)),1,'last');
%     bb=excPFCpsths(i,maxind:end);
%     riseend=maxind+find(bb<mean(excPFCpsths(i,:))+std(excPFCpsths(i,:)),1,'first');
%     curdur=riseend-risestart;
%     risedurPFCexc=[risedurPFCexc curdur];
% end
% 
% for i=1:size(inhPFCpsths,1)
%     [a maxind]=min(inhPFCpsths(i,251:750));
%     maxind=maxind+250;
%     aa=inhPFCpsths(i,1:maxind);
%     risestart=find(aa>mean(inhPFCpsths(i,:))-std(inhPFCpsths(i,:)),1,'last');
%     bb=inhPFCpsths(i,maxind:end);
%     riseend=maxind+find(bb>mean(inhPFCpsths(i,:))-std(inhPFCpsths(i,:)),1,'first');
%     curdur=riseend-risestart;
%     risedurPFCinh=[risedurPFCinh curdur];
% end
% 
% for i=1:size(excCA1psths,1)
%     [a maxind]=max(excCA1psths(i,251:750));
%     maxind=maxind+250;
%     aa=excCA1psths(i,1:maxind);
%     risestart=find(aa<mean(excCA1psths(i,:))+std(excCA1psths(i,:)),1,'last');
%     bb=excCA1psths(i,maxind:end);
%     riseend=maxind+find(bb<mean(excCA1psths(i,:))+std(excCA1psths(i,:)),1,'first');
%     curdur=riseend-risestart;
%     risedurCA1=[risedurCA1 curdur];
% end
% 
% bins2=0:10:250;
% 
% %figure;plot(risestarts);hold on;plot(risestartsCA1,'r')
% [risePFCexc t]=hist(risedurPFCexc,bins2)
% [risePFCinh t]=hist(risedurPFCinh,bins2)
% [riseCA1 t]=hist(risedurCA1,bins2)
% 
% figure
% plot(bins2,riseCA1/sum(riseCA1),'k','linewidth',2)
% hold on;
% plot(bins2,risePFCexc/sum(risePFCexc),'r','linewidth',2);
% plot(bins2,risePFCinh/sum(risePFCinh),'linewidth',2)
% ylabel('fraction of units')
% xlabel('SWR modulation duration (ms)')
% legend('CA1','PFC swr-exc','PFC swr-inh')
% 
% figfile = [figdir,'RipTrig_ModlnDurn',filestr]
% if savefig1==1,
%     print('-dpdf', figfile); print('-dpng', figfile, '-r300'); saveas(gcf,figfile,'fig'); print('-depsc2', figfile); print('-djpeg', figfile);
% end
% 
% 
% 
% mean(risedurPFCexc-500), sem(risedurPFCexc-500),
% mean(risedurPFCinh-500), sem(risedurPFCinh-500),
% mean(risedurCA1-500), sem(risedurCA1-500),
% 
% p1 = ranksum(risedurCA1,risedurPFCexc) %p=
% p2 = ranksum(risedurCA1,risedurPFCinh) %p=
% p3 = ranksum(risedurPFCexc,risedurPFCinh) %p=
