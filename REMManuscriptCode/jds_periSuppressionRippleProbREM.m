
day = 1;
suppHistREM = [];
suppHistREM_s = [];

animalprefixlist = {'ZT2','JS17','JS15','JS14','JS12','JS13','JS34','BG1','JS21','KL8'};

% for m = 1:size(anims,1)
    for a = 1:length(animalprefixlist)
%     for a = anims(m,:)
        animalprefix = char(animalprefixlist(a));
        dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
        load(sprintf('%s%sPFCREMmuasuppression0%d.mat',dir,animalprefix,day));
        supp{a} = suppressiontimes;
        clear suppressiontimes
    end
    for s = 1:1000
        tmpShufProb = [];
        for a = 1:length(animalprefixlist)
%         for a = anims(m,:)

            animalprefix = char(animalprefixlist(a));

            dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);
            load(sprintf('%s%sctxrippletime_chainREM0%d.mat',dir,animalprefix,day));
            %         load(sprintf('%s%sPFCNREMmuasuppression0%d.mat',dir,animalprefix,day));

            load(sprintf('%s%sremeps0%d.mat',dir,animalprefix,day));
            epochs = remeps;
            for ep=1:length(epochs)

                epoch = epochs(ep);

                ctx = ctxripple{day}{epoch};
                suppr = supp{a}{day}{epoch};
                %             suppr = suppressiontimes{day}{epoch};
                supprTimes = suppr.troughtime{1};
                supprTimes_s = suppr.shuf{1}{s};
                currRips = ctx.starttimeC;
                if (~isempty(ctx))
                    if s == 1
                        allHist = [];
                        for r = 1:size(supprTimes,1)
                            currSup = supprTimes(r,1);
                            currRipsTmp =  currRips(find( (currRips>=(currSup-3))...
                                & (currRips<=(currSup+3)) ));
                            currRipsTmp = currRipsTmp-(currSup);
                            histRips = histc(currRipsTmp,[-3:0.1:3]);
                            allHist = [allHist; histRips(:).'];
                        end
                        suppHistREM = [suppHistREM; allHist];
                    end
                    allHist_s = [];
                    for r = 1:size(supprTimes_s,1)
                        currSup = supprTimes_s(r,1);
                        currRipsTmp =  currRips(find( (currRips>=(currSup-3))...
                            & (currRips<=(currSup+3)) ));
                        currRipsTmp = currRipsTmp-(currSup);
                        histRips = histc(currRipsTmp,[-3:0.1:3]);
                        allHist_s = [allHist_s; histRips(:).'];
                    end
                    tmpShufProb = [tmpShufProb; allHist_s];
                end
            end
        end
        suppHistREM_s = [suppHistREM_s; sum(tmpShufProb)./size(tmpShufProb,1)];
    end
    clear supp
% end

%%
upper = prctile(suppHistREM_s,99.5);
bottom = prctile(suppHistREM_s,0.5);

figure; hold on
bar(-30:30, sum(suppHistREM)./size(suppHistREM,1))
y = [0 0.1];
x = [0 0];
plot(x,y,'--k')
plot(-30:30, upper, '--r')
plot(-30:30, bottom, '--r')
% plot(-50:50, conv(sum(allHist)./size(allHist,1),g1,'same'))
xlim([-20 20])
xticks([-20:10:20])
xticklabels({'-2','-1','0','1','2'})
xlabel('Time from suppression trough (s)')
ylabel('PFC ripple probability')
set(gcf, 'renderer', 'painters')

keyboard
