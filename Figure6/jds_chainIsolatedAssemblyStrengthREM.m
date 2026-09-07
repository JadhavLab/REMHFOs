function jds_chainIsolatedAssemblyStrengthREM(animalprefixlist,area,state)
%JDS_CHAINISOLATEDASSEMBLYSTRENGTHREM Assembly reactivation in chain vs isolated ripples.
%   jds_chainIsolatedAssemblyStrengthREM(animalprefixlist,area,state)
%   compares peak z-scored assembly reactivation strength within chain
%   cortical ripples vs isolated cortical ripples during REM (sign-rank,
%   boxplot).
%
%   animalprefixlist - cell array of animal prefix strings
%   area             - brain area of the reactivation file, 'CA1' or 'PFC'
%   state            - sleep-state label in the reactivation filename

day = 1;

isolated_react = [];
chain_react = [];

for a = 1:length(animalprefixlist)

    animalprefix = char(animalprefixlist(a));
    dir = sprintf('/Volumes/JUSTIN/SingleDay/%s_direct/',animalprefix);

    %Load reactivation strength file for all assemblies and epochs
    load(sprintf('%s%s%s_RTimeStrength%sNewSpk_20_%02d.mat',dir,animalprefix,area,state,day));
    %Load ripples
    load(sprintf('%s%sctxrippletime_chainREM%02d.mat',dir,animalprefix,day));
    load(sprintf('%s%sremeps%02d.mat',dir,animalprefix,day));

    epochs = remeps;

    for e = 1:length(epochs)
        ep = epochs(e);
        if ep == 1
            continue
        end
        assemblytmp = RtimeStrength{ep}.reactivationStrength;
        chaintimes = [ctxripple{day}{ep}.starttimeC ctxripple{day}{ep}.endtimeC];

        isotimes = [ctxripple{day}{ep}.starttimeNC ctxripple{day}{ep}.endtimeNC];

        %Chain ripples
        if ~isempty(assemblytmp)
            for ii = 1:size(assemblytmp,2)
                react_idx = [];
                for t = 1:size(chaintimes,1)
                    sttmp = lookup(chaintimes(t,1), assemblytmp{ii}(:,1));
                    endtmp = lookup(chaintimes(t,2), assemblytmp{ii}(:,1));
                    react_idx = [react_idx; [sttmp endtmp]];
                end
                atmp = [];
                strengthstmp = zscore(assemblytmp{ii}(:,2));
                for r = 1:size(react_idx,1)
                    tmp = strengthstmp(react_idx(r,1):react_idx(r,2)); %get vector of reactivation strenths for specified time period
                    atmp = [atmp; max(tmp)];
                end
                chain_react = [chain_react; mean(atmp)];
            end
            %Isolated ripples
            for ii = 1:size(assemblytmp,2)
                react_idx = [];
                for t = 1:size(isotimes,1)
                    sttmp = lookup(isotimes(t,1), assemblytmp{ii}(:,1));
                    endtmp = lookup(isotimes(t,2), assemblytmp{ii}(:,1));
                    react_idx = [react_idx; [sttmp endtmp]];
                end
                atmp = [];
                strengthstmp = zscore(assemblytmp{ii}(:,2));
                for r = 1:size(react_idx,1)
                    tmp = strengthstmp(react_idx(r,1):react_idx(r,2)); %get vector of reactivation strenths for specified time period
                    atmp = [atmp; max(tmp)];
                end
                isolated_react = [isolated_react; mean(atmp)];
            end
        end
    end
end
[p h] = signrank(chain_react,isolated_react)

datacombinedReact = [chain_react; isolated_react];
g1 = repmat({'Chain'},length(chain_react),1);
g2 = repmat({'Isolated'},length(isolated_react),1);
g = [g1;g2];

figure;
h = boxplot(datacombinedReact,g,'OutlierSize',7,'Symbol','k+'); set(h(7,:),'Visible','off');
title(['Reactivation strength-p = ' num2str(p)])
ylabel('Reactivation strength (AU)')
set(gcf, 'renderer', 'painters')
