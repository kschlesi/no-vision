% script for doing flexibility / time-dependent analysis on the synths.

% for each saveString
strList = {'s8nrun_g1o1',...
           'arun20_g1o1',...
           's8nrun_g1o1_rem8_1',...
           's8nrun_g1o1_rem8',...
           'aruns8_g1o1',...
           's8nrun_g1o1_rem20_1',...
           's8nrun_g1o1_rem20_2',...
           's8nrun_g1o1_rem20',...
           };

for i=1:numel(strList)
    % set paramString
    paramString = char(strList{i});
    disp(['Calculating for ' num2str(i) ': ' paramString]);

    % load all structural params, dynamic params, CD params
    % leave out A_ens, Cens, Rass, Rg
    load(['Results/' paramString '.mat'],'-regexp','^((?!A_ens|Cens).)*$');
    
    % calculate flexibilities over time windows
    partnS = squeeze(mode(Rg,1));
    flexS = zeros(nR,sims); 
    for s=1:sims
        flexS(:,s) = flexibility(partnS(:,:,s)',0);
    end
    
    % histogram flexibilities and plot averages
    figure; hist(flexS); title(insertBefore(paramString,'_','\'));
    figure; plot(mean(flexS,2),'o','MarkerFaceColor',next_color); hold on;
            flexPlot = NaN.*ones(size(flexS));
            flexPlot(~~flexS) = flexS(~~flexS);
            plot(flexPlot,'x'); 
            title(insertBefore(paramString,'_','\'));
            legend([{'mean'},cellfun(@num2str,num2cell(1:sims),'UniformOutput',false)]);
            hold off;
            
    % plot ROC curves for each tw
    figure; hold on;
        %ROC_figure(mean(Rass(:,:,t)),N,M,m,toRemove,sims,'ExtFigureCmd');
        tpList = zeros(sims+1,T);
        fpList = zeros(sims+1,T);
        for t=1:T
            %subplot(2,ceil(T/2),t)
            [tp,fp] = ROC_figure(Rass(:,:,t),N,M,m,toRemove,sims,...
                      'ExtFigureCmd','plotargs',{':o'});
            tpList(:,t) = tp;
            fpList(:,t) = fp;
        end
        plot(mean(fpList,2),mean(tpList,2),'o',...
                                           'MarkerFaceColor','k',...
                                           'MarkerEdgeColor','None');
        title(insertBefore(paramString,'_','\'));
        legend(cellfun(@num2str,num2cell(1:T),'UniformOutput',false),...
               'location','southeast');
        hold off;
end