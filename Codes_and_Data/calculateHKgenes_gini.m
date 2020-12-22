function [GINI,hkGenes,validation,f] = calculateHKgenes_gini(expressionData,remGenesFlag,gcCut,knownGenes)

if nargin < 3
    gcCut = [];
    knownGenes = [];
end
if isnumeric(expressionData.gene)
    expressionData.gene = cellstr(char(num2str(expressionData.gene)));
    expressionData.gene = strrep(expressionData.gene,' ','');
end
if isfield(expressionData,'value')
    expressionData.valuebyTissue = expressionData.value;
    clear expressionData.value
end
% expressionData.valuebyTissue(expressionData.valuebyTissue==0) = 1e-1;
gFields = getFieldsSameSize(expressionData,'gene');
gFields = gFields.sameRows;
if remGenesFlag==1
    rem_genes = sum(expressionData.valuebyTissue==0,2)==length(expressionData.Tissue);
    expressionData.gene(rem_genes) = [];
    for i=1:length(gFields)
        eval(strcat('expressionData.',gFields{i},'(rem_genes,:) = [];'));
    end
elseif remGenesFlag==2
    rem_genes = median(expressionData.valuebyTissue,2)==0;
    expressionData.gene(rem_genes) = [];
    for i=1:length(gFields)
        eval(strcat('expressionData.',gFields{i},'(rem_genes,:) = [];'));
    end
elseif remGenesFlag==3
    rem_genes = sum(expressionData.valuebyTissue==0,2)==length(expressionData.Tissue) | ...
        median(expressionData.valuebyTissue,2)==0;
    expressionData.gene(rem_genes) = [];
    for i=1:length(gFields)
        eval(strcat('expressionData.',gFields{i},'(rem_genes,:) = [];'));
    end
elseif remGenesFlag==0
    rem_genes = [];
end
expressionData.valuebyTissue(expressionData.valuebyTissue==0) = min(expressionData.valuebyTissue(expressionData.valuebyTissue~=0));
fprintf('No. of genes removed from analysis = %d\n',sum(rem_genes));
my_expr = expressionData.valuebyTissue;
f = sum(my_expr,1)'*100/sum(sum(my_expr));
my_genes = expressionData.gene;

gc = ones(length(my_genes),1);
for i=1:length(my_genes)
    %     keep_samp = my_expr(i,:)~=0;
    %     geneExpr = my_expr(i,keep_samp);
    %     fi = f(keep_samp);
    %     gc(i,1) = gini(fi,geneExpr,false);
    gc(i,1) = gini(f,my_expr(i,:),false);
end
% gc = (gc-min(gc))/max(gc-min(gc)); % normalize gini coefficient
if ~isempty(knownGenes) || ~isempty(gcCut)
    figure;
    hold on
    median_expr = median(my_expr,2);
    if ~isempty(knownGenes)
        scatter(log10(median_expr(~ismember(my_genes,knownGenes))),gc(~ismember(my_genes,knownGenes)),30,'marker','o','markerfacecolor',[0.5 0.5 0.5],...
            'markeredgecolor','none','markerfacealpha',0.4);
        scatter(log10(median_expr(ismember(my_genes,knownGenes))),gc(ismember(my_genes,knownGenes)),5,'marker','o','markerfacecolor','k',...
            'markeredgecolor','k');
    else
        scatter(log10(median_expr),gc,30,'marker','o','markerfacecolor',[0.5 0.5 0.5],...
            'markeredgecolor','none','markerfacealpha',0.4);
        if ~isempty(gcCut)
            scatter(log10(median_expr(gc <= prctile(gc,gcCut))),gc(gc <= prctile(gc,gcCut)),30,'marker','o','markerfacecolor',[1 0.8 0],...
                'markeredgecolor','none');
        end
    end
    ft = fittype('m*x+C');
    [L,gof] = fit(log10(median_expr),gc,ft,'StartPoint',[-1 1]);
    plot(L,'--r'); ylim([0 1]);
    if ~isempty(knownGenes)
        legend('Data','known housekeeping genes',strcat(num2str(L.m),'x+(',num2str(L.C),');R^2 = ',num2str(gof.rsquare)));
    else
        legend('Data',strcat(num2str(L.m),'x+(',num2str(L.C),');R^2 = ',num2str(gof.rsquare)));
    end
    xlabel(['Median expression (log_{10}' expressionData.unit ')']);
    ylabel('Gini coeefficient');
    hold off
    % figure;
    % histogram(gc,0.005:0.01:1.005);
    knownGenes = intersect(knownGenes,expressionData.gene);
    if ~isempty(gcCut)
        [hkGenes,gene_gc,isknown,~,~,p] = findHK(my_genes,my_expr,gc,gcCut,knownGenes,true,expressionData.unit);
        GINI.threshold = p;
        % fix this section
        for i=1:length(hkGenes)
            C{i,1} = find(ismember(expressionData.gene,hkGenes{i,1}));
        end
        hkGenes(cellfun(@length,C)==2) = [];
        % % %
        fprintf('No. of genes known and also predicted: %d\n',length(intersect(hkGenes,knownGenes)));
        fprintf('No. of new genes predicted: %d\n',length(setdiff(hkGenes,knownGenes)));
        validation.hkGenes = hkGenes;
        validation.giniCoeff = gene_gc; validation.isknown = isknown;
    else
        validation = [];
        gcCut = 2:2:100;
        common_hk = zeros(length(gcCut),1); both = common_hk; onlyPred = common_hk; accuracy = common_hk;
        pred_hk = common_hk; data_hk = common_hk;
        for i=1:length(gcCut)
            [hkGenes,~,~,M(i),C(i)] = findHK(my_genes,my_expr,gc,gcCut(i),knownGenes,false,expressionData.unit);
            common_hk(i) = length(intersect(hkGenes,knownGenes));
            pred_hk(i) = length(setdiff(hkGenes,knownGenes));
            data_hk(i) = length(setdiff(knownGenes,hkGenes));
            both(i) = length(intersect(hkGenes,knownGenes))/length(knownGenes);
            onlyPred(i) = length(setdiff(hkGenes,knownGenes))/length(hkGenes);
            accuracy(i) = length(intersect(hkGenes,knownGenes))/length(hkGenes);
        end
        if ~isempty(knownGenes)
            fprintf('GC to cover at least 90%% of known genes = %0.2f\n',gcCut(find(both>=0.9,1)));
            figure;
            subplot(1,2,1)
            plot(gcCut,[common_hk pred_hk data_hk]);
            grid on
            legend('both','only GINI','only known');
            xlabel('Gini coefficient cutoff');
            ylabel('number of genes');
            subplot(1,2,2)
            plot(gcCut,[both onlyPred accuracy]);
            grid on
            legend('fraction captured','fraction new (FP)','accuracy x 10^{-2} (TP)')
            xlabel('Gini coefficient cutoff');
            figure;
            plot(gcCut(2:end),[accuracy(1:end-1)-accuracy(2:end) both(2:end)-both(1:end-1)]);
            legend('rate of change in accuracy','rate of inclusion in both');
            figure;
            plot(gcCut,[M; C]');
            legend('slope of fit','intercept of fit');
        end
    end
end
GINI.gc = gc; GINI.genes = my_genes;
for i=1:length(gFields)
    eval(strcat('GINI.',gFields{i},'=expressionData.',gFields{i},';'));
end
if isfield(expressionData,'source')
    GINI.source = expressionData.source;
end
if isfield(expressionData,'geneName')
    GINI.gene_name = GINI.geneName;
end

function [hkGenes,gene_gc,isknown,M,C,p] = findHK(my_genes,my_expr,gc,gcCut,knownGenes,plotFlag,unit)

median_expr = median(my_expr,2);
max_expr = max(my_expr,[],2);
if gcCut==0
    ihk = gc == min(gc);
elseif gcCut < 1
    ihk = gc <= gcCut;
else
    p = prctile(gc,gcCut);
    ihk = gc <= prctile(gc,gcCut);
    fprintf('%0.1f percentile Gini coefficient = %0.4f\n',gcCut,prctile(gc,gcCut));
end
if ~isempty(knownGenes)
    iknown = ismember(my_genes,knownGenes);
    ft = fittype('m*x+c');
    [L,gof] = fit(log10(max_expr(ihk)),log10(median_expr(ihk)),ft,'StartPoint',[1 0.5]);
    M = L.m;
    C = L.c;
else
    ft = fittype('m*x+c');
    [L,gof] = fit(log10(max_expr(ihk)),log10(median_expr(ihk)),ft,'StartPoint',[1 0.5]);
    M = L.m;
    C = L.c;
    iknown = false(length(my_genes),1);
end
hkGenes = my_genes(ihk);
gene_gc = [];
isknown = [];
if plotFlag
    gene_gc = gc(ihk);
    isknown = ismember(hkGenes,knownGenes);
    figure; hold on;
    scatter(log10(max_expr),log10(median_expr),10,'marker','o','markerfacecolor',[0.5 0.5 0.5],...
        'markeredgecolor','none','markerfacealpha',0.4);
    plot(L,'--g');
    if ~isempty(knownGenes)
        scatter(log10(max_expr(iknown & ~ihk)),log10(median_expr(iknown & ~ihk)),10,'marker','o','markerfacecolor',[1 0 0],...
            'markeredgecolor','none','markerfacealpha',1);
        ssd1 = sum((log10(median_expr(iknown & ~ihk))-L(log10(max_expr(iknown & ~ihk)))).^2);
        scatter(log10(max_expr(iknown & ihk)),log10(median_expr(iknown & ihk)),10,'marker','o','markerfacecolor',[0 1 0],...
            'markeredgecolor','none','markerfacealpha',1);
        ssd0 = sum((log10(median_expr(iknown & ihk))-L(log10(max_expr(iknown & ihk)))).^2);
    end
    scatter(log10(max_expr(ihk & ~iknown)),log10(median_expr(ihk & ~iknown)),10,'marker','o','markerfacecolor',[1 0.8 0],...
        'markeredgecolor','none');
    ssd2 = sum((log10(median_expr(~iknown & ihk))-L(log10(max_expr(~iknown & ihk)))).^2);
    plot([0 5],[0 5],'--k');
    set(gca,'dataaspectratiomode','manual');
    hold off;
    if ~isempty(knownGenes)
        legend('other genes',strcat(num2str(L.m),'x+(',num2str(L.c),');R^2 = ',num2str(gof.rsquare)),...
            strcat('known but not predicted genes (SSE=',num2str(ssd1),')'),...
            strcat('known and predicted genes (',num2str(roundn(sum(ihk & iknown)*100/sum(ihk),-2)),'%, SSE=',num2str(gof.sse),')'),...
            strcat('not known but predicted genes (SSE=',num2str(ssd2),')'),'X = Y','location','northwest');
    else
        legend('other genes',strcat(num2str(L.m),'x+(',num2str(L.c),');R^2 = ',num2str(gof.rsquare)),...
            'predicted genes','X = Y','location','northwest');
    end
    ylabel(['Median expression (log_{10}' unit ')']);
    xlabel(['Max expression (log_{10}' unit ')']);
    grid on
    if ~isempty(knownGenes)
        figure;
        hold on;
        fprintf('SSE for known housekeeping genes which were also predicted = %0.8f\n',ssd0);
        plot(log10(max_expr(iknown & ~ihk)),log10(median_expr(iknown & ~ihk)),'^','color',[0.8 0 0])
        fprintf('SSE for known housekeeping genes which were not predicted = %0.8f\n',ssd1);
        plot(log10(max_expr(ihk & ~iknown)),log10(median_expr(ihk & ~iknown)),'o','color',[0.1 0.5 0.8])
        fprintf('SSE for predicted housekeeping genes which were not known = %0.8f\n',ssd2);
        fprintf('SE for known housekeeping genes which were not predicted = %0.8f\n',...
            sum((log10(median_expr(iknown & ~ihk))-L(log10(max_expr(iknown & ~ihk))))));
        plot([-1 5],L([-1 5]),'marker','none','color',[0.5 0.8 0]);
        plot([-1 5],[-1 5],'--k');
        ylabel(['Median expression (log_{10}' unit ')']);
        xlabel(['Max expression (log_{10}' unit ')']);
        hold off;
        legend('known but not predicted','not known but predicted',strcat(num2str(L.m),'x+(',num2str(L.c),');R^2 = ',num2str(gof.rsquare)),'X = Y','location','northwest');
    end
end