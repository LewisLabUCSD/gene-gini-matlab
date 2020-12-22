function [sigTable,pvalues,nums,ugos,x] = hypergeo_GOtable(baseGO,GO,ichip,ihit,pCutoff,rCutoff)

chipGO = createGenesetGO(baseGO,ichip);
hitGO = createGenesetGO(baseGO,ihit);
ugos = unique(chipGO.GOid);

pvalues = ones(length(ugos),1);
chipnums = zeros(length(ugos),1); hitnums = zeros(length(ugos),1);
for i=1:length(ugos)
    if ~isempty(find(chipGO.GOid==ugos(i)))
        chipnums(i,1) = length(find(chipGO.GOid==ugos(i)));
    end
    if ~isempty(find(hitGO.GOid==ugos(i)))
        hitnums(i,1) = length(find(hitGO.GOid==ugos(i)));
    end
end
nums = [hitnums chipnums];
MNhits = hygestat(length(chipGO.GOid),chipnums,length(hitGO.GOid));
FC = hitnums./MNhits;
enrich = FC >= 1;

pvalues(enrich) = 1 - hygecdf(hitnums(enrich),length(chipGO.GOid),chipnums(enrich),length(hitGO.GOid));
pvalues(~enrich) = hygecdf(hitnums(~enrich),length(chipGO.GOid),chipnums(~enrich),length(hitGO.GOid));
pvaluesBH = mafdr(pvalues,'bhfdr',true);

% pvalues = 1 - hygecdf(hitnums,length(chipGO.GOid),chipnums,length(hitGO.GOid));
% pvalues(~enrich) = 1 - pvalues(~enrich);
% pvaluesBH = mafdr(pvalues,'bhfdr',true);

kpind = hitnums~=0 & chipnums~=0 & pvalues~=0;
hitnums = hitnums(kpind); chipnums = chipnums(kpind);
ugos = ugos(kpind); 
pvalues = pvalues(kpind); pvaluesBH = pvaluesBH(kpind);
MNhits = MNhits(kpind); FC = FC(kpind); 

x = [hitnums chipnums MNhits pvalues pvaluesBH];
indicesEN = (pvaluesBH <= pCutoff & log2(FC)>rCutoff);
indicesDP = (pvaluesBH <= pCutoff & log2(FC)<=-rCutoff); 
indices = indicesEN | indicesDP; % my list of GO terms
figure;
scatter(log2(FC(~indices)),-log10(pvaluesBH(~indices)),4,...
    'Marker','o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','none');
hold on;
scatter(log2(FC(indicesEN)),-log10(pvaluesBH(indicesEN)),hitnums(indicesEN),...
    'Marker','o','MarkerFaceColor',[0.1 0.8 0],'MarkerEdgeColor','none','MarkerFaceAlpha',0.4);
scatter(log2(FC(indicesDP)),-log10(pvaluesBH(indicesDP)),hitnums(indicesDP),...
    'Marker','o','MarkerFaceColor',[0.8 0 0.1],'MarkerEdgeColor','none','MarkerFaceAlpha',0.4);
yl = ylim; xl = xlim;
plot([rCutoff rCutoff],[-log10(pCutoff) yl(2)],'--r');
plot([-rCutoff -rCutoff],[-log10(pCutoff) yl(2)],'--r');
plot([xl(1) -rCutoff],[-log10(pCutoff) -log10(pCutoff)],'--r');
plot([rCutoff xl(2)],[-log10(pCutoff) -log10(pCutoff)],'--r');
hold off;
xlabel('log_{2}(Ratio)');
ylabel('-log_{10}(P_{adj})');
set(gca,'DataAspectRatioMode','manual');
fprintf('%d terms are enriched and pass coverage criterion.\n',sum(indicesEN));
fprintf('%d terms are depleted and pass coverage criterion.\n',sum(indicesDP));

sigTable.enriched.GOid = ugos(indicesEN);
sigTable.enriched.pvalue = pvaluesBH(indicesEN);
sigTable.enriched.hitNums = hitnums(indicesEN);
sigTable.enriched.chipNums = chipnums(indicesEN);
sigTable.enriched.ratio= FC(indicesEN);

sigTable.depleted.GOid = ugos(indicesDP);
sigTable.depleted.pvalue = pvaluesBH(indicesDP);
sigTable.depleted.hitNums = hitnums(indicesDP);
sigTable.depleted.chipNums = chipnums(indicesDP);
sigTable.depleted.ratio= FC(indicesDP);

[~,ixEN] = sort(sigTable.enriched.hitNums,'descend');
[~,ixDP] = sort(sigTable.depleted.hitNums,'descend');

enriched = sigTable.enriched;
e = min([sum(indicesEN) 15]);
fprintf('GO Term         GO Name      p-val   counts\n');
for i=1:e
    fprintf('%d\t%s\t%1.4f\t%d / %d\n',...
        enriched.GOid(ixEN(i)),get(GO(enriched.GOid(ixEN(i))).Terms,'name'),enriched.pvalue(ixEN(i)),...
        enriched.hitNums(ixEN(i)),enriched.chipNums(ixEN(i)));
printGOnames{i,1} = get(GO(enriched.GOid(ixEN(i))).Terms,'name');
end
hold on;
% text(log2(enriched.ratio(ixEN(1:e))),-log10(enriched.pvalue(ixEN(1:e))),...
%     strcat(num2str(enriched.GOid(ixEN(1:e))),repmat({':'},e,1),printGOnames),'HorizontalAlignment','left');

printGOnames = [];
depleted = sigTable.depleted;
d = min([sum(indicesDP) 7]);
if d~=0
    fprintf('GO Term         GO Name      p-val   counts\n');
    for i=1:d
        fprintf('%d\t%s\t%1.4f\t%d / %d\n',...
            depleted.GOid(ixDP(i)),get(GO(depleted.GOid(ixDP(i))).Terms,'name'),depleted.pvalue(ixDP(i)),...
            depleted.hitNums(ixDP(i)),depleted.chipNums(ixDP(i)));
        printGOnames{i,1} = get(GO(depleted.GOid(ixDP(i))).Terms,'name');
    end
%     text(log2(depleted.ratio(ixDP(1:d))),-log10(depleted.pvalue(ixDP(1:d))),...
%         strcat(num2str(depleted.GOid(ixDP(1:d))),repmat({':'},d,1),printGOnames),'HorizontalAlignment','right');
end
hold off;


function genesGO = createGenesetGO(baseGO,indices)

genesGO = baseGO(indices,:);