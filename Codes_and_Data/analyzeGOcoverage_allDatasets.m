function [enriched,coverage] = analyzeGOcoverage_allDatasets

% load transcriptomics data
load(strcat(pwd,'\Data\transcriptomics\data_gtex.mat'));
load(strcat(pwd,'\Data\transcriptomics\data_hpa.mat'));
load(strcat(pwd,'\Data\transcriptomics\data_cm.mat'));
load(strcat(pwd,'\Data\transcriptomics\data_klijn.mat'));
% calculate gini index
gini_gtex = calculateHKgenes_gini(data_gtex,3,17.5,[]);
gini_hpa = calculateHKgenes_gini(data_hpa,3,17.5,[]);
gini_cm = calculateHKgenes_gini(data_cm,3,17.5,[]);
gini_klijn = calculateHKgenes_gini(data_klijn,3,17.5,[]);
load(strcat(pwd,'\Data\hkEisen.mat')); % load housekeeping list in Eisenberg 

% perform GO enrichement analysis
GO = geneont('live',true);
[sigTable_gtex,GOcover_gtex] = performGOenrichment(gini_gtex,GO,length(hkEisen));
[sigTable_hpa,GOcover_hpa] = performGOenrichment(gini_hpa,GO,length(hkEisen));
[sigTable_cm,GOcover_cm] = performGOenrichment(gini_cm,GO,length(hkEisen));
[sigTable_klijn,GOcover_klijn] = performGOenrichment(gini_klijn,GO,length(hkEisen));

% close all

[enrichedMat,~,enrGOids,~] = makeSigTableMatrix({sigTable_gtex,sigTable_hpa,sigTable_cm,sigTable_klijn});

[GOids,numMat,covMat] = makeGOcoverageMatrix({GOcover_gtex;GOcover_hpa;GOcover_cm;GOcover_klijn});

[~,ia,ib] = intersect(enrGOids,GOids);
enrGOids = enrGOids(ia); enrichedMat = enrichedMat(ia,:);
covMat = covMat(ib,:);
numMat = numMat(ib,:);

% [~,i] = sort(sum(numMat,2),'descend');
% enrGOids = enrGOids(i);
% covMat = covMat(i,:);
% numMat = numMat(i,:);

enrichFlag(:,1) = ismember(enrGOids,sigTable_gtex.enriched.GOid);
enrichFlag(:,2) = ismember(enrGOids,sigTable_hpa.enriched.GOid);
enrichFlag(:,3) = ismember(enrGOids,sigTable_cm.enriched.GOid);
enrichFlag(:,4) = ismember(enrGOids,sigTable_klijn.enriched.GOid);

O = {'GTEx','HPA','CellMiner','Klijn et al.'};
enriched.enrichedMat = enrichedMat; enriched.GOids = enrGOids; enriched.organism = O;
coverage.GOids = enrGOids; coverage.numMat = numMat; coverage.covMat = covMat;
[rho,pval] = corr(numMat);
coverage.rho = rho; coverage.pval = pval; coverage.enriched = enrichFlag;

figure;
hold on;
c = {'r','g','b','k'};
for j=1:size(numMat,2)
    histogram(numMat(:,j),0:10:300,'displaystyle','stairs','edgecolor',c{j});
end
legend('GTEx','HPA','CellMiner','Klijn et al.');

figure;
hold on;
c = {'r','g','b','k'};
for j=1:size(numMat,2)
    histogram(covMat(:,j),0:0.02:1,'displaystyle','stairs','edgecolor',c{j});
end
legend(O);