function [enriched, coverage,iAnimal,L,Nhk] = multi_organism_analysis

% cd celegans\tissue_GFP\resolveExpressionData
% organism = {'Human','Mouse','Bonobo','Chicken','Gorilla','Chimpanzee','Orangutan','Macaque','Opossum','Platypus'};
T = readtable(strcat(pwd,'\Data\Ortho_1to1_AllSpecies.txt'));
O = T.Properties.VariableNames;
hkGenes = getBrawandHousekeeping(O);
Nhk = cellfun(@length,hkGenes);
for i=1:length(O)
    L(:,i) = eval(['ismember(T.',O{i},',hkGenes{strcmp(O,''',O{i},''')});']);
end

N = sum(L,2);

for i=1:length(O)
    x(i) = sum(N==i);
end
bar(1:1:9,x);

calcJaccardSimilarity(L,O,'matrix',true,0);
unzip(strcat(pwd,'\Data\goa_human.zip'));
humanGO = goannotread_ch(strcat(pwd,'\Data\goa_human.gaf'),'aspect','p');
humanGO = struct2table(humanGO);
GO = geneont('live',true);

% get 1:1 orthologs
load(strcat(pwd,'\Data\orthologs_1_1_AllSpecies_GeneSymbol.mat'),'TgeneInfo');

ichip = ismember(humanGO.DB_Object_Symbol,TgeneInfo.To);
for i=1:length(O)
    iAnimal(:,i) = ismember(humanGO.DB_Object_Symbol,TgeneInfo.To(ismember(TgeneInfo.From,T.Human(L(:,i)))));
    sigTable{i} = hypergeo_GOtable(humanGO,GO,ichip,iAnimal(:,i),0.05,0);
    GOcover{i} = getGOtermCoverage(humanGO,ichip,iAnimal(:,i));
    close
end
iAnimal = [iAnimal ichip];

[enrichedMat,~,enrGOids,~] = makeSigTableMatrix(sigTable);
for i=1:length(O)
    enrichFlag(:,i) = ismember(enrGOids,sigTable{i}.enriched.GOid);
end

[GOids,numMat,covMat] = makeGOcoverageMatrix(GOcover);
[rho,pval] = corr(numMat);

[~,ia,ib] = intersect(enrGOids,GOids);
enrGOids = enrGOids(ia); enrichedMat = enrichedMat(ia,:); enrichFlag = enrichFlag(ia,:);
covMat = covMat(ib,:);
numMat = numMat(ib,:);

clusterTissues (enrichedMat,O,0,'euclidean','complete','tissues');
calcJaccardSimilarity(any(enrichedMat,3),O,'matrix',true);

enriched.enrichedMat = enrichedMat; enriched.GOids = enrGOids; enriched.organism = O;
coverage.GOids = enrGOids; coverage.numMat = numMat; coverage.covMat = covMat;
coverage.rho = rho; coverage.pval = pval; coverage.enrich = enrichFlag;

function hkGenes = getBrawandHousekeeping(organisms)

for i=1:length(organisms)
    load([pwd '\Data\transcriptomics\' organisms{i} '.mat'],'rnaData');
    giniData = calculateHKgenes_gini(rnaData,3,17.5,[]);
    hkGenes{i} = giniData.genes(giniData.gc <= giniData.threshold);
    fprintf('Threshold at %0.2f\n',length(hkGenes{i})*100/length(giniData.genes));
    close all
end
