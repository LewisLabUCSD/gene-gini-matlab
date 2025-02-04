function [enriched,coverage] = compareCancerAllEssentials(cancerData)

% cancerData is either from CellMiner or Klijn et al.
[cell_lines,crisprFiles] = existCRISPRfileCellLine(cancerData.Tissue);
% get the guide RNA depletion data
depl_p = getCRISPRscores(crisprFiles);
% get essential genes
geneConv = getBiomartGeneConversion;
[~,ia,ib] = intersect(depl_p.genes,geneConv.NCBIGene_formerlyEntrezgene_ID);
depl_p.values = depl_p.values(ia,:); depl_p.genes = depl_p.genes(ia);
depl_p.gene_names = geneConv.GeneName(ib);
iess_gene = sum(depl_p.values < 0,2)>=20;
fprintf('%d genes were essential\n',sum(iess_gene));

% get GO terms
unzip(strcat(pwd,'\Data\goa_human.zip'),strcat(pwd,'\Data'));
humanGO = goannotread_ch(strcat(pwd,'\Data\goa_human.gaf'),'aspect','p');
humanGO = struct2table(humanGO);

% get ontology
GO = geneont('live',true);
mask = strcmp(get(GO.terms,'ontology'),'biological process');
bp = GO.terms(mask);
BP = GO(bp);

% prepare hits and chip
ihits = ismember(humanGO.DB_Object_Symbol,depl_p.gene_names(iess_gene));
ichip = ismember(humanGO.DB_Object_Symbol,depl_p.gene_names);

% perform GO enrichment tests
sigTableEss = hypergeo_GOtable(humanGO,BP,ichip,ihits,0.05,0);

% get coverage information
GOcoverEss = getGOtermCoverage(humanGO,ichip,ihits);

% limit samples to the cell-lines for which we have depletion ratios
cancerData.valuebyTissue(:,~ismember(cancerData.Tissue,cell_lines)) = [];
cancerData.Tissue(~ismember(cancerData.Tissue,cell_lines)) = [];

% calculate gini coefficients and HK genes
cancerGini = calculateHKgenes_gini(cancerData,3,17.5,depl_p.genes(iess_gene));
cancerHK_names = cancerGini.gene_name(cancerGini.gc < ...
    prctile(cancerGini.gc,...
    sum(iess_gene)*100/length(cancerGini.genes)));

% prepare hits and chip
ihits = ismember(humanGO.DB_Object_Symbol,cancerHK_names);
ichip = ismember(humanGO.DB_Object_Symbol,cancerGini.gene_name);

% perform GO enrichment tests
sigTableHK = hypergeo_GOtable(humanGO,BP,ichip,ihits,0.05,0);

% get coverage information
GOcoverHK = getGOtermCoverage(humanGO,ichip,ihits);

[enrichedMat,~,enrGOids,~] = makeSigTableMatrix({sigTableEss,sigTableHK});

[GOids,numMat,covMat] = makeGOcoverageMatrix({GOcoverEss;GOcoverHK});

[~,ia,ib] = intersect(enrGOids,GOids);
enrGOids = GOids(ib);
covMat = covMat(ib,:);
numMat = numMat(ib,:);

enrichFlag(:,1) = ismember(enrGOids,sigTableEss.enriched.GOid);
enrichFlag(:,2) = ismember(enrGOids,sigTableHK.enriched.GOid);

O = {'Essential genes','Cancer data'};
enriched.enrichedMat = enrichedMat; enriched.GOids = enrGOids; enriched.organism = O;
coverage.GOids = enrGOids; coverage.numMat = numMat; coverage.covMat = covMat;
[rho,pval] = corr(numMat);
coverage.rho = rho; coverage.pval = pval; coverage.enriched = enrichFlag;
