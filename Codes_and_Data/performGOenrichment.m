function [sigTable,GOcover] = performGOenrichment(giniData,GO,n)


% % load GO table
unzip(strcat(pwd,'\Data\goa_human.zip'));
humanGO = goannotread_ch(strcat(pwd,'\Data\goa_human.gaf'),'aspect','p');
humanGO = struct2table(humanGO);

if isfield(giniData,'gc') || nargin > 3
    % n: number of gini genes to be chosen
    geneList = giniData.gene_name(giniData.gc < prctile(giniData.gc,n*100/length(giniData.genes)));
    % % map genes to human GO terms
    ihits = ismember(humanGO.DB_Object_Symbol,geneList);
    ichip = ismember(humanGO.DB_Object_Symbol,giniData.gene_name);
else
    ihits = ismember(humanGO.DB_Object_Symbol,giniData.hitGenes);
    ichip = ismember(humanGO.DB_Object_Symbol,giniData.chipGenes);
end

% perform GO enrichment tests
sigTable = hypergeo_GOtable(humanGO,GO,ichip,ihits,0.05,0);

% get coverage information
[GOcover] = getGOtermCoverage(humanGO,ichip,ihits);
