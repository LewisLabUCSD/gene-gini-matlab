function [cellLines,crisprFile] = existCRISPRfileCellLine(cellLines)

f = false(length(cellLines),1);
c = 0;
for i=1:length(cellLines)
    if exist(strcat(pwd,'\Data\geneEssentiality\CRISPR_Avana\',cellLines{i},'_Avana.xlsx'),'file')~=0
        c = c + 1;
        f(i) = true;
        crisprFile{c,1} = strcat(pwd,'\Data\geneEssentiality\CRISPR_Avana\',cellLines{i},'_Avana.xlsx');
    elseif exist(strcat(pwd,'\Data\geneEssentiality\CRISPR_Avana\',cellLines{i},'_Gecko.xlsx'),'file')~=0
        c = c + 1;
        f(i) = true;
        crisprFile{c,1} = strcat(pwd,'\Data\geneEssentiality\CRISPR_Avana\',cellLines{i},'_Gecko.xlsx');
    else
        f(i) = false;
    end
end
cellLines = cellLines(f)';
