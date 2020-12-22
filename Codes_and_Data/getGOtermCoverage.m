function [GOcover] = getGOtermCoverage(baseGO,ichip,ihit)

chipGO = createGenesetGO(baseGO,ichip);
hitGO = createGenesetGO(baseGO,ihit);
ugos = unique(chipGO.GOid);

chipnums = zeros(length(ugos),1); hitnums = zeros(length(ugos),1);
for i=1:length(ugos)
    if ~isempty(find(chipGO.GOid==ugos(i)))
        chipnums(i,1) = length(find(chipGO.GOid==ugos(i)));
    end
    if ~isempty(find(hitGO.GOid==ugos(i)))
        hitnums(i,1) = length(find(hitGO.GOid==ugos(i)));
    end 
end
[hitnums,ix] = sort(hitnums,'descend');
chipnums = chipnums(ix); ugos = ugos(ix);
nums = [hitnums chipnums];
GOcover.GOterms = ugos; GOcover.nums = nums;

function genesGO = createGenesetGO(baseGO,indices)

genesGO = baseGO(indices,:);