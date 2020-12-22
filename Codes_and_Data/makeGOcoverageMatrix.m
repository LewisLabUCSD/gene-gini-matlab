function [GOids,covMat,xMat] = makeGOcoverageMatrix(GOcover)

GOids = GOcover{1}.GOterms;
covMat(:,1) = GOcover{1}.nums(:,1);
xMat(:,1) = GOcover{1}.nums(:,1)./GOcover{1}.nums(:,2);
for i=2:length(GOcover)
    for j=1:length(GOcover{i}.GOterms)
        if sum(GOids==GOcover{i}.GOterms(j))~=0
            covMat(GOids==GOcover{i}.GOterms(j),i) = GOcover{i}.nums(j,1);
            xMat(GOids==GOcover{i}.GOterms(j),i) = GOcover{i}.nums(j,1)/GOcover{i}.nums(j,2);
        else
            covMat(end+1,i) = GOcover{i}.nums(j,1);
            xMat(end+1,i) = GOcover{i}.nums(j,1)/GOcover{i}.nums(j,2);
            GOids(end+1) = GOcover{i}.GOterms(j);
        end
    end
end