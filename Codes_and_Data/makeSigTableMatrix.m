function [enrichedMat,depletedMat,enrGOids,depGOids] = makeSigTableMatrix(sigTable)

enrGOids = sigTable{1}.enriched.GOid;
enrichedMat(:,1) = sigTable{1}.enriched.hitNums;
for i=2:length(sigTable)
    for j=1:length(sigTable{i}.enriched.GOid)
        if sum(enrGOids==sigTable{i}.enriched.GOid(j))~=0
            enrichedMat(enrGOids==sigTable{i}.enriched.GOid(j),i) = sigTable{i}.enriched.hitNums(j);
        else
            enrichedMat(end+1,i) = sigTable{i}.enriched.hitNums(j);
            enrGOids(end+1) = sigTable{i}.enriched.GOid(j);
        end
    end
end

depGOids = sigTable{1}.depleted.GOid;
depletedMat(:,1) = sigTable{1}.depleted.hitNums;
for i=2:length(sigTable)
    for j=1:length(sigTable{i}.depleted.GOid)
        if sum(depGOids==sigTable{i}.depleted.GOid(j))~=0
            depletedMat(depGOids==sigTable{i}.depleted.GOid(j),i) = sigTable{i}.depleted.hitNums(j);
        else
            depletedMat(end+1,i) = sigTable{i}.depleted.hitNums(j);
            depGOids(end+1) = sigTable{i}.depleted.GOid(j);
        end
    end
end