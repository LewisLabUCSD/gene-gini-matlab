function crScores = getCRISPRscores(filenames)

for i=1:length(filenames)
    depl_p = getCellCRISPRdata(filenames{i});
    if i > 1
        g = length(crScores.genes);
        [~,ia,ib] = intersect(depl_p.genes,crScores.genes);
%         [ib,ix] = sort(ib,'ascend'); ia = ia(ix);
        [~,ic] = setdiff(depl_p.genes,crScores.genes);
        crScores.genes(g+1:g+length(ic)) = crScores.genes(ic);
        % pre-assign nan to the new column
        crScores.values(:,i) = nan(g,1);
        % assign values to genes which are already present in the list
        crScores.values(ib,i) = depl_p.values(ia);
        % pre-assign nan to the new rows, if any
        if ~isempty(ic)
            crScores.values(end+1:end+length(ic),:) = nan(length(ic),i);
            % assign values to genes which are not present in the list
            crScores.values(g+1:g+length(ic),i) = depl_p.values(ic);
        end
    else
        crScores = depl_p;
    end
end

function genes_ratios = getCellCRISPRdata(filename)

% Load depletion (CRISPR) and ID data
[gGeck,gEntr] = textread(strcat(pwd,'\Data\gecko_id_to_entrez.txt'),'%s %s');
[~,~,raw] = xlsread(filename);
gDepl = raw(:,1);
for i = 1:length(gDepl)
    if isnumeric(gDepl{i})
        gDepl{i} = num2str(gDepl{i});
    end
end

dRatio = raw(:,2);
% Only get those genes which can be translated from gecko to entrez
[~, entrInds] = ismember(gDepl,gGeck);
genes_ratios.genes = gEntr(entrInds(entrInds~=0));
genes_ratios.values = cell2mat(dRatio(find(entrInds)));