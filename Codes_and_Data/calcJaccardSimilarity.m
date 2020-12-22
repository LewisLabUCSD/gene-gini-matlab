function [J,outperm,T] = calcJaccardSimilarity(A,tisNames,typeFlag,figFlag,k)

nTis = length(tisNames)
J = zeros(nTis,nTis);
for i=1:nTis
    for j=i:nTis
        if strcmp(typeFlag,'matrix')
            J(j,i) = sum(A(:,i) & A(:,j))/sum(A(:,i) | A(:,j));
            J(i,j) = sum(A(:,i) & A(:,j))/sum(A(:,i) | A(:,j));
        elseif strcmp(typeFlag,'struct-rxns')
            J(j,i) = length(intersect(A{i}.rxns,A{j}.rxns))/length(union(A{i}.rxns,A{j}.rxns));
            J(i,j) = length(intersect(A{i}.rxns,A{j}.rxns))/length(union(A{i}.rxns,A{j}.rxns));
        elseif strcmp(typeFlag,'struct-genes')
            J(j,i) = length(intersect(A{i}.genes,A{j}.genes))/length(union(A{i}.genes,A{j}.genes));
            J(i,j) = length(intersect(A{i}.genes,A{j}.genes))/length(union(A{i}.genes,A{j}.genes));
        end
    end
end

Z = linkage(J,'complete');
Y = pdist(J,'euclidean');
leafOrderg = optimalleaforder(Z,Y);
fprintf('Cophenetic correlation coeffcient using %s linkage and %s distance = %0.4f\n',...
    'complete','euclidean',cophenet(Z,Y));
if figFlag
    figure;
    subplot(1,2,1);
    if iscell(tisNames)
        tisNames = strrep(tisNames,'_','-');
    end
    if length(tisNames)<=100
        [~,~,outperm] = dendrogram(Z,0,'Orientation','left','labels',tisNames,'Reorder',leafOrderg);
    else
        [~,T,outperm] = dendrogram(Z,k,'Orientation','left','Reorder',leafOrderg);
    end
    if strcmp(typeFlag,'matrix') && length(tisNames) <= 100
        title('Clustered based on metabolic similarity');
        hold on
        A = A(:,outperm);
        barh(-sum(A,1)/max(sum(A,1)),'facecolor',[0.2 0.4 0.8]);
        xlim([-1.2 1.2]);
        set(gca,'yaxislocation','origin');
        yticklabels(strrep(tisNames(outperm),'_',':'));
        T = [300 600 900 1200 1500];%[1000 2000 3000 4000 5000];
        t = T./max(sum(A,1));
        xticks([-sort(t,'descend') 0:0.4:1.2]);
        xticklabels([sort(T,'descend') 0:0.4:1.2]);
        xlabel('linkage distance    number of reactions');
        set(gca,'layer','top');
        hold off
    elseif strcmp(typeFlag,'struct-rxns')
        title('Clustered based on metabolic similarity of reactions');
    elseif strcmp(typeFlag,'struct-genes')
        title('Clustered based on metabolic similarity of genes');
    end
    
    subplot(1,2,2);
%     figure;
    if length(tisNames) <= 100
        outpermg = outperm(sort(1:1:nTis,'descend'));
        imagesc(tril(J(outpermg,outpermg)));
        ax = gca; ax.TickDir = 'out'; ax.DataAspectRatioMode = 'manual';
        xticks(1:1:nTis); yticks(1:1:nTis);
        xticklabels(strrep(tisNames(outpermg),'_',':'));
        yticklabels(strrep(tisNames(outpermg),'_',':'));
        xtickangle(45);
    else
        [~,ix] = sort(T,'descend');
        imagesc(tril(J(ix,ix)));
        ax = gca; ax.TickDir = 'out'; ax.DataAspectRatioMode = 'manual';
    end
    cbi = colorbar;
    cbi.TickDirection = 'out';
    colormap jet
end
