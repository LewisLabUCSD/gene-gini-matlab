function [Hg,Tg,outpermg] = clusterTissues(exprMat,Tissue,k,distMethod,linkageMethod,elementFlag)

exprMat(sum(exprMat,2)==0,:) = [];
% exprMat = log10(exprMat+1);
if strcmp(elementFlag,'genes')
    % hierarchical clustering
    Y = pdist(exprMat,distMethod);
    % y = squareform(Y);
    Z = linkage(Y,linkageMethod);
    % Y = pdist(Yhist,distMethod);
    leafOrderg = optimalleaforder(Z,Y);
    fprintf('Cophenetic correlation coeffcient using %s linkage and %s distance = %0.4f\n',...
        linkageMethod,distMethod,cophenet(Z,Y));
    % idx = cluster(Z,'maxclust',k,'criterion','distance');
    % cObj = clustergram(Yhist);
    figure;
    [Hg,Tg,outpermg] = dendrogram(Z,k,'Orientation','left','Reorder',leafOrderg);
else
    % hierarchical clustering
    Y = pdist(exprMat',distMethod);
    y = squareform(Y);
    Z = linkage(Y,linkageMethod);
    % Y = pdist(Yhist,distMethod);
    leafOrderg = optimalleaforder(Z,Y);
    fprintf('Cophenetic correlation coeffcient using %s linkage and %s distance = %0.4f\n',...
        linkageMethod,distMethod,cophenet(Z,Y));
    % idx = cluster(Z,'maxclust',k,'criterion','distance');
    % cObj = clustergram(Yhist);
    figure;
    subplot(1,2,1);
    [Hg,Tg,outpermg] = dendrogram(Z,0,'Orientation','left','Reorder',leafOrderg,'labels',strrep(Tissue,'_','\_'));
    xtickangle(45);
end