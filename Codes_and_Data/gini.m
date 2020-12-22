function [g,l,a] = gini(P,v,figFlag)

% GINI computes the Gini coefficient and the Lorentz curve.

% P: population, the distribution of total transcripts in different samples
% v: value, the distribution of a given transcript in different samples
% figFlag: false (default), if using within a loop and don't need a plot

if nargin < 3
    figFlag = false; % will only be true if needed for one gene
end

P = [0;P(:)]; v = [0;v(:)]; 

% filter out NaNs need to do this because of preprocessing the transcriptomics
isok = all(~isnan([P,v]'))';
P = P(isok); v = v(isok);


% process input
val_pop_tot = v .* P;
[~,ix] = sort(v); % sort the transcript levels
P = P(ix); val_pop_tot    = val_pop_tot(ix);
% cumulative sums
P = cumsum(P);  val_pop_tot    = cumsum(val_pop_tot);
% normalize
rel_pop = P/P(end); rel_val = val_pop_tot/val_pop_tot(end);

% Gini coefficient calculation
g = 1 - sum((rel_val(1:end-1) + rel_val(2:end)) .* diff(rel_pop));

% Lorentz curve
l = [rel_pop, rel_val];
a = [P, val_pop_tot];

if figFlag
    % draw the lorenz curve
    hold on
    area(rel_pop,rel_val,'FaceColor',[0.5,0.5,1.0]); % the Lorentz curve
    plot([0,1],[0,1],'--k'); % X = Y line for visual comparison
    set(gca,'dataaspectratiomode','manual'); 
    set(gca,'XTick',get(gca,'YTick'))   
    set(gca,'Layer','top');            
    title(['\bfGini coefficient = ',num2str(g)]);
    xlabel('share of population');
    ylabel('share of value');
    hold off
end
