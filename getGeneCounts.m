function [gCts] = getGeneCounts(q,nGenes)
% getGeneCounts transforms a list of transcript calls into a vector of gene
% abundances
%Inputs:
%   q: list of transcript calls. Column 1: spatial coordinate 1, 
% column 2: spatial coordinate 2, column 3: gene identity by codebook row
%   nGenes: Number of genes in the codebook.
%Outputs:
%   gCts: vector giving gene expression abundances

gCts = zeros(nGenes,1);
for x = 1:nGenes
    gCts(x) = size(find(q(:,3)==x),1);
end
end

