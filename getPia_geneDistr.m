function [distr_gene] = getPia_geneDistr(gg,pDist,ginds,nBins)
% getPia_geneDistr gives the distribution of gene expression by their
% distance from a set of spatial positions
%Inputs:
%   gg: Cells-by-genes count table
%   pDist: Distance of each cell from a set of spatial positions (in our
%   pipeline, from the pia).
%   ginds: Gene indices for which to compute distribution
%   nBins: Number of bins in which to divide distances to the set of
%   spatial positions
%Outputs:
%   distr_gene: Binned distribution of gene expression by distance from a
%   set of spatial positions.

% Separate pDist values into bins
pDisc = discretize(pDist,nBins);
nGenes = length(ginds);
distr_gene = zeros(nBins,nGenes);

% Get expression levels of genes within each distance bin.
for y = 1:nGenes
for x = 1:nBins
    distr_gene(x,y) = mean(gg(ginds(y),pDisc==x));
end
end

end

