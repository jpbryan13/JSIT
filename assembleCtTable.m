function [ctTable] = assembleCtTable(q,nGenes)
% assembleCtTable transforms a list of transcript calls, paired with
% segmented cells, into a genes-by-cells count table
%Inputs:
%   q: list of transcript calls. Column 1: spatial coordinate 1, 
% column 2: spatial coordinate 2, column 3: spatial coordinate 3,
% column 4: gene identity by codebook row, column 5: segmented cell
% identity
%   nGenes: Number of genes in the codebook.
%Outputs:
%   ctTable: genes-by-cells gene expression count table

segInds = unique(q(:,end));
nCells = length(segInds);
ctTable = zeros(nGenes,nCells);
for x = 1:nCells
    ctTable(:,x) = getGeneCounts(q(q(:,end)==segInds(x),[1,2,4]),nGenes);
end

end

