function [ggn,cDist] = normalizeByRadius(gg,cc,nBins,doLog)
% normalizeByRadius normalizes a count table, adjusting the distribution 
% of each gene based on the distance of each cell to the center of the FOV
%
%Inputs:
%   gg: cells-by-genes count table
%   cc: table of cell locations, within FOV
%   nBins: number of bins into which to divide the FOV
%   doLog: boolean, normalize logarithm of count table if true.
%Outputs:
%   ggn: normalized count table
%   cDist: distance of each cell from center of FOV.

%Initialize
nGenes = size(gg,1);
cDist = zeros(size(cc,1),1);
for x = 1:size(cc,1)
    cDist(x) = norm(cc(x,:)-[1024,1024]);
end

% Split cells into bins by radial distance from center.
[yDist,yEdges] = discretize(cDist,nBins);

% Compute percentage of all cells in each radial distance bin
distPcts = zeros(nBins,1);
for x = 1:10
    distPcts(x) = length(find(yDist==x))./length(cDist);
end

if(doLog) %Take logarithm of count table
    gg = log(gg+1);
end


% Normalize by radial distance bin
ggn = zeros(size(gg));
for x = 1:nGenes
    binCts = zeros(nBins,1);
    radNorm = zeros(nBins,1);
    for y = 1:nBins
        binCts(y) = sum(gg(x,yDist==y))./distPcts(y);
    end
    for y = 1:nBins
        if(binCts(y)~=0)    
            radNorm(y) = binCts(2)./binCts(y);
            ggn(x,yDist==y) = gg(x,yDist==y).*radNorm(y);
        end
    end
end

% Take exponential if logarithm was previously taken
if(doLog)
    ggn = exp(ggn)-1;
end
end

