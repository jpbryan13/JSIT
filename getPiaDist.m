function [pDist] = getPiaDist(subs,cc)
% getPiaDist computes the distance of a list of coordinates from a set of 
% spatial positions (used in our pipeline to compute the distance of cells
% from the pia).
%Inputs:
%   subs: set of spatial positions from which to compute distance
%   cc: list of coordinates, compute distance from each coordinate to subs
%Outputs:
%   pDist: list of distances of coordinates in cc from the set of positions
%   in subs

pDist = zeros(size(cc,1),1);
for x = 1:size(cc,1)
    [~,pDist(x)] = knnsearch(subs,cc(x,:));
end
end
