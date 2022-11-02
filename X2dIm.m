function [ dIm,iIm ] = X2dIm( X,s1,s2)
% X2dIm transforms the X matrix into a map of mRNA transcript calls, and a
% map of the magnitude of the elements of the X matrix
%Inputs:
%   X: Input X matrix
%   s1, s2: dimensions of the full FOV
%Outputs:
%   dIm: map of mRNA transcript calls. Value of each element gives gene
%   identity of detected transcript (if element is equal to zero, no
%   transcript is present).
%   iIm: map of X matrix element magnitudes.

%Initialize
nGenes = size(X,2);
dIm = zeros(s1,s2);
iIm = zeros(s1,s2);

%Create maps from X matrix
for a = 1:nGenes
    geneA = reshape(X(:,a),[s1,s2]);
    dIm(geneA~=0) = a;
    iIm(geneA~=0) = geneA(geneA~=0);
end

end

