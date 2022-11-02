function [q,dIm_c] = dIm2q_ex(dIm,iIm,distIm,cls,codebook)
% dIm2q_ex transforms a set of maps giving transcript call information into
% a list of transcript calls.
%Inputs:
%   dIm: spatial map of transcript calls in FOV, with element values indicating
%   gene identity
%   iIm: spatial map of averaged measured image intensity at each location
%   in FOV
%   distIm: spatial map of magnitude of X matrix at each location in FOV
%   cls: minimum cluster size: calls in dIm smaller than cls are filtered
%   out.
%   codebook: codebook matrix
%Outputs:
%   q: list of transcript calls
%   dIm_c: dIm map, filtered by cluster size

% Initialize
q = [];
dIm_c = zeros(size(dIm));

% Filter dIm by cluster size
for x = 1:size(codebook,1)
    A = dIm==x; %Pick out calls of gene x
    B = bwconncomp(A,4); % Separate map by connected components (aka clusters)
    for y = 1:B.NumObjects
        if(numel(B.PixelIdxList{y})>=cls) %filter connected components by size
            dIm_c(B.PixelIdxList{y})=dIm(B.PixelIdxList{y});
        end
    end
end

%Transform maps into list of transcript calls
for x = 1:size(codebook,1)
    A = dIm_c==x;
    B = bwconncomp(A,4);
    codeInd = x.*ones(B.NumObjects,1);
    c = regionprops(B,'Centroid');
    mIn = zeros(length(c),1);
    mDist = zeros(length(c),1);
    cArea = zeros(length(c),1);
    for y = 1:length(c)
        mIn(y) = mean(iIm(B.PixelIdxList{y})); %Average Image intensity over cluster
        mDist(y) = mean(distIm(B.PixelIdxList{y})); %Average distIm value over cluster
        cArea(y) = length(B.PixelIdxList{y}); %Cluster size
    end
    qarr = [cat(1,c.Centroid),codeInd,cArea,mIn,mDist];
    if(isempty(qarr))
        qarr = double.empty(0,6);
    end
    q = [q;qarr];
end

end

