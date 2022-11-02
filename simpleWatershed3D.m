function [Lcl] = simpleWatershed3D(I)
% simpleWatershed3D performs the watershed segmentation algorithm on an
% image
%Inputs:
%   I: Image or stack of images to be segmented
%Outputs:
%   Lcl: Segmentation image
%Written by Loïc Binan

%Preprocessing steps
myimage=imbinarize(I,'adaptive','sensitivity',0.25);
myimage=imdilate(myimage,strel('disk',4));
myimage=bwareaopen(myimage,7000);
BW=imdilate(myimage, strel('disk',8));
D = zeros(size(BW));
for x = 1:size(BW,3)
D(:,:,x) = bwdist(~BW(:,:,x));
end
D = -D;
D(~BW) = Inf;
DsansMinimas = imhmin(D,4);

%Perform watershed algorithm
L = watershed(DsansMinimas);

%Post-processing
L(~BW) = 0;
LCleaned = zeros(size(L));
for x = 1:size(L,3)
    LCleaned(:,:,x) = imclearborder(L(:,:,x));
end
l=bwareaopen(LCleaned,12000);
Lcl = LCleaned;
Lcl(~l) = 0;
end

