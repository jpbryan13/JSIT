function [segObjects] = getSegObjects_3D(I,sizeThresh)
% getSegObjects_3D turns an image of segmented cells into a set of
% segmented cell objects, in the "cell" datatype.
%Inputs:
%   I: Image depicting cell segmentation. Elements with the same value
%   belong to the same cell.
%   sizeThresh: Cells smaller than sizeThresh are filtered out.
%Outputs:
%   segObjects: "cell" datatype, each element is a list of indices giving
%   the spatial coordinates of a segmented cell.

A = unique(I); %Get all values corresponding to a segmented cell
segObjects = cell(0,1); %Initialize
for x = 1:size(A,1)
    segPix = find(I==A(x)); %get indices corresponding to cell x
    if(size(segPix,1)>=sizeThresh)
        segObjects{end+1} = segPix;
    end
end
end

