function [locs] = getCellCtrs(seg_im,segThresh)
% getCellCtrs transforms a cell segmentation image into a list of cell 
% centroid locations 
%Inputs:
%   seg_im: cell segmentation image
%   segThresh: cells with area below threshold are filtered out
%Outputs:
%   locs: list of cell centroid locations

uCells = unique(seg_im);
uCells = uCells(2:end);
locs = [];
if(~isempty(uCells))
    for x = 1:size(uCells,1)
        BW = seg_im==uCells(x);
        cc = bwconncomp(BW,26);
        props = regionprops(cc);
        ctrs = props.Centroid;
        pas = zeros(size(props,1),1);
        for y = 1:size(props,1)
            pas(y) = props(y).Area;
        end
        if(max(pas)>segThresh)
            locs = [locs;ctrs(1:2)];
        end
    end
end
end

