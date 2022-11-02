function [qSorted,inds_s] = sortMolecules_3D(q,I,segThresh)
% sortMolecules_3D matches transcript calls to segmented cells
%Inputs:
%   q: list of transcript calls. Column 1: spatial coordinate 1, 
% column 2: spatial coordinate 2, column 3: spatial coordinate 3,
% column 4: gene identity by codebook row.
%   I: Image depicting cell segmentation. Elements with the same value
%   belong to the same cell.
%   segThresh: Cells smaller than segThresh are filtered out.
%Outputs:
%   qSorted: list of transcript calls, with cell assignment
%   inds_s: gives indices of the segmented cells which contain nonzero
%   number of transcripts.

% Create segmented cell objects based on image I
segObjects = getSegObjects_3D(I,segThresh);


% Sort transcript calls by segmented cell
qr = q;
if(size(qr,2)==4)
    qr(:,1:3) = round(qr(:,1:3));
    qSorted = zeros(size(qr,1),5);
    qSorted(:,1:4) = qr(:,1:4);
elseif(size(qr,2)==3)
    qr(:,1:2) = round(qr(:,1:2));
    qSorted = zeros(size(qr,1),4);
    qSorted(:,1:3) = qr(:,1:3);
end
inds_s = [];
for x = 2:size(segObjects,2) % Object 1 will be background, so start with 2
    if(size(qr,2)==4)
        [r,c,z] = ind2sub(size(I),segObjects{x});
        sOsubs = [c,r,z];
    elseif(size(qr,2)==3)
        [r,c] = ind2sub(size(I),segObjects{x});
        sOsubs = [c,r];
    end
    inds = ismember(qr(:,1:end-1),sOsubs,'rows');
    if(~isempty(inds))
        qSorted(inds,end) = x;
        inds_s = [inds_s,x];
    end
end

end

