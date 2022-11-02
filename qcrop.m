function [ q2,inds ] = qcrop( q,xmin,xmax,ymin,ymax )
% qcrop filters a list of transcript calls by spatial position
%Inputs:
%   q: list of transcript calls. Column 1: spatial coordinate 1, 
% column 2: spatial coordinate 2, column 3: spatial coordinate 3,
% column 4: gene identity by codebook row.
%   xmin: calls with first spatial coordinate < xmin are filtered out
%   xmax: calls with first spatial coordinate > xmax filtered out
%   ymin: calls with second spatial coordinate < ymin filtered out
%   ymax: calls with second spatial coordinate > ymax filtered out
%Outputs:
%   q2: filtered list of transcript calls
%   inds: boolean vector identifying transcript calls not filtered out.


% Identify calls falling within specified spatial bounds
inds_1 = q(:,1)>=ymin;
inds_2 = q(:,1)<=ymax;
inds_3 = q(:,2)>=xmin;
inds_4 = q(:,2)<=xmax;
inds = find(inds_1.*inds_2.*inds_3.*inds_4);

% Filter list
q2 = q(inds,:);
q2(:,1) = q2(:,1)-(ymin-1);
q2(:,2) = q2(:,2)-(xmin-1);

end

