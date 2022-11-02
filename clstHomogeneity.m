function [h] = clstHomogeneity(lbl,clst)
% clstHomogeneity computes the homogeneity between a set of class labels
% and an unsupervised clustering
%Inputs:
%   lbl: class labels to be compared against
%   clst: cluster labels from unsupervised clustering
%Outputs:
%   h: cluster homogeneity

%Initialize
if(length(lbl)~=length(clst))
    error('Labeled sets must have same length');
end
if(min(lbl)<=0)
    lbl = lbl-min(lbl)+1;
end
if(min(clst)<=0)
    clst = clst-min(clst)+1;
end
nl = max(lbl);
nc = max(clst);
cl_ent = 0;
l_ent = 0;
n = length(lbl);

%Compute homogeneity
for x = 1:nl
    nx = sum(lbl==x); 
    if(nx~=0)
    for y = 1:nc
        ncl = sum((lbl==x).*(clst==y));
        ny = sum(clst==y);
        if(ncl~=0)
        cl_ent = cl_ent+((ncl./n).*log(ncl./ny));
        end
    end
    l_ent = l_ent+((nx./n).*log(nx./n));
    end
end
h = 1-(cl_ent./l_ent);
end
