function [ xp ] = proxSGL( alpha,t,x)
%This function calculates the proximal function for the Sparse Group LASSO
%regularization function.
%Inputs:
%   alpha - gives relative importance of mixed l1/l2 row norm versus overall
%   l1 norm
%   x - argument of the proximal function
%   t - threshold for soft-thresholding function. Related to the desired level of sparsity
%Outputs:
%   xp - x, with proximal function applied.

xt = softThresh(x,alpha.*t); %First do soft-thresholding

%Initialize variables
nRows = size(xt,1);
xp = zeros(size(xt));
tmat = ones(size(xt)).*t;

%Then apply the proximal function of the mixed norm
for r = 1:nRows 
    if(norm(xt(r,:),2)>=((1-alpha)*t))
        xp(r,:) = (1-(((1-alpha).*tmat(r,:))./norm(xt(r,:),2))).*xt(r,:);
    end
end
end

