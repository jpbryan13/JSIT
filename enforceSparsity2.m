function [Xf] = enforceSparsity2(X,k,t)
% enforceSparsity2 enforces k-sparsity on the rows of X, keeping only the
% highest k values of each row.
%Inputs:
%   X: the original matrix
%   k: k-sparsity to be enforced on rows of X
%   t: hard threshold to be applied to X
%Outputs:
%   Xf: X, with k-sparsity applied to rows.

Xf = zeros(size(X)); %Initialize
X(X<t) = 0; %Hard-threshold X

%Enforce k-sparsity
for x = 1:size(X,1)
    Xvec = X(x,:);
    for ks = 1:k
        [a,ind] = max(abs(Xvec));
        if(a>0)
            Xf(x,ind) = Xvec(ind);
        else
            break;
        end
        Xvec(ind) = 0;
    end
end


end

