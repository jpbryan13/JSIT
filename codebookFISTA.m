function [ X ] = codebookFISTA( A,C,Y,lambda,kmax,proxOp,alpha,beta,eK,eM )
%codebookFISTA performs the FISTA optimization algorithm on a problem
%of the form ||AXC-Y||_{F}^{2}. This function performs the algorithm
%directly, without taking advantage of any tricks or structure of the
%matrices.
%Inputs:
%   A: PSF matrix
%   C: Codebook
%   Y: iST Measurement Matrix
%   lambda: regularization parameter 1
%   kmax: number of iterations
%   proxOp: string giving regularization mode for which to compute proximal
%   operator. Options are "L1", "L0", and "SGL" (Sparse Group Lasso).
%   alpha: regularization parameter 2, used in SGL proximal operation
%   beta: regularization parameter 3, used in L0 proximal operation
%   eK: largest eigenvalue of K=C*C'
%   eM: largest eigenvalue of M=A*A'
%Outputs:
%   X: Estimate of spatial mRNA transcript distribution

% Compute Lipschitz constant
Lf = eK*eM;

% Initialize variables
s1 = size(A,2);
s2 = size(C,1);
sy = sqrt(size(A,1));
z = zeros(s1,s2);
X = z;
t = 1;
Zs = zeros(s1,s2,kmax);
Ydiff = [];
Ymin = [];

% Run FISTA algorithm
for k = 1:kmax
    gfz = A'*(A*z*C-Y)*C';
    Xp = X;
    if(strcmp(proxOp,'L1'))
        X = softThresh(X-(gfz/Lf),lambda/Lf);
    elseif(strcmp(proxOp,'L0'))
        X = proxL0(lambda/Lf,beta,X-(gfz/Lf));
    elseif(strcmp(proxOp,'SGL'))
        X = proxSGL(alpha,lambda/Lf,X-(gfz/Lf));
    end
    tp = t;
    t = 0.5.*(1+sqrt(1+4*(t^2)));
    z = X+((t/tp).*(X-Xp));
    Zs(:,:,k) = z;
end

end

