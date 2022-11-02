function [ X ] = processXstack( Xstack,s1,s2,h1,h2 )
% ProcessXstack takes a stack of recovered X matrices, and turns it into a
% single X matrix.
%Inputs:
%   Xstack: stack of X matrices
%   s1, s2: the "block" dimensions of the original image
%   h1, h2: the dimensions of the original image patches
%Outputs:
%   X: Assembled X matrix

%Initialize
nGenes = size(Xstack,3);
nXs = size(Xstack,1);
Xstack_r = zeros(nXs,h1,h2,nGenes);

%Reshape all X stack members into 3D tensors, so they can be put together
for a = 1:nXs
    Xstack_r(a,:,:,:) = reshape(Xstack(a,:,:),[h1,h2,nGenes]);
end

%Initialization step 2
Xstack_r2 = [];
Xstack_r2x = [];

%Assemble all X stack members into a single 3D tensor
for a = 1:nXs
    Xstack_r2x = cat(3,Xstack_r2x,Xstack_r(a,:,:,:));
    if(mod(a,s2)==0)
        Xstack_r2 = cat(2,Xstack_r2,Xstack_r2x);
        Xstack_r2x = [];
    end
end

%Reshape into a single 2D matrix
X = [];
for p = 1:nGenes
    Xstack_r3x = reshape(Xstack_r2(1,:,:,p),[s1*s2*h1*h2,1]);
    X = [X,Xstack_r3x];
end
    
end

