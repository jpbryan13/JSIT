function [ Ystack ] = makeYstack( Y,s1,s2 )
% makeYstack breaks up the measured data matrix Y into patches 
%Inputs:
%   Y: the measured data matrix
%   s1, s2: size, along axis 1 and 2, of patches.
%Outputs:
%   Ystack: stack of patches

%Initialize
a1 = size(Y,1);
a2 = size(Y,2);
d1 = round(a1/s1);
d2 = round(a2/s2);
Ystack = [];

%Split up Y into patches
for x = 1:d1
    for y = 1:d2
        Yxy = Y(1+(x-1)*s1:x*s1,1+(y-1)*s2:y*s2,:);
        Ystack = cat(4,Ystack,Yxy);
    end
end

Ystack = reshape(Ystack,[s1*s2,size(Y,3),size(Ystack,4)]);

end

