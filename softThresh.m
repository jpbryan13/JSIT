function [ yt ] = softThresh(y,t)
%softThresh performs the soft-thresholding operation.
%Inputs:
%   y: the original signal
%   t: the threshold
%Outputs:
%   yt: the thresholded signal
tmat = ones(size(y)).*t;
yt = max(abs(y) - tmat,0).*sign(y);
end

