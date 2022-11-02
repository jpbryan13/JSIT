function [ xp ] = proxL0(alpha,beta,X)
%proxL0 performs the a relaxed version of the hard threshold, the proximal 
% function of the L0 norm. 
%Inputs:
%   alpha: the position of the threshold transition
%   beta: controls the steepness of the threshold transition
%   X: the original signal
%Outputs:
%   xp: the signal, with the proximal function applied.

rx = X;
rx(X<0)=0;
xp = rx./(1+exp(-beta.*(abs(X) - alpha)));

end

