function [ psfMat ] = getPsfMat2( sz,sf,sig )
% The purpose of this function is to generate the matrix form of a
% gaussian point-spread function, such that element (i,j) of psfMat is the
% proportion of photons emitted at location j which are detected at
% location i.
% Inputs:
%   sz: the dimensions of the FoV in the high-resolution grid
%   sf: the scale factor by which the high-resolution grid is downsampled.
%   Must be an integer value, and sz should be evenly divisible by sf.
%   sig: the standard deviation of the gaussian.
% Outputs:
%   psfMat: the matrix representing the operation of PSF.
sz_sf = sz/sf;
if(sz_sf~=round(sz_sf))
    error('Size must be an even multiple of the scale factor.'); 
elseif(sf~=round(sf))
    error('Scale factor must be an integer.');
end
psf = fspecial('gaussian',2*sz,sig); %generate psf kernel
psf_v = psf(:); %vectorize the PSF
N = sz^2;
psfMat = [];

% Create PSF matrix from vectorized PSF kernel
for y = 1:N
    x = y-1;
    s1 = mod(x,sz);
    s2 = floor(x/sz);
    psf_r = imtranslate(psf,[s2,s1]);
    psf_s = psf_r(sz+1:end,sz+1:end);
    psf_sf = zeros(sz_sf,sz_sf);
    for z = 1:sz_sf
        for w = 1:sz_sf
            psf_sf(z,w) = sum(sum(psf_s((z-1)*sf+1:z*sf,(w-1)*sf+1:w*sf)));
        end
    end
    psfMat = [psfMat,psf_sf(:)];
end

end

