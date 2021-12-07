function psf = wfmpsf(lambdaEx, lambdaEm, numAper, magObj, rindexObj, ...
    ccdSize, dz, xysize, nslices, varargin)
% Generate the point-spread function for a widefield microscope using a 
% scalar diffraction-limited model (Stokseth refer to [1] and [3] below)
% Input parameters:
% ----------------
% lambdaEx: Excitation Wavelength (in nm)
% lambdaEm: Emission wavelegnth (in nm)
% numAper: Numerical aperture of the objective
% magObj: Objective total magnification
% rindexObj: Refractive index of the objective immersion medium
% ccdSize: Pixel dimension of the CCD (in the plane of the camera)
% rindex_sp: Refractive index of the specimen medium
% xysize: Size of the desired image (specimen view size/pixel dimension)
% nslices: Number of slices desired (Depth view/Z axis sampling)
% depth: depth of the specimen under the cover-slip in nm
% dxy: CCD pixel size (in nm)
% dz: Optical axis Z sampling or defocusing (in nm)
% nor: Normalization on the PSF (default: no normalization)
%      0: l-infinity normalization
%      1: l-1 normalization
%
% References:
% ----------
% [1] P. A. Stokseth (1969). `Properties of a defocused optical system?. 
% J. Opt. Soc. Am. A 59:1314?1321. 
% [2] P. Pankajakshan, et al. (2010). `Point-spread function model for 
% fluorescence macroscopy imaging?. In Proc. of Asilomar Conference on 
% Signals, Systems and Computers.
% [3] P. Pankajakshan (2009). Blind Deconvolution for Confocal Laser 
% Scanning Microscopy. Ph.D. thesis, Universit? de Nice Sophia-Antipolis.
%
% Author: Praveen Pankajakshan
% Email: praveen.pankaj@gmail.com
% Last modified: June 24, 2011 17:21
% Copyright (C) 2011 Praveen Pankajakshan
% All rights reserved.
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, 
% this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, 
% this list of conditions and the following disclaimer in the documentation 
% and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY <COPYRIGHT HOLDER> ``AS IS'' AND ANY EXPRESS OR 
% IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
% EVENT SHALL <COPYRIGHT HOLDER> OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
% INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
% BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
% ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% The views and conclusions contained in the software and documentation are those 
% of the authors and should not be interpreted as representing official policies, 
% either expressed or implied, of Praveen Pankajakshan.
%---------------
%%Example usage:
%%Note that for large z sections, the algorithms speed is dependent on the
%%2D fFTs
%%Without spherical aberration, and normalized to l-infinity:
%%Microscope objective: 63X/1.4, oil immersion
%%CCD camera pixel size: 6.45 microns
%%Fluorophore is GFP and excited in blue
% psf = wfmpsf(488, 520, 1.4, 63, 1.518, 6450, 300, 128, 64, 1.518, 0, 0);
% 
%%With spherical aberration:
%%Microscope objective: 63X/1.4, oil immersion
%%CCD camera pixel size: 6.45 microns
%%Specimen immersed has a buffer index of 1.45 
%%and imaged at 15 microns in depth
%%Fluorophore is GFP and excited in blue
% psf = wfmpsf(488, 520, 1.2, 40, 1.33, 6450, 300, 256, 128, 1.45, 15000, 0);
%
%%To plot the PSFs:
% load cmap_fire;
%%Remove the errors due to the Fourier transform 
% psf(psf<5.7e-06) = 5.7e-06;
% figure, imagesc(maxintensityproj(log(psfl), 1)'); 
% colormap(cmap_fire); colorbar;
%---------------
nor = 0; % default no normalization on the PSF
% calculation of pixel size in the plane of the specimen
dxy = ccdSize/magObj;
if nargin>9
    switch nargin
        case 10
            display(['Refractive index of specimen or the imaging', ...
                ' depth is missing. Assuming no spherical aberration.']);
            rindexSp = rindexObj;
            depth = 0;
        case 11
            rindexSp = varargin{1};
            depth = varargin{2};
        otherwise
            rindexSp = varargin{1};
            depth = varargin{2};
            nor = varargin{3};
    end
       
    psfEx = aberratedpsf(lambdaEx, numAper, rindexObj, dxy, dz, ...
        xysize, nslices, rindexSp, depth, 'illumination');
    psfEm = aberratedpsf(lambdaEm, numAper, rindexObj, dxy, dz, ...
        xysize, nslices, rindexSp, depth, 'emission');    
    
end
psf = psfEx.*psfEm;
psf(isnan(psf)) = eps;
if nor == 0
    psf = psf./sum(psf(:));
elseif nor == 1
    psf = psf./max(psf(:));
end
function psf = aberratedpsf(lambda, numAper, rindexObj, dxy, dz, ...
    xysize, nslices, rindexSp, depth, varargin)
% Generate the unaberrated and the aberrated PSF using P.A.Stokseth model
% lambda: Excitation/Emission Wavelength (in nm)
% numAper: Numerical aperture of the microscope
% rindexObj: Refractive index of the objective immersion medium
% dxy: Radial pixel dimension (in nm)
% dz: Axial Optical axis Z sampling or slicing width (in nm)
% xysize: Size of the desired image (specimen view size/pixel dimension)
% nslices: Number of slices desired (Depth view/Z axis sampling)
% rindexSp: Refractive index of the specimen immersion medium
% depth: Depth under the coverslip
% initializing
pupil = zeros(xysize, xysize);
psf = pupil;
N = xysize/2;
n = nslices/2;
% Pupil space pixel dimensions dkx, dky
dkxy = (2*pi)/(xysize*dxy);
% Calculate the defocus
defocus = (-n:(n-1)).*dz;
%Calculated the wavelength of light inside the objective lens and specimen
lambdaObj = lambda/rindexObj;
lambdaSp = lambda/rindexSp;
% Calculate the wave vectors in vaccuum, objective and specimens
k0 = 2*pi/lambda;
kObj = 2*pi/lambdaObj;
kSp = 2*pi/lambdaSp;
% Radius of the pupil function disk
%kMax = 4*xysize*((dxy*numAper)/lambda)^2;
kMax = (2*pi*numAper)/(lambda*dkxy);
% Generate the pupil function amplitude
kxcord = (1:xysize)'-N-1;%Setting N+1 as the center of the pupil
kycord = kxcord;
[kx, ky] = meshgrid(kxcord, kycord);
k = sqrt(kx.^2+ky.^2);
pupil = (k< kMax);
% Calculate the sine of the semi-aperture angle in the objective lens
sinthetaObj = (k.*(dkxy))/kObj;
sinthetaObj(sinthetaObj>1) = 1;
% Calculate the cosine of the semi-aperture angle in the objective lens
costhetaObj = eps+sqrt(1-(sinthetaObj.^2));
% Calculate the sine of the semi-aperture angle in the specimen
sinthetaSp = (k.*(dkxy))/kSp;
sinthetaSp(sinthetaSp>1) = 1;
% Calculate the cosine of the semi-aperture angle in the specimen
costhetaSp = eps+sqrt(1-(sinthetaSp.^2));
% Defocus Phase calculation
phid = (sqrt(-1)*kObj).*costhetaObj;
% Spherical aberration phase calculation
phisa = (sqrt(-1)*k0*depth).*((rindexSp.*costhetaSp)-(rindexObj.*costhetaObj));
% Calculate the optical path difference due to spherical aberrations
OPDSA = exp(phisa);
if nargin>9
    switch varargin{1};
        % Apodizing the pupils for excitation and emission
        case 'illumination'
            pupil = (pupil.*sqrt(costhetaObj));%for illumination
        case 'emission'
            pupil = (pupil./sqrt(costhetaObj));% for emission
    end
end
for zk = 1:nslices
    OPDDefocus = exp(defocus(1, zk).*phid);
    pupilDefocus = pupil.*OPDDefocus;
    pupilSA = pupilDefocus.*OPDSA;
    % Calculate the coherent PSF by using inverse Fourier Transform
    psf(:, :, zk) = ifft2(pupilSA);
    % Calculate the incoherent PSF from the coherent PSF
    psf(:, :, zk) = fftshift(abs(psf(:, :, zk)).^2);
end
psf = sqrt(psf);