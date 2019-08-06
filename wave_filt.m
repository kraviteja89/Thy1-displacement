function filtered_data = wave_filt(raw_data, max_lvl, wtype)
% function filtered_data = wave_filt(raw_data, max_lvl, wtype)
% Filter(de-noise) single channel data using wavelet decomposition
% Data can contain 'nan' values.
% The data is median filtered(filter width 3) to remove nan values.
% Any remaining nan values are replaced with the last valid value
% 1. The data is decomposed into wavelet coefficients at max_lvl levels
% 2. A threshold for the wavelet coefficients is calculated using the
%   visushrink method (Donoho, David L., et al. "Wavelet shrinkage:
%   asymptopia?." Journal of the Royal Statistical Society: Series B
%   (Methodological) 57.2 (1995): 301-337.)
% 3. The values below threshold are set to zero(hard threhsolding)
% 4. The thresholded wavelet coefficients are transformed back into the
%   signal domain.
% Written by Ravi Kedarasetti (2019), department of Engineering Science and
% Mechanics, Pennsylvania State University 
%
% Inputs
% raw_data  : row or column vector containing single channel noisy data
% max_lvl   : Maximum level of wavelet decomposition. If unspecified, a value of 
%             floor(log2(length(raw_data)) is used
% wtype     : family of wavelet to be used. If unspecified 'bior3.3' will
%             be used
%
% Outputs
% filtered_data     : denoised data

% If the wavelet family to be used is not specified, use the default
% 'bior3.3'
if nargin<3
    wtype = 'bior3.3';
end

if nargin<2
    max_lvl = floor(log2(length(raw_data)));
end

raw_data = medfilt1(raw_data,3);
nan_inds = find(isnan(raw_data));

try
    if nan_inds(1)==1
        raw_data(nan_inds(1)) = 0;
        nan_inds = nan_inds(2:end);
    end
catch
    fprintf('No frames with movement more than 2.5 micron \n')
end
for n3 = 1:length(nan_inds)
    raw_data(nan_inds(n3)) = raw_data(nan_inds(n3)-1);
end

[wcoeffs,l] = wavedec(raw_data, max_lvl, wtype);
sig = median(abs(wcoeffs))/0.6745;
lambda = sig*sqrt(2*log(length(raw_data)));
wcoeffs_hard = wcoeffs;
wcoeffs_hard(abs(wcoeffs_hard)<lambda) = 0;
filtered_data = waverec(wcoeffs_hard,l,wtype);