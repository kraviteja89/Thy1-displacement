function filtered_data = wave_filt(raw_data, max_lvl, wtype)
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