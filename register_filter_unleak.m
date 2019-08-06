function register_filter_unleak(animal_ID, day_ID, file_num,reg_green,pixel_scale)
% function register_filter_unleak(animal_ID, day_ID, file_num,reg_green,pixel_scale)
% Register, remove crosstalk from green channel into red channel and median filter the data
% creates a structure Thy1 and saves it in the same folder as the one
% specified by animal_ID
% The data structure contains:
% Fs                  : Sampling frequency for imaging
% green_data.original : raw data from the green channel
% red_data.original   : raw data from the red channel
% pixel_shift         : movement in pixels of each frame from frame 1
% no_shift_inds       : the frames where movement is less than 2.5 microns
% filt_size           : filter size for median filter ([3 3 5])
% green_data.processed: processed data from the green channel
% red_data.processed  : processed data from the red channel
% The processed data is also saved as Tiff files for post processing in
% ImageJ or Avizo
% Image registration can be based on the green or the red channel.
% default for image registration is red channel
% 
% Written by Ravi Kedarasetti (2019), department of Engineering Science and
% Mechanics, Pennsylvania State University 
% 
% Input
% animal_ID     : directory containing all the data from the target animal
% day_ID        : date of image collection in format YYMMDD (string)
% file_num      : 3 digit file number for the day (string)
% reg_green     : Logical value (0,1) to indicate green channel for image
%                 registration. Default is 0
% pixel_scale   : scale of the image in microns/pixel. Default is 0.82/3
%                 (16x 0.8NA Nikon objective at 3x digital zoom)


if nargin<5
    pixel_scale = 0.82/3;
end

if nargin<4
    reg_green=0;
end

file_name = [animal_ID '/' day_ID '_' file_num];
mkdir(file_name)
Thy1.fname = file_name;

[Fs, green_data] = read_tiff_downsample([ animal_ID '/Green/' day_ID '_' file_num], 1);
[~, red_data] = read_tiff_downsample([animal_ID '/Red/' day_ID '_' file_num], 1);
Thy1.Fs = Fs;
Thy1.reg_green = reg_green;
Thy1.green_data.original = uint16(green_data);
Thy1.red_data.original = uint16(red_data);
%% Image registration using dft and the green channel
iter = 5;
[ht, wd, n_frames] = size(green_data);
pixel_shift = zeros(2,size(green_data,2));


if reg_green
    first_frame = green_data(:,:,1);
else
    first_frame = red_data(:,:,1);
end
fft_first_frame = fft2(first_frame);

parfor n = 2:n_frames
    if reg_green
        curr_frame = green_data(:,:,n);
    else
        curr_frame = red_data(:,:,n);
    end
    u_curr = zeros(2, iter);
     for it = 1:iter
         D1 = zeros([size(curr_frame) 2]);
         
         D1(:,:,1) = sum(u_curr(1,1:it-1));
         D1(:,:,2) = sum(u_curr(2,1:it-1));
         
         im_curr = imwarp(curr_frame,D1);
         
         fft_curr_frame = fft2(im_curr);
         [curr_disp, ~] = dftregistration(fft_first_frame, fft_curr_frame,4);
         u_curr(1,it) = -curr_disp(4);
         u_curr(2,it) = -curr_disp(3);
     end
    pixel_shift(:,n) = -sum(u_curr,2);
end
Thy1.pixel_shift = pixel_shift;
%% Remove frames where there is excess movement
shift_mag = sqrt(diff(pixel_shift(1,:)).^2 + diff(pixel_shift(2,:)).^2);
shift_thresh = 2.5/pixel_scale;%
shift_bin = double(shift_mag>shift_thresh);
no_shift_inds1 = find(shift_bin<1);
shift_bin = conv(shift_bin,[1 1 1],'full');
shift_bin = shift_bin(2:end-1);
no_shift_inds = find(shift_bin<1);
no_shift_inds = [1 no_shift_inds+1];
Thy1.no_shift_inds = no_shift_inds;
%% Create registered images
reg_wd = ceil(wd+ max(pixel_shift(1,:)) - min(pixel_shift(1,:)));
reg_ht = ceil(ht+ max(pixel_shift(2,:)) - min(pixel_shift(2,:)));
registered_green = zeros(reg_ht, reg_wd, n_frames);
registered_red = zeros(reg_ht, reg_wd, n_frames);

for n = 1:n_frames
    XY_shift = floor(pixel_shift(:,n) - min(pixel_shift,[],2));
    registered_green((1:ht)+XY_shift(2), (1:wd) + XY_shift(1), n) = green_data(:,:,n);
    registered_red((1:ht)+XY_shift(2), (1:wd) + XY_shift(1), n) = red_data(:,:,n);
end
%% Median filter
filt_size = [3 3 5];
filt_iter = 2;
for n = 1:filt_iter
    thresh_green = medfilt3(registered_green, filt_size);
    thresh_red = medfilt3(registered_red,filt_size);
end

thresh_green_no_shift = medfilt3(registered_green(:,:,no_shift_inds), filt_size);
thresh_red_no_shift = medfilt3(registered_red(:,:,no_shift_inds),filt_size);
Thy1.filt_size = filt_size;

%% Remove green data from red using a linear model
fun = @(x) norm(thresh_red(:) - x*thresh_green(:));
options = optimset('Display','iter','MaxIter',100);
x1 = fminsearch(fun,1,options);

thresh_red = thresh_red - x1*thresh_green;
thresh_red = medfilt3(thresh_red,[3 3 5]);
thresh_red(thresh_red<0) = 0;
%% Threshold images
lim = 65535; %16 bit images
thresh_green = thresh_green/prctile(thresh_green(:),99);
thresh_green(thresh_green>1) = 1;
thresh_green = thresh_green*lim;

thresh_red = thresh_red/prctile(thresh_red(:),99);
thresh_red(thresh_red>1) = 1;
thresh_red = thresh_red*lim;
%% Save the processed image files - for post-processing in ImageJ
mkdir([animal_ID '/Green/Processed/'])
mkdir([ animal_ID '/Red/Processed/'])


green_name = [animal_ID '/Green/Processed/' day_ID '_' file_num '.tif'];
red_name = [animal_ID '/Red/Processed/' day_ID '_' file_num '.tif'];

n = 1;
imwrite(uint16(thresh_green(:,:,n)),green_name);
imwrite(uint16(thresh_red(:,:,n)), red_name);


for n = 2:n_frames
    imwrite(uint16(thresh_green(:,:,n)),green_name, 'WriteMode','append');
    imwrite(uint16(thresh_red(:,:,n)), red_name, 'WriteMode','append');
end

%% Save the data structure
Thy1.green_data.processed = uint16(thresh_green);
Thy1.red_data.processed = uint16(thresh_red);
save([file_name '.mat'], 'Thy1', '-v7.3')