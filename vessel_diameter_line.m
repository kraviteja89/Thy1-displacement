function vessel_diameter_line(animal_ID,day_ID,file_num, pixel_scale)
% Calculate vessel diameter using two methods
% 1. Thresholding in radon space - Gao and Drew, JCBFM 2014
% 2. Manually drawn line representing the diameter of the vessel
% The second method could be necessary when intravenous texas red dye is
% used. Texas red can be taken up by macrophages and other cells near blood
% vessels and can affect the calculations when an area of cross-section
% based approach is used.
%
% This function reads the "Thy1" data structure created by 
% the function register_filter_unleak and adds the following variables 
% walking_data          : the velocity data collected using a rotary
%                         encoder(collected at 1000 Hz) in [time ; veloctiy] format
% walking_bin              : binarized walk/rest vector downsampled to imaging
%                         frequency
% centroid              : The centroid of the vessel
% pixel_scale           : scale of the image in microns/pixel
% diameter_radon.raw    : diameter values calculated by method 1
% dimater_line.raw      : diameter values calculated by method 2
% diameter_radon.filt    : diameter values from method 1, filtered
% dimater_line.filt      : diameter values from method 2, filtered
% This function also saves figures for the diameter calculation
% Written by Ravi Kedarasetti (2019), department of Engineering Science and
% Mechanics, Pennsylvania State University 
% 
% Input
% animal_ID     : directory containing all the data from the target animal
% day_ID        : date of image collection in format YYMMDD (string)
% file_num      : 3 digit file number for the day (string)
% pixel_scale   : scale of the image in microns/pixel. Default is 0.82/3
%                 (16x 0.8NA Nikon objective at 3x digital zoom)

if nargin<4
    pixel_scale = 0.82/3;
end
%% Read the Thy1 data structure
file_name = [animal_ID '/' day_ID '_' file_num];
load([file_name '.mat']);
vessel_data = double(Thy1.red_data.processed);
%% Read the walking data from optical sensor
Fs = Thy1.Fs;
walking_data = dlmread([animal_ID '/' day_ID '_' file_num '.txt']);
global acc_cutoff
acc_cutoff= 1e-4;
walking_bin = velocity_proc2(walking_data(:,2),1000,Fs);
n_frames = size(Thy1.green_data.processed,3);
try
    walking_bin = walking_bin(1:n_frames);
catch
    walking_bin = [walking_bin zeros(1,n_frames - length(walking_bin))];
end
Thy1.walking_data = walking_data;
Thy1.walking_bin = walking_bin;
%% Draw ROI for vessel
figure(100)
imagesc(vessel_data(:,:,1))
axis equal
title('Draw rectangular ROI for vessel')
h = imrect(gca);
h.wait();
lims = floor(h.getPosition);

title('Draw Ellipse to find centroid of vessel')
h = imellipse(gca);
h.wait();
vessel_edges = h.getPosition;
Thy1.centroid = [vessel_edges(1)+0.5*vessel_edges(3) vessel_edges(2)+0.5*vessel_edges(4)];

title('Draw Line for radius measurement')
h = imline(gca);
line_ends = h.getPosition;
%% Rotate image so that the selected line is horizontal
hold_image=vessel_data(lims(2):lims(2)+lims(4),lims(1):lims(1)+lims(3),1);
center = [lims(2) lims(1)] + floor((size(hold_image)+1)/2);
a = (line_ends(2,2) - line_ends(1,2))/(line_ends(2,1) - line_ends(1,1));
b = -1;
c = -a*line_ends(1,1) + line_ends(1,2);
dist = (a*center(2) + b*center(1) + c)/sqrt(a^2 + b^2);
ang = rad2deg(atan2((line_ends(2,2) - line_ends(1,2)),(line_ends(2,1) - line_ends(1,1))));
hold_image1 = imrotate(hold_image,ang);
c1 = floor((size(hold_image1)+1)/2);
line_ind = floor(c1(1) +sign(line_ends(2,1) - line_ends(1,1))*dist);
%% Calculate the radius using both methods
pixel_area = zeros(1,n_frames);
pixel_diameter_line = zeros(1,n_frames);
the_angles=1:1:180;
rtd_threhold=.5;%0.5 'half max'
irtd_threhold=.2;%0.2 gets rid of the noise from the inverse radon

sigma = 2;
sz = 12;    % length of gaussFilter vector
x = linspace(-sz / 2, sz / 2, sz);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize

for n = 1:n_frames
     hold_image=vessel_data(lims(2):lims(2)+lims(4),lims(1):lims(1)+lims(3),n);
     
     % measure the diameter along the line
    hold_image1 = imrotate(hold_image,ang);
    mean_line = conv(mean(hold_image1(line_ind-5:line_ind+5,:)),gaussFilter,'same');
    pixel_diameter_line(n) =  length(find(mean_line>0.5*max(mean_line)));
     
     % measure the diameter using radon transform
     hold_image=hold_image-mean(hold_image(:));
     radon_hold_image=radon(hold_image,the_angles);
     
     for k=1:length(the_angles)
         % normalize so that for each angle, the transfomed image is between
         % 0 and 1
         radon_hold_image(:,k)=radon_hold_image(:,k)-min(radon_hold_image(:,k));
         radon_hold_image(:,k)=radon_hold_image(:,k)/max(radon_hold_image(:,k));
         % find the peak of the projection at this angle
         [maxpoint(k),maxpointlocation(k)]=max(radon_hold_image(:,k));
         % threshold at half maximum
         try
             [~,min_edge(k)]=max(find(radon_hold_image(1:maxpointlocation(k),k)<rtd_threhold));
             [~,max_edge(k)]=max(find(radon_hold_image(maxpointlocation(k)+1:end,k)>rtd_threhold));
         catch
             min_edge(k)=0;
             max_edge(k)=0;
         end
         
         radon_hold_image(1:min_edge(k),k)=0;
         radon_hold_image((maxpointlocation(k)+max_edge(k)):end,k)=0;
     end
     
     irtd_norm=iradon(double(radon_hold_image>rtd_threhold*max(radon_hold_image(:))),(the_angles),'linear','Hamming');
    
    
    [cc,l] = bwboundaries(irtd_norm>irtd_threhold*max(irtd_norm(:)),4);
    numPixels = cellfun(@length,cc);
    [~,idx] = max(numPixels);
    area_filled=regionprops(l,'FilledArea','Image','FilledImage');
    pixel_area(n) = length(find(area_filled(idx).FilledImage));
end
%% Wavelet based filtering
pixel_area1 = nan(size(pixel_area));
pixel_area1(Thy1.no_shift_inds) = pixel_area(Thy1.no_shift_inds);

pixel_diameter_line1 = nan(size(pixel_diameter_line));
pixel_diameter_line1(Thy1.no_shift_inds) = pixel_diameter_line(Thy1.no_shift_inds);


pixel_area_filt = wave_filt(pixel_area1, 12,'bior3.3');
pixel_diameter_line_filt = wave_filt(pixel_diameter_line1, 12,'bior3.3');
%% Save the calculated diameters
Thy1.pixel_scale = pixel_scale;
Thy1.diameter_radon.raw = pixel_scale*2*sqrt(pixel_area/pi);
Thy1.diameter_radon.filt = pixel_scale*2*sqrt(pixel_area_filt/pi);

Thy1.diameter_line.raw = pixel_scale*pixel_diameter_line;
Thy1.diameter_line.filt = pixel_scale*pixel_diameter_line_filt;
%% Plot the radon calculated diameters
walking_bin(~walking_bin) = nan;

radoncontours=contourc(double(irtd_norm(1:end,1:end)),double([irtd_threhold irtd_threhold]*max(irtd_norm(:))));
okpoints=find(radoncontours(1,:)>1);
use_points=NaN*zeros(size(radoncontours));
use_points(1,okpoints)=radoncontours(1,okpoints);
use_points(2,okpoints)=radoncontours(2,okpoints);
fig = figure();
set(fig, 'Position', [1 1 1200 800])
subplot(2,2,1)
hold off
imagesc(hold_image)
axis image
axis xy
colormap gray
title('Last Frame')

subplot(2,2,2)
imagesc(irtd_norm)
hold on
plot(use_points(1,:),use_points(2,:))
hold off
axis image
axis xy
title('After Radon Cleanup')

t = (1:n_frames)/Fs;
subplot(2,1,2)
plot(t,Thy1.diameter_radon.raw)
hold on
plot(t,Thy1.diameter_radon.filt)
hold off
xlabel('time(s)')
ylabel('Vessel Diameter(\mum)')
yyaxis right
plot(t, walking_bin, 'ro', 'MarkerFaceColor','r')
set(gca, 'YTick',[])
saveas(fig, [file_name '_Diameters.png'])
% saveas(fig, [file_name '_Diameters.fig'])
%% Plot the diameters calculated from the line
fig = figure();
set(fig, 'Position', [1 1 1200 800])
subplot(2,2,1)
imagesc(hold_image1)
hold on
axis image
axis xy
title('Box for diameter measurements')
plot([1, size(hold_image1,2)],[line_ind-5,line_ind-5 ], 'r','LineWidth',2)
plot([1, size(hold_image1,2)],[line_ind+5,line_ind+5 ], 'r','LineWidth',2)
hold off

subplot(2,2,2)
plot(mean_line)
title('Mean pixel values in the box')


subplot(2,1,2)
plot(t, Thy1.diameter_line.raw)
hold on
plot(t, Thy1.diameter_line.filt)
hold off
xlabel('time(s)')
ylabel('Vessel Diameter(\mum)')
yyaxis right
plot(t, walking_bin, 'ro', 'MarkerFaceColor','r')
set(gca, 'YTick',[])
saveas(fig, [file_name '_Diameters_line.png'])
% saveas(fig, [file_name '_Diameters_line.fig'])


save([file_name '.mat'], 'Thy1', '-v7.3')
close all
