function thy1_displacement_fields_iterative(animal_ID,day_ID,file_num, pixel_scale, grid_size, step_size)
% Calculate the displacement fields from two-photon data(tiff files)
% 1. Fit the vessel diameter calculated to a gamma function
%    If goodness of fit R_squared is < 0.5, skip the displacementfield
%    calculation.
% 2. Find a reference frame - based on the binarized walking data
% 3. Divide the reference frame into blocks and select bright blocks
% 4. Calculate displacement values for the selected blocks
% 5. Calculate direction of displacement
% 6. Plot and save displacements that are radially outward, and have a good
%    correlation with the diameter change
%
% Written by Ravi Kedarasetti (2019), department of Engineering Science and
% Mechanics, Pennsylvania State University 
% 
% Input
% animal_ID     : directory containing all the data from the target animal
% day_ID        : date of image collection in format YYMMDD (string)
% file_num      : 3 digit file number for the day (string)
% pixel_scale   : scale of the image in microns/pixel. Default is 0.82/3
% grid_size     : size of each block in pixels(default 64)
% step_size     : distance between consecutive blocks(default is 16)

if nargin<5
    grid_size = 64;
    step_size = 16;
end

if nargin<4
    pixel_scale = 0.82/3;
end


%% Load data
file_name = [animal_ID '/' day_ID '_' file_num];
load([file_name '.mat']);
%% calculate the walking data
Fs = Thy1.Fs;
[ht, wd, n_frames] = size(Thy1.green_data.processed);
walking_bin = Thy1.walking_bin;
walk_perc = sum(walking_bin)/n_frames;
%% HRF parameters for vessel diameter
global wave_len;
wave_len = floor(Fs*10);%10 second maximum
wave_t = 0:wave_len;

x0 = [1,5,5,2];
err_fun = @(x)gamma_err(x,Thy1.diameter_radon.filt, walking_bin);
options = optimset('Display','iter','MaxIter',1000, 'TolFun',1e-2,'TolX',1e-4,'MaxFunEvals',5000);
dia_params = fminsearch(err_fun,x0,options);

err_fun_line = @(x)gamma_err(x,Thy1.diameter_line.filt, walking_bin);
dia_params_line = fminsearch(err_fun_line,x0,options);

Rsquared = 1-gamma_err(dia_params,Thy1.diameter_radon.filt, walking_bin);
Rsquared_line = 1-gamma_err(dia_params_line,Thy1.diameter_line.filt, walking_bin);

if Rsquared>Rsquared_line
    Thy1.dia_fit.Rsquared = Rsquared;
    Thy1.dia_fit.data = Thy1.diameter_radon.filt;
    Thy1.dia_fit.params = dia_params;
else
    Thy1.dia_fit.Rsquared = Rsquared_line;
    Thy1.dia_fit.data = Thy1.diameter_line.filt;
    Thy1.dia_fit.params = dia_params_line;
end

%% Plot data for vessel_diamater fit
the_data = Thy1.dia_fit.data;
the_data = the_data - mean(the_data);
HRF = gamma_fun(Thy1.dia_fit.params);
pred = conv(walking_bin,HRF,'same');
pred = pred - mean(pred);
walking_bin(~walking_bin) = nan;

fig1 = figure(80);
set(fig1, 'Position', [100 100 1600 800])
t = (1:n_frames)/Fs;

plot(t,the_data,'b-' ,'LineWidth',2)
hold on
plot(t,pred, 'r-','LineWidth',2)
plot(t,walking_bin*max(the_data), '^b')
hold off
set(gca, 'FontSize',12)
xlabel('time(s)', 'FontSize',15)
legend('Actual data', 'HRF Prediction', 'Binarized Running', 'Location','southwest');
ylim([min(the_data) max(the_data)]*1.2)
xlim([t(1) t(end)+50])
ylabel('Vessel Diameter change(\mum)', 'FontSize',15)
title(['Vessel diamater fit to HRF, R^2 = ' num2str(round(Thy1.dia_fit.Rsquared,2))], 'FontSize',20)


fig2 = figure(90);
set(fig2, 'Position', [100 100 1600 800])
plot((1:length(HRF))/Fs, HRF)
title('HRF')
xlabel('time(s)')
set(gca, 'FontSize',14)

fig = inset(100,fig1,fig2);
saveas(fig, [file_name '_Diameter_fit.png'])
saveas(fig, [file_name '_Diameter_fit.fig'])


if Thy1.dia_fit.Rsquared < 0.5
    save([file_name '.mat'], 'Thy1', '-v7.3')
    return
end
%% Find reference frame - the average frame for a long resting period
green_data = Thy1.green_data.processed;
rest_length = 20; % seconds
run_frames = find(~isnan(walking_bin));
while rest_length>1
    long_rest = run_frames(diff(run_frames)>rest_length*Thy1.Fs);
    if ~isempty(long_rest)
        st_frame = floor(long_rest(1) + Thy1.Fs);
        last_frame = floor(long_rest(1) + rest_length*Thy1.Fs);
        base_dia = median(Thy1.dia_fit.data(st_frame:last_frame));
        reference_frame = median(green_data(:,:,st_frame:last_frame),3);
        break
    else
        rest_length = rest_length-1;
        if rest_length==1
            fprintf('No rest period found for %s...taking first frame as reference \n',file_name)
            reference_frame = green_data(:,:,1);
            base_dia = median(Thy1.dia_fit.data);
            st_frame = 1;
            last_frame = n_frames;
        end
    end
end
Thy1.base_dia = base_dia;
Thy1.rest_length = rest_length;
Thy1.rest_frames = [st_frame last_frame];
Thy1.reference_frame = reference_frame;
%% roi mask based on brightness
iter = 5;
Thy1.grid_size = grid_size;
Thy1.step_size = step_size;
yboxes = floor((ht-grid_size)/step_size)+1;
xboxes = floor((wd-grid_size)/step_size)+1;

ux = nan(yboxes, xboxes, n_frames);
uy = nan(yboxes, xboxes, n_frames);
err_data = ones(yboxes, xboxes, n_frames);

bright_mask = zeros(yboxes,xboxes);
for n2 = 1:yboxes
    for n1 = 1:xboxes
        first_box = reference_frame((n2-1)*step_size+1:(n2-1)*step_size+grid_size,(n1-1)*step_size+1:(n1-1)*step_size+grid_size);
        bright_mask(n2,n1) = sum(first_box(:));
    end
end

[~,bright_inds] = sort(bright_mask(:));
bright_inds = bright_inds(floor(length(bright_inds)*85/100):end);
%% Calculate displacement fields
no_shift_inds = Thy1.no_shift_inds;
parfor f = 1:n_frames
    if ~ismember(f, no_shift_inds)
        continue
    end

    curr_frame = green_data(:,:,f);    
    for n1 =  1:xboxes
        for n2 = 1:yboxes
              if ~ismember((n1-1)*yboxes + n2,bright_inds)
                  continue
              end
            
            first_box = reference_frame((n2-1)*step_size+1:(n2-1)*step_size+grid_size,(n1-1)*step_size+1:(n1-1)*step_size+grid_size);
            fft_first_box = fft2(first_box);
            u_curr = zeros(2, iter);
            
            for it = 1:iter
                D1 = zeros([size(reference_frame) 2]);
                
                D1(:,:,1) = sum(u_curr(1,1:it-1));
                D1(:,:,2) = sum(u_curr(2,1:it-1));
                
                im_curr = imwarp(curr_frame,D1);
                curr_box = im_curr((n2-1)*step_size+1:(n2-1)*step_size+grid_size,(n1-1)*step_size+1:(n1-1)*step_size+grid_size);
                
                
                fft_curr_box = fft2(curr_box);
                
                [curr_disp, ~] = dftregistration(fft_first_box, fft_curr_box,4);
                u_curr(1,it) = -curr_disp(4); 
                u_curr(2,it) = -curr_disp(3);
            end
            u_mag = sqrt(sum(u_curr(1,:))^2 +  sum(u_curr(2,:))^2);
            u_last = sqrt(u_curr(1,iter)^2 +  u_curr(2,iter)^2);
            if u_last/u_mag < 0.01 || u_mag<0.25
                ux(n2,n1,f) = pixel_scale*sum(u_curr(1,:));
                uy(n2,n1,f) = pixel_scale*sum(u_curr(2,:));
                err_data(n2,n1,f) = curr_disp(1);
            end
        end
    end
end
%% Pick displacement values based on error
err_thresh = 0.7;
Thy1.ux = ux;
Thy1.uy = uy;
Thy1.err_data = err_data;
Thy1.err_thresh = err_thresh;

err_mask = err_data<err_thresh;
ux1 = nan(size(ux));
ux1(err_mask) = ux(err_mask);

uy1 = nan(size(uy));
uy1(err_mask) = uy(err_mask);
%% Make a mask based on expected direction
[~,top_def_inds] = sort(Thy1.dia_fit.data(Thy1.no_shift_inds));
top_def_inds = top_def_inds(floor((1-walk_perc)*length(top_def_inds)):end);
utheta_med = nan(yboxes,xboxes);

for n1 = 1:xboxes
    for n2 = 1:yboxes
        
        if ~ismember((n1-1)*yboxes + n2,bright_inds)
            continue
        end
        the_line_x = squeeze(ux1(n2,n1,Thy1.no_shift_inds));
        the_line_x = medfilt1(the_line_x,3,'omitnan');
        
        the_line_y = squeeze(uy1(n2,n1,Thy1.no_shift_inds));
        the_line_y = medfilt1(the_line_y,3,'omitnan');

        utheta = atan2(the_line_y, the_line_x);
        utheta_med(n2,n1) = median(utheta(top_def_inds), 'omitnan');
    end
end

[centers_x, centers_y] = meshgrid(1:xboxes, 1:yboxes);
centers_x = (centers_x-1)*step_size + grid_size/2 - Thy1.centroid(1);
centers_y = (centers_y-1)*step_size + grid_size/2 - Thy1.centroid(2);
thetas_expected = atan2(centers_y, centers_x);

theta_mask = abs(thetas_expected - utheta_med);
theta_mask(theta_mask>pi) = theta_mask(theta_mask>pi) - 2*pi; 
theta_mask = theta_mask < 30/180*pi;
%%
Thy1.utheta = utheta_med;
Thy1.theta_mask = theta_mask;
%%
reference_frame = Thy1.reference_frame;
green_image = zeros([size(reference_frame),3]);
green_image(:,:,2) = reference_frame/max(reference_frame(:));
red_frame = double(Thy1.red_data.processed(:,:,1));
red_image = zeros([size(reference_frame),3]);
red_image(:,:,1) = red_frame/max(red_frame(:));
red_image(:,:,3) = red_frame/max(red_frame(:));

red_thresh = prctile(red_frame(:), 95);
red_mask = red_frame>red_thresh;
%%
dias = Thy1.dia_fit.data - Thy1.base_dia;
dias_norm = dias - mean(dias);
dias_norm = dias_norm/norm(dias_norm);

corr_vals = zeros(yboxes,xboxes);
trusted_n1 = [];
trusted_n2 = [];

for n1 = 1:xboxes
    for n2 = 1:yboxes
        
        if ~theta_mask(n2,n1)
            continue
        end
        
        the_line_x = medfilt1(squeeze(ux1(n2,n1,:)),3,'omitnan');
        the_line_y = medfilt1(squeeze(uy1(n2,n1,:)),3,'omitnan');
        
        
        ur = cos(thetas_expected(n2,n1))*the_line_x + sin(thetas_expected(n2,n1))*the_line_y;
        nan_inds = sum(isnan(ur));
        
        if nan_inds > 0.3*length(Thy1.no_shift_inds)
            continue
        end
        
        
       ur = wave_filt(ur,12,'bior3.3');
       ur1 = ur - mean(ur);
       ur1 = ur1/norm(ur1);
       corr_vals(n2,n1) = sum(dias_norm.*ur1');
        
       if corr_vals(n2,n1)<0.8
           continue
       end
       
       trusted_n1 = [trusted_n1 n1];
       trusted_n2 = [trusted_n2 n2];
       
        fig = figure();
    set(fig, 'Position', [1 1 2000 700])
    subplot(1,3,1)
    imshow(green_image*255)
    hold on
    h = imshow(255*red_image);
    set(h, 'AlphaData',red_mask)
    
    x1 = (n1-1)*Thy1.step_size+(Thy1.grid_size+1)/2;
    y1 = (n2-1)*Thy1.step_size+(Thy1.grid_size+1)/2;
    dx = 50*cos(thetas_expected(n2,n1));
    dy = 50*sin(thetas_expected(n2,n1));
    quiver(x1,y1,dx,dy,'Color', 'y', 'LineWidth',3, 'MaxHeadSize',5)
    
    axis equal
    set(gca, 'YDir', 'normal')
    axis off
    
    t = (1:n_frames)/Fs;
    t2 = (1:length(walking_bin))/Fs;
    subplot(1,3,[2 3])
    plot(t,dias/2, 'LineWidth',2)
    hold on
    plot(t,ur, 'k' , 'LineWidth',2)
    xlabel('Time (s)')
    
    
    ylabel('Radial displacement(\mum)')
    yyaxis right
    plot(t2, walking_bin, 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
    set(gca, 'YTick', [])
    ylim([0 1.2])
    legend( 'Vessel wall', 'Brain tissue', 'Binarized Running')
    set(gca, 'FontSize', 16)
    
    saveas(fig, [animal_ID '/refined_displacements/' day_ID '_' file_num '_radial_displacement' num2str(length(trusted_n1)) '.png'])
    saveas(fig, [animal_ID '/refined_displacements/' day_ID '_' file_num '_radial_displacement' num2str(length(trusted_n1)) '.pdf'])
       
        
    end
end
%%
Thy1.corr_vals = corr_vals;
Thy1.trusted_cols = trusted_n1;
Thy1.trusted_rows = trusted_n2;
save([file_name '.mat'], 'Thy1', '-v7.3')
