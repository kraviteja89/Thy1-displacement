% Calculate the average normalized arterial wall and brain tissue displacements
% Select the excel file that has the following colums (separated by %%) 
% Animal ID	%% Day ID %%File Number%%Pixel Scale%%Location Number
% The excel file should be in the same folder as the Data
%
% Written by Ravi Kedarasetti (2019), department of Engineering Science and
% Mechanics, Pennsylvania State University 
%%
clear all
close all
clc
%% Read the excel file
[file_name, path_name] = uigetfile('E:/*.xlsx', 'Select summary excel file');
data_IDS = readtable([path_name file_name]);
Fs = 3;
IIR_size = 21; % size of the impulse response function (based on an imaging frequency of 3Hz)
vessel_impulse = zeros(IIR_size, size(data_IDS,1));
tissue_impulse = zeros(IIR_size, size(data_IDS,1));
%% Calculate the impulse response for each dataset
for n = 1:size(data_IDS,1)
    % read the mat file
    animal_ID = [path_name data_IDS.AnimalID{n}];
    day_ID = data_IDS.DayID{n};
    file_num = data_IDS.FileNumber{n};
    load([animal_ID '/' day_ID '_' file_num '.mat'])
    
    % Recalculate the brain tissue displacement at the selected location
    thetas_expected = Thy1.thetas_expected;
    [ht, wd, n_frames] = size(Thy1.green_data.processed);
    yboxes = floor((ht-Thy1.grid_size)/Thy1.step_size)+1;
    xboxes = floor((wd-Thy1.grid_size)/Thy1.step_size)+1;
    
    err_mask = Thy1.err_data<Thy1.err_thresh;
    ux1 = nan(size(Thy1.ux));
    ux1(err_mask) = Thy1.ux(err_mask);
    
    uy1 = nan(size(Thy1.uy));
    uy1(err_mask) = Thy1.uy(err_mask);
    
    trusted_n1 = Thy1.trusted_cols ;
    trusted_n2 = Thy1.trusted_rows ;
    
    n_disp = data_IDS.LocationNumber(n);
    n1 = trusted_n1(n_disp);
    n2 = trusted_n2(n_disp);
    
    % filter the displacement using wavelet based filter
    the_line_x = medfilt1(squeeze(ux1(n2,n1,:)),3,'omitnan');
    the_line_y = medfilt1(squeeze(uy1(n2,n1,:)),3,'omitnan');
    ur = cos(thetas_expected(n2,n1))*the_line_x + sin(thetas_expected(n2,n1))*the_line_y;
    ur = wave_filt(ur,12,'bior3.3');
    
    %calculate the impulse response using deconvolution
    dias = Thy1.dia_fit.data;
    toffset = ceil(Thy1.dia_fit.params(4))+3;
    dias = dias';
    walk_bin = Thy1.walking_bin;
    
    IR_vessel=IR_analytic(walk_bin - mean(walk_bin),dias - mean(dias),toffset,IIR_size-1);
    IR_disp=IR_analytic(walk_bin - mean(walk_bin),ur - mean(ur),toffset,IIR_size-1);
    
    vessel_impulse(:,n) = IR_vessel.IR/2;
    tissue_impulse(:,n) = IR_disp.IR;
end

%% Normalize the impulse response functions
[max_vessel_disp, vessel_ids] = max(vessel_impulse(5:15,:));
max_tissue_disp = max(tissue_impulse(5:15,:));
vessel_impulse_norm = zeros(IIR_size+ max(vessel_ids) - min(vessel_ids), size(data_IDS,1));
tissue_impulse_norm = zeros(IIR_size+ max(vessel_ids) - min(vessel_ids), size(data_IDS,1));
for n = 1:size(data_IDS,1)
    vessel_impulse_norm(max(vessel_ids) - vessel_ids(n) + 1: IIR_size +max(vessel_ids)- vessel_ids(n), n) = vessel_impulse(:,n)/max_vessel_disp(n);
    tissue_impulse_norm(max(vessel_ids) - vessel_ids(n) + 1: IIR_size +max(vessel_ids)- vessel_ids(n), n) = tissue_impulse(:,n)/max_tissue_disp(n);
end

%% Calculate and plot the mean and standard deviation of impulse response
t_vec = (1:size(vessel_impulse_norm))/Fs;
fig = figure();
set(fig, 'Position', [100 100 1000 800])
hold on
plot(t_vec, mean(vessel_impulse_norm,2),'m','LineWidth', 2)
plot(t_vec, mean(tissue_impulse_norm,2),'g','LineWidth', 2)
h_tissue = fill([t_vec fliplr(t_vec)],  [mean(tissue_impulse_norm') + std(tissue_impulse_norm'),...
    fliplr(mean(tissue_impulse_norm') - std(tissue_impulse_norm'))], 'g', 'EdgeColor', 'none');
set(h_tissue, 'FaceAlpha', 0.1)

xlabel('Time(s)')
legend('Normalized vessel response', 'Normalized tissue response', '\pm std')
set(gca, 'FontSize', 16)
saveas(fig, 'Impulse_response.png')
saveas(fig, 'Impulse_response.pdf')
