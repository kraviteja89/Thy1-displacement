% Calculate the arterial wall and brain tissue displacements for each dataset
% 
% This file should be run as the first part of data processing. The
% resulting displacements calculated from running this file are
% manually verified using ImageJ and the locations are added to the excel
% file. 
% The displacement_summary.m file is run to give the final summary
% Select the excel file that has the following colums (separated by %%) 
% Animal ID	%% Day ID %%File Number%%Pixel Scale
% The excel file should be in the same folder as the Data
%
% Written by Ravi Kedarasetti (2019), department of Engineering Science and
% Mechanics, Pennsylvania State University 
%% Read the excel file
[file_name, path_name] = uigetfile('E:/*.xlsx', 'Select summary excel file');
data_IDS = readtable([path_name file_name]);

reg_green = 0;
grid_size = 64;
step_size = 16;
%% Preprocessing
for n= 1:size(data_IDS,1)
    register_filter_unleak([path_name data_IDS.AnimalID{n}], data_IDS.DayID{n}, data_IDS.FileNumber{n} ,reg_green, data_IDS.PixelScale(n));
end    
%% Diameter calculation
for n= 1:size(data_IDS,1)
    vessel_diameter_line([path_name data_IDS.AnimalID{n}], data_IDS.DayID{n}, data_IDS.FileNumber{n}, data_IDS.PixelScale(n));
end
%% Tissue displacement calculation
for n= 1:size(data_IDS,1)
    thy1_displacement_fields_iterative([path_name data_IDS.AnimalID{n}], data_IDS.DayID{n}, data_IDS.FileNumber{n},data_IDS.PixelScale(n), grid_size, step_size)
    close all
end