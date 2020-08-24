clear all
clc
clf

%FLOCIMAGEANALYSIS: Imports a folder of photos, and exports the matrix to an output
%txt file
%Written by: ssz
%Last edited: June 10, 2019

%%%%%%%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get current directory
directoryFolder = cd();

% Set up folder import information for photos, and type of experiment for
% export
photoFolder = '/Users/sarahzeichner/Documents/MATLAB/Flocs/FlocJarPhotos/1000mL90Mesh_v3_080519/';

% Note: the experiment type here should match the one in "Floc Data
% Analysis" so the right file gets imported
experimentType = "1000mL90Mesh_v3_4.0g"; %change this so that the correct name gets added to the output file

%Choose specific crop dimensions to propogate across photos
croppedDim = [1.481510000000000e+03,8.315100000000000e+02,1.979800000000000e+02,7.419800000000000e+02];
%heightCalCrop = [1.321510000000000e+03,10.515100000000000e+02,1.979800000000000e+02,7.419800000000000e+02];

%%%%%%%%%%%%%%%%%%%%%%%%%%  MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Import and calibrate photos

%crop photos and import to evaluate average pixel values
croppedImage = CropPhotoForFocusedAnalysis(char(photoFolder), croppedDim);
[x1, y1, z1] = size(croppedImage);

%get matrix to look at settling contour per experiment over time
importedPhotoMatrix = ImportPhotosWithMatrix(char(strcat(photoFolder)), ...
    croppedDim,...
    x1);

dlmwrite(strcat(experimentType,"_output.txt"), importedPhotoMatrix);