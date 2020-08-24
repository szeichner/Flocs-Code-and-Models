function [croppedImage, croppedDim] = CropPhotoForFocusedAnalysis(photoFolderLocation, croppedDim)
%CHOOSECROPDIMENSIONS Based on first photo in folder series, user chooses
%where to crop the photos for analysis

directoryInfo = dir(photoFolderLocation);
fileNames = {directoryInfo.name};
fileNames = fileNames(~strncmp(fileNames, '.', 1));

firstImage = imread(char(strcat(photoFolderLocation, '/', fileNames(1))));
% figure(1)
% imshow(firstImage);

[croppedImage] = imcrop(firstImage, croppedDim);
% figure(2)
% imshow(croppedImage);
end

