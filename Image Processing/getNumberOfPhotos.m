function [numberOfPhotos] = getNumberOfPhotos(photoFolderLocation)
%GETNUMBEROFPHOTOS Get number of photos at a specific folder location
directoryInfo = dir(photoFolderLocation);
fileNames = {directoryInfo.name};
fileNames = fileNames(~strncmp(fileNames, '.', 1));
numberOfPhotos = length(fileNames);
end

