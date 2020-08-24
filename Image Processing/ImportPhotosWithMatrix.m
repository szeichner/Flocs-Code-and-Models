function [averageVectOverTime] = ImportPhotosWithMatrix(photoDirectory,...
    pixelCropArray, ...
    croppedPhotoPixelHeight)
%IMPORTPHOTOS Import photos from a directory and output matrix of average
%values over time

%get a list of file names in directory
dinfo = dir(photoDirectory);
fileNames = {dinfo.name};
fileNames = fileNames(~strncmp(fileNames, '.', 1));
averageVectOverTime = zeros(croppedPhotoPixelHeight, 1);

%iterate over all the photos in the folder, and add average pixel values to
%matrix
for i = 1:length(fileNames)
    photo = Photo(photoDirectory, fileNames(i), pixelCropArray);
    avgVect =  photo.PixelMatrixAvgVector;
    averageVectOverTime = [averageVectOverTime avgVect]; %#ok<*AGROW>
end

%remove first column of zeros, to return cleaned up matrix
averageVectOverTime(:,1) = [];

end
