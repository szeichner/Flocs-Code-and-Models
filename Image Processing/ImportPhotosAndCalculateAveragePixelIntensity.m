function [rtnConcentrationVersusIntensityArray, concentration] = ImportPhotosAndCalculateAveragePixelIntensity(calibrationPhotoDirectory,...
    pixelCropArray)
%IMPORTPHOTOANDCALCULATEAVERAGEPIXELINTENSITY Import photos from folder 
%and calculate average pixel intensity array
r = '\d\.\d+';
beakerVolume = 0.5;

%get a list of file names in directory
dinfo = dir(char(calibrationPhotoDirectory));
fileNames = {dinfo.name};
fileNames = fileNames(~strncmp(fileNames, '.', 1));
height = 742; %todo: fix this hardcoded number
rtnConcentrationVersusIntensityArray = zeros(height+1, length(fileNames)); 

%iterate over each photo in folder, get the name and add the concentration
%with the intensity to a rtn array
for i = 1:length(fileNames)
    name = char(fileNames(i));
    
    %get concentration from fileName
    concentration = str2num(char(regexp(name, r, 'match'))) / beakerVolume; %#ok<ST2NM>
    rtnConcentrationVersusIntensityArray(i,2) = concentration;
    
    %get average intensity from phtos
    photo = Photo(calibrationPhotoDirectory, name, pixelCropArray);
    avgVect =  photo.PixelMatrixAvgVector;
    rtnConcentrationVersusIntensityArray(1,i) = concentration;
    rtnConcentrationVersusIntensityArray(2:height+1,i) = avgVect';
end

end

