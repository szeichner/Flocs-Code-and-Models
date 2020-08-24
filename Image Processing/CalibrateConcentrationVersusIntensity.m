function [rtnCalibrationPolynomial] = CalibrateConcentrationVersusIntensity(experimentType, ...
    calibrationPhotoFolder, ...
    croppedDim, ...
    polyFitOrder,...
    box_split)
%CALIBRATECONCENTRATIONVERSUSINTENSITY Parse the name of the experiment and
%get concentration to intensity calibration
rtnCalibrationPolynomial = zeros(box_split, 2);

%parse folder name to get the type of experiment to calibrate for
splitExpTypeName = strsplit(char(experimentType),'_');

%concatenate name to get the correct calibration relationship
calibrateFolderName = strcat(calibrationPhotoFolder, '/', splitExpTypeName(1));

%get matrix to look at settling contour per experiment over time
%(concentration, intensity)
concentrationVersusIntensityArray = ImportPhotosAndCalculateAveragePixelIntensity(calibrateFolderName, croppedDim);
concentrationVector = concentrationVersusIntensityArray(1, :); %extract header row to reference later

%divide into subboxes for analysis
subBoxLength = round(size(concentrationVersusIntensityArray,1) / box_split) - 2; %subtract 2 to account for header row
startHeight = 2; %skip header row

for j=1:box_split
    %divide up submatrices
    d_subBox = concentrationVersusIntensityArray(startHeight:startHeight+subBoxLength, :);
    %average each subbox 
    meanMatrix = mean(d_subBox(:, :));
    
    %plot the subbox versus concentration header
    ft=fittype('exp1');
    expFt = fit(meanMatrix',concentrationVector',ft);
    rtnCalibrationPolynomial(j, 1) = expFt.a;
    rtnCalibrationPolynomial(j, 2) = expFt.b;
    
    figure(2)
    hold on
    plot(expFt, meanMatrix, concentrationVector)
    title("Experiment Concentration versus Intensity Calibration")
    legend;
    hold off
    
    startHeight = startHeight + subBoxLength;
end
end

