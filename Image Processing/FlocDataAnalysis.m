clear all
clc
clf

%FLOCDATAANALYSIS: Import matrices, convert intensities to concentrations and
%calculate settling velocities
%Written by: ssz
%Last edited: June 10, 2019

%%%%%%%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up folders, plotting settings, etc.
% Get current directory
directoryFolder = cd();

%Set up photo calibration for experiment of interest
photoFolder = 'CaMont'; %use 'CaMont_0.2g' or 'Kaolinite_0.2g' or '90Sand_2.0g' for folder naming
calibrationPhotoFolder = '/Users/sarahzeichner/Documents/MATLAB/Flocs/FlocJarPhotos/calibration';
calibrationPolynomialOrder = 1; %specify the order polynomial wanted to fit the calibration curve
intialConcentration = .2  / 0.5; %g/L
croppedDim = [1.481510000000000e+03,8.315100000000000e+02,1.979800000000000e+02,7.419800000000000e+02];
delta_Z = 0.059; %meters, based on height of cropped dimensions

%Set up data culling and processing
dataCullingParam = 0.2;
n_box_split  = 1;
slopeThreshhold = 0.005; %threshhold above which it counts as settling

%Define files to import to analyze; separate folder for each
%exopolymer/clay combination
importDataFolder = '/Users/sarahzeichner/Documents/MATLAB/Flocs/code_data_input';
dinfo = dir(char(importDataFolder));
fileNames = {dinfo.name};
fileNames = fileNames(~strncmp(fileNames, '.', 1));

%Set up output information for csv export
outputColumnNames = {'startHeight';'endHeight';'velocity'};
velocityOutputMatrix = zeros(n_box_split*length(fileNames), length(outputColumnNames));
rowLabels = { };
csvName = 'settlingVelocity.csv';

%%%%%%%%%%%%%%%%%%%%%%%%%%  MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calibrate concentration versus intensity for each experiment
concentrationCalibrationPolyMatrix = CalibrateConcentrationVersusIntensity(photoFolder, ...
    calibrationPhotoFolder, ...
    croppedDim, ...
    calibrationPolynomialOrder, ...
    n_box_split);

% Process data for experiments of interest
outputTableRow = 1;

figure(1)
hold on
for i=1:length(fileNames)
    %read in matrix for analysis
    M = dlmread(char(strcat(importDataFolder,'/',fileNames(i))));
    mean_M = mean(M);
    
    %get length of time of experiment based on matrix size
    maxTime = length(mean_M);
    timeAxisVector = 0:1:maxTime-1;
    timeAxisVector = timeAxisVector';
    startHeight = 1;
    
    for j=1:n_box_split
        
        %split up into a sub box
        subBoxLength = round(size(M,1) / n_box_split) - 1;
        subM = M(startHeight:startHeight+subBoxLength, :);
        
        %convert matrix intensity pixel values to concentration
        % a * exp(b*x) exponential form of concentration relationship
        concentrationM = concentrationCalibrationPolyMatrix(j,1) * exp(concentrationCalibrationPolyMatrix(j,2) * subM);
        %concentrationM2 = mean((11.5 * exp(-6 * subM)) / intialConcentration);
        %concentrationM3 = mean((11.5 * exp(-6.7 * subM)) / intialConcentration);
        
        %divide each concentration by intial concentration to normalize
        normConcentrationM = concentrationM / intialConcentration;
        normConcentrationVector = mean(normConcentrationM);
        timeVersusConcentrationVector = [timeAxisVector, normConcentrationVector'];
        
        % Cull the data and return the culled vector
        culledTimeVersusConcentrationVector = CullSettlingData(timeVersusConcentrationVector, slopeThreshhold); 
        
        %Plot all the data (unculled)
        scatter(timeVersusConcentrationVector(:,1), timeVersusConcentrationVector(:,2));
        %scatter(timeAxisVector, concentrationM2', "+")
        %scatter(timeAxisVector,concentrationM3', "o")
        legend(fileNames, 'Location', 'best');
        
        %Perform least squares regression on culled data
        [intercept, slope] = CalculateLeastSquares(culledTimeVersusConcentrationVector(2:5,1), culledTimeVersusConcentrationVector(2:5,2));
        
        %Calculate settling velocity based on slope and convert to m/s
        %(velocity = concentration / sec units)
        velocity = -1 * slope * (delta_Z/n_box_split); %the height of the box is constant even if the concentration settles out at different points
        
        %Write to vector to output
        endHeight = startHeight + subBoxLength;
        
        %output the heights in meters
        velocityOutputMatrix(outputTableRow,1) = startHeight / delta_Z;
        velocityOutputMatrix(outputTableRow,2) = endHeight / delta_Z;
        velocityOutputMatrix(outputTableRow,3) = velocity;

        %todo: fix output code to csv file
        rowLabels = [rowLabels; char(strcat(char(fileNames(i)), '_', num2str(outputTableRow)))]; %#ok<AGROW>
        startHeight = endHeight;
        outputTableRow = outputTableRow + 1;
    end
end
hold off

T = array2table(velocityOutputMatrix,'RowNames',rowLabels,'VariableNames',outputColumnNames);
writetable(T,csvName,'WriteRowNames',true);
