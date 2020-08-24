function [rtnCulledConcentrationVector] = CullSettlingData(data, slopeThreshhold)
%CULLSETTLINGDATA Use moving averages to cull data not part of meaningful
%settling regime
%Written by: ssz
%Last edited: June 11, 2019

%Set sub analysis size to analyze change in slopes
subAnalysisSize = 10;
startPoint = 1;

allSlopes = zeros(size(data,1) - 1, 2);
endSize = size(data,1) - 1;

%Split into 10 data point segments and analyze change in slope over this
%period of time
for i=1:endSize
    if((endSize - subAnalysisSize) <= i && i <= endSize)
        subMatrix = data(startPoint:size(data,1),:);
    else
        subMatrix = data(startPoint:startPoint+subAnalysisSize, :);
    end
    
    [~, slope] = CalculateLeastSquares(subMatrix(:,1), subMatrix(:,2));
    allSlopes(i, 1) = i;
    allSlopes(i, 2) = -1 * slope; %convert to positive
    startPoint = startPoint + 1;
end

[index,~,~] = find(allSlopes(:,2)>slopeThreshhold);

if(isempty(index))
     rtnCulledConcentrationVector = data(:, :);    
else
    startPoint = index(1);
    endPoint= index(end);

%if the start and end point are close to each other or there isn't a lot of 
%slope happening, then just return the whole thing
if((abs(endPoint - startPoint) < 10))
    %return all rows and columns
    rtnCulledConcentrationVector = data(:,:); 
else
    rtnCulledConcentrationVector = data(startPoint:endPoint, :); 
end

end

