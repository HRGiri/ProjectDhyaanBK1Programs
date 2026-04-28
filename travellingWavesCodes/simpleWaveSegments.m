% taken from "https://github.com/gouthambhargava/TravellingWaveProject/blob/main/commonAnalysisCodes/analysisCodes/simpleWaveSegments.m"

function [waveVector,boundries,allUniqueDirs] = simpleWaveSegments(data,lengthLimit)
    if nargin<2
        lengthLimit = 0;
    end
    % preservedData = data;
    preservedData = reshape(data,[length(data),1]);%convert to column vector 
    waveVector = nan(1,length(data));
    data = preservedData;
    data(isnan(data)) = 0;
    data(data~=0) = 1;
    data = find(data==0);
    dataEpochs = find(diff(data)>1);
    boundries = [data(dataEpochs)+1,data(dataEpochs+1)-1]';
    boundries(:,diff(boundries)<lengthLimit) = []; 
    allUniqueDirs = [];
    for k = 1:size(boundries,2)
        waveVector(boundries(1,k):boundries(2,k)) = preservedData(boundries(1,k):boundries(2,k));
        uniqueDirs = circMeanNan(preservedData(boundries(1,k):boundries(2,k)));
        allUniqueDirs = cat(1,allUniqueDirs,uniqueDirs);
    end
end