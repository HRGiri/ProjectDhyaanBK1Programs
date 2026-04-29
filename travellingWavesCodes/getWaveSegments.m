% taken from "https://github.com/gouthambhargava/TravellingWaveProject/blob/main/commonAnalysisCodes/analysisCodes/getWaveSegments.m"

function [waveVector,uniqueDirs,waveBounds,uniqueWavelengths,uniqueSpeeds,uniqueClusterSizes] = getWaveSegments(outputsTW,timeVals,wobbleLim,segOption,boundryLims, lengthLimit)
% set boundryLims as a vector of times outside of which direction vectors will not be considered
% Inputs: -outputsTW - output from getTWCircParams which contains the PGD, direction, speed, tempFreq and clusters for a single trial
%         -timeVals - vector of time values for the given trial
%         -wobbleLim - sets a threshold of an angle(in deg). Consecutive
%         values in direction under this threshold are considered part of a
%         single wave. Leave this option empty if segmentation according to
%         direction values is not required
%         -segOption - 1 - waves are considered whererever significant PGD values (>0) are observed. Waves <25ms are ignored.
%                      2 - Based on the 'wave wobble', successive time points with significant PGD(>0) are said to belong to the
%                      same wave only the wave direction varies by less than the 'wave wobble'.
%                      3 - A time segment is considered a wave if the variation between any two time points within
%                       the wave is less than the 'wave wobble'.
%                      4 - Determines the presence or absence of a wave at a given time point based on both the PGD and direction
%                          Details specified in Das, A., Zabeh, E., Jacobs, J. (2023). How to Detect and Analyze Traveling Waves in Human Intracranial EEG Oscillations?. In: Axmacher, N. (eds) Intracranial EEG. Studies in Neuroscience, Psychology and Behavioral Economics. Springer, Cham. https://doi.org/10.1007/978-3-031-20910-9_30
%         -boundryLims - give time limits, ex [0.25 0.75]. Waves will be identified only within this limit
%         -lengthLimit - set a length value (in ms). Only waves above this
%         limit be considered.
% Outputs
%          waveVector - vector of time indices which contain waves
%          uniqueDirs - mean direction values for each wave segment
%          waveBounds - start and end time indices for each wave in waveVector

% define some parameters
boundryLims = dsearchn(timeVals',boundryLims(1)):dsearchn(timeVals',boundryLims(2)); %indices of the specified time limit
Fs = 1/(timeVals(2)-timeVals(1));
lengthLimit = lengthLimit*Fs*10^-3;
elecNum = 64;  
if segOption == 4 % wave strength detection requires both pgd and direction
    pgd = outputsTW.pgd;
    direction = outputsTW.direction;
    if min(size(pgd))==1 % if this method is to be applied to TW estimations from global estimations of PGD, PGD and direction are repeated over all electrodes.
        pgd = repmat(pgd',elecNum,1);
        direction = repmat(direction',elecNum,1);
    end
    direction(:,setdiff(1:length(timeVals),boundryLims)) = nan;
    pgd(:,setdiff(1:length(timeVals),boundryLims)) = nan;
else
    % pgd = mean(outputsTW.pgd,'omitnan');
    % direction = mean(outputsTW.direction,'omitnan');
    direction = outputsTW.direction;
    waveLength = outputsTW.Wavelength;
    speed = outputsTW.speed;
    speed = [NaN;speed];
    clusters = outputsTW.clusters;
    
    pgd = outputsTW.pgd;
    pgdThreshold = 0.5;
    
    direction(isnan(pgd) | (pgd<pgdThreshold)) = nan;
    
    waveLength(isnan(pgd) | (pgd<pgdThreshold)) = nan;
    waveLength(abs(waveLength)==inf) = nan;  
    waveLength = waveLength*360/rad2deg(1); % Convert to mm/deg

    speed(isnan(pgd) | (pgd<pgdThreshold)) = nan;       
    speed(abs(speed)==inf) = nan;
    speed = speed/1000;   % Convert to m/s
    
    clusters(isnan(pgd) & (pgd<pgdThreshold)) = [];
    clusterSizes = double(size(clusters));
    for i = 1:length(pgd)
        clusterSizes(i) = length(clusters{i});
    end


    % if size(direction,1)>1
    %     direction = direction(1,:);
    % end
    % direction(setdiff(1:length(timeVals),boundryLims)) = nan;
end
%% diverge into type of segmentation required
if segOption==1 % get wave segments where segments are purely seperated by significant PGD. Segmentation is very lenient
    [waveVector,waveBounds] = getsimpleWaveSegments(direction,lengthLimit);
    uniqueDirs = [];
    uniqueWavelengths = [];
    uniqueSpeeds = [];
    uniqueClusterSizes = [];
    for j = 1:size(waveBounds,2)
        uniqueDirs = cat(2,uniqueDirs,circ_mean(direction(waveBounds(1,j):waveBounds(2,j))));
        uniqueWavelengths = cat(2,uniqueWavelengths,nanmean(waveLength(waveBounds(1,j):waveBounds(2,j))));
        uniqueSpeeds = cat(2,uniqueSpeeds,nanmean(speed(waveBounds(1,j):waveBounds(2,j))));
        uniqueClusterSizes = cat(2,uniqueClusterSizes,mean(clusterSizes(waveBounds(1,j):waveBounds(2,j))));
    end

elseif segOption==2 % Imposes an additional criteria on wave segmentation, based on wobbleLim.
    if min(size(direction))>1
        direction = direction(1,:);
    end
    % waves within this earlier segmentation are further degraded based on
    % the angular difference between successive time points. Each time
    % point can vary only within the wobbleLim. Variation beyond this limit
    % marks the start of a new wave.
    % [newDirection,boundries] = getsimpleWaveSegments(direction,lengthLimit);
    % if ~isempty(boundries)
    %     for k = 1:size(boundries,2)
    %         tempBounds = newDirection(boundries(1,k):boundries(2,k));
    %         tempBoundsInd = boundries(1,k):boundries(2,k);
    %         tempDiff = abs(diff(tempBounds));
    %         revisedBounds = [];
    %         if max(tempDiff)<deg2rad(wobbleLim)
    %             revisedBounds = cat(2,revisedBounds,[boundries(1,k);boundries(2,k)]);
    %         else
    %             stpPts = find(tempDiff>deg2rad(wobbleLim))+1;
    %             stpPts = [1 stpPts;stpPts-1 numel(tempBounds)];
    %             revisedBounds = cat(2,revisedBounds,[tempBoundsInd(stpPts(1,:));tempBoundsInd(stpPts(2,:))]);
    %         end
    %         revisedBounds(:,diff(revisedBounds)<lengthLimit) = [];
    %     end
    % else
    %     revisedBounds = [];
    % end
    % waveVector = nan(1,length(direction));
    % uniqueDirs = [];
    % if ~isempty(revisedBounds)
    %     for j = 1:size(revisedBounds,2)
    %         waveVector(revisedBounds(1,j):revisedBounds(2,j)) = direction(revisedBounds(1,j):revisedBounds(2,j));
    %         uniqueDirs = cat(2,uniqueDirs,circ_mean(direction(revisedBounds(1,j):revisedBounds(2,j))'));
    %     end
    % end
    % waveBounds = revisedBounds;
    directionNew = wrapTo2Pi(reshape(direction,[1,length(direction)]));
    direcDiff = diff([directionNew,0]);
    sigPts = find(abs(direcDiff)<=deg2rad(wobbleLim));
    directionNew(setdiff(1:length(directionNew),sigPts)) = nan;
    % waves within this earlier segmentation are further degraded based on
    % the angular difference between successive time points. Each time
    % point can vary only within the wobbleLim. Variation beyond this limit
    % marks the start of a new wave.

    [waveVector,waveBounds] = getsimpleWaveSegments(directionNew,lengthLimit);

    % waveVector = nan(1,length(direction));
    uniqueDirs = [];
    if ~isempty(waveBounds)
        for j = 1:size(waveBounds,2)
            % waveVector(revisedBounds(1,j):revisedBounds(2,j)) = direction(revisedBounds(1,j):revisedBounds(2,j));
            uniqueDirs = cat(1,uniqueDirs,circ_mean(direction(waveBounds(1,j):waveBounds(2,j))));
        end
    end

elseif segOption==3 % ensures that wave segments from option 1 are further fragmented so that the total amount of variation between
    % the start and end of a wave is within the wobbleLim
    [directionTemp,boundries] = getsimpleWaveSegments(direction,lengthLimit);
    directionTemp = wrapTo2Pi(directionTemp);
    revisedBounds = [];
    for k = 1:size(boundries,2)
        tempBounds = directionTemp(boundries(1,k):boundries(2,k));
        boundaryIndices = boundries(1,k):boundries(2,k);
        waveIndices = getSegmentedWaves(tempBounds,boundaryIndices,wobbleLim);
        revisedBounds = cat(2,revisedBounds,waveIndices);
    end
    revisedBounds(:,diff(revisedBounds)<lengthLimit) = [];

    waveVector = nan(1,length(direction));
    uniqueDirs = [];
    waveBounds = [];
    if ~isempty(revisedBounds)
        for j = 1:size(revisedBounds,2)
            waveVector(revisedBounds(1,j):revisedBounds(2,j)) = direction(revisedBounds(1,j):revisedBounds(2,j));
            uniqueDirs = cat(2,uniqueDirs,circ_mean(direction(revisedBounds(1,j):revisedBounds(2,j))));
        end
        waveBounds = revisedBounds;
    end
else % the pgd and directions are used to calculate wave strength. Based on Das, A., Zabeh, E., Jacobs, J. (2023). How to Detect and Analyze Traveling Waves in Human Intracranial EEG Oscillations?. In: Axmacher, N. (eds) Intracranial EEG. Studies in Neuroscience, Psychology and Behavioral Economics. Springer, Cham. https://doi.org/10.1007/978-3-031-20910-9_30

    numElecs = numel(find(nansum(pgd,2))); %#ok<NANSUM>
    stability = nan(1,length(pgd));
    uniqueDirs = [];
    for i = 1:size(pgd,1)
        for j = 1:size(pgd,2)-1
            waveStab = abs(pgd(i,j+1)*exp(sqrt(-1)*direction(i,j+1)) - pgd(i,j)*exp(sqrt(-1)*direction(i,j)));
            stability(j+1) = -sum(waveStab)/numElecs;
        end
    end
    stability = zscoreNan(stability);
    stability(stability<0) = nan;
    % stability(stability>0) = 1;
    
    [waveVector,waveBounds] = getsimpleWaveSegments(stability,lengthLimit); 

end
end

%% additional functions
function [waveVector,boundries] = getsimpleWaveSegments(dataOrig,lengthLimit)
if nargin<2
    lengthLimit = 0;
end
dataOrig = reshape(dataOrig,[length(dataOrig),1]);%convert to column vector
waveVector = nan(1,length(dataOrig));
data = dataOrig;
data(isnan(data)) = 0;
data(data~=0) = 1;
data = find(data==0);
dataEpochs = find(diff(data)>1);
boundries = [data(dataEpochs)+1,data(dataEpochs+1)-1]';
boundries(:,diff(boundries)<lengthLimit) = [];
for k = 1:size(boundries,2)
    waveVector(boundries(1,k):boundries(2,k)) = dataOrig(boundries(1,k):boundries(2,k));
end
end

function [waveIndices] = getSegmentedWaves(data,boundaryIndices,devAngle)
data = reshape(data,[1,length(data)]);
dataFull = [nan data nan];
uniqueAngles = unique(dataFull);
uniqueAngles(isnan(uniqueAngles)) = [];

segWaveBounds = [];
inMat = [];
for i = 1:numel(uniqueAngles)
    data = dataFull;
    initAngle = uniqueAngles(i);
    deviation = deg2rad(devAngle);
    dataDiff = abs(data-initAngle);
    data(abs(dataDiff)>deviation) = 0;
    data(data==0) = nan;

    [~,boundsInt] = getsimpleWaveSegments(data);
    segWaveBounds = cat(2,segWaveBounds,boundsInt);
end

segWaveBounds = segWaveBounds-1;
[~,b] = sort(diff(segWaveBounds),'descend');
sortedBounds = segWaveBounds(:,b);
durations = diff(sortedBounds);

% finalBounds = zeros(length(sortedBounds),length(sortedBounds)-1);
if size(sortedBounds,2)>1
    for i = 1:length(sortedBounds)
        wave1 = sortedBounds(1,i):sortedBounds(2,i);
        for j = 1:length(sortedBounds)
            wave2 = sortedBounds(1,j):sortedBounds(2,j);
            inMat(i,j) = numel(intersect(wave1,wave2));
        end
    end
    inMat(find(diag(diag(inMat)))) = 1;
    inMat(inMat==0) = 0;

    for i = 1:length(sortedBounds)
        selection = inMat(i,:);
        overlapWaves = find(selection>0);
        if ~isempty(overlapWaves)
            for j = 1:numel(overlapWaves)
                if durations(i)== durations(overlapWaves(j))
                    inMat(i,overlapWaves(j)) = 1;
                    inMat(overlapWaves(j),i) = 1;
                elseif durations(i)>durations(overlapWaves(j))
                    inMat(overlapWaves(j),:) = nan;
                    inMat(:,overlapWaves(j)) = nan;
                end
            end
        end
    end
    nonOverlappingWaves = sortedBounds(:,find(diag(inMat)==1));
else
    nonOverlappingWaves = sortedBounds;
end
waveIndices = zeros(2,size(nonOverlappingWaves,2));
for i = 1:size(nonOverlappingWaves,2)
    waveIndices(:,i) = [boundaryIndices(nonOverlappingWaves(1,i));boundaryIndices(nonOverlappingWaves(2,i))];
end
end

function [nanScored] = zscoreNan(data) % zscoring when data has nans
meanData = mean(data(~isnan(data)));
stdData = std(data(~isnan(data)));
nanScored = (data-meanData)/stdData;
end