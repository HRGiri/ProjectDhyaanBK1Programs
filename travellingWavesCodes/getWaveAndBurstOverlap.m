% modified from 'https://github.com/gouthambhargava/TravellingWaveProject/blob/main/commonAnalysisCodes/analysisCodes/getWaveAndBurstOverlap.m'

function [slowGammaOverlap,fastGammaOverlap,burstTS] = getWaveAndBurstOverlap(burstTS,outputs,timeVals,burstLengthLimit,waveLengthLimit,waveWobble,binEdges,goodElectrodes,thresh,segOption)
% function only for generating figure 5. Identifies overlapping gamma
% bursts and subsequently the variation of wave locations within these
% overlapping bursts
% Inputs
% burstTS - burst matrix for all trials
% outputs - cell with output structures of both SG and FG
% timeVals
% burstLengthLimit - (in ms) bursts under this limit will be discarded
% waveLengthLimit - (in ms) bursts under this limit will be discarded
% waveWobble - angle (in deg) that is considered for wave segmentation
% binEdges - for plotting

if nargin<10
    segOption = 1;    %3;
end
% generate the binned data
% some parameters
% minBurstSize = 100;
lims = [0.25 1.25];
srate = 1/(timeVals(2)-timeVals(1));
% burstLimit = 100;
% wobbleLim = 0;
burstLimit = burstLengthLimit*srate*10^-3; % set wavelimit in samples
% lengthLimit = 100;
lengthLimit = waveLengthLimit*srate*10^-3; % set wavelimit in samples

burstTS(isnan(burstTS)) = 0;
elecFrac = numel(goodElectrodes)*thresh;
burstFrac = squeeze(sum(burstTS)); % sum across all electrodes for all trials to get a trails x time x frequencies matrix
burstFrac(burstFrac<elecFrac) = 0;
burstFrac(burstFrac==0) = nan;

[allBoundsSG,~] = getSegmentedBursts(burstFrac(:,:,1),lims,timeVals,burstLimit);
[allBoundsFG,~] = getSegmentedBursts(burstFrac(:,:,2),lims,timeVals,burstLimit);

waveVector = zeros(size(burstFrac));
for i = 1:numel(outputs(1,:))
    waveVector(i,:,1) = getWaveSegments(outputs{1,i},timeVals,waveWobble,segOption,lims, lengthLimit);
    waveVector(i,:,2) = getWaveSegments(outputs{2,i},timeVals,waveWobble,segOption,lims, lengthLimit);
end

waveVector(~isnan(waveVector)) = 1;
waveVector(isnan(waveVector)) = 0;


slowGammaOverlap = [];
fastGammaOverlap = [];

% get dir values corrosponding to the bursts
for i = 1:numel(allBoundsSG)
    for j = 1:size(allBoundsSG{i},2)
        int = waveVector(i,allBoundsSG{i}(1,j):allBoundsSG{i}(2,j),1);
        int1 = allBoundsSG{i}(1,j):allBoundsSG{i}(2,j);
        int1 = (int1-min(int1))/(max(int1)-min(int1));
        int1(int==0) = [];
        y = histcounts(int1,binEdges);
        slowGammaOverlap = cat(1,slowGammaOverlap,y);
    end
    for j = 1:size(allBoundsFG{i},2)
        int = waveVector(i,allBoundsFG{i}(1,j):allBoundsFG{i}(2,j),2);
        int1 = allBoundsFG{i}(1,j):allBoundsFG{i}(2,j);
        int1 = (int1-min(int1))/(max(int1)-min(int1));
        int1(int==0) = [];
        y = histcounts(int1,binEdges);
        fastGammaOverlap = cat(1,fastGammaOverlap,y);
    end
end
clear int int1 i j y
end


% additional functions
function [allBounds,allBursts] = getSegmentedBursts(bursts,Lims,timeVals,burstLimit)
% find the boundries of the bursts and remove the bursts of length less
% than burstLimit

boundryLims = [dsearchn(timeVals',Lims(1)),dsearchn(timeVals',Lims(2))];
bursts(:,setdiff(1:length(bursts),boundryLims(1):boundryLims(2))) = nan;

% find the burst indices for all trials and plot
allBounds = cell(1,size(bursts,1));
allBursts = nan(size(bursts));
for i = 1:size(bursts,1)
    boundries = [];
    burstInt = zeros(1,length(timeVals));
    burstInt(~isnan(bursts(i,:))) = 1;
    burstData = find(burstInt==0);
    burstEpochs = find(diff(burstData)>1);
    boundries = cat(2,boundries,[burstData(burstEpochs)+1;burstData(burstEpochs+1)-1]);
    boundries(:,diff(boundries)<burstLimit) = [];
    for k = 1:size(boundries,2)
        allBursts(i,boundries(1,k):boundries(2,k)) = 1;
    end
    allBounds{i} = boundries;
end
end
