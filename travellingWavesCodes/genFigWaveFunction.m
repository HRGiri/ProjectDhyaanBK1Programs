%% 
clear
clc

%% load metaData
subjectName = '013AR';
date = '280122';
protocol = 'M2';
folderSource = '/Users/aniketmandal/Documents/MATLAB/NSP_Sem3/grantData/segmentedData';

segmentedDataFolder = fullfile(folderSource,subjectName,'EEG',date,protocol,'segmentedData','LFP');
lfpInfo = load(fullfile(segmentedDataFolder, 'lfpInfo.mat'));
timeVals = lfpInfo.timeVals;

badTrialsStruct = load(fullfile(segmentedDataFolder,'..','badTrials_wo_v8'));
badTrials = badTrialsStruct.badTrials;
badElecs = badTrialsStruct.badElecs.noisyElecs;

montageFolder = '/Users/aniketmandal/Documents/MATLAB/NSP_Sem3/connectivityHelpers/Montages/Layouts/actiCap64_UOL';
labelFilePath = fullfile(montageFolder,'actiCap64_UOLLabels.mat');
montageLabels = load(labelFilePath).montageLabels;
chanlocFilePath = fullfile(montageFolder,'actiCap64_UOL.mat');
chanlocs = load(chanlocFilePath).chanlocs;

%% set up some parameters
minBurstSize = 25; % in ms
wobble = 0; % in degww
thresh = 0.5;
binEdges = 0:0.05:1;
electrodeFraction = 0.5;
electrodeChoice = 'selected';
waveDetectionMethod = 1;

freqRangeList{1} = [20 40]; freqRangeList{2} = [40 60];
gammaBands = {'fast', 'slow'};
numFrequencyRanges = numel(freqRangeList);
[X,Y] = meshgrid(1:1:9);
waveLengthLimit = 10;

thresholdFactor = 3;
% stimulusDurationS = [0 0.8]; % Stimulus duration to be highlighted
baselinePeriodS = [-1 0];
stimulusPeriodS = [0.25 1.25];
% analysisPeriodS = [-0.5 1];
filterOrder = 4;
stimPeriod = [0.25 1.25];

%% Load electrode data
numTrials = length(lfpInfo.goodStimPos);
numElectrodes = length(lfpInfo.electrodesStored);
numTimePoints = length(timeVals);
allData = zeros(numTrials,numElectrodes,numTimePoints);

for i=1:numElectrodes
    allData(:,i,:) = load(fullfile(segmentedDataFolder,['elec' num2str(i) '.mat'])).analogData;
end

% Only select the Good EEG electrodes and good trials
goodElecs = 1:64;%setdiff(1:64,badElecs);
goodTrials = setdiff(1:numTrials,badTrials);
allData = allData(:,1:64,:);  % 64 EEG electrodes
allData = allData(goodTrials,goodElecs,:);

%% load output data
numTrials = size(allData, 1);
TWDataFolder = fullfile('/Users/aniketmandal/Documents/MATLAB/NSP_Sem3/grantData', 'TWData',subjectName);

outputs = cell(numFrequencyRanges, numTrials);
for gammaBandIndex=1:numFrequencyRanges
    subjectName = '013AR';
    TWDataFile = fullfile(TWDataFolder, sprintf('%s_%sGamma.mat', protocol, gammaBands{gammaBandIndex}));
    temp = load(TWDataFile);
    outputTemp = temp.outputs;

    outputs(gammaBandIndex, :) = outputTemp;
end

%%
numGoodElectrodes = length(goodElecs);
numTrials = size(allData,1);
burstTS = nan(numGoodElectrodes,numTrials,length(timeVals),numFrequencyRanges);
req = 1;
segOption = 1;

 for iFreq=1:numFrequencyRanges
    for iElec=1:numGoodElectrodes
        [~,~,~,burstTS(iElec,:,:,iFreq),~,~] = getHilbertBurst(squeeze(allData(:,iElec,:)),timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,...
            freqRangeList{iFreq},filterOrder,req);
    end
 end
 
[slowGammaOverlap1,fastGammaOverlap1,burstTS1] = getWaveAndBurstOverlap(burstTS,outputs,timeVals,minBurstSize,waveLengthLimit,wobble,binEdges,goodElecs,thresh,segOption);

 %% Plotting
colorVals = cat(1,[52 148 186]./255,[236 112 22]./255);

burstData = sum(slowGammaOverlap1);
figure;
plot(binEdges(2:end),burstData/max(burstData),'-o','LineWidth',1.2,'Color',colorVals(1,:))
hold on
burstData = sum(fastGammaOverlap1);
plot(binEdges(2:end),burstData/max(burstData),'-o','LineWidth',1.2,'Color',colorVals(2,:))
% ylabel(sprintf('Protocol: %s', protocol))
% xlim([0.1 1])
title('TW distribution along \gamma bursts')

% figure;
% burstData = sum(slowGammaOverlap1);
% burstData = burstData./max(burstData,[],2);
% burstDataErr = std(burstData);
% plot(binEdges(2:end),mean(burstData),'-o','LineWidth',1.2,'Color',colorVals(1,:))
% hold on
% errorbar(binEdges(2:end),mean(burstData),burstDataErr,'LineWidth',1.2,'Color',colorVals(1,:))
% 
% burstData = sum(slowGammaOverlap1);
% burstData = burstData./max(burstData,[],2);
% burstDataErr = std(burstData);
% plot(binEdges(2:end),mean(burstData),'-o','LineWidth',1.2,'Color',colorVals(2,:))
% hold on
% errorbar(binEdges(2:end),mean(burstData),burstDataErr,'LineWidth',1.2,'Color',colorVals(2,:))
% xlim([0.1 1])
% xlabel('Gamma Burst Bins')
% ylabel('All Ori')
% title('TW distribution along \gamma bursts-all ori')



