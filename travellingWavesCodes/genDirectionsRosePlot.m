%% plots the roseplot of directions of the travelling waves
% run displayTWSim.m 1st to get output, output would be 1*numTrials or
% numFreqRanges*numTrials based on freqRanges to consider; use one of two
% below variations accordingly

% for M2/M1/G1
outputsM2 = outputs;

% For M2
numTrials = numel(outputsM2(1,:)); % same number of trials in both output structures

% define some parameters for wave detection
wobbleLim1 = 0; %degree
segOption1 = 1;
wobbleLim2 = 5; %degree
segOption2 = 2;
wobbleLim3 = 10; %degree
segOption3 = 3;
lengthLimit = 1;   %10; %ms
boundryLims = [0.25 1.25];
numFrequencyRanges = size(outputsM2,1);

%%
% when outputs has size 1*numTrials (only 1 freq band)

%for 0 deg
%initialize outputs
waveVector = nan(numTrials,length(timeVals));
uniqueDirs = cell(numFrequencyRanges,numTrials);
waveBounds = cell(numFrequencyRanges,numTrials);

% % for ease of working
% outputsM1 = outputsG1;

for i = 1:numTrials
    [waveVector(i,:),uniqueDirs{1,i},waveBounds{1,i}] = getWaveSegments(outputsM1{1,i},timeVals,wobbleLim3,segOption2,boundryLims, lengthLimit);
end

% make so that each element of the uniqueDirs cells are row vectors
if size(uniqueDirs{1, 1}, 2)==1   %if column vector, do this
    for i = 1:numTrials
        uniqueDirs{1, i} = uniqueDirs{1, i}';
    end
end

m1_0deg = [];
for i=1:numTrials
    m1_0deg = [m1_0deg, uniqueDirs{1, i}];
end


%%
% when outputs has size numFreqBands*numTrials (multiple freq bands)

% for 0 deg
%initialize outputs
waveVector = nan(numTrials,length(timeVals),numFrequencyRanges);
uniqueDirs = cell(numFrequencyRanges,numTrials);
waveBounds = cell(numFrequencyRanges,numTrials);

for i = 1:numTrials
    for j = 1:numFrequencyRanges
    [waveVector(i,:,j),uniqueDirs{j,i},waveBounds{j,i}] = getWaveSegments(outputsM1{j,i},timeVals,wobbleLim2,segOption2,boundryLims, lengthLimit);
    end
end

% find overlapping waves 
allUniqueDirs = [];
dirSG = nan(numTrials,length(timeVals));
dirFG = nan(numTrials,length(timeVals));
waveBoundsOv = cell(1,numTrials);
emptyCells = zeros(1,numTrials);
overlap = 0.5;
for i = 1:numTrials
    [waveBoundsOv{1,i},dirSG(i,:),dirFG(i,:),uniqueDirsTemp,emptyCells(i)] = getOverlappingWaves(waveVector(i,:,1),waveBounds{1,i},waveVector(i,:,2),waveBounds{2,i},overlap);   
    allUniqueDirs = cat(2,allUniqueDirs,uniqueDirsTemp);
end
allUniqueDirs(:,isnan(allUniqueDirs(1,:))) = [];
%%

% get all wave angles
waveBoundsOv(emptyCells==1) = [];
ovTrials = find(emptyCells==0);

allSGWaves = [];
allFGWaves = [];
for i = 1:length(waveBoundsOv)
    waveTemp = waveBoundsOv{i};
    for j = 1:size(waveTemp{1},2)
        overlappingPts = intersect(waveTemp{1}(1,j):waveTemp{1}(2,j),waveTemp{2}(1,j):waveTemp{2}(2,j));
        sgWaves = {dirSG(ovTrials(i),overlappingPts)};
        allSGWaves = cat(1,allSGWaves,sgWaves);
        fgWaves = {dirFG(ovTrials(i),overlappingPts)};
        allFGWaves = cat(1,allFGWaves,fgWaves);
    end
end 
allOvWaves = cat(2,allSGWaves,allFGWaves);

m1_5deg = {uniqueDirs,allUniqueDirs,allOvWaves,ovTrials};    
% get circular correlation for M1
[rho1(2), pval1(2)] = circ_corrcc(allUniqueDirs(1,:),allUniqueDirs(2,:));

%% % 
% plotting
figure; 
tiledlayout(1, 2, 'TileSpacing', 'loose');

nexttile;
title("40-60 Hz gamma")
polarhistogram(m1_0deg)

nexttile
% ........


sgtitle("Condition M2, rose plot of directions of travelling waves across all trials for 013AR")

