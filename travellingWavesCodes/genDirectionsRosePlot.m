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
    [waveVector(i,:),uniqueDirs{1,i},waveBounds{1,i}] = getWaveSegments(outputsM1{1,i},timeVals,wobbleLim3,segOption3,boundryLims, lengthLimit);
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

figure; title("Condition M2, rose plot of directions of travelling waves across all trials for 013AR")
polarhistogram(m1_0deg)