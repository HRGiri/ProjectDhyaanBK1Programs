clear;
%% Load Meta Data
subjectName = '013AR';
date = '280122';
protocol = 'M2';
gammaBand = "fast";
folderSource = 'E:\NeoLabData\projectDhyaan\segmentedData';

segmentedDataFolder = fullfile(folderSource,subjectName,'EEG',date,protocol,'segmentedData','LFP');
lfpInfo = load(fullfile(segmentedDataFolder, 'lfpInfo.mat'));
timeVals = lfpInfo.timeVals;

badTrialsStruct = load(fullfile(segmentedDataFolder,'..','badTrials_wo_v8'));
badTrials = badTrialsStruct.badTrials;
badElecs = badTrialsStruct.badElecs.noisyElecs;

montageFolder = 'E:\OneDrive - Indian Institute of Science\NeoLabData\Programs\Montages\Layouts\actiCap64_UOL';
labelFilePath = fullfile(montageFolder,'actiCap64_UOLLabels.mat');
montageLabels = load(labelFilePath).montageLabels;
chanlocFilePath = fullfile(montageFolder,'actiCap64_UOL.mat');
chanlocs = load(chanlocFilePath).chanlocs;

%% Load electrode data
numTrials = length(lfpInfo.goodStimPos);
numElectrodes = length(lfpInfo.electrodesStored);
numTimePoints = length(timeVals);
data = zeros(numTrials,numElectrodes,numTimePoints);

for i=1:numElectrodes
    data(:,i,:) = load(fullfile(segmentedDataFolder,['elec' num2str(i) '.mat'])).analogData;
end

% Only select the Good EEG electrodes and good trials
goodElecs = 1:64;%setdiff(1:64,badElecs);
goodTrials = setdiff(1:numTrials,badTrials);
data = data(:,1:64,:);  % 64 EEG electrodes
data = data(goodTrials,goodElecs,:);

%% Select Electrodes to analyse
% occipitalElectrodeLabels = {"O1","Oz","O2"};
% parietoOccipitalElectrodeLabels = {"PO7","PO3","POz","PO4","PO8"};
% parietalElectrodeLabels = {"P7","P5","P3","P1","Pz","P2","P4","P6","P8"};
% centroParietalElectrodeLabels = {"CP5","CP3","CP1","CPz","CP2","CP4","CP6"};
% centralElectrodeLabels = {"C5","C3","C1","Cz","C2","C4","C6"};
% frontoCentralElectrodeLabels = {"FC5","FC3","FC1","FC2","FC4","FC6"};
% electrodesToAnalyseLabels = [occipitalElectrodeLabels,...
%  parietoOccipitalElectrodeLabels,...
%   parietalElectrodeLabels,...
%    ...centroParietalElectrodeLabels,...
%     ...centralElectrodeLabels,...
%     ...frontoCentralElectrodeLabels...
%     ];
% labels = montageLabels(:,2);
% electrodesToAnalyse = zeros(1,length(electrodesToAnalyseLabels));
% for i=1:length(electrodesToAnalyse)
%     electrodesToAnalyse(i) = find(labels==electrodesToAnalyseLabels{i});
% end

% electrodesToAnalyse = 1:64; % Uncomment to use the whole electrode grid
% Find clusters of electrodes
X = projectTo2DPlane(chanlocs);
% X = [chanlocs.Y; chanlocs.X; chanlocs.Z]';
% [clusters, clusterPeakFreqs] = getClusters(data,timeVals,X,[0 80]);
clusters{1} = 1:64;
if strcmp(gammaBand,'slow')
    freqRange = [20 40];    %#ok<UNRCH>
else
    freqRange = [40 60]; %#ok<UNRCH>
end
outputs = getTWCircParams(squeeze(data(1,:,:)),timeVals,...
        clusters{1},freqRange,1,[],chanlocs);

%% Plot Simulation
% Coarse simulation
% simulationPeriod = [0.25 1.25];
% simulationSpeed = 0.01;

% Fine Simulation
simulationPeriod = [0.345 0.35];
simulationSpeed = 1;

% Single Frame snapshot
% timePoint = 0.718;
% simulationPeriod = [timePoint timePoint];
% simulationSpeed = 1; % Placeholder
clusterIndex = 1;
runSimulation(outputs,chanlocs,clusters{clusterIndex},timeVals,simulationPeriod,simulationSpeed,0);    % Comment to just get the outputs
pgd = outputs.pgd;
pgd(abs(pgd)==inf)=0;
pgd(isnan(pgd))=0;
trapz(timeVals,pgd)
% figure;
% plot(timeVals,pgd);
% xlim([0.2 1.25])