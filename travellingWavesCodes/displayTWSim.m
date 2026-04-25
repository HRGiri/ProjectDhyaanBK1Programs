clear;
%% Load Meta Data
subjectName = '013AR';
date = '280122';
protocol = 'G1';
folderSource = 'E:\NeoLabData\projectDhyaan\segmentedData';

segmentedDataFolder = fullfile(folderSource,subjectName,'EEG',date,protocol,'segmentedData','LFP');
lfpInfo = load(fullfile(segmentedDataFolder, 'lfpInfo.mat'));
timeVals = lfpInfo.timeVals;
load('E:\OneDrive - Indian Institute of Science\Coursework\Neural Signal Processing\Assignments\Assignment2\actiCap64_UOL.mat')

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

%% Select Electrodes to analyse
occipitalElectrodeLabels = {"O1","Oz","O2"};
parietoOccipitalElectrodeLabels = {"PO7","PO3","POz","PO4","PO8"};
parietalElectrodeLabels = {"P7","P5","P3","P1","Pz","P2","P4","P6","P8"};
centroParietalElectrodeLabels = {"CP5","CP3","CP1","CPz","CP2","CP4","CP6"};
centralElectrodeLabels = {"C5","C3","C1","Cz","C2","C4","C6"};
frontoCentralElectrodeLabels = {"FC5","FC3","FC1","FC2","FC4","FC6"};
electrodesToAnalyseLabels = [occipitalElectrodeLabels,...
 parietoOccipitalElectrodeLabels,...
  parietalElectrodeLabels,...
   ...centroParietalElectrodeLabels,...
    ...centralElectrodeLabels,...
    ...frontoCentralElectrodeLabels...
    ];
labels = montageLabels(:,2);
electrodesToAnalyse = zeros(1,length(electrodesToAnalyseLabels));
for i=1:length(electrodesToAnalyse)
    electrodesToAnalyse(i) = find(labels==electrodesToAnalyseLabels{i});
end

% electrodesToAnalyse = 1:64; % Uncomment to use the whole electrode grid
[outputs] = getTWCircParams(squeeze(data(1,1:64,:)),timeVals,electrodesToAnalyse,[40 60],1,[]);

%% Plot Simulation
% Coarse simulation
% simulationPeriod = [0.5 0.75];
% simulationSpeed = 0.01;

% Fine Simulation
simulationPeriod = [0.716 0.724];
simulationSpeed = 1;

% Single Frame snapshot
% timePoint = 0.718;
% simulationPeriod = [timePoint timePoint];
% simulationSpeed = 1; % Placeholder

runSimulation(outputs,chanlocs,electrodesToAnalyse,timeVals,simulationPeriod,simulationSpeed,0);    % Comment to just get the outputs