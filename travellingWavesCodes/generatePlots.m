%% Parameters
% run displayTWSim.m 1st to get output, output would be 1*numTrials or
% numFreqRanges*numTrials based on freqRanges to consider; use one of two
% below variations accordingly
folderSource = 'E:\NeoLabData\projectDhyaan';
subjectName = '013AR';
tWDataFolder = fullfile(folderSource,"TWData",subjectName);

protocols = {"G1","G2","M2"};
gammaBands = {"slow","fast"};

% define some parameters for wave detection
wobbleLim = 0; %degree
segOption = 1;
% wobbleLim = 5; %degree
% segOption = 2;
% wobbleLim = 10; %degree
% segOption = 3;
lengthLimit = 10;   %10; %ms
boundryLims = [0.25 1.25];

timeVals = -1.249:0.001:1.25;

%% Load Data
% protocol = "G1";
for protocolIndex = 1:length(protocols)
    protocol = protocols{protocolIndex};
outputs = {};
for gammaBandIndex = 1:length(gammaBands)
    gammaBand = gammaBands{gammaBandIndex};
    tWDataFile = fullfile(tWDataFolder,sprintf("%s_%sGamma.mat",protocol,gammaBand));
    out = load(tWDataFile).outputs;    
    outputs = [outputs; out];
end
outputAllProtocols.(protocol) = outputs;
fprintf("Loaded data for %s\n",protocol);
end

%% Plot
% waveProperty = 'Wave Propagation Speed';
% waveProperty = 'Wavelength';
waveProperty = 'Wave Propagation Direction';
% waveProperty = 'Electrode Cluster Sizes';
% waveProperty = 'Wave Strength AUC';

% Get properties
for protocolIndex = 1:length(protocols)
    protocol = protocols{protocolIndex};
    outputs = outputAllProtocols.(protocol);
    
    numTrials = numel(outputs(1,:)); % same number of trials in both output structures
    numFrequencyRanges = size(outputs,1);
%initialize outputs
    waveVector = nan(numTrials,length(timeVals),numFrequencyRanges);
    dirs = cell(numFrequencyRanges,numTrials);
    waveLengths = cell(numFrequencyRanges,numTrials);
    speeds = cell(numFrequencyRanges,numTrials);
    clusterSizes = cell(numFrequencyRanges,numTrials);
    auc = cell(numFrequencyRanges,numTrials);
    waveBounds = cell(numFrequencyRanges,numTrials);

    for i = 1:numTrials
        for j = 1:numFrequencyRanges            
            [waveVector(i,:,j),dirs{j,i},waveBounds{j,i},waveLengths{j,i},speeds{j,i},clusterSizes{j,i}] = getWaveSegments(outputs{j,i},timeVals,wobbleLim,segOption,boundryLims, lengthLimit);
            % AUC
            pgd = outputs{j,i}.pgd;
            pgd(isnan(pgd) | (abs(pgd)==inf)) = 0;
            auc{j,i} = trapz(timeVals,pgd);
        end
    end
    
    if strcmp(waveProperty,'Wavelength')
        uniqueVals = waveLengths;
        xAxisLabel = 'Wavelength (mm)';
    elseif strcmp(waveProperty,'Wave Propagation Direction')
        uniqueVals = dirs;        
    elseif strcmp(waveProperty,'Wave Propagation Speed')
        uniqueVals = speeds;
        xAxisLabel = 'Speed (m/s)';
    elseif strcmp(waveProperty,'Electrode Cluster Sizes')
        uniqueVals = clusterSizes;
        xAxisLabel = '# Electrodes in Cluster';
    elseif strcmp(waveProperty,'Wave Strength AUC')
        uniqueVals = auc;
        xAxisLabel = 'Area Under Curve (Wave Strength)';
    end
    % make so that each element of the uniqueDirs cells are row vectors
    
    for i = 1:numTrials
        for j = 1:numFrequencyRanges
            if size(uniqueVals{j, i}, 2)==1   %if column vector, do this
                uniqueVals{j, i} = uniqueVals{j, i}';
            end
        end
    end

    values = cell(1,numFrequencyRanges);
    for j=1:numFrequencyRanges
        values{j} = [];
    for i=1:numTrials        
        values{j} = [values{j}, uniqueVals{j, i}];
    end
    end
    allValues.(protocol) = values;
end

allData = [];
for protocolIndex=1:length(protocols)
    protocol = protocols{protocolIndex};
    values = allValues.(protocol);    
    allData = [allData values{1,1} values{1,2}]; %#ok<*AGROW>
end
edges = linspace(min(allData), max(allData), 21);

% plotting
figure; 
tiledlayout(1, length(protocols), 'TileSpacing', 'loose');
gammaBandLabels = {"Slow Gamma", "Fast Gamma"};
for protocolIndex=1:length(protocols)
    protocol = protocols{protocolIndex};
    values = allValues.(protocol);
    nexttile; 
    if strcmp(waveProperty,'Wave Propagation Direction')
        polarhistogram(values{1},20,'Normalization','probability') %#ok<*UNRCH>
        hold on;
        polarhistogram(values{2},20,'Normalization','probability')        
    else
        histogram(values{1},edges,'Normalization','probability');        
        hold on;
        histogram(values{2},edges,'Normalization','probability');        
        
        xlabel(xAxisLabel);
        xlim([edges(1) edges(end)]);
        legend({gammaBandLabels{1},gammaBandLabels{2}})
    end
    title(protocol)
end
if ~strcmp(waveProperty,'Wave Propagation Direction')
    ax = findall(gcf,'type','axes');
    linkaxes(ax, 'xy');
else
    legend({gammaBandLabels{1},gammaBandLabels{2}})
end
set(findall(gcf,'-property','FontSize'), 'FontSize', 14);
% set(findall(gcf, 'Type', 'text', 'Tag', 'Title'), 'FontSize',18);
% set(findall(gcf, 'Type', 'text', 'Tag', 'XLabel'), 'FontSize',16);
sgtitle(waveProperty,'FontWeight','bold','FontSize',20)