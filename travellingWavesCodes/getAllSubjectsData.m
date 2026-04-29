%% Parameters
subjectNames = {"013AR"};
runDisplayDemographicDetailsPairedSubjects;
folderSource = 'E:\NeoLabData\projectDhyaan\TWData';
d = dir(folderSource);
names = {d([d.isdir]).name};
subjectNames = names(~ismember(names, {'.','..'}));
medCount=0;
for i=1:length(subjectNames)
    [row, col] = find(strcmp(pairedSubjectNameList,subjectNames{i}));
    if col == 1
        medCount = medCount + 1;
    end
end
%% Load Data
for subjectIndex=1:length(subjectNames)
    subjectName = subjectNames{subjectIndex};
    tWDataFolder = fullfile(folderSource,subjectName);
    if isfile(fullfile(tWDataFolder,"properties.mat"))
        fprintf("Data for subject %s already exists!\n",subjectName);
        continue
    end

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

    % Load Data
    % protocol = "G1";
    clear outputs;
    clear outputAllProtocols;
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
    waveProperties = {"Wave Propagation Speed", "Wavelength", "Wave Propagation Direction", "Electrode Cluster Sizes", "Wave Strength AUC"};

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

        for propertyIndex=1:length(waveProperties)
            waveProperty = waveProperties{propertyIndex};
            if strcmp(waveProperty,'Wavelength')
                uniqueVals = waveLengths;
                xAxisLabel = 'Wavelength (mm)';
                fieldName = 'Wavelength';
            elseif strcmp(waveProperty,'Wave Propagation Direction')
                uniqueVals = dirs;
                fieldName = 'Direction';
            elseif strcmp(waveProperty,'Wave Propagation Speed')
                uniqueVals = speeds;
                xAxisLabel = 'Speed (m/s)';
                fieldName = 'Speed';
            elseif strcmp(waveProperty,'Electrode Cluster Sizes')
                uniqueVals = clusterSizes;
                xAxisLabel = '# Electrodes in Cluster';
                fieldName = 'ClusterSize';
            elseif strcmp(waveProperty,'Wave Strength AUC')
                uniqueVals = auc;
                xAxisLabel = 'Area Under Curve (Wave Strength)';
                fieldName = 'WaveStrength';
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
                subjectProperties.(protocol).(gammaBands{j}).(fieldName) = values{j};
            end
        end
    end
    save(fullfile(tWDataFolder,"properties.mat"),"subjectProperties");
    fprintf("Saved Properties for subject: %s\n",subjectName);
end