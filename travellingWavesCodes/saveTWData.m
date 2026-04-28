%% Params
subjectNames = {"013AR"};
protocols = {"G1","G2","M2"};
gammaBands = {"slow", "fast"};

folderSource = 'E:\NeoLabData\projectDhyaan\segmentedData';
montageFolder = 'E:\OneDrive - Indian Institute of Science\NeoLabData\Programs\Montages\Layouts\actiCap64_UOL';
labelFilePath = fullfile(montageFolder,'actiCap64_UOLLabels.mat');
montageLabels = load(labelFilePath).montageLabels;
chanlocFilePath = fullfile(montageFolder,'actiCap64_UOL.mat');
chanlocs = load(chanlocFilePath).chanlocs;

for subjectIndex = 1:length(subjectNames)
    subjectName = subjectNames{subjectIndex};
    for protocolIndex = 1:length(protocols)
        protocol = protocols{protocolIndex};

        eegFolder = fullfile(folderSource,subjectName,'EEG');
        d = dir(eegFolder);
        names = {d([d.isdir]).name};
        names = names(~ismember(names, {'.','..'}));
        date = names{1};
        %% Load Meta Data
        segmentedDataFolder = fullfile(folderSource,subjectName,'EEG',date,protocol,'segmentedData','LFP');
        lfpInfo = load(fullfile(segmentedDataFolder, 'lfpInfo.mat'));
        timeVals = lfpInfo.timeVals;

        badTrialsStruct = load(fullfile(segmentedDataFolder,'..','badTrials_wo_v8'));
        badTrials = badTrialsStruct.badTrials;
        badElecs = badTrialsStruct.badElecs.noisyElecs;

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
        numGoodTrials = length(goodTrials);

        
        %% Get TW Properties
        for gammaBandIndex = 1:length(gammaBands)
            gammaBand = gammaBands{gammaBandIndex};
            electrodesToAnalyse = 1:64;
            if strcmp(gammaBand,'slow')
                freqRange = [20 40];    %#ok<UNRCH>
            else
                freqRange = [40 60]; %#ok<UNRCH>
            end

            fprintf("Processing subject: %s, Protocol: %s, Gamma Band: %s\n", subjectName, protocol, gammaBand)
            outputs = cell(1,numGoodTrials);
            for trial = 1:numGoodTrials
                outputs{1,trial} = getTWCircParams(squeeze(data(trial,:,:)),timeVals,...
                    electrodesToAnalyse,freqRange,1,[],chanlocs);
                fprintf("Trial %d/%d\n", trial, numGoodTrials)
            end
            targetFolder = fullfile(folderSource,"..","TWData",subjectName);
            makeDirectory(targetFolder);
            save(fullfile(targetFolder,sprintf("%s_%sGamma.mat", protocol,gammaBand)),"outputs")
            disp("Done!")
        end
    end
end