%% Parameters
runDisplayDemographicDetailsPairedSubjects;
folderSource = 'E:\NeoLabData\projectDhyaan\TWData';
d = dir(folderSource);
names = {d([d.isdir]).name};
subjectNames = names(~ismember(names, {'.','..'}));
waveProperties = {"Wave Propagation Speed", "Wavelength", "Wave Propagation Direction", "Electrode Cluster Sizes", "Wave Strength AUC"};
protocols = {"G1","G2","M2"};
gammaBands = {"slow","fast"};

%% Load Data
meditators = struct();
controls = struct();
for protocolIndex = 1:length(protocols)
    protocol = protocols{protocolIndex};
    for gammaBandIndex = 1:length(gammaBands)
        gammaBand = gammaBands{gammaBandIndex};
        for propertyIndex=1:length(waveProperties)
                waveProperty = waveProperties{propertyIndex};
                if strcmp(waveProperty,'Wavelength')                    
                    xAxisLabel = 'Wavelength (mm)';
                    fieldName = 'Wavelength';
                elseif strcmp(waveProperty,'Wave Propagation Direction')                    
                    fieldName = 'Direction';
                elseif strcmp(waveProperty,'Wave Propagation Speed')                    
                    xAxisLabel = 'Speed (m/s)';
                    fieldName = 'Speed';
                elseif strcmp(waveProperty,'Electrode Cluster Sizes')                    
                    xAxisLabel = '# Electrodes in Cluster';
                    fieldName = 'ClusterSize';
                elseif strcmp(waveProperty,'Wave Strength AUC')                    
                    xAxisLabel = 'Area Under Curve (Wave Strength)';
                    fieldName = 'WaveStrength';
                end

                controls.(protocol).(gammaBand).(sprintf("mean%s",fieldName)) = [];
                meditators.(protocol).(gammaBand).(sprintf("mean%s",fieldName)) = [];
        end
    end
end

for subjectIndex=1:length(subjectNames)
    subjectName = subjectNames{subjectIndex};    
    tWDataFile = fullfile(folderSource,subjectName,'properties.mat');
    subjectProperties = load(tWDataFile).subjectProperties;    
    [~,isControl] = find(strcmp(pairedSubjectNameList,subjectName));    
    isControl = isControl - 1;  

    % Load Data
    for protocolIndex = 1:length(protocols)
        protocol = protocols{protocolIndex};        
        for gammaBandIndex = 1:length(gammaBands)
            gammaBand = gammaBands{gammaBandIndex};
            for propertyIndex=1:length(waveProperties)
                waveProperty = waveProperties{propertyIndex};
                if strcmp(waveProperty,'Wavelength')                    
                    xAxisLabel = 'Wavelength (mm)';
                    fieldName = 'Wavelength';
                elseif strcmp(waveProperty,'Wave Propagation Direction')                    
                    fieldName = 'Direction';
                elseif strcmp(waveProperty,'Wave Propagation Speed')                    
                    xAxisLabel = 'Speed (m/s)';
                    fieldName = 'Speed';
                elseif strcmp(waveProperty,'Electrode Cluster Sizes')                    
                    xAxisLabel = '# Electrodes in Cluster';
                    fieldName = 'ClusterSize';
                elseif strcmp(waveProperty,'Wave Strength AUC')                    
                    xAxisLabel = 'Area Under Curve (Wave Strength)';
                    fieldName = 'WaveStrength';
                end
            if isControl
                controls.(protocol).(gammaBand).(sprintf("mean%s",fieldName)) = [controls.(protocol).(gammaBand).(sprintf("mean%s",fieldName)) mean(subjectProperties.(protocol).(gammaBand).(fieldName))];                
            else
                meditators.(protocol).(gammaBand).(sprintf("mean%s",fieldName)) = [meditators.(protocol).(gammaBand).(sprintf("mean%s",fieldName)) mean(subjectProperties.(protocol).(gammaBand).(fieldName))];            
            end
            end
        end
    end
    fprintf("Done loading for subject: %s\n",subjectName);
end

%% Plot
    % waveProperty = 'Wave Propagation Speed';
    % waveProperty = 'Wavelength';
    % % waveProperty = 'Wave Propagation Direction';
    % waveProperty = 'Electrode Cluster Sizes';
    % waveProperty = 'Wave Strength AUC';
    % Get properties
    figure
    for protocolIndex = 1:length(protocols)
        protocol = protocols{protocolIndex};
         for gammaBandIndex = 1:length(gammaBands)
            gammaBand = gammaBands{gammaBandIndex};
        if strcmp(waveProperty,'Wavelength')                
            xAxisLabel = 'Wavelength (mm)';
            fieldName = 'Wavelength';
        elseif strcmp(waveProperty,'Wave Propagation Direction')                
            fieldName = 'Direction';
        elseif strcmp(waveProperty,'Wave Propagation Speed')                
            xAxisLabel = 'Speed (m/s)';
            fieldName = 'Speed';
        elseif strcmp(waveProperty,'Electrode Cluster Sizes')                
            xAxisLabel = '# Electrodes in Cluster';
            fieldName = 'ClusterSize';
        elseif strcmp(waveProperty,'Wave Strength AUC')                
            xAxisLabel = 'Area Under Curve (Wave Strength)';
            fieldName = 'WaveStrength';
        end

        data = [meditators.(protocol).(gammaBand).(sprintf("mean%s",fieldName))' controls.(protocol).(gammaBand).(sprintf("mean%s",fieldName))'];

        subplot(length(gammaBands),length(protocols),3*(gammaBandIndex-1)+protocolIndex);
        violinplot([1 2],data);
        xticks([1 2])
        xticklabels({'Meditators','Controls'})
        if protocolIndex == 1
            ylabel(xAxisLabel)        
        end
        if protocolIndex == 3
            if strcmp(gammaBand,'slow')
                text(2.75,0.25,"Slow Gamma", 'FontWeight','bold','Rotation',90)
            else
                text(2.75,0.75,"Fast Gamma",'FontWeight','bold','Rotation',90)
            end
        end
        if gammaBandIndex == 1
            title(protocol)
        end
        end

    end
ax = findall(gcf,'type','axes');
linkaxes(ax, 'xy');         
set(findall(gcf,'-property','FontSize'), 'FontSize', 14);
% set(findall(gcf, 'Type', 'text', 'Tag', 'Title'), 'FontSize',18);
% set(findall(gcf, 'Type', 'text', 'Tag', 'XLabel'), 'FontSize',16);
sgtitle(waveProperty,'FontWeight','bold','FontSize',20)