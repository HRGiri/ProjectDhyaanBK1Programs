function [outputs] = getTWCircParams(data,timeVals,goodElectrodes,freqs,req,nPerm,chanlocs)
%% Inputs
%%Inputs
% data - single trial data in the electrodes x time points format
% goodElectrodes - list of good electrodes
% timeVals - vector of time values
% freqs - vector of frequency limits (ex. [30 60])
% req - 0 (no filtering required, if MP is being used) or 1 (butterworth filter being used)
%% 
if nargin<5
    nPerm = [];
end
%% Parameters    
    fs = 1000; %sampling frequency
    numPCs = 5;
%% filter and extract instant phase of the lfp data
    data = data(goodElectrodes,:);
    burstTS = zeros(size(data,1),size(data,2));
    filterOrder = 4; 
    thresholdFactor = 3;
    baselinePeriodS = [-1 0]; 
    stimulusPeriodS = [0.25 1.25];

    for i = 1:size(data,1)
        % Get the burst time series and the filtered data
        [~,~,~,burstTS(i,:),data(i,:),~] = getHilbertBurst(data(i,:),timeVals,thresholdFactor,0,stimulusPeriodS,baselinePeriodS,freqs,filterOrder,req);
    end
    
    burstTS(isnan(burstTS)) = 0;      
    data = data - mean(data, 2);        % DC subtraction
    phiMat = angle(hilbert(data'))';       

    % Load Channel Locations
    X = projectTo2DPlane(chanlocs);
    X = X(goodElectrodes,:);
    X = X * 90;   % Head Radius in mm
    % deltaSpace = 90*getHighestSpacing(X);    
    
    timePoints = size(data,2);
    
    %initialize results
    pgd = zeros(timePoints,1);
    direction = zeros(timePoints,1);    
    cluster = cell(timePoints,1);
    circVmean = zeros(timePoints,1);
    sFreq = zeros(timePoints,1);
    if ~isempty(nPerm)
    pgdPerm = zeros(timePoints,nPerm);
    end

    % Clustering
    % PCA    
    [coeff, ~] = pca(data');
    % coeff = coeff(:, 1);    % Take the first PC
    % weights = abs(coeff);
    coeff = coeff(:,numPCs);    % Take the first n PCs
    weights = abs(sum(coeff,2));
    pcElecs = weights >= median(weights);   % Select the top 50% electrodes contributing to the PCs
    % pcElecs = true(size(pcElecs));  % Uncomment to not use PCA
    % run circ reg after getting significant electrodes from burst
    % detection
        for timei = 1:size(phiMat,2)
            phiGrid = phiMat(:,timei);
            elecs = find(burstTS(:,timei) & pcElecs);            
            %skip to next time point if cluster has <4 electrodes
            if numel(elecs) < 4
                % pgd(:,timei) = 0;
                % direction(:,timei) = 0;
                % cluster(:,timei) = 0;
                pgd(timei) = 0;
                direction(timei) = 0;
                cluster{timei} = [];
            else
                cluster{timei} = goodElectrodes(elecs);
                % elecs = getContiguousElectrodes(X,elecs);
                clusters{1} = elecs;
                for clusterIndex = 1:length(clusters)                
                    elecs = clusters{clusterIndex};
                    cluster{timei,clusterIndex} = goodElectrodes(elecs);
                circularCord = phiGrid(elecs);
                circVmean(timei,1) = circ_mean(circularCord);                
                linearCord = X(elecs,:);                
                % do regression analysis on the polar and linear coordinates
                [direction(timei,1),sFreq(timei,1),~,~,pgd(timei,1)] = circRegMod(circularCord,linearCord);
                if ~isempty(nPerm) 
                    for perm = 1:nPerm 
                        permVar = zeros(length(circularCord),nPerm);
                        permVar(:,perm) = circularCord(randperm(length(circularCord)));
                    end
                    for perm = 1:nPerm
                        [~,~,~,~,pgdPerm(timei,perm)] = circRegMod(permVar(:,perm),linearCord);
                    end
                end
                end
            end
%               display(['Done with timepoint:',num2str(timei)])
        end
%%  get additional params
    direction(isnan(direction)) = 0;
    pgd(isnan(pgd)) = 0;
    outputs.direction = direction+pi;
    outputs.Wavelength = (1./sFreq);  % in mm/rad
    deltaT = 1/fs; % in s
    phi_dot= abs(diff(circVmean))/deltaT; 
    outputs.tempFreq = phi_dot/(2*pi); % convert to Hz
    outputs.speed = outputs.tempFreq./sFreq(2:end); % convert to mm/s
    outputs.pgd = pgd;
    outputs.clusters = cluster; %the electrodes involved in the TW 
    outputs.phi = phiMat;
    if ~isempty(nPerm)
    outputs.pgdPerm = prctile(pgdPerm',0.99);
    end
end