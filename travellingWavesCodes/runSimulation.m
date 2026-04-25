function runSimulation(outputs,chanlocs,electrodesToAnalyse,timeVals,simulationPeriod,simulationSpeed,pauseAfterEveryFrame)
%RUNSIMULATION Summary of this function goes here
%   Detailed explanation goes here
electrodeLocs = projectTo2DPlane(chanlocs);
electrodeLocs = electrodeLocs(electrodesToAnalyse, :);
% Precompute grid (do this once)
gridRes = 100;
xlin = linspace(min(electrodeLocs(:,1)), max(electrodeLocs(:,1)), gridRes);
ylin = linspace(min(electrodeLocs(:,2)), max(electrodeLocs(:,2)), gridRes);
[Xq, Yq] = meshgrid(xlin, ylin);
% Compute convex hull of electrode positions
K = convhull(electrodeLocs(:,1), electrodeLocs(:,2));

% Create mask for valid region
inMask = inpolygon(Xq, Yq, electrodeLocs(K,1), electrodeLocs(K,2));

figure;
% start_t = find(timeVals>=0.716,1);
% end_t = find(timeVals>=0.724,1);
start_t = find(timeVals>=simulationPeriod(1),1);
end_t = find(timeVals>=simulationPeriod(2),1);
for t = start_t:end_t

    % --- Extract phase ---
    phi_t = outputs.phi(:, t);

    % --- Complex interpolation (correct for phase) ---
    complexPhase = exp(1i * phi_t);

    Zq_real = griddata(electrodeLocs(:,1), electrodeLocs(:,2), real(complexPhase), Xq, Yq, 'cubic');
    Zq_imag = griddata(electrodeLocs(:,1), electrodeLocs(:,2), imag(complexPhase), Xq, Yq, 'cubic');
    Zq = angle(Zq_real + 1i * Zq_imag);
    % center = mean(X);
    % radius = max(vecnorm(X - center, 2, 2));
    % 
    % circleMask = (Xq - center(1)).^2 + (Yq - center(2)).^2 <= radius^2;
    % Zq(~circleMask) = NaN;
    Zq(~inMask) = NaN;  % mask outside region

    % --- Plot ---
    clf; % clear figure each frame

    h = imagesc(xlin, ylin, Zq);
    set(gca, 'YDir', 'normal');
    axis equal;
    ylim([-1.1 1.1]);
    xlim([-1.1 1.1]);
    % Make NaNs transparent
    set(h, 'AlphaData', ~isnan(Zq));
    set(gca, 'Color', 'w'); % background white instead of red
    hold on;

    colormap(hsv);
    clim([-pi pi]);
    colorbar;

    % Electrodes
    scatter(electrodeLocs(:,1), electrodeLocs(:,2), 20, 'k', 'filled');

    % --- Cluster + direction ---
    clusterElectrodes = outputs.clusters{t};
    if ~isempty(clusterElectrodes)
        clusterElectrodesIndices = zeros(1,length(clusterElectrodes));
        for i=1:length(clusterElectrodesIndices)
            clusterElectrodesIndices(i) = find(electrodesToAnalyse == clusterElectrodes(i));
        end
    else
        clusterElectrodesIndices = [];
    end
    if ~isempty(clusterElectrodesIndices)
        theta = outputs.direction(t);

        % If needed, uncomment:
        % theta = deg2rad(theta);

        u = cos(theta);
        v = sin(theta);

        arrowScale = 0.05 * mean(range(electrodeLocs));
        
        if clusterElectrodesIndices ~= 0
            quiver( ...
                electrodeLocs(clusterElectrodesIndices,1), ...
                electrodeLocs(clusterElectrodesIndices,2), ...
                arrowScale * u * ones(length(clusterElectrodesIndices),1), ...
                arrowScale * v * ones(length(clusterElectrodesIndices),1), ...
                0, ...
                'k', 'LineWidth', 1.5, 'MaxHeadSize', 2 ...
            );
        end
    end

    title(sprintf('t = %.3f ms', 1000*timeVals(t)));
    xlabel('X'); ylabel('Y');

    drawnow;
    if pauseAfterEveryFrame
        pause();
    else
        pause(simulationSpeed);  
    end
end
end

