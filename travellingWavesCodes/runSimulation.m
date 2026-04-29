function runSimulation(outputs,chanlocs,electrodesToAnalyse,timeVals,simulationPeriod,simulationSpeed,pauseAfterEveryFrame)
%RUNSIMULATION Summary of this function goes here
%   Detailed explanation goes here
electrodeLocs = projectTo2DPlane(chanlocs);
electrodeLocs = electrodeLocs(electrodesToAnalyse, :);
% Precompute grid (do this once)
gridRes = 100;
xlin = linspace(min(electrodeLocs(:,1))-0.1, max(electrodeLocs(:,1))+0.1, gridRes);
ylin = linspace(min(electrodeLocs(:,2))-0.1, max(electrodeLocs(:,2))+0.1, gridRes);
[Xq, Yq] = meshgrid(xlin, ylin);

% Construct an artificial circle
center = mean(electrodeLocs);
radius = max(vecnorm(electrodeLocs - center, 2, 2)) + 0.05;
theta = linspace(0, 2*pi, 100);

xc = center(1) + radius * cos(theta);
yc = center(2) + radius * sin(theta);

% Write into a video
% videoWriter = VideoWriter('travellingWave.mp4', 'MPEG-4');
% videoWriter.FrameRate = 1/simulationSpeed; % adjust speed
% open(videoWriter);

figure;
start_t = find(timeVals>=simulationPeriod(1),1);
end_t = find(timeVals>=simulationPeriod(2),1);
for t = start_t:end_t

    % --- Extract phase ---
    phi_t = outputs.phi(:, t);

    % Extrapolation of phase to the boundary
    phi_boundary = zeros(size(xc));

    for k = 1:length(xc)
        d = vecnorm(electrodeLocs - [xc(k), yc(k)], 2, 2);
        [~, idx] = min(d);
        phi_boundary(k) = phi_t(idx);
    end
    
    % --- Complex interpolation (correct for phase) ---
    complexPhase = exp(1i * phi_t);    
    complex_boundary = exp(1i * phi_boundary);

    X_aug = [electrodeLocs(:,1); xc(:)];
    Y_aug = [electrodeLocs(:,2); yc(:)];
    
    Z_real_aug = [real(complexPhase); real(complex_boundary(:))];
    Z_imag_aug = [imag(complexPhase); imag(complex_boundary(:))];

    Zq_real = griddata(X_aug, Y_aug, Z_real_aug, Xq, Yq, 'natural');
    Zq_imag = griddata(X_aug, Y_aug, Z_imag_aug, Xq, Yq, 'natural');

    Zq = angle(Zq_real + 1i * Zq_imag);
    
    % --- Plot ---
    % clf; % clear figure each frame    
    subplot(2,3,t-start_t+1)
    h = imagesc(xlin, ylin, Zq);
    set(gca, 'YDir', 'normal');
    axis equal;
    ylim([-1.1 1.1]);
    xlim([-1.1 1.1]);
    axis off;
    % Make NaNs transparent
    set(h, 'AlphaData', ~isnan(Zq));
    set(gca, 'Color', 'w'); % background white instead of red
    hold on;

    colormap(hsv);
    clim([-pi pi]);
    colorbar;
    
    % Plot circle  
    plot(xc, yc, 'k', 'LineWidth', 1);

    % Electrodes
    scatter(electrodeLocs(:,1), electrodeLocs(:,2), 30, 'k', 'filled');

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

        arrowScale = 0.1 * mean(range(electrodeLocs));
        
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

    % title(sprintf('Phase Map\nt = %.0f ms', 1000*timeVals(t)));
    title(sprintf('t = %.0f ms', 1000*timeVals(t)));
    % xlabel('X'); ylabel('Y');

    drawnow;
    if pauseAfterEveryFrame
        pause();
    else
        pause(simulationSpeed);  
    end
    
    % --- Capture frame ---
    % frame = getframe(gcf);
    % writeVideo(videoWriter, frame);
end
% close(videoWriter);
end

