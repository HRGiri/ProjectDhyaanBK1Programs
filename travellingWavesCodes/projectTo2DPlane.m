function X = projectTo2DPlane(chanlocs)
%PROJECTTO2DPLANE Summary of this function goes here
%   Detailed explanation goes here
    R = [chanlocs.Y; chanlocs.X; chanlocs.Z]';
    R_c = R-mean(R);
    [~, score, ~] = pca(R_c);
    X = score(:, 1:2);
end

