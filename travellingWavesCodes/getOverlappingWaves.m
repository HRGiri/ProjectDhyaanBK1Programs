% taken from "https://github.com/gouthambhargava/TravellingWaveProject/blob/main/commonAnalysisCodes/analysisCodes/getOverlappingWaves.m"

function [newBounds,allDirSg,allDirFg,allUniqueDirs,emptyCell,intPts] = getOverlappingWaves(dirSg,sgBounds,dirFg,fgBounds,overlap,metFlag)
%Find the time indices of overlapping waves between slow and fast gamma
%waves.
% Inputs - the vector of direction values for waves identified in slow and fast gamma are given by dirSg and dirFg
% and the start of end points of each wave is given by boundriesTimesSg and boundriesTimesFg
% Overlap - this defines the amount of overlap between the waves that is to be considered significant. Set to 50% by default.
% metFlag = set to 1 to get mean of overlapping waves, set to 2 to get mean of only overlapping timepoints between waves
%
if nargin<5
    overlap = 0.5;
end

if nargin<6
    metFlag = 2;
end


fgCounts = zeros(1,size(fgBounds,2));
intPts = nan(1,length(dirSg));
overlapTW = [];
for j = 1:size(sgBounds,2)
    singleSg = sgBounds(1,j):sgBounds(2,j);
    sgLength = length(singleSg);
    for k = 1:length(fgCounts)
        singleFg = fgBounds(1,k):fgBounds(2,k);
        fgLength = length(singleFg);
        intersectLength = numel(intersect(singleSg,singleFg));
        if intersectLength>=overlap*sgLength || intersectLength>=overlap*fgLength
            overlapTW = cat(1,overlapTW,[j,k]);
            intPts(intersect(singleSg,singleFg)) = 1;
        end
    end
end

%% get new boundries, only the TW's which overlap
if ~isempty(overlapTW)
    newBounds = cell(2,1);
    for j = 1:size(overlapTW,1)
        newBounds{1}(:,j) = sgBounds(:,overlapTW(j,1));
        newBounds{2}(:,j) = fgBounds(:,overlapTW(j,2));
    end

    % get unique directions of overlapping waves
    allDirSg = nan(size(dirSg));
    allDirFg = nan(size(dirFg));

    if metFlag==1
        uniqueDirs1 = [];
        uniqueDirs2 = [];
        for j = 1:size(newBounds{1},2)
            bounds1 = newBounds{1}(:,j);
            bounds2 = newBounds{2}(:,j);
            allDirSg(bounds1(1):bounds1(2)) = dirSg(bounds1(1):bounds1(2));
            allDirFg(bounds2(1):bounds2(2)) = dirFg(bounds2(1):bounds2(2));
            vals1 = circ_mean(dirSg(bounds1(1):bounds1(2))');
            vals2 = circ_mean(dirFg(bounds2(1):bounds2(2))');
            vals1(isnan(vals1)) = [];
            vals2(isnan(vals2)) = [];
            uniqueDirs1 = cat(2,uniqueDirs1,vals1);
            uniqueDirs2 = cat(2,uniqueDirs2,vals2);
        end
        allUniqueDirs = cat(1,uniqueDirs1,uniqueDirs2);
    else

        for j = 1:size(newBounds{1},2)
            bounds1 = newBounds{1}(:,j);
            bounds2 = newBounds{2}(:,j);
            allDirSg(bounds1(1):bounds1(2)) = dirSg(bounds1(1):bounds1(2));
            allDirFg(bounds2(1):bounds2(2)) = dirFg(bounds2(1):bounds2(2));
        end
        ov1 = allDirSg;
        ov1(isnan(allDirFg)) = nan;
        ov2 = allDirFg;
        ov2(isnan(allDirSg)) = nan;
        [~,~,uniqueDirs1] = simpleWaveSegments(ov1);
        [~,~,uniqueDirs2] = simpleWaveSegments(ov2);
        allUniqueDirs = cat(1,uniqueDirs1',uniqueDirs2');
    end
    emptyCell = 0;
else
    allUniqueDirs = [nan;nan];
    newBounds = nan;
    allDirSg = nan(size(dirSg));
    allDirFg = nan(size(dirFg));
    emptyCell = 1;
end
end