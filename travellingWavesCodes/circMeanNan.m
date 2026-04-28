% taken from "https://github.com/gouthambhargava/TravellingWaveProject/blob/main/commonAnalysisCodes/analysisCodes/circMeanNan.m"

function [mu] = circMeanNan(alpha, dim)
%
% [mu ul ll] = circ_mean(alpha, w, dim)
%   Computes the mean direction for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		weightings in case of binned angle data]
%     [dim  compute along this dimension, default: 1st non-singular dimension]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_mean(alpha, [], dim)
%
%   Output:
%     mu		mean direction

%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al. 
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 2
  dim = find(size(alpha) > 1, 1, 'first');
  if isempty(dim)
    dim = 1;
  end
end

% compute weighted sum of cos and sin of angles
mu = angle(nansum(exp(1i*alpha),dim));
% mu(mu==0) = nan;
end