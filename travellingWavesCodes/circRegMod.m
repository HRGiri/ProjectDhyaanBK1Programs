function [direction,spatial_frequency,sl,Rsquared,PGD] = circRegMod(circularV,linearV)
% circ_lin_regress_2D is a function used to fit a 2D linear variables to a
% circular variable. One exmaple of the application is to fit phase traveling wave
% from linear coordinates.
%
% For method details, see Zabeh et al. 2023
% as per Theta and Alpha Oscillations Are Traveling Waves in the Human
% Neocortex -Honghui Zhang, Andrew J. Watrous, Ansh Patel, Joshua Jacobs
%
% linearV is a N by 2 matrix convey the linear variables. circularV (in radian) is a N
% by 1 array, which we are trying to fit. The algorithm handles wrapping-up cases on its own.
% visulization is a logic variable. It generlizes 2 visualization plot. One is the residual.
%
% direction of circular variable ascending with linear variables. Please
% note the direction of wave propagation will be 180 degree from this
% number as wave propagates to phase descending direction.
% spatial_frequency is the change rate of circular variable to linear
% variables at the ascending direction. The unit is radian per linear unit
% sl is the slopes of phase change relative to each column in linearV.
% Rsquared denotes how much variance of the circular
% variable is explained by the regression model. To access the statitical
% significace, please perform a permutaion procedure.
%% from https://github.com/erfanzabeh/WaveMonk
pos_x = linearV(:,1);
pos_y = linearV(:,2);
phase = mod(circularV,2*pi);
%helper function to calculate Residual resultant vector length between
%fit and actual phase
% here 'phase' is the actual phase and 'pos_x-slope2*pos_y' is the predicted phase
myfun = @(slope1,slope2)sqrt((sum(cos(phase-slope1*pos_x-slope2*pos_y)/length(phase)).^2 + (sum(sin(phase-(slope1*pos_x)-slope2*pos_y))/length(phase)).^2));


% arrange range and steplength for parameter space. angle ranges 2pi
% and spatial frequency range from 0 to 18 degree per unit. the upper
% limit of spatial frequency is depend on the spatial nyquist
% frequency. You may change the steplength to anything but watch out
% for computation time.

% angle_range=pi*(0:5:360)/180;
% spatial_frequency_range=(0:1:18)*pi/180;

angle_range=pi*(0:10:360)/180;
spatial_frequency_range=(0:1:18)*pi/180;

[angleMatrix,spatial_frequency_Matrix] = meshgrid(angle_range,spatial_frequency_range); % make it to a matrix for arrayfun

% transfer angle and spatial_frequency into 2d linear slopes
slope1_Matrix=cos(angleMatrix).*spatial_frequency_Matrix; % slope for pos_x
slope2_Matrix=sin(angleMatrix).*spatial_frequency_Matrix; % slope for pos_y
% this is a range of mean vector lengths r.
Residual_resultant_vector_length=arrayfun(myfun,slope1_Matrix,slope2_Matrix); % calculate the resultant_vector_length for each possible slope combos
[row,column]=find(Residual_resultant_vector_length==max(max(Residual_resultant_vector_length))); % find the largest resultant_vector_length

% Erfan edit 
row = max(row);
column = max(column);

%get the direction and spatial_frequency. If running traveling
%wave analysis, the direction should be flipped 180 degrees as waves
%propagates to the phase descending directions.
direction=angleMatrix(row,column); % find the best fit propagtion direction
spatial_frequency=spatial_frequency_Matrix(row,column); % find the best fit spatial frequency

sl=spatial_frequency*[cos(direction) sin(direction)]; %transform to linear slopes

% calculate offset for fitted values
offs = atan2(sum(sin(phase-sl(1)*pos_x-sl(2)*pos_y)),sum(cos(phase-sl(1)*pos_x-sl(2)*pos_y))) ;

% circular-linear correlation:
pos_circ = mod(sl(1)*pos_x+sl(2)*pos_y+offs, 2*pi); % circular variable derived from the position
phase_mean = mod(angle(sum(exp(1i*phase))/length(phase)),2*pi); % circular mean of the theta phase
pos_circ_mean = mod(angle(sum(exp(1i*pos_circ))/length(phase)),2*pi); % circular mean of the circular position variable
% cc is the circular correlation between predicted and actual phases
cc = sum(sin(phase - phase_mean) .* sin(pos_circ - pos_circ_mean)) / sqrt( sum(sin(phase - phase_mean).^2) * sum(sin(pos_circ - pos_circ_mean).^2) );
% calculating the circ-correlation coefficient between measured phase and fitted
Rsquared=cc^2;
% pos_circ=mod(pos_circ,2*pi);

% to calculate pgd
PGD = 1- (((1-Rsquared)*(length(circularV)-1))/(length(circularV)-size(linearV,2)-1));
end



