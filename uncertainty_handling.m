function U_Handling = uncertainty_handling(std_mean,volume)
% function out = uncertainty_handling(std_mean_ac, std_mean_non, volume)

% Calculate the contribution of uncertainty of the total bias due to sample handling to the total POC uncertainty
% Inputs:

% std_mean : the standard error of the mean of the three corresponding estimates of acidified and non-acidified filter blanks
% volume   : the volume of seawater filtered; in litres


U_Handling = sqrt(std_mean.^2 + std_mean.^2)./ volume; % [ug/l]

end

