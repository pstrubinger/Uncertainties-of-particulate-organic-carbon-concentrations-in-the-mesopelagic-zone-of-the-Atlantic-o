function out = total_uncertainty_u_POC(Volume, PI_uPOC, PI_cap, PI_acidified, PI_nonacidified, M_uPOC, M_cap, M_acidified, M_nonacidified , sigma_volume)
%function out = total_uncertainty_u_POC(Volume, PI_uPOC, PI_cap, PI_acidified, PI_nonacidified, M_uPOC, M_cap, M_acidified, M_nonacidified sigma_volume)

%   Calculate the total uncertainty of the uPOC concentration after applying the standard law of propagation of uncertainty
%   to the measurement equiation POC = (M_uPOC - M_aDOC)/ V

% Inputs:
% Volume       : the volume of seawater filtered; in litres
% PI           : Prediction intervals, uncertainties of the carbon masses predicted by the calibration equations Altman(2000)
% M            : Carbon mass; in ug
% Sigma_volume : The combined uncertainty of volume in each sample

out = sqrt(Volume.^2 .* (PI_uPOC.^2 + PI_cap.^2 + PI_acidified.^2 + PI_nonacidified.^2) + (M_uPOC - M_cap - M_acidified + M_nonacidified).^2 .* sigma_volume.^2) ./ Volume.^2; 

end