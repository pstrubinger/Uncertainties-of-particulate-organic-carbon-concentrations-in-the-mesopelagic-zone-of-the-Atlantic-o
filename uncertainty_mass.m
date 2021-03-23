function U_Mp = uncertainty_mass(PI_uPOC, PI_aDOC, r)
%function U_Mp = uncertainty_mass(PI_uPOC, PI_aDOC, r)

%   Calculate the uncertainty in mass estimates predicted by the calibration equation after applying the standard law of propagation of uncertainty
%   to the measurement equiation POC = (M_uPOC - M_aDOC)

% Inputs:
% PI           : Prediction intervals, uncertainties of the carbon masses predicted by the calibration equations Altman(2000)
% r            : Correlation coefficient between M_uPOC and M_aDOC


U_Mp = sqrt(PI_uPOC.^2 + PI_aDOC.^2 - 2 * PI_uPOC .* PI_aDOC * r); %[ug]

end
