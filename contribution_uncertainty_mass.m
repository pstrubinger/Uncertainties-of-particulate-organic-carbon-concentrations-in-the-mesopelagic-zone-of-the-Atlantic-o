function U_Mass = contribution_uncertainty_mass(U_Mp, volume)
%function U_Mass = contribution_uncertainty_mass(sigma_mass, volume)

%Calculate the contributions of mass estimates predicted by the calibration equation uncertainty to the combined uncertainty of POC

% Inputs:
% sigma_mass       : The uncertainties of carbon masses estimated from uPOC and aDOC filters by the regression models
% volume           : the volume of seawater filtered; in litres

U_Mass = U_Mp./volume;    % [ug/l]

end