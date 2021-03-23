function U_POC = uncertainty_total(PI_uPOC, PI_aDOC, r, M_POC, volume, sigma_volume, U_Handling)
%function U_POC = uncertainty_total(PI_uPOC, PI_aDOC, r, M_POC, volume, sigma_volume, U_Handling)

%   Calculate the total uncertainty of the POC concentration after applying the standard law of propagation of uncertainty
%   to the measurement equiation POC = (M_uPOC - M_aDOC)/ V

% Inputs:
% M_POC            : POC Carbon mass ; in ug
% sigma_mass       : The uncertainties of carbon masses estimated from uPOC and aDOC filters by the regression models
% Volume           : the volume of seawater filtered; in litres
% Sigma_volume     : The combined uncertainty of volume in each sample
% U_Handling      : The contribution of uncertainty of the total bias due to sample handling to the total POC uncertainty in [ug/l]
% keyboard

U_Mp = sqrt(PI_uPOC.^2 + PI_aDOC.^2 - 2 * PI_uPOC .* PI_aDOC * r); % [ug]

U_Mass = U_Mp./volume; % [ug/l]
 
U_Volume = M_POC.* sigma_volume ./ volume.^2; % [ug * l / l^2] = [ug/l]

U_POC = sqrt(U_Mass.^2 + U_Volume.^2 + U_Handling.^2); % in [ug/l]

end

