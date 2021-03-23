function U_Volume = uncertainty_volume( M_POC, Volume, sigma_volume)
%function U_Volume = uncertainty_volume( M_POC, Volume, sigma_volume)

%   Calculate the uncertainty in volumen after applying the standard law of propagation of uncertainty
%   to the measurement equiation POC = M_POC / V

% Inputs:
% Volume       : the volume of seawater filtered; in litres
% M            : POC Carbon mass ; in ug
% Sigma_volume : The combined uncertainty of volume in each sample; in litres

U_Volume = M_POC.* sigma_volume ./ Volume.^2; %[ug * l / l^2] = [ug/l]

end