function out = relative_uncertainty_contribution(U_source, U_POC)
%function out = relative_uncertainty_contribution(U_source, U_POC)

% Calculate the relative contributions of each source of uncertainty to the combined uncertainty of POC

%Inputs:

% U_source : The contribution of uncertainty of each source to the total POC uncertainty in [ug/l]
% U_POC    : Total POC uncertainty in [ug/l]

out = U_source./ U_POC; % dimensionless

end

