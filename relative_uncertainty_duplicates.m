function RUD = relative_uncertainty_duplicates(Std_RD, mean_D, D1, D2)
%function RUD = relative_uncertainty_duplicates(Std_RD, mean_D, D1, D2)
% 
% calculate the relative uncertainty of duplicate estimates derived from...
% the standard deviation of all relative duplicate differences(Eq.4 in Hyslop and white(2009)..
% that is equivalent to the relative uncertainty of the difference in...
% duplicate estimates. Eq.10 in the manuscript.
% Inputs:
% Std_RD  : standard deviation of all the relative duplicate difference
% mean_D  : mean POC concentration of the two duplicate measurements
% D1      : First duplicate (Sample)
% D2      : Second duplicate (Replicate)

RUD = Std_RD * mean_D * sqrt(2) / sqrt(D1.^2 + D2.^2);

end

