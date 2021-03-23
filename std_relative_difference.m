function Std_RD = std_relative_difference(N,sum_relative_difference)
%function Std_RD = std_relative_difference(N,sum_relative_difference)
% 
% calculate the standard deviation of all the relative duplicate differences following Eq 4 in Hyslop and white (2009)
%
% Inputs:
% N                        : number of duplicate pairs
% sum_relative_difference  : sum of all the relative duplicate absolute differences

Std_RD = sqrt(pi/2) * (sum_relative_difference/N); % standard deviation of all the relative duplicate difference (Hyslop and white(2009)

end

