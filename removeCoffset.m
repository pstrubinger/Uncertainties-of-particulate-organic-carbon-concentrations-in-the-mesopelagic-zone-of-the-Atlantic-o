function out = removeCoffset(sample, capsule, acidified, non_acidified)
% function out = removeCoffset(sample,capsule,non_acidified,acidified)
%
%   Remove the mass contains in the capsule, non-acidified, and acidified
%   filter to the mass of the sample and the blanks 

%Inputs: 
% sample        : carbon mass of the sample
% capsule       : carbon mass of the capsule
% non_acidified : carbon mass of the  non-acidified filter
% acidified     : carbon mass of the acidified filter

out = sample - capsule - (acidified - non_acidified);

end

