function out = prcrng(x)
% function out = prcrng(x)
%

% Calculate the robust standard deviation of the residuals 
% It use the percentile 84 and the percentile 16 of the residual from the calibration equation.

% Inputs:

% x : regression residuals  

out = (prctile(x,84)- prctile(x,16))/2;
end

