function PI = prediction_interval(x, interval, sigma_res, S, mean_x, std_x, df)
%function PI = prediction_interval(x, interval, sigma_res, S, mean_x, std_x)
% 
% calculate the prediction interval for a given regression following
% Altman(2000)
%
% Inputs:
% x         : independent variable
% interval  : width of the prediction interval (0 to 1)
% sigma_res : std of the regression residuals
% S         : number of points used for the regression
% mean_x    : mean of the independent variable used for the regression
% std_x     : std of the independent variable used for the regression
% df        : degree of freedom


alpha = 1-interval; 

t = tinv(1-alpha/2,df);


PI = t*sigma_res*sqrt(1 +1/S +(x - mean_x).^2/((S-1)*std_x^2)); % prediction interval (Altman,2000)


end

