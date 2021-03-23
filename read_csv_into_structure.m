clear all
%% Import data from csv file.
%% select current directory
% directories
fileDir = "D:\Measuring POC_Paper\poc\POC\CHN_RUNs\";
dataDir = "D:\Measuring POC_Paper\poc\POC\Data\";
figDir = "D:\Measuring POC_Paper\poc\POC\Figures\";
fn = dir(strcat(fileDir, "*csv"));

for ifn = 1:length(fn)
%% Initialise variables.
filename = fn(ifn).name;
delimiter = ',';
startRow = 2;

%% Format used in each column:
%   column1: text (%s.%s.%f.%f) %C.n00.000.01
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)

fmt = '%s%f%f%f%f%f%[^\n\r]'; % fmt = formatSpec for textscan

%% Open the text file.
FID = fopen(strcat(fileDir, filename),'r');

%% Taking just the firt rows and saving it into a large string
labels = fgetl(FID); 

%% Extracting data using textscan into a new variable called dataArray
dataArray = textscan(FID, fmt, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-2, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose (FID); %% closing the filedata
labels_name = strsplit(labels,',');% slipting the string using the delimiters comma into a cell array

for ifield = 1:length (labels_name ) %% for loop to create and fill the structure called Results 
      Results.(labels_name{ifield} ) = dataArray{ifield};
end
Runame = string(filename(:,1:10));
Results.Run = [];
Results.Run = [Results.Run;repmat(Runame,length(Results.ID),1)];

%% Clear temporary variables
clearvars delimiter startRow fmt FID dataArray ans ifield labels labels_name

%% Working with the standards
igood = []; %creating a vector to save the index of the standards

%%% contains(Results_01.ID(1),'S') checking if the position 1 in....
...Results_01.ID contains/begins with the string 'S'%%%%

for i = 1:length(Results.ID)
    if contains (Results.ID(i),'S')
        igood = [igood, i];
    end
end

%% Creating structure called Standards for the standard calculations

field1 = 'RUN'; 
value1 = repmat(Runame,length(igood(:)),1);
field2 = 'ID'; 
value2 = Results.ID(igood(:));
field3 = 'Weight'; 
value3 = Results.Weight(igood(:));
field4 = 'CarbonPeak'; 
value4 = Results.Carbon(igood(:));
CarbonInStand = 0.7109;
field5 = 'CarbonMass'; % it is in ug; JROSS measured in mg.
value5 = (Results.Weight(igood)*1000*CarbonInStand);
Standards = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5);

%% Calibration Equation %%%
%%%Create a linear model to get the slope, y-intercept, residuals, 
mdl = fitlm (Standards.CarbonPeak(:), Standards.CarbonMass(:), 'linear', 'RobustOpts','on');%,'Intercept',false) %%Linear model 
IpValue = mdl.Coefficients{1,4};
if IpValue > 0.05
    mdl = fitlm (Standards.CarbonPeak(:), Standards.CarbonMass(:), 'linear','Intercept',false, 'RobustOpts','on');
end

%% Creating structure called Samples for the sample calculations
clearvars igood field1 value1 field2 value2 field3 value3 field4 value4 field5 value5 CarbonInStand i
igood = []; %creating a vector to save the index of the standards

for i = 1:length(Results.ID)
    if contains (Results.ID(i),'A')
        igood = [igood, i];
    end
end

field1 = 'RUN';
value1 = repmat(Runame,length(igood(:)),1);
field2 = 'ID';
value2 = Results.ID(igood(:));
field3 = 'Depth';
value3 = Results.Depth(igood(:));
field4 = 'Volumen'; 
value4 = Results.Volume(igood(:));
field5 = 'CarbonPeak'; 
value5 = Results.Carbon(igood(:));
field6 = 'CarbonMass';
value6 = predict(mdl,Results.Carbon(igood(:)));
Samples = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6);


%% %% Creating structure called Blanks for the sample calculations %%
clearvars igood field1 value1 field2 value2 field3 value3 field4 value4 field5 value5 field6 value6 i
igood = []; %creating a vector to save the index of the standards

for i = 1:length(Results.ID)
    if contains (Results.ID(i),'B')
        igood = [igood, i];
    end
end

field1 = 'RUN';
value1 = repmat(Runame,length(igood(:)),1);
field2 = 'ID';
value2 = Results.ID(igood(:));
field3 = 'Depth';
value3 = Results.Depth(igood(:));
field4 = 'Volumen'; 
value4 = Results.Volume(igood(:));
field5 = 'CarbonPeak'; 
value5 = Results.Carbon(igood(:));
field6 = 'CarbonMass';
value6 = predict(mdl,Results.Carbon(igood(:)));
Blanks = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6);

%% %% %% Creating structure called replicate of samples for the sample calculations %%
clearvars igood field1 value1 field2 value2 field3 value3 field4 value4 field5 value5 field6 value6 i
igood = []; %creating a vector to save the index of the standards

for i = 1:length(Results.ID)
    if contains (Results.ID(i),'R')
        igood = [igood, i];
    end
end


field1 = 'RUN';
value1 = repmat(Runame,length(igood(:)),1);
field2 = 'ID';
value2 = Results.ID(igood(:));
field3 = 'Depth';
value3 = Results.Depth(igood(:));
field4 = 'Volumen'; 
value4 = Results.Volume(igood(:));
field5 = 'CarbonPeak'; 
value5 = Results.Carbon(igood(:));
field6 = 'CarbonMass';
value6 = predict(mdl,Results.Carbon(igood(:)));
R_samples = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6);

%% %% Creating structure called replicate of blanks for the sample calculations %%
clearvars igood field1 value1 field2 value2 field3 value3 field4 value4 field5 value5 field6 value6 i
igood = []; %creating a vector to save the index of the standards

for i = 1:length(Results.ID)
    if contains (Results.ID(i),'P')
        igood = [igood, i];
    end
end

field1 = 'RUN';
value1 = repmat(Runame,length(igood(:)),1);
field2 = 'ID';
value2 = Results.ID(igood(:));
field3 = 'Depth';
value3 = Results.Depth(igood(:));
field4 = 'Volumen'; 
value4 = Results.Volume(igood(:));
field5 = 'CarbonPeak'; 
value5 = Results.Carbon(igood(:));
field6 = 'CarbonMass';
value6 = predict(mdl,Results.Carbon(igood(:)));
R_blanks = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6);

%% %% %% Creating structure called capsules for the sample calculations %%
clearvars igood field1 value1 field2 value2 field3 value3 field4 value4 field5 value5 field6 value6 i
igood = []; %creating a vector to save the index of the standards

for i = 1:length(Results.ID)
    if contains (Results.ID(i),'C')
        igood = [igood, i];
    end
end

field1 = 'RUN';
value1 = repmat(Runame,length(igood(:)),1);
field2 = 'ID';
value2 = Results.ID(igood(:));
field3 = 'CarbonPeak'; 
value3 = Results.Carbon(igood(:));
field4 = 'CarbonMass';
value4 = predict(mdl,Results.Carbon(igood(:)));
Capsules = struct(field1,value1,field2,value2,field3,value3,field4,value4);

%% %% %% %% Creating structure called filter_Acidified for the sample calculations %%
clearvars igood field1 value1 field2 value2 field3 value3 field4 value4 i
igood = []; %creating a vector to save the index of the standards

for i = 1:length(Results.ID)
    if contains (Results.ID(i),'T')
        igood = [igood, i];
    end
end

field1 = 'RUN';
value1 = repmat(Runame,length(igood(:)),1);
field2 = 'ID';
value2 = Results.ID(igood(:));
field3 = 'CarbonPeak'; 
value3 = Results.Carbon(igood(:));
field4 = 'CarbonMass';
value4 = predict(mdl,Results.Carbon(igood(:)));
Filter_acidified = struct(field1,value1,field2,value2,field3,value3,field4,value4);

%% %% %% %% %% Creating structure called filter_non_acidified for the sample calculations %%
clearvars igood field1 value1 field2 value2 field3 value3 field4 value4 i
igood = []; %creating a vector to save the index of the standards

for i = 1:length(Results.ID)
    if contains (Results.ID(i),'F')
        igood = [igood, i];
    end
end

field1 = 'RUN';
value1 = repmat(Runame,length(igood(:)),1);
field2 = 'ID';
value2 = Results.ID(igood(:));
field3 = 'CarbonPeak'; 
value3 = Results.Carbon(igood(:));
field4 = 'CarbonMass';
value4 = predict(mdl,Results.Carbon(igood(:)));
Filter_non_acidified = struct(field1,value1,field2,value2,field3,value3,field4,value4);

%% %% Clear temporary variables
clearvars igood field1 value1 field2 value2 field3 value3 field4 value4 i Runame

%% preparation input to calculate prediction interval
S = length(Standards.ID(:,:)); % number of standards used to fit model.
df = S-2; % degree of freedom
residuals = mdl.Residuals{:,1}; %
median_residuals = median(residuals);
Sigma_res = prcrng(residuals); % Robust standard deviation of the residuals
%% Calculation of PI for the samples
mean_area = mean(Samples.CarbonMass(:));
std_area = std(Samples.CarbonMass(:));
PI = prediction_interval(Samples.CarbonMass(:), 0.68, Sigma_res, S, mean_area, std_area, df);
Samples.PI = PI;
clearvars mean_area std_area PI
%% %% Calculation of PI for the blanks
mean_area = mean(Blanks.CarbonMass(:));
std_area = std(Blanks.CarbonMass(:));
PI = prediction_interval(Blanks.CarbonMass(:), 0.68, Sigma_res, S, mean_area, std_area, df);
Blanks.PI = PI;
clearvars mean_area std_area PI
%% %% Calculation of PI for the replicate of samples
mean_area = mean(R_samples.CarbonMass(:));
std_area = std(R_samples.CarbonMass(:));
PI = prediction_interval(R_samples.CarbonMass(:), 0.68, Sigma_res, S, mean_area, std_area, df);
R_samples.PI = PI;
clearvars mean_area std_area PI
%% %% %% Calculation of PI for the Replicate of blanks
mean_area = mean(R_blanks.CarbonMass(:));
std_area = std(R_blanks.CarbonMass(:));
PI = prediction_interval(R_blanks.CarbonMass(:), 0.68, Sigma_res, S, mean_area, std_area, df);
R_blanks.PI = PI;
clearvars mean_area std_area PI
%% %% %% %% Calculation of PI for the capsules
mean_area = mean(Capsules.CarbonMass(:));
std_area = std(Capsules.CarbonMass(:));
PI = prediction_interval(Capsules.CarbonMass(:), 0.68, Sigma_res, S, mean_area, std_area, df);
Capsules.PI = PI;
clearvars mean_area std_area PI
%% %% %% %% Calculation of PI for the acidified filters
mean_area = mean(Filter_acidified.CarbonMass(:));
std_area = std(Filter_acidified.CarbonMass(:));
PI = prediction_interval(Filter_acidified.CarbonMass(:), 0.68, Sigma_res, S, mean_area, std_area, df);
Filter_acidified.PI = PI;
clearvars mean_area std_area PI
%% %% %% %% Calculation of PI for the acidified filters
mean_area = mean(Filter_non_acidified.CarbonMass(:));
std_area = std(Filter_non_acidified.CarbonMass(:));
PI = prediction_interval(Filter_non_acidified.CarbonMass(:), 0.68, Sigma_res, S, mean_area, std_area, df);
Filter_non_acidified.PI = PI;
clearvars mean_area std_area PI

%% Saving the data
fnmat = strrep(filename, 'csv', 'mat')
fnmat = strcat (dataDir, fnmat);
save(fnmat, 'Standards', 'Samples', 'Blanks', 'R_samples', 'R_blanks', 'Capsules', 'Filter_acidified', 'Filter_non_acidified', 'mdl')
%% 
% predict C for everything that is not a standard%%%
% calculate prediction intervals for everything that is not a standard%%%%
% save carbon prediction intervals, mdl, inputs for PI

%% figures%%%%
figure(1) 
clf
subplot(1,2,1);
h = plot(mdl);
title('Regression model');
xlabel('Carbon Peak Area');
ylabel('Carbon mass ($\mu g$)', 'Interpreter', 'latex');
ylim ([0, 350]);
legend ( 'Calibration Standards', 'Linear Regression', 'Confidence bounds', 'Location', 'Southeast'); 
str1 = num2str(mdl.Coefficients.Estimate)
text(35000,300,str1,'FontSize', 8, 'Color', 'k')%,'FontWeight','bold');

subplot (1,2,2);
hold on
    plot(mdl.Variables.x1, mdl.Residuals.Raw, 'o');
    plot(mdl.Variables.x1, zeros(size(mdl.Variables.x1)) , 'k-');
    title ('Residuals');
    xlabel ('Carbon Peak Area');
    ylabel ('Residuals (\mug)');
    ylim([ -5, 5]);
    set (gcf, 'paperposition', [3.0917    9.2937   21.8167   11.1125]);
    box('on');
    w = num2str(median_residuals);
    str1 = '$\bar{\chi}$ =  ';
    str2 = [str1, w, ' $\mu g$'];
    str3 = ['\sigma_{res} = ' num2str(Sigma_res), ' \mug'];
    text(35000, 4.5, str2, 'Interpreter','latex', 'FontSize', 8, 'Color', 'k')%, 'FontWeight', 'bold');
    text(35000, 4, str3, 'FontSize', 8, 'Color', 'k')%, 'FontWeight', 'bold');


%% saving figure %%%
fnout = strrep(filename, 'csv', 'png');
print ('-dpng', strcat(figDir, 'Calibration_', fnout), '-r350');
%% 
end