%Script for convert all individual structures of Samples, Blanks and Replicates into one total structure
clear all;
dataDir = "D:\Measuring POC_Paper\poc\POC\Data\";
fn = dir(strcat(dataDir, "*mat"));
etopoFile = dir(strcat(dataDir, "*tif"));
duplicateDir = "D:\Measuring POC_Paper\poc\POC\Data\Duplicates\";
JossDir = 'D:\Measuring POC_Paper\poc\POC\POC_BLANKS\';
fn1 = dir(strcat(JossDir, '*mat'));
fname = dir(strcat(duplicateDir, "*mat"));
fname = fname.name;
figDir = "D:\Measuring POC_Paper\poc\POC\Figures\";

flds = {'RUN', 'ID', 'Depth', 'Volumen','CarbonMass','PI','uPOCuncertainty','Comb_Unc_Volume','TrueCarbonMass','MassPOC','POC', 'uncertainty_handling'};
flds1 = {'RUN', 'ID', 'Depth', 'Volumen','CarbonMass','PI','aDOCuncertainty','Comb_Unc_Volume','TrueCarbonMass','MassPOC','POC', 'uncertainty_handling'};

for iflds = 1:length(flds)
    Total_POC.(flds{iflds}) = []; %Creating a empty structure to store data.
    Total_uPOC.(flds{iflds}) = [];
end
for iflds1 = 1:length(flds1)
    Total_aDOC.(flds1{iflds1}) = [];
end
for ifn = 1:length (fn) % calling each file .mat and loading it.
    filename = fn(ifn).name;
    load(strcat(dataDir, filename));
  
    for iflds = 1:length(flds) % vertical concatenation of each field.
        Total_POC.(flds{iflds}) = [Total_POC.(flds{iflds}); Samples.(flds{iflds})];
        Total_uPOC.(flds{iflds}) = [Total_uPOC.(flds{iflds}); Samples.(flds{iflds})];
    end
    for iflds1 = 1:length(flds1)
        Total_aDOC.(flds1{iflds1}) = [Total_aDOC.(flds1{iflds1}); Blanks.(flds1{iflds1})];
    end
end

clearvars 'Standards' 'Samples' 'Blanks' 'R_samples' 'R_blanks' 'Capsules' 'Filter_acidified' 'Filter_non_acidified' 'mdl' 'ifn' 'iflds' 'flds' 'iflds1' 'flds1';


%% Calculating the Correlation coefficient between uPOC and aDOC
Tmp = corrcoef(Total_uPOC.TrueCarbonMass,Total_aDOC.TrueCarbonMass);
Total_POC.corr = Tmp(2,1);
clearvars 'Tmp'
%% Calculating the Total Uncertainty of POC concentration

Total_POC.U_POC = uncertainty_total(Total_uPOC.PI, Total_aDOC.PI, Total_POC.corr, Total_POC.MassPOC, Total_POC.Volumen, Total_POC.Comb_Unc_Volume, Total_POC.uncertainty_handling);
% U2_POC = uncertainty_total(PI_uPOC, PI_aDOC, r, M_POC, volume, sigma_volume, U2_Handling)

%% %% Dividing the samples into zone, Productive and Mesopelagic.

for idepth = 1:length(Total_POC.ID)
       if Total_POC.Depth(idepth) <= 200 %Zones : 1 for Productive zone; 2 for mesopelagic zone
          Total_POC.Zone(idepth) = 1;
       else
          Total_POC.Zone(idepth)= 2;
       end
end
Total_POC.Zone = Total_POC.Zone';
iP = find (Total_POC.Zone == 1);
iM = find (Total_POC.Zone == 2);

%% Calculating experimental POC uncertaity using percentile precision

load(strcat(duplicateDir, fname)); %Loading data from duplicate analysis. 

Total_POC.Expe_POCuncer_PP = nan(size(Total_POC.ID)); 
Total_POC.Expe_POCuncer_PP(iP) = RU_P_PP*Total_POC.POC(iP);
Total_POC.Expe_POCuncer_PP(iM) = RU_M_PP*Total_POC.POC(iM);

%% Calculating relative contribution of uncertainty from calibration equation to the combined uncertainty of POC

%   Calculate the uncertainty in mass estimates predicted by the calibration equation after applying the standard law of propagation of uncertainty
%   to the measurement equiation POC = (M_uPOC - M_aDOC)
%   U_Mp = uncertainty_mass(PI_uPOC, PI_aDOC, r)
Total_POC.U_Mp = uncertainty_mass(Total_uPOC.PI, Total_aDOC.PI, Total_POC.corr); 

%   Calculate the contributions of mass estimates predicted by the calibration equation uncertainty to the combined uncertainty of POC
%   U2_Mass = contribution_uncertainty_mass(sigma_mass, volume)
Total_POC.U_Mass = contribution_uncertainty_mass(Total_POC.U_Mp, Total_POC.Volumen);

%   Calculate the relative contributions of mass estimates predicted by the calibration equation uncertainty to the combined uncertainty of POC
Total_POC.relative_U2_Mass = Total_POC.U_Mass.^2./ Total_POC.U_POC.^2;
relative_U2_Mass_productive = median(Total_POC.relative_U2_Mass(iP));
relative_U2_Mass_mesopelagic = median(Total_POC.relative_U2_Mass(iM));

%out = relative_uncertainty_contribution(U_source, U_POC)
Total_POC.relative_U2_Mass_PP = relative_uncertainty_contribution(Total_POC.U_Mass, Total_POC.Expe_POCuncer_PP);
relative_U2_Mass_experimental = median(Total_POC.relative_U2_Mass_PP)*100;
relative_U2_Mass_experimental_std = mad(Total_POC.relative_U2_Mass_PP)*100;
relative_U2_Mass_productive_experimental = median(Total_POC.relative_U2_Mass_PP(iP))*100;
relative_U2_Mass_productive_experimental_std = mad(Total_POC.relative_U2_Mass_PP(iP))*100;
relative_U2_Mass_mesopelagic_experimental = median(Total_POC.relative_U2_Mass_PP(iM))*100;
relative_U2_Mass_mesopelagic_experimental_std = mad(Total_POC.relative_U2_Mass_PP(iM))*100;

%% %% Calculating relative contribution of uncertainty in volume to the combined uncertainty of POC
% U2_Volume = uncertainty_volume( M_POC, Volume, sigma_volume)
Total_POC.U_volume = uncertainty_volume(Total_POC.MassPOC, Total_POC.Volumen, Total_POC.Comb_Unc_Volume);

Total_POC.relative_U2_volume = Total_POC.U_volume.^2 ./ Total_POC.U_POC.^2;
relative_U2_volume = median(Total_POC.relative_U2_volume)*100; % percent
relative_U2_volume_std = mad(Total_POC.relative_U2_volume)*100;
relative_U2_volume_productive = median(Total_POC.relative_U2_volume(iP))*100;
relative_U2_volume_productive_std = mad(Total_POC.relative_U2_volume(iP))*100;
relative_U2_volume_mesopelagic = median(Total_POC.relative_U2_volume(iM))*100;
relative_U2_volume_mesopelagic_std = mad(Total_POC.relative_U2_volume(iM))*100;

Total_POC.relative_U2_volume_PP = relative_uncertainty_contribution (Total_POC.U_volume, Total_POC.Expe_POCuncer_PP);
relative_U2_volume_experimental = median(Total_POC.relative_U2_volume_PP)*100;
relative_U2_volume_experimental_std = mad(Total_POC.relative_U2_volume_PP)*100;
relative_U2_volume_productive_experimental = median(Total_POC.relative_U2_volume_PP(iP))*100;
relative_U2_volume_productive_experimental_std = mad(Total_POC.relative_U2_volume_PP(iP))*100;
relative_U2_volume_mesopelagic_experimental = median(Total_POC.relative_U2_volume_PP(iM))*100;
relative_U2_volume_mesopelagic_experimental_std = mad(Total_POC.relative_U2_volume_PP(iM))*100;

%% Calculating relative contribution of uncertainty from handling to the combined uncertainty of POC

Total_POC.relative_U2_Handling = Total_POC.uncertainty_handling.^2 ./ Total_POC.U_POC.^2;
Total_POC.sigma_handling = Total_POC.uncertainty_handling .* Total_POC.Volumen; 
relative_U2_Handling_productive = median(Total_POC.relative_U2_Handling(iP));
relative_U2_Handling_mesopelagic = median(Total_POC.relative_U2_Handling(iM));

Total_POC.relative_U2_Handling_PP = relative_uncertainty_contribution(Total_POC.uncertainty_handling, Total_POC.Expe_POCuncer_PP);
relative_U2_Handling_experimental = median(Total_POC.relative_U2_Handling_PP)*100;
relative_U2_Handling_experimental_std = mad(Total_POC.relative_U2_Handling_PP)*100;
relative_U2_Handling_productive_experimental = median(Total_POC.relative_U2_Handling_PP(iP))*100;
relative_U2_Handling_productive_experimental_std = mad(Total_POC.relative_U2_Handling_PP(iP))*100;
relative_U2_Handling_mesopelagic_experimental = median(Total_POC.relative_U2_Handling_PP(iM))*100;
relative_U2_Handling_mesopelagic_experimental_std = mad(Total_POC.relative_U2_Handling_PP(iM))*100;

%% sum of relative uncertainties
Total_POC.sum_relative_uncertainties = (Total_POC.relative_U2_Mass + Total_POC.relative_U2_volume + Total_POC.relative_U2_Handling);
Total_POC.sum_relative_uncertainties_PP = (Total_POC.relative_U2_Mass_PP + Total_POC.relative_U2_volume_PP + Total_POC.relative_U2_Handling_PP);
