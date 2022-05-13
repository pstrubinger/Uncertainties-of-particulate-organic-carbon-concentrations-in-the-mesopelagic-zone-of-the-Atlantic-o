%Script for convert all individual structures of Samples, Blanks and Replicates into one total structure
clear all;
dataDir = "D:\Measuring POC_Paper\poc\POC\Data\";
TotalDir = "D:\Measuring POC_Paper\poc\POC\Data\Total\";
fn = dir(strcat(dataDir, "*mat"));
etopoFile = dir(strcat(dataDir, "*tif"));
duplicateDir = "D:\Measuring POC_Paper\poc\POC\Data\Duplicates\";
JossDir = 'D:\Measuring POC_Paper\poc\POC\POC_BLANKS\';
fn1 = dir(strcat(JossDir, '*mat'));
%fname = dir(strcat(duplicateDir, "*mat"));'New_Duplicates.mat'
fname = dir(strcat(duplicateDir, "New_Duplicates.mat"));
fname = fname.name;
figDir = "D:\Measuring POC_Paper\poc\POC\Figures\New\";
fname1 = dir(strcat(TotalDir, "Cap_filters_total.mat"));
fname1 = fname1.name;

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
Total_POC.corr = Tmp(2,1); %Value = 0.5111 (17/03/22)
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
iP = find (Total_POC.Zone == 1);%Value = 225 
iM = find (Total_POC.Zone == 2);%Value = 167 
Bad_runs = ["Results_11"; "Results_12"];
ibad = contains(Total_POC.RUN, Bad_runs);
igood = ~ibad;
iP_good = find (Total_POC.Zone == 1 & igood == 1);%Value = 190 
iP_bad = find (Total_POC.Zone == 1 & igood == 0);%Value = 35 
iM_good = find (Total_POC.Zone == 2 & igood == 1);%Value = 134 
iM_bad = find (Total_POC.Zone == 2 & igood == 0);%Value = 33 

%% Calculating experimental POC uncertaity using percentile precision

load(strcat(duplicateDir, fname)); %Loading data from duplicate analysis. 

% Total_POC.Expe_POCuncer_PP = nan(size(Total_POC.ID)); 
% Total_POC.Expe_POCuncer_PP(iP) = Sr_P*Total_POC.POC(iP);
% Total_POC.Expe_POCuncer_PP(iM) = Sr_M*Total_POC.POC(iM);
% 
% Total_POC.Expe_POCuncer_PP_relative = (Total_POC.Expe_POCuncer_PP * 100) ./ Total_POC.POC;
% Total_POC.Expe_POCuncer_PP_relative_median_Productive= median(Total_POC.Expe_POCuncer_PP_relative(iP_good)); %Value = 12.3841 (17/03/22)
% Total_POC.Expe_POCuncer_PP_relative_std_Productive= mad(Total_POC.Expe_POCuncer_PP_relative(iP_good),1); %Value = 0 (17/03/22)
% 
% Total_POC.Expe_POCuncer_PP_relative_median_Mesopelagic= median(Total_POC.Expe_POCuncer_PP_relative(iM_good)); %Value = 35.2351 (17/03/22)
% Total_POC.Expe_POCuncer_PP_relative_std_Mesopelagic= mad(Total_POC.Expe_POCuncer_PP_relative(iM_good),1); %Value = 0.5677 (17/03/22)
% 
% Total_POC.Expe_POCuncer_PP_median_Productive= median(Total_POC.Expe_POCuncer_PP(iP_good)); %Value = 2.2960 (17/03/22)
% Total_POC.Expe_POCuncer_PP_std_Productive= mad(Total_POC.Expe_POCuncer_PP(iP_good),1); %Value = 1.0748 (17/03/22)
% 
% Total_POC.Expe_POCuncer_PP_median_Mesopelagic= median(Total_POC.Expe_POCuncer_PP(iM_good)); %Value = 2.5916 (17/03/22)
% Total_POC.Expe_POCuncer_PP_std_Mesopelagic= mad(Total_POC.Expe_POCuncer_PP(iM_good),1); %Value = 0.5677 (17/03/22)

Total_POC.Expe_POCuncer_PP = nan(size(Total_POC.ID)); 
Total_POC.Expe_POCuncer_PP(iP) = Sr_P*Total_POC.POC(iP);
Total_POC.Expe_POCuncer_PP(iM) = Sr_M*Total_POC.POC(iM);

Total_POC.Expe_POCuncer_PP_relative = (Total_POC.Expe_POCuncer_PP * 100) ./ Total_POC.POC;
Total_POC.Expe_POCuncer_PP_relative_median_Productive= median(Total_POC.Expe_POCuncer_PP_relative(iP_good)); %Value = 12.3841 (13/05/22)
Total_POC.Expe_POCuncer_PP_relative_std_Productive= prcrng(Total_POC.Expe_POCuncer_PP_relative(iP_good)); %Value = 8.88 e-16 (13/05/22)

Total_POC.Expe_POCuncer_PP_relative_median_Mesopelagic= median(Total_POC.Expe_POCuncer_PP_relative(iM_good)); %Value = 35.2351 (13/05/22)
Total_POC.Expe_POCuncer_PP_relative_std_Mesopelagic= prcrng(Total_POC.Expe_POCuncer_PP_relative(iM_good)); %Value = 3.55e-15 (13/05/22)

Total_POC.Expe_POCuncer_PP_median_Productive= median(Total_POC.Expe_POCuncer_PP(iP_good)); %Value = 2.2960 (13/05/22)
Total_POC.Expe_POCuncer_PP_std_Productive= prcrng(Total_POC.Expe_POCuncer_PP(iP_good)); %Value = 2.4273 (13/05/22)

Total_POC.Expe_POCuncer_PP_median_Mesopelagic= median(Total_POC.Expe_POCuncer_PP(iM_good)); %Value = 2.5916 (13/05/22)
Total_POC.Expe_POCuncer_PP_std_Mesopelagic= prcrng(Total_POC.Expe_POCuncer_PP(iM_good)); %Value = 0.8989 (13/05/22)

%% Calculating relative contribution of uncertainty from calibration equation to the combined uncertainty of POC

%   Calculate the uncertainty in mass estimates predicted by the calibration equation after applying the standard law of propagation of uncertainty
%   to the measurement equiation POC = (M_uPOC - M_aDOC)
%   U_Mp = uncertainty_mass(PI_uPOC, PI_aDOC, r)
%   U_Mp = sqrt(PI_uPOC.^2 + PI_aDOC.^2 - 2 * PI_uPOC .* PI_aDOC * r); %[ug]

Total_POC.U_Mp = uncertainty_mass(Total_uPOC.PI, Total_aDOC.PI, Total_POC.corr); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculate the contributions of mass estimates predicted by the calibration equation uncertainty to the combined uncertainty of POC
%   U2_Mass = contribution_uncertainty_mass(sigma_mass, volume) 
%   U_Mass = U_Mp./volume; 
%   Sigma C(M) Eq.15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total_POC.U_Mass = contribution_uncertainty_mass(Total_POC.U_Mp, Total_POC.Volumen); 
% Total_POC.U_Mass_median= median(Total_POC.U_Mass); %Value = 0.4295 (17/03/22)
% Total_POC.U_Mass_std= mad(Total_POC.U_Mass,1); %Value = 0.2060 (17/03/22)
% 
% Total_POC.U_Mass_median_good= median(Total_POC.U_Mass(igood)); %Value = 0.3279 (17/03/22)
% Total_POC.U_Mass_std_good= mad(Total_POC.U_Mass(igood),1); %Value = 0.1175 (17/03/22)
% 
% Total_POC.U_Mass_median_good_Productive= median(Total_POC.U_Mass(iP_good)); %Value = 0.3679 (17/03/22)
% Total_POC.U_Mass_std_good_Productive= mad(Total_POC.U_Mass(iP_good),1); %Value = 0.1263 (17/03/22)
% 
% Total_POC.U_Mass_median_good_Mesopelagic= median(Total_POC.U_Mass(iM_good)); %Value = 0.2427 (17/03/22)
% Total_POC.U_Mass_std_good_Mesopelagic= mad(Total_POC.U_Mass(iM_good),1); %Value = 0.0799 (17/03/22)
% 
% Total_POC.U_Mass_median_bad= median(Total_POC.U_Mass(ibad)); %Value = 1.4477 (17/03/22)
% Total_POC.U_Mass_std_bad= mad(Total_POC.U_Mass(ibad),1); %Value = 0.4056 (17/03/22)

Total_POC.U_Mass = contribution_uncertainty_mass(Total_POC.U_Mp, Total_POC.Volumen); 
Total_POC.U_Mass_median= median(Total_POC.U_Mass); %Value = 0.4295 (17/03/22)
Total_POC.U_Mass_std= prcrng(Total_POC.U_Mass); %Value = 0.4134 (17/03/22)

Total_POC.U_Mass_median_good= median(Total_POC.U_Mass(igood)); %Value = 0.3279 (17/03/22)
Total_POC.U_Mass_std_good= prcrng(Total_POC.U_Mass(igood)); %Value = 0.2240 (17/03/22)

Total_POC.U_Mass_median_good_Productive= median(Total_POC.U_Mass(iP_good)); %Value = 0.3679 (17/03/22)
Total_POC.U_Mass_std_good_Productive= prcrng(Total_POC.U_Mass(iP_good)); %Value = 0.2415 (17/03/22)

Total_POC.U_Mass_median_good_Mesopelagic= median(Total_POC.U_Mass(iM_good)); %Value = 0.2427 (17/03/22)
Total_POC.U_Mass_std_good_Mesopelagic= prcrng(Total_POC.U_Mass(iM_good)); %Value = 0.1506 (17/03/22)

Total_POC.U_Mass_median_bad= median(Total_POC.U_Mass(ibad)); %Value = 1.4477 (17/03/22)
Total_POC.U_Mass_std_bad= prcrng(Total_POC.U_Mass(ibad)); %Value = 0.7294 (17/03/22)

%   Calculate the relative contributions of mass estimates predicted by the calibration equation uncertainty to the combined uncertainty of POC
% Total_POC.relative_U2_Mass = Total_POC.U_Mass.^2./ Total_POC.U_POC.^2; %U_POC is theoretical or calculated squared due to only account for our uncertainties
% Total_POC.relative_U2_Mass_median_good_productive = median(Total_POC.relative_U2_Mass(iP_good)); %Value = 0.9311 (17/03/22)
% Total_POC.relative_U2_Mass_std_good_productive = mad(Total_POC.relative_U2_Mass(iP_good),1); %Value = 0.0346 (17/03/22)
% Total_POC.relative_U2_Mass_median_good_mesopelagic = median(Total_POC.relative_U2_Mass(iM_good)); %Value = 0.9540 (17/03/22)
% Total_POC.relative_U2_Mass_std_good_mesopelagic = mad(Total_POC.relative_U2_Mass(iM_good),1); %Value = 0.0258 (17/03/22)

Total_POC.relative_U2_Mass = Total_POC.U_Mass.^2./ Total_POC.U_POC.^2; %U_POC is theoretical or calculated squared due to only account for our uncertainties
Total_POC.relative_U2_Mass_median_good_productive = median(Total_POC.relative_U2_Mass(iP_good)); %Value = 0.9311 (17/03/22)
Total_POC.relative_U2_Mass_std_good_productive = prcrng(Total_POC.relative_U2_Mass(iP_good)); %Value = 0.1030 (17/03/22)
Total_POC.relative_U2_Mass_median_good_mesopelagic = median(Total_POC.relative_U2_Mass(iM_good)); %Value = 0.9540 (17/03/22)
Total_POC.relative_U2_Mass_std_good_mesopelagic = prcrng(Total_POC.relative_U2_Mass(iM_good)); %Value = 0.0600 (17/03/22)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table.4 new manuscript (17/03/22)
%out = relative_uncertainty_contribution(U_source, U_POC)
% Total_POC.relative_U2_Mass_PP = relative_uncertainty_contribution(Total_POC.U_Mass, Total_POC.Expe_POCuncer_PP);
% Total_POC.relative_U2_Mass_experimental_median = median(Total_POC.relative_U2_Mass_PP);%Value = 0.1676 (17/03/22)
% Total_POC.relative_U2_Mass_experimental_std = mad(Total_POC.relative_U2_Mass_PP,1);%Value = 0.0846 (17/03/22)
% Total_POC.relative_U2_Mass_experimental_median_good = median(Total_POC.relative_U2_Mass_PP(igood));%Value = 0.1460 (17/03/22)
% Total_POC.relative_U2_Mass_experimental_std_good = mad(Total_POC.relative_U2_Mass_PP(igood),1);%Value = 0.0627 (17/03/22)
% Total_POC.relative_U2_Mass_experimental_median_bad = median(Total_POC.relative_U2_Mass_PP(ibad));%Value = 0.5752 (17/03/22)
% Total_POC.relative_U2_Mass_experimental_std_bad = mad(Total_POC.relative_U2_Mass_PP(ibad),1);%Value = 0.3205 (17/03/22)
% 
% Total_POC.relative_U2_Mass_experimental_median_good_productive = median(Total_POC.relative_U2_Mass_PP(iP_good));%Value = 0.1686 (17/03/22)
% Total_POC.relative_U2_Mass_experimental_std_good_productive = mad(Total_POC.relative_U2_Mass_PP(iP_good),1);%Value = 0.0641 (17/03/22)
% Total_POC.relative_U2_Mass_experimental_median_good_mesopelagic = median(Total_POC.relative_U2_Mass_PP(iM_good));%Value = 0.1002(17/03/22)
% Total_POC.relative_U2_Mass_experimental_std_good_mesopelagic = mad(Total_POC.relative_U2_Mass_PP(iM_good),1);%Value = 0.0411 (17/03/22)

Total_POC.relative_U2_Mass_PP = relative_uncertainty_contribution(Total_POC.U_Mass, Total_POC.Expe_POCuncer_PP);
Total_POC.relative_U2_Mass_experimental_median = median(Total_POC.relative_U2_Mass_PP);%Value = 0.1676 (17/03/22)
Total_POC.relative_U2_Mass_experimental_std = prcrng(Total_POC.relative_U2_Mass_PP);%Value = 0.1637 (17/03/22)
Total_POC.relative_U2_Mass_experimental_median_good = median(Total_POC.relative_U2_Mass_PP(igood));%Value = 0.1460 (17/03/22)
Total_POC.relative_U2_Mass_experimental_std_good = prcrng(Total_POC.relative_U2_Mass_PP(igood));%Value = 0.1079 (17/03/22)
Total_POC.relative_U2_Mass_experimental_median_bad = median(Total_POC.relative_U2_Mass_PP(ibad));%Value = 0.5752 (17/03/22)
Total_POC.relative_U2_Mass_experimental_std_bad = prcrng(Total_POC.relative_U2_Mass_PP(ibad));%Value = 0.4153 (17/03/22)

Total_POC.relative_U2_Mass_experimental_median_good_productive = median(Total_POC.relative_U2_Mass_PP(iP_good));%Value = 0.1686 (17/03/22)
Total_POC.relative_U2_Mass_experimental_std_good_productive = prcrng(Total_POC.relative_U2_Mass_PP(iP_good));%Value = 0.1294 (17/03/22)
Total_POC.relative_U2_Mass_experimental_median_good_mesopelagic = median(Total_POC.relative_U2_Mass_PP(iM_good));%Value = 0.1002(17/03/22)
Total_POC.relative_U2_Mass_experimental_std_good_mesopelagic = prcrng(Total_POC.relative_U2_Mass_PP(iM_good));%Value = 0.0781 (17/03/22)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table.3 new manuscript (21/03/22)
% Total_uPOC.PI_median_good_productive = median(Total_uPOC.PI(iP_good));%Value = 1.7247 (21/03/22)
% Total_uPOC.PI_std_good_productive = mad(Total_uPOC.PI(iP_good),1);%Value = 0.4455 (21/03/22)
% Total_uPOC.PI_median_good_mesopelagic = median(Total_uPOC.PI(iM_good));%Value = 1.6236 (21/03/22)
% Total_uPOC.PI_std_good_mesopelagic = mad(Total_uPOC.PI(iM_good),1);%Value = 0.5392 (21/03/22)
% 
% Total_aDOC.PI_median_good_productive = median(Total_aDOC.PI(iP_good));%Value = 1.6878 (21/03/22)
% Total_aDOC.PI_std_good_productive = mad(Total_aDOC.PI(iP_good),1);%Value = 0.4990 (21/03/22)
% Total_aDOC.PI_median_good_mesopelagic = median(Total_aDOC.PI(iM_good));%Value = 1.6115 (21/03/22)
% Total_aDOC.PI_std_good_mesopelagic = mad(Total_aDOC.PI(iM_good),1);%Value = 0.5346 (21/03/22)

Total_uPOC.PI_median_good_productive = median(Total_uPOC.PI(iP_good));%Value = 1.7247 (21/03/22)
Total_uPOC.PI_std_good_productive = prcrng(Total_uPOC.PI(iP_good));%Value = 0.9978 (21/03/22)
Total_uPOC.PI_median_good_mesopelagic = median(Total_uPOC.PI(iM_good));%Value = 1.6236 (21/03/22)
Total_uPOC.PI_std_good_mesopelagic = prcrng(Total_uPOC.PI(iM_good));%Value = 1.0752 (21/03/22)

Total_aDOC.PI_median_good_productive = median(Total_aDOC.PI(iP_good));%Value = 1.6878 (21/03/22)
Total_aDOC.PI_std_good_productive = prcrng(Total_aDOC.PI(iP_good));%Value = 0.9966 (21/03/22)
Total_aDOC.PI_median_good_mesopelagic = median(Total_aDOC.PI(iM_good));%Value = 1.6115 (21/03/22)
Total_aDOC.PI_std_good_mesopelagic = prcrng(Total_aDOC.PI(iM_good));%Value = 1.0723 (21/03/22)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% Calculating relative contribution of uncertainty in volume to the combined uncertainty of POC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table.3 new manuscript (17/03/22)
% Total_POC.Comb_Unc_Volume_median_good_productive = median(Total_POC.Comb_Unc_Volume(iP_good));%Value = 0.0141 (17/03/22)
% Total_POC.Comb_Unc_Volume_std_good_productive = mad(Total_POC.Comb_Unc_Volume(iP_good),1);%Value = 0.0032 (17/03/22)
% Total_POC.Comb_Unc_Volume_median_good_mesopelagic = median(Total_POC.Comb_Unc_Volume(iM_good));%Value = 0.0173 (17/03/22)
% Total_POC.Comb_Unc_Volume_std_good_mesopelagic = mad(Total_POC.Comb_Unc_Volume(iM_good),1);%Value = 0.0016 (17/03/22)

Total_POC.Comb_Unc_Volume_median_good_productive = median(Total_POC.Comb_Unc_Volume(iP_good));%Value = 0.0141 (17/03/22)
Total_POC.Comb_Unc_Volume_std_good_productive = prcrng(Total_POC.Comb_Unc_Volume(iP_good));%Value = 0.0037 (17/03/22)
Total_POC.Comb_Unc_Volume_median_good_mesopelagic = median(Total_POC.Comb_Unc_Volume(iM_good));%Value = 0.0173 (17/03/22)
Total_POC.Comb_Unc_Volume_std_good_mesopelagic = prcrng(Total_POC.Comb_Unc_Volume(iM_good));%Value = 0.0016 (17/03/22)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculate the uncertainty in volumen after applying the standard law of propagation of uncertainty
%   to the measurement equiation POC = M_POC / V
% U2_Volume = uncertainty_volume( M_POC, Volume, sigma_volume)
%U_Volume = M_POC.* sigma_volume ./ Volume.^2; %[ug * l / l^2] = [ug/l]
%   Sigma C(V) Eq.13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total_POC.U_volume = uncertainty_volume(Total_POC.MassPOC, Total_POC.Volumen, Total_POC.Comb_Unc_Volume);
% 
% Total_POC.U_volume_median= median(Total_POC.U_volume); %Value = 0.0279 (17/03/22)
% Total_POC.U_volume_std= mad(Total_POC.U_volume,1); %Value = 0.0136 (17/03/22)
% 
% Total_POC.U_volume_median_good= median(Total_POC.U_volume(igood)); %Value = 0.0285 (17/03/22)
% Total_POC.U_volume_std_good= mad(Total_POC.U_volume(igood),1); %Value = 0.0140 (17/03/22)
% 
% Total_POC.U_volume_median_good_Productive= median(Total_POC.U_volume(iP_good)); %Value = 0.0543 (17/03/22)
% Total_POC.U_volume_std_good_Productive= mad(Total_POC.U_volume(iP_good),1); %Value = 0.0294 (17/03/22)
% 
% Total_POC.U_volume_median_good_Mesopelagic= median(Total_POC.U_volume(iM_good)); %Value = 0.0184 (17/03/22)
% Total_POC.U_volume_std_good_Mesopelagic= mad(Total_POC.U_volume(iM_good),1); %Value = 0.0042 (17/03/22)
% 
% Total_POC.U_volume_median_bad= median(Total_POC.U_volume(ibad)); %Value = 0.0264 (17/03/22)
% Total_POC.U_volume_std_bad= mad(Total_POC.U_volume(ibad),1); %Value = 0.0142 (17/03/22)

Total_POC.U_volume = uncertainty_volume(Total_POC.MassPOC, Total_POC.Volumen, Total_POC.Comb_Unc_Volume);

Total_POC.U_volume_median= median(Total_POC.U_volume); %Value = 0.0279 (17/03/22)
Total_POC.U_volume_std= prcrng(Total_POC.U_volume); %Value = 0.0369 (17/03/22)

Total_POC.U_volume_median_good= median(Total_POC.U_volume(igood)); %Value = 0.0285 (17/03/22)
Total_POC.U_volume_std_good= prcrng(Total_POC.U_volume(igood)); %Value = 0.0388 (17/03/22)

Total_POC.U_volume_median_good_Productive= median(Total_POC.U_volume(iP_good)); %Value = 0.0543 (17/03/22)
Total_POC.U_volume_std_good_Productive= prcrng(Total_POC.U_volume(iP_good)); %Value = 0.0613 (17/03/22)

Total_POC.U_volume_median_good_Mesopelagic= median(Total_POC.U_volume(iM_good)); %Value = 0.0184 (17/03/22)
Total_POC.U_volume_std_good_Mesopelagic= prcrng(Total_POC.U_volume(iM_good)); %Value = 0.0064 (17/03/22)

Total_POC.U_volume_median_bad= median(Total_POC.U_volume(ibad)); %Value = 0.0264 (17/03/22)
Total_POC.U_volume_std_bad= prcrng(Total_POC.U_volume(ibad)); %Value = 0.0275 (17/03/22)

%   Calculate the relative contributions of volume to the combined uncertainty of POC
% Total_POC.relative_U2_volume = Total_POC.U_volume.^2./ Total_POC.U_POC.^2; %U_POC is theoretical or calculated squared duet to only account for our uncertainties
% Total_POC.relative_U2_volume_median_good_productive = median(Total_POC.relative_U2_volume(iP_good)); %Value = 0.0169 (17/03/22)
% Total_POC.relative_U2_volume_std_good_productive = mad(Total_POC.relative_U2_volume(iP_good),1); %Value = 0.0127 (17/03/22)
% Total_POC.relative_U2_volume_median_good_mesopelagic = median(Total_POC.relative_U2_volume(iM_good)); %Value = 0.0044 (17/03/22)
% Total_POC.relative_U2_volume_std_good_mesopelagic = mad(Total_POC.relative_U2_volume(iM_good),1); %Value = 0.0031 (17/03/22)

Total_POC.relative_U2_volume = Total_POC.U_volume.^2./ Total_POC.U_POC.^2; %U_POC is theoretical or calculated squared duet to only account for our uncertainties
Total_POC.relative_U2_volume_median_good_productive = median(Total_POC.relative_U2_volume(iP_good)); %Value = 0.0169 (17/03/22)
Total_POC.relative_U2_volume_std_good_productive = prcrng(Total_POC.relative_U2_volume(iP_good)); %Value = 0.0508 (17/03/22)
Total_POC.relative_U2_volume_median_good_mesopelagic = median(Total_POC.relative_U2_volume(iM_good)); %Value = 0.0044 (17/03/22)
Total_POC.relative_U2_volume_std_good_mesopelagic = prcrng(Total_POC.relative_U2_volume(iM_good)); %Value = 0.0067 (17/03/22)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table.4 new manuscript (17/03/22)

% Total_POC.relative_U2_volume_PP = relative_uncertainty_contribution(Total_POC.U_volume, Total_POC.Expe_POCuncer_PP);
% Total_POC.relative_U2_volume_experimental_median = median(Total_POC.relative_U2_volume_PP);%Value = 0.0210 (17/03/22)
% Total_POC.relative_U2_volume_experimental_std = mad(Total_POC.relative_U2_volume_PP,1);%Value = 0.0049 (17/03/22)
% Total_POC.relative_U2_volume_experimental_median_good = median(Total_POC.relative_U2_volume_PP(igood));%Value = 0.0211 (17/03/22)
% Total_POC.relative_U2_volume_experimental_std_good = mad(Total_POC.relative_U2_volume_PP(igood),1);%Value = 0.0049 (17/03/22)
% Total_POC.relative_U2_volume_experimental_median_bad = median(Total_POC.relative_U2_volume_PP(ibad));%Value = 0.0135 (17/03/22)
% Total_POC.relative_U2_volume_experimental_std_bad = mad(Total_POC.relative_U2_volume_PP(ibad),1);%Value = 0.0064 (17/03/22)
% 
% Total_POC.relative_U2_volume_experimental_median_good_productive = median(Total_POC.relative_U2_volume_PP(iP_good));%Value = 0.0257 (17/03/22)
% Total_POC.relative_U2_volume_experimental_std_good_productive = mad(Total_POC.relative_U2_volume_PP(iP_good),1);%Value = 0.0046 (17/03/22)
% Total_POC.relative_U2_volume_experimental_median_good_mesopelagic = median(Total_POC.relative_U2_volume_PP(iM_good));%Value = 0.0074(17/03/22)
% Total_POC.relative_U2_volume_experimental_std_good_mesopelagic = mad(Total_POC.relative_U2_volume_PP(iM_good),1);%Value = 0.0004 (17/03/22)

Total_POC.relative_U2_volume_PP = relative_uncertainty_contribution(Total_POC.U_volume, Total_POC.Expe_POCuncer_PP);
Total_POC.relative_U2_volume_experimental_median = median(Total_POC.relative_U2_volume_PP);%Value = 0.0210 (17/03/22)
Total_POC.relative_U2_volume_experimental_std = prcrng(Total_POC.relative_U2_volume_PP);%Value = 0.0093 (17/03/22)
Total_POC.relative_U2_volume_experimental_median_good = median(Total_POC.relative_U2_volume_PP(igood));%Value = 0.0211 (17/03/22)
Total_POC.relative_U2_volume_experimental_std_good = prcrng(Total_POC.relative_U2_volume_PP(igood));%Value = 0.0093 (17/03/22)
Total_POC.relative_U2_volume_experimental_median_bad = median(Total_POC.relative_U2_volume_PP(ibad));%Value = 0.0135 (17/03/22)
Total_POC.relative_U2_volume_experimental_std_bad = prcrng(Total_POC.relative_U2_volume_PP(ibad));%Value = 0.0093 (17/03/22)

Total_POC.relative_U2_volume_experimental_median_good_productive = median(Total_POC.relative_U2_volume_PP(iP_good));%Value = 0.0257 (17/03/22)
Total_POC.relative_U2_volume_experimental_std_good_productive = prcrng(Total_POC.relative_U2_volume_PP(iP_good));%Value = 0.0024 (17/03/22)
Total_POC.relative_U2_volume_experimental_median_good_mesopelagic = median(Total_POC.relative_U2_volume_PP(iM_good));%Value = 0.0074(17/03/22)
Total_POC.relative_U2_volume_experimental_std_good_mesopelagic = prcrng(Total_POC.relative_U2_volume_PP(iM_good));%Value = 6.66 e-4 (17/03/22)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Calculating relative contribution of uncertainty from handling to the combined uncertainty of POC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table.3 new manuscript (21/03/22)

%Total_POC.sigma_handling = Total_POC.uncertainty_handling .* Total_POC.Volumen; old equation
% Total_POC.sigma_handling = (Total_POC.uncertainty_handling .* Total_POC.Volumen)./sqrt(2); %New equation (21/03/22)
% Total_POC.sigma_handling_median_good_productive = median(Total_POC.sigma_handling(iP_good));%Value = 0.2729 (21/03/22)
% Total_POC.sigma_handling_std_good_productive = mad(Total_POC.sigma_handling(iP_good),1);%Value = 0.1336 (21/03/22)
% Total_POC.sigma_handling_median_good_mesopelagic = median(Total_POC.sigma_handling(iM_good));%Value = 0.2398 (21/03/22)
% Total_POC.sigma_handling_std_good_mesopelagic = mad(Total_POC.sigma_handling(iM_good),1);%Value = 0.1031 (21/03/22)
Total_POC.sigma_handling = (Total_POC.uncertainty_handling .* Total_POC.Volumen)./sqrt(2); %New equation (21/03/22)
Total_POC.sigma_handling_median_good_productive = median(Total_POC.sigma_handling(iP_good));%Value = 0.2729 (21/03/22)
Total_POC.sigma_handling_std_good_productive = prcrng(Total_POC.sigma_handling(iP_good));%Value = 0.2043 (21/03/22)
Total_POC.sigma_handling_median_good_mesopelagic = median(Total_POC.sigma_handling(iM_good));%Value = 0.2398 (21/03/22)
Total_POC.sigma_handling_std_good_mesopelagic = prcrng(Total_POC.sigma_handling(iM_good));%Value = 0.1773 (21/03/22)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculate the uncertainty in handling after applying the standard law of propagation of uncertainty
%   to the measurement equiation POC = M_POC / V
% std_mean : the standard error of the mean of the three corresponding estimates of acidified and non-acidified filter blanks
% volume   : the volume of seawater filtered; in litres
% U_Handling = sqrt(std_mean.^2 + std_mean.^2)./ volume; % [ug/l]
%   Sigma C(n) Eq.16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Total_POC.uncertainty_handling_median= median(Total_POC.uncertainty_handling); %Value = 0.0676 (21/03/22)
% Total_POC.uncertainty_handling_std= mad(Total_POC.uncertainty_handling,1); %Value = 0.0347 (21/03/22)
% 
% Total_POC.uncertainty_handling_median_good= median(Total_POC.uncertainty_handling(igood)); %Value = 0.0668 (21/03/22)
% Total_POC.uncertainty_handling_std_good= mad(Total_POC.uncertainty_handling(igood),1); %Value = 0.0331 (21/03/22)
% 
% Total_POC.uncertainty_handling_median_good_Productive= median(Total_POC.uncertainty_handling(iP_good)); %Value = 0.0800 (21/03/22)
% Total_POC.uncertainty_handling_std_good_Productive= mad(Total_POC.uncertainty_handling(iP_good),1); %Value = 0.0418 (21/03/22)
% 
% Total_POC.uncertainty_handling_median_good_Mesopelagic= median(Total_POC.uncertainty_handling(iM_good)); %Value = 0.0512 (21/03/22)
% Total_POC.uncertainty_handling_std_good_Mesopelagic= mad(Total_POC.uncertainty_handling(iM_good),1); %Value = 0.0223 (21/03/22)
% 
% Total_POC.uncertainty_handling_median_bad= median(Total_POC.uncertainty_handling(ibad)); %Value = 0.0749 (21/03/22)
% Total_POC.uncertainty_handling_std_bad= mad(Total_POC.uncertainty_handling(ibad),1); %Value = 0.0477 (21/03/22)

Total_POC.uncertainty_handling_median= median(Total_POC.uncertainty_handling); %Value = 0.0676 (21/03/22)
Total_POC.uncertainty_handling_std= prcrng(Total_POC.uncertainty_handling); %Value = 0.0597 (21/03/22)

Total_POC.uncertainty_handling_median_good= median(Total_POC.uncertainty_handling(igood)); %Value = 0.0668 (21/03/22)
Total_POC.uncertainty_handling_std_good= prcrng(Total_POC.uncertainty_handling(igood)); %Value = 0.0581 (21/03/22)

Total_POC.uncertainty_handling_median_good_Productive= median(Total_POC.uncertainty_handling(iP_good)); %Value = 0.0800 (21/03/22)
Total_POC.uncertainty_handling_std_good_Productive= prcrng(Total_POC.uncertainty_handling(iP_good)); %Value = 0.0746 (21/03/22)

Total_POC.uncertainty_handling_median_good_Mesopelagic= median(Total_POC.uncertainty_handling(iM_good)); %Value = 0.0512 (21/03/22)
Total_POC.uncertainty_handling_std_good_Mesopelagic= prcrng(Total_POC.uncertainty_handling(iM_good)); %Value = 0.0438 (21/03/22)

Total_POC.uncertainty_handling_median_bad= median(Total_POC.uncertainty_handling(ibad)); %Value = 0.0749 (21/03/22)
Total_POC.uncertainty_handling_std_bad= prcrng(Total_POC.uncertainty_handling(ibad)); %Value = 0.0742 (21/03/22)

%   Calculate the relative contributions of volume to the combined uncertainty of POC
% Total_POC.relative_U2_Handling = Total_POC.uncertainty_handling.^2./ Total_POC.U_POC.^2; %U_POC is theoretical or calculated squared duet to only account for our uncertainties
% Total_POC.relative_U2_Handling_median_good_productive = median(Total_POC.relative_U2_Handling(iP_good)); %Value = 0.0413 (21/03/22)
% Total_POC.relative_U2_Handling_std_good_productive = mad(Total_POC.relative_U2_Handling(iP_good),1); %Value = 0.0238 (21/03/22)
% Total_POC.relative_U2_Handling_median_good_mesopelagic = median(Total_POC.relative_U2_Handling(iM_good)); %Value = 0.0376 (21/03/22)
% Total_POC.relative_U2_Handling_std_good_mesopelagic = mad(Total_POC.relative_U2_Handling(iM_good),1); %Value = 0.0223 (21/03/22)

Total_POC.relative_U2_Handling = Total_POC.uncertainty_handling.^2./ Total_POC.U_POC.^2; %U_POC is theoretical or calculated squared duet to only account for our uncertainties
Total_POC.relative_U2_Handling_median_good_productive = median(Total_POC.relative_U2_Handling(iP_good)); %Value = 0.0413 (21/03/22)
Total_POC.relative_U2_Handling_std_good_productive = prcrng(Total_POC.relative_U2_Handling(iP_good)); %Value = 0.0405 (21/03/22)
Total_POC.relative_U2_Handling_median_good_mesopelagic = median(Total_POC.relative_U2_Handling(iM_good)); %Value = 0.0376 (21/03/22)
Total_POC.relative_U2_Handling_std_good_mesopelagic = prcrng(Total_POC.relative_U2_Handling(iM_good)); %Value = 0.0493 (21/03/22)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table.4 new manuscript (21/03/22)

% Total_POC.relative_U2_Handling_PP = relative_uncertainty_contribution(Total_POC.uncertainty_handling, Total_POC.Expe_POCuncer_PP);
% Total_POC.relative_U2_Handling_experimental_median = median(Total_POC.relative_U2_Handling_PP);%Value = 0.0283 (21/03/22)
% Total_POC.relative_U2_Handling_experimental_std = mad(Total_POC.relative_U2_Handling_PP,1);%Value = 0.0148 (21/03/22)
% Total_POC.relative_U2_Handling_experimental_median_good = median(Total_POC.relative_U2_Handling_PP(igood));%Value = 0.0272 (21/03/22)
% Total_POC.relative_U2_Handling_experimental_std_good = mad(Total_POC.relative_U2_Handling_PP(igood),1);%Value = 0.0142 (21/03/22)
% Total_POC.relative_U2_Handling_experimental_median_bad = median(Total_POC.relative_U2_Handling_PP(ibad));%Value = 0.0335 (21/03/22)
% Total_POC.relative_U2_Handling_experimental_std_bad = mad(Total_POC.relative_U2_Handling_PP(ibad),1);%Value = 0.0181 (21/03/22)
% 
% Total_POC.relative_U2_Handling_experimental_median_good_productive = median(Total_POC.relative_U2_Handling_PP(iP_good));%Value = 0.0333 (21/03/22)
% Total_POC.relative_U2_Handling_experimental_std_good_productive = mad(Total_POC.relative_U2_Handling_PP(iP_good),1);%Value = 0.0146 (21/03/22)
% Total_POC.relative_U2_Handling_experimental_median_good_mesopelagic = median(Total_POC.relative_U2_Handling_PP(iM_good));%Value = 0.0200(21/03/22)
% Total_POC.relative_U2_Handling_experimental_std_good_mesopelagic = mad(Total_POC.relative_U2_Handling_PP(iM_good),1);%Value = 0.0103 (21/03/22)
Total_POC.relative_U2_Handling_PP = relative_uncertainty_contribution(Total_POC.uncertainty_handling, Total_POC.Expe_POCuncer_PP);
Total_POC.relative_U2_Handling_experimental_median = median(Total_POC.relative_U2_Handling_PP);%Value = 0.0283 (21/03/22)
Total_POC.relative_U2_Handling_experimental_std = prcrng(Total_POC.relative_U2_Handling_PP);%Value = 0.0279 (21/03/22)
Total_POC.relative_U2_Handling_experimental_median_good = median(Total_POC.relative_U2_Handling_PP(igood));%Value = 0.0272 (21/03/22)
Total_POC.relative_U2_Handling_experimental_std_good = prcrng(Total_POC.relative_U2_Handling_PP(igood));%Value = 0.0267 (21/03/22)
Total_POC.relative_U2_Handling_experimental_median_bad = median(Total_POC.relative_U2_Handling_PP(ibad));%Value = 0.0335 (21/03/22)
Total_POC.relative_U2_Handling_experimental_std_bad = prcrng(Total_POC.relative_U2_Handling_PP(ibad));%Value = 0.0295 (21/03/22)

Total_POC.relative_U2_Handling_experimental_median_good_productive = median(Total_POC.relative_U2_Handling_PP(iP_good));%Value = 0.0333 (21/03/22)
Total_POC.relative_U2_Handling_experimental_std_good_productive = prcrng(Total_POC.relative_U2_Handling_PP(iP_good));%Value = 0.0291 (21/03/22)
Total_POC.relative_U2_Handling_experimental_median_good_mesopelagic = median(Total_POC.relative_U2_Handling_PP(iM_good));%Value = 0.0200(21/03/22)
Total_POC.relative_U2_Handling_experimental_std_good_mesopelagic = prcrng(Total_POC.relative_U2_Handling_PP(iM_good));%Value = 0.0204 (21/03/22)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% sum of relative uncertainties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table.4 new manuscript (21/03/22)

% Total_POC.sum_relative_uncertainties = (Total_POC.relative_U2_Mass + Total_POC.relative_U2_volume + Total_POC.relative_U2_Handling);
% Total_POC.sum_relative_uncertainties_PP = (Total_POC.relative_U2_Mass_PP + Total_POC.relative_U2_volume_PP + Total_POC.relative_U2_Handling_PP);
% Total_POC.sum_relative_uncertainties_PP_median = median(Total_POC.sum_relative_uncertainties_PP(igood));%Value = 0.1945 (21/03/22)
% Total_POC.sum_relative_uncertainties_PP_std = mad(Total_POC.sum_relative_uncertainties_PP(igood),1);%Value = 0.0793 (21/03/22)
% 
% Total_POC.sum_relative_uncertainties_PP_median_good_productive = median(Total_POC.sum_relative_uncertainties_PP(iP_good));%Value = 0.2318 (21/03/22)
% Total_POC.sum_relative_uncertainties_PP_std_good_productive = mad(Total_POC.sum_relative_uncertainties_PP(iP_good),1);%Value = 0.0676 (21/03/22)
% Total_POC.sum_relative_uncertainties_PP_median_good_mesopelagic = median(Total_POC.sum_relative_uncertainties_PP(iM_good));%Value = 0.1224 (21/03/22)
% Total_POC.sum_relative_uncertainties_PP_std_good_mesopelagic = mad(Total_POC.sum_relative_uncertainties_PP(iM_good),1);%Value = 0.0462 (21/03/22)
% 
% Total_POC.Unquantified_relative_uncertainties_PP = 1 - Total_POC.sum_relative_uncertainties_PP;
% Total_POC.Unquantified_relative_uncertainties_PP_median_good_productive = median(Total_POC.Unquantified_relative_uncertainties_PP(iP_good));%Value = 0.7682 (21/03/22)
% Total_POC.Unquantified_relative_uncertainties_PP_std_good_productive = mad(Total_POC.Unquantified_relative_uncertainties_PP(iP_good),1);%Value = 0.0676 (21/03/22)
% Total_POC.Unquantified_relative_uncertainties_PP_median_good_mesopelagic = median(Total_POC.Unquantified_relative_uncertainties_PP(iM_good));%Value = 0.8746 (21/03/22)
% Total_POC.Unquantified_relative_uncertainties_PP_std_good_mesopelagic = mad(Total_POC.Unquantified_relative_uncertainties_PP(iM_good),1);%Value = 0.0462 (21/03/22)

Total_POC.sum_relative_uncertainties = (Total_POC.relative_U2_Mass + Total_POC.relative_U2_volume + Total_POC.relative_U2_Handling);
Total_POC.sum_relative_uncertainties_PP = (Total_POC.relative_U2_Mass_PP + Total_POC.relative_U2_volume_PP + Total_POC.relative_U2_Handling_PP);
Total_POC.sum_relative_uncertainties_PP_median = median(Total_POC.sum_relative_uncertainties_PP(igood));%Value = 0.1945 (21/03/22)
Total_POC.sum_relative_uncertainties_PP_std = prcrng(Total_POC.sum_relative_uncertainties_PP(igood));%Value = 0.1296 (21/03/22)

Total_POC.sum_relative_uncertainties_PP_median_good_productive = median(Total_POC.sum_relative_uncertainties_PP(iP_good));%Value = 0.2318 (21/03/22)
Total_POC.sum_relative_uncertainties_PP_std_good_productive = prcrng(Total_POC.sum_relative_uncertainties_PP(iP_good));%Value = 0.1336 (21/03/22)
Total_POC.sum_relative_uncertainties_PP_median_good_mesopelagic = median(Total_POC.sum_relative_uncertainties_PP(iM_good));%Value = 0.1224 (21/03/22)
Total_POC.sum_relative_uncertainties_PP_std_good_mesopelagic = prcrng(Total_POC.sum_relative_uncertainties_PP(iM_good));%Value = 0.0869 (21/03/22)

Total_POC.Unquantified_relative_uncertainties_PP = 1 - Total_POC.sum_relative_uncertainties_PP;
Total_POC.Unquantified_relative_uncertainties_PP_median_good_productive = median(Total_POC.Unquantified_relative_uncertainties_PP(iP_good));%Value = 0.7682 (21/03/22)
Total_POC.Unquantified_relative_uncertainties_PP_std_good_productive = prcrng(Total_POC.Unquantified_relative_uncertainties_PP(iP_good));%Value = 0.1336 (21/03/22)
Total_POC.Unquantified_relative_uncertainties_PP_median_good_mesopelagic = median(Total_POC.Unquantified_relative_uncertainties_PP(iM_good));%Value = 0.8746 (21/03/22)
Total_POC.Unquantified_relative_uncertainties_PP_std_good_mesopelagic = prcrng(Total_POC.Unquantified_relative_uncertainties_PP(iM_good));%Value = 0.0869 (21/03/22)

% Table.4 new manuscript (21/03/22) IF THE RATIO IS SQUARED!!!!!!!!!!!!!!!
% A = Total_POC.U_Mass.^2./ Total_POC.Expe_POCuncer_PP.^2; 
% B = Total_POC.U_volume.^2./ Total_POC.Expe_POCuncer_PP.^2; 
% C = Total_POC.uncertainty_handling.^2./ Total_POC.Expe_POCuncer_PP.^2; 
% 
% D = (A + B + C);
% D_median = median(D);%Value = 0.0304 (3.04%)(21/03/22)
% D_std = mad(D);%Value = 0.2710 (27.10%)(21/03/22)
% 
% D_median_good_productive = median(D(iP_good));%Value = 0.0308 (3.1%) (21/03/22)
% D_std_good_productive = mad(D(iP_good),1);%Value = 0.0191 (19.1%) (21/03/22)
% D_median_good_mesopelagic = median(D(iM_good));%Value = 0.0105 (1.1%) (21/03/22)
% D_std_good_mesopelagic = mad(D(iM_good),1);%Value = 0.0077 (0.8%) (21/03/22)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Creating a new structure called Paul_Results 

filename1 = fn1.name;
load(strcat(JossDir, filename1));

field1 = 'RUN';
value1 = Total_POC.RUN;
field2 = 'ID';
value2 = Total_POC.ID;
field3 = 'Date';
value3 = NaT(392,1,'Format','dd/MM/yyyy');
field4 = 'Time';
value4 = NaT(392,1,'Format','HH:mm');
field5 = 'CTD';
value5 = NaN(392,1);
field6 = 'Lat'; 
value6 = NaN(392,1);
field7 = 'Lon'; 
value7 = NaN(392,1);
field8 = 'Depth';
value8 = Total_POC.Depth;
field9 = 'Volume'; 
value9 = Total_POC.Volumen;
field10 = 'POC'; 
value10 = Total_POC.POC;
field11 = 'U_POC';
value11 = Total_POC.U_POC;
field12 = 'U_POC_Experimental_PP';
value12 = Total_POC.Expe_POCuncer_PP;
field13 = 'Sum_U_POC_Experimental_PP';
value13 = Total_POC.sum_relative_uncertainties_PP;
field14 = 'Mass_aDOC';
value14 = Total_aDOC.TrueCarbonMass;
field15 = 'Mass_uPOC';
value15 = Total_uPOC.TrueCarbonMass;


Paul_Results = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,...
    field7,value7,field8,value8,field9,value9,field10,value10,field11,value11,field12,value12,field13,value13,field14,value14,field15,value15);

clearvars field1 value1 field2 value2 field3 value3 field4 value4 field5 value5 field6 value6...
    field7 value7 field8 value8 field9 value9 field10 value10 field11 value11 field12 value12 field13 value13 field14 value14 field15 value15

Paul_Results.Date(2:end) = JROSS_results.Date(1:391);
Paul_Results.Time(2:end) = JROSS_results.Time(1:391);
Paul_Results.CTD(2:end) = JROSS_results.CTD(1:391);
Paul_Results.Lon(2:end) = JROSS_results.Lon(1:391);
Paul_Results.Lat(2:end) = JROSS_results.Lat(1:391);
Paul_Results.Zone = Total_POC.Zone;

%% Plotting latitudinal section of uPOC, aDOC, POC, U_POC_Experimental_PP, and Sum_U_POC_Experimental_PP

Paul_Results.Sum_U_POC_Experimental_PP_percent = Paul_Results.Sum_U_POC_Experimental_PP * 100;  

f2 = figure (2); % uPOC mass
clf
section_plot(Paul_Results.Lat(igood), Paul_Results.Depth(igood), Paul_Results.Mass_uPOC(igood), 'uPOC (\mug)', [0 250]) % plotting using the function section_plot 
% c.Label.String = ('uPOC (\mug)');

f3= figure (3); % aDOC mass
clf
section_plot(Paul_Results.Lat(igood), Paul_Results.Depth(igood), Paul_Results.Mass_aDOC(igood), 'aDOC (\mug)', [0 30]) % plotting using the function section_plot 
% c.Label.String = ('aDOC (\mug)'); 

f4 = figure (4); % POC
clf
section_plot(Paul_Results.Lat(igood), Paul_Results.Depth(igood), Paul_Results.POC(igood), 'POC (mg m^{-3})', [0 50]) % plotting using the function section_plot 
% c.Label.String = ('POC (mg m^-3)'); 

f5 = figure (5); %U_POC_Experimental_PP
clc
section_plot(Paul_Results.Lat(igood), Paul_Results.Depth(igood), Paul_Results.U_POC_Experimental_PP(igood), 'POC Total Uncertainty (mg m^{-3})', [0 10]);
% c.Label.String = ('POC Total Uncertainty (mg m^-3)'); 

f6 = figure (6); %Sum_U_POC_Experimental_PP for good
clc
section_plot(Paul_Results.Lat(igood), Paul_Results.Depth(igood), Paul_Results.Sum_U_POC_Experimental_PP_percent(igood), 'Total Uncertainty Explained (%)', [0 50]);
% c.Label.String = ('Total Uncertainty Explained (%)');

% f6 = figure (6); %Sum_U_POC_Experimental_PP for all including bad
% clc
% section_plot(Paul_Results.Lat, Paul_Results.Depth , Paul_Results.Sum_U_POC_Experimental_PP_percent, 'Total Uncertainty Explained (%)', [0 50]);
% % c.Label.String = ('Total Uncertainty Explained (%)');

f7 = figure (7); %relative_U2_Mass_PP_percent
clc
section_plot(Paul_Results.Lat(igood), Paul_Results.Depth(igood), (Total_POC.relative_U2_Mass_PP(igood) * 100), 'Contribution to Total Uncertainty (%)',[0 40]);
% c.Label.String = ('Contribution to Total Uncertainty (%)');

f8 = figure (8); % POC in aDOC
clf
section_plot(Paul_Results.Lat(igood), Paul_Results.Depth(igood), Total_aDOC.POC(igood), 'aDOC (mg m^{-3})', [0 10]) % plotting using the function section_plot 
% c.Label.String = ('aDOC (mg m^-3)'); 

% f8 = figure (8); % POC
% clf
% section_plot(Paul_Results.Lat, Paul_Results.Depth , Total_aDOC.aDOC_concentration, 'aDOC (mg m^{-3})', [0 10]) % plotting using the function section_plot 
% % c.Label.String = ('aDOC (mg m^-3)'); 

f9 = figure (9); % POC in uPOC
clf
section_plot(Paul_Results.Lat(igood), Paul_Results.Depth(igood), Total_uPOC.POC(igood), 'uPOC (mg m^{-3})', [0 50]) % plotting using the function section_plot 
% c.Label.String = ('POC (mg m^-3)'); 

% f9 = figure (9); % POC in uPOC
% clf
% section_plot(Paul_Results.Lat(igood), Paul_Results.Depth(igood), Total_uPOC.uPOC_concentration, 'uPOC (mg m^{-3})', [0 50]) % plotting using the function section_plot 
% % c.Label.String = ('POC (mg m^-3)'); 


%% Saving figures
f2name = 'New_uPOC_mass_Section';
print (f2, '-dpng', strcat(figDir,f2name), '-r350');
print (f2, '-depsc', strcat(figDir,f2name));

f3name = 'New_aDOC_mass_Section';
print (f3, '-dpng', strcat(figDir,f3name), '-r350');
print (f3, '-depsc', strcat(figDir,f3name));

f4name = 'New_POC_concentration_Section';
print (f4, '-dpng', strcat(figDir,f4name), '-r350');
print (f4, '-depsc', strcat(figDir,f4name));

f5name = 'New_U_POC_Experimental_PP_Section';
print (f5, '-dpng', strcat(figDir,f5name), '-r350');
print (f5, '-depsc', strcat(figDir,f5name));

f6name = 'New_Sum_U_POC_Experimental_PP_Section';
print (f6, '-dpng', strcat(figDir,f6name), '-r350');
print (f6, '-depsc', strcat(figDir,f6name));

f7name = 'New_relative_U2_Mass_PP_percent_Section';
print (f7, '-dpng', strcat(figDir,f7name), '-r350');
print (f7, '-depsc', strcat(figDir,f7name));

% f6name = 'Sum_U_POC_Experimental_PP_Section_all';
% print (f6, '-dpng', strcat(figDir,f6name), '-r350');
% print (f6, '-depsc', strcat(figDir,f6name));
% f7name = 'relative_U2_Mass_PP_percent_Section_all';
% print (f7, '-dpng', strcat(figDir,f7name), '-r350');
% print (f7, '-depsc', strcat(figDir,f7name));

f8name = 'New_aDOC_concentration_Section';
print (f8, '-dpng', strcat(figDir,f8name), '-r350');
print (f8, '-depsc', strcat(figDir,f8name));

f9name = 'New_uPOC_concentration_Section';
print (f9, '-dpng', strcat(figDir,f9name), '-r350');
print (f9, '-depsc', strcat(figDir,f9name));

%% POC distribution across the sampled biogeographical provinces

NADR = 43.51;   %province = 1          
NAST= 25.50;    %province = 2
NATL = 10.0;    %province = 3
WTRA = -5.65;   %province = 4
SATL = -35.0;   %province = 5
% SSTC > -35    %province = 6

for iprovince = 1:length(Paul_Results.ID)
       if Paul_Results.Lat(iprovince) >= NADR 
          Paul_Results.Province(iprovince) = 1;
       elseif Paul_Results.Lat(iprovince) < NADR && Paul_Results.Lat(iprovince)>= NAST
          Paul_Results.Province(iprovince) = 2;
       elseif Paul_Results.Lat(iprovince) < NAST && Paul_Results.Lat(iprovince)>= NATL
          Paul_Results.Province(iprovince) = 3;
       elseif Paul_Results.Lat(iprovince) < NATL && Paul_Results.Lat(iprovince)>= WTRA
          Paul_Results.Province(iprovince) = 4;
       elseif Paul_Results.Lat(iprovince) < WTRA && Paul_Results.Lat(iprovince)>= SATL
          Paul_Results.Province(iprovince) = 5;   
       else
          Paul_Results.Province(iprovince) = 6;
       end
end

Paul_Results.Province = Paul_Results.Province';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table.2 new manuscript (24/03/22) / Results 3.1
iNADRP = find (Paul_Results.Zone == 1 & Paul_Results.Province == 1 & igood == 1);%Value = 19 (24/03/22)
iNADRM = find (Paul_Results.Zone == 2 & Paul_Results.Province == 1 & igood == 1);%Value = 10 (24/03/22)

iNASTP = find (Paul_Results.Zone == 1 & Paul_Results.Province == 2 & igood == 1);%Value = 7 (24/03/22)
iNASTM = find (Paul_Results.Zone == 2 & Paul_Results.Province == 2 & igood == 1);%Value = 6 (24/03/22)

iNATLP = find (Paul_Results.Zone == 1 & Paul_Results.Province == 3 & igood == 1);%Value = 36 (24/03/22)
iNATLM = find (Paul_Results.Zone == 2 & Paul_Results.Province == 3 & igood == 1);%Value = 18 (24/03/22)

iWTRAP = find (Paul_Results.Zone == 1 & Paul_Results.Province == 4 & igood == 1);%Value = 32 (24/03/22)
iWTRAM = find (Paul_Results.Zone == 2 & Paul_Results.Province == 4 & igood == 1);%Value = 16 (24/03/22)

iSATLP = find (Paul_Results.Zone == 1 & Paul_Results.Province == 5 & igood == 1);%Value = 57 (24/03/22)
iSATLM = find (Paul_Results.Zone == 2 & Paul_Results.Province == 5 & igood == 1);%Value = 44 (24/03/22)

iSSTCP = find (Paul_Results.Zone == 1 & Paul_Results.Province == 6 & igood == 1);%Value = 39 (24/03/22)
iSSTCM = find (Paul_Results.Zone == 2 & Paul_Results.Province == 6 & igood == 1);%Value = 40 (24/03/22)

% Paul_Results.Median_NADRP = median(Paul_Results.POC(iNADRP),'omitnan');%Value = 48.9451 (24/03/22)
% Paul_Results.RSTD_NADRP = mad(Paul_Results.POC(iNADRP),1);%Value = 19.1521 (24/03/22)
% Paul_Results.Median_NADRM = median(Paul_Results.POC(iNADRM));%Value = 7.2218 (24/03/22)
% Paul_Results.RSTD_NADRM = mad(Paul_Results.POC(iNADRM),1);%Value = 2.0674 (24/03/22)
% 
% Paul_Results.Median_NASTP = median(Paul_Results.POC(iNASTP));%Value = 26.7537 (24/03/22)
% Paul_Results.RSTD_NASTP = mad(Paul_Results.POC(iNASTP,1));%Value = 8.9462 (24/03/22)
% Paul_Results.Median_NASTM = median(Paul_Results.POC(iNASTM));%Value = 5.6144 (24/03/22)
% Paul_Results.RSTD_NASTM = mad(Paul_Results.POC(iNASTM),1);%Value = 0.2148 (24/03/22)
% 
% Paul_Results.Median_NATLP = median(Paul_Results.POC(iNATLP));%Value = 14.5174 (24/03/22)
% Paul_Results.RSTD_NATLP = mad(Paul_Results.POC(iNATLP),1);%Value =5.9573 (24/03/22)
% Paul_Results.Median_NATLM = median(Paul_Results.POC(iNATLM),'omitnan');%Value = 7.8594 (24/03/22)
% Paul_Results.RSTD_NATLM = mad(Paul_Results.POC(iNATLM),1);%Value = 1.0324 (24/03/22)
% 
% Paul_Results.Median_WTRAP = median(Paul_Results.POC(iWTRAP));%Value = 19.4882 (24/03/22)
% Paul_Results.RSTD_WTRAP = mad(Paul_Results.POC(iWTRAP),1);%Value = 8.5700 (24/03/22)
% Paul_Results.Median_WTRAM = median(Paul_Results.POC(iWTRAM));%Value = 8.7503 (24/03/22)
% Paul_Results.RSTD_WTRAM = mad(Paul_Results.POC(iWTRAM),1);%Value = 2.1121 (24/03/22)
% 
% Paul_Results.Median_SATLP = median(Paul_Results.POC(iSATLP));%Value = 13.9178 (24/03/22)
% Paul_Results.RSTD_SATLP = mad(Paul_Results.POC(iSATLP),1);%Value = 4.3155 (24/03/22)
% Paul_Results.Median_SATLM = median(Paul_Results.POC(iSATLM));%Value = 6.7016 (24/03/22)
% Paul_Results.RSTD_SATLM = mad(Paul_Results.POC(iSATLM),1);%Value = 1.1330 (24/03/22)
% 
% Paul_Results.Median_SSTCP = median(Paul_Results.POC(iSSTCP));%Value = 46.0122 (24/03/22)
% Paul_Results.RSTD_SSTCP = mad(Paul_Results.POC(iSSTCP),1);%Value = 14.5583 (24/03/22)
% Paul_Results.Median_SSTCM = median(Paul_Results.POC(iSSTCM));%Value = 8.6542 (24/03/22)
% Paul_Results.RSTD_SSTCM = mad(Paul_Results.POC(iSSTCM),1);%Value = 1.8256 (24/03/22)
% 
% Paul_Results.Median_POC_good_productive = median(Paul_Results.POC(iP_good),'omitnan');%Value = 18.5402 (24/03/22)
% Paul_Results.RSTD_POC__good_productive = mad(Paul_Results.POC(iP_good),1);%Value = 8.6791 (24/03/22)
% Paul_Results.Median_POC_good_mesopelagic = median(Paul_Results.POC(iM_good),'omitnan'); %Value = 7.3550 (24/03/22)
% Paul_Results.RSTD_POC_good_mesopelagic = mad(Paul_Results.POC(iM_good),1);%Value = 1.6111 (24/03/22)

Paul_Results.Median_NADRP = median(Paul_Results.POC(iNADRP),'omitnan');%Value = 48.9451 (24/03/22)
Paul_Results.RSTD_NADRP = prcrng(Paul_Results.POC(iNADRP));%Value = 28.3069 (24/03/22)
Paul_Results.Median_NADRM = median(Paul_Results.POC(iNADRM));%Value = 7.2218 (24/03/22)
Paul_Results.RSTD_NADRM = prcrng(Paul_Results.POC(iNADRM));%Value = 3.5092 (24/03/22)

Paul_Results.Median_NASTP = median(Paul_Results.POC(iNASTP));%Value = 26.7537 (24/03/22)
Paul_Results.RSTD_NASTP = prcrng(Paul_Results.POC(iNASTP));%Value = 11.0975 (24/03/22)
Paul_Results.Median_NASTM = median(Paul_Results.POC(iNASTM));%Value = 5.6144 (24/03/22)
Paul_Results.RSTD_NASTM = prcrng(Paul_Results.POC(iNASTM));%Value = 1.1197 (24/03/22)

Paul_Results.Median_NATLP = median(Paul_Results.POC(iNATLP));%Value = 14.5174 (24/03/22)
Paul_Results.RSTD_NATLP = prcrng(Paul_Results.POC(iNATLP));%Value = 9.7728 (24/03/22)
Paul_Results.Median_NATLM = median(Paul_Results.POC(iNATLM),'omitnan');%Value = 7.8594 (24/03/22)
Paul_Results.RSTD_NATLM = prcrng(Paul_Results.POC(iNATLM));%Value = 2.0787 (24/03/22)

Paul_Results.Median_WTRAP = median(Paul_Results.POC(iWTRAP));%Value = 19.4882 (24/03/22)
Paul_Results.RSTD_WTRAP = prcrng(Paul_Results.POC(iWTRAP));%Value = 10.4352 (24/03/22)
Paul_Results.Median_WTRAM = median(Paul_Results.POC(iWTRAM));%Value = 8.7503 (24/03/22)
Paul_Results.RSTD_WTRAM = prcrng(Paul_Results.POC(iWTRAM));%Value = 4.3424 (24/03/22)

Paul_Results.Median_SATLP = median(Paul_Results.POC(iSATLP));%Value = 13.9178 (24/03/22)
Paul_Results.RSTD_SATLP = prcrng(Paul_Results.POC(iSATLP));%Value = 6.3870 (24/03/22)
Paul_Results.Median_SATLM = median(Paul_Results.POC(iSATLM));%Value = 6.7016 (24/03/22)
Paul_Results.RSTD_SATLM = prcrng(Paul_Results.POC(iSATLM));%Value = 1.7981 (24/03/22)

Paul_Results.Median_SSTCP = median(Paul_Results.POC(iSSTCP));%Value = 46.0122 (24/03/22)
Paul_Results.RSTD_SSTCP = prcrng(Paul_Results.POC(iSSTCP));%Value = 22.7466 (24/03/22)
Paul_Results.Median_SSTCM = median(Paul_Results.POC(iSSTCM));%Value = 8.6542 (24/03/22)
Paul_Results.RSTD_SSTCM = prcrng(Paul_Results.POC(iSSTCM));%Value = 3.2662 (24/03/22)

Paul_Results.Median_POC_good_productive = median(Paul_Results.POC(iP_good),'omitnan');%Value = 18.5402 (24/03/22)
Paul_Results.RSTD_POC__good_productive = prcrng(Paul_Results.POC(iP_good));%Value = 19.6006 (24/03/22)
Paul_Results.Median_POC_good_mesopelagic = median(Paul_Results.POC(iM_good),'omitnan'); %Value = 7.3550 (24/03/22)
Paul_Results.RSTD_POC_good_mesopelagic = prcrng(Paul_Results.POC(iM_good));%Value = 2.5511 (24/03/22)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%% plotting cruise track
% f8 = figure(8);
% clc
% h = worldmap([-55 55],[-60 10]);
% set(findall(h,'Tag','PLabel'),'visible','on','FontSize',14);
% set(findall(h,'Tag','MLabel'),'visible','off');
% latlim = [-55 55];
% lonlim = [-60 10];
% % [Z, refvec] = etopo('etopo1_ice_c_f4.flt', 1, latlim, lonlim); 
% [Z, refvec] = etopo('etopo1_bed_c_f4.flt', 1, latlim, lonlim); 
% load coastlines
% geoshow(Z, refvec, 'DisplayType', 'surface');
% axis off
% geoshow('landareas.shp', 'FaceColor', [0.15 0.5 0.15], 'EdgeColor', 'black');
% hold on
% geoshow (Paul_Results.Lat(Paul_Results.Province == 1), Paul_Results.Lon(Paul_Results.Province == 1), 'DisplayType','point','Marker', 'o', 'MarkerEdgeColor','k','MarkerFaceColor','y','MarkerSize', 5);
% geoshow (Paul_Results.Lat(Paul_Results.Province == 2), Paul_Results.Lon(Paul_Results.Province == 2), 'DisplayType','point','Marker', 'o', 'MarkerEdgeColor','k','MarkerFaceColor','m','MarkerSize', 5);
% geoshow (Paul_Results.Lat(Paul_Results.Province == 3), Paul_Results.Lon(Paul_Results.Province == 3), 'DisplayType','point','Marker', 'o', 'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize', 5);
% geoshow (Paul_Results.Lat(Paul_Results.Province == 4), Paul_Results.Lon(Paul_Results.Province == 4), 'DisplayType','point','Marker', 'o', 'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize', 5);
% geoshow (Paul_Results.Lat(Paul_Results.Province == 5), Paul_Results.Lon(Paul_Results.Province == 5), 'DisplayType','point','Marker', 'o', 'MarkerEdgeColor','k','MarkerFaceColor','c','MarkerSize', 5);
% geoshow (Paul_Results.Lat(Paul_Results.Province == 6), Paul_Results.Lon(Paul_Results.Province == 6), 'DisplayType','point','Marker', 'o', 'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize', 5);
% textm(48.5, -25, 'NADR','Color','y','FontSize', 14);
% textm(38, -40, 'NAST','Color','m','FontSize', 14);
% textm(18, -42, 'NATL','Color','b','FontSize', 14);
% textm(0, -21, 'WTRA','Color','w','FontSize', 14);
% textm(-20, -21, 'SATL','Color','c','FontSize', 14);
% textm(-46, -35, 'SSTC','Color','r','FontSize', 14);
% tightmap
% %% Saving figure
% f8name = 'AMT24_Track';
% print (f8, '-dpng', strcat(figDir,f8name), '-r350');
% print (f8, '-depsc', strcat(figDir,f8name));
%% Correction for biases / Results 3.3
% Total_aDOC.MedianTrueCarbonMassP = median(Total_aDOC.TrueCarbonMass(iP_good));%Value = 9.4473 (24/03/22)
% Total_aDOC.RSDTrueCarbonMassP = mad(Total_aDOC.TrueCarbonMass(iP_good),1);%Value = 2.6944 (24/03/22)
% Total_aDOC.MedianTrueCarbonMassM = median(Total_aDOC.TrueCarbonMass(iM_good));%Value = 6.3858 (24/03/22)
% Total_aDOC.RSDTrueCarbonMassM = mad(Total_aDOC.TrueCarbonMass(iM_good),1);%Value = 1.4057 (24/03/22)
% 
% Total_aDOC.MedianCarbonMassP = median(Total_aDOC.CarbonMass(iP_good));%Value = 12.0529 (24/03/22)
% Total_aDOC.RSDCarbonMassP = mad(Total_aDOC.CarbonMass(iP_good),1);%Value = 2.7477 (24/03/22)
% Total_aDOC.MedianCarbonMassM = median(Total_aDOC.CarbonMass(iM_good));%Value = 9.0872 (24/03/22)
% Total_aDOC.RSDCarbonMassM = mad(Total_aDOC.CarbonMass(iM_good),1);%Value = 1.7183 (24/03/22)

Total_aDOC.MedianTrueCarbonMassP = median(Total_aDOC.TrueCarbonMass(iP_good));%Value = 9.4473 (24/03/22)
Total_aDOC.RSDTrueCarbonMassP = prcrng(Total_aDOC.TrueCarbonMass(iP_good));%Value = 4.3636 (24/03/22)
Total_aDOC.MedianTrueCarbonMassM = median(Total_aDOC.TrueCarbonMass(iM_good));%Value = 6.3858 (24/03/22)
Total_aDOC.RSDTrueCarbonMassM = prcrng(Total_aDOC.TrueCarbonMass(iM_good));%Value = 2.4857 (24/03/22)

Total_aDOC.MedianCarbonMassP = median(Total_aDOC.CarbonMass(iP_good));%Value = 12.0529 (24/03/22)
Total_aDOC.RSDCarbonMassP = prcrng(Total_aDOC.CarbonMass(iP_good));%Value = 4.3762 (24/03/22)
Total_aDOC.MedianCarbonMassM = median(Total_aDOC.CarbonMass(iM_good));%Value = 9.0872 (24/03/22)
Total_aDOC.RSDCarbonMassM = prcrng(Total_aDOC.CarbonMass(iM_good));%Value = 2.9559 (24/03/22)

[Tmp1,Tmp2] = corrcoef(Total_uPOC.TrueCarbonMass(Paul_Results.Province == 1 & igood == 1),Total_aDOC.TrueCarbonMass(Paul_Results.Province == 1 & igood == 1));
Tmp1= Tmp1(2,1) %Value = 0.8272 (24/03/22)
Tmp2= Tmp2(2,1) %Value = 3.1654e-08 (24/03/22) Reported as P<.001

clearvars 'Tmp1' 'Tmp2'

[Tmp1,Tmp2] = corrcoef(Total_uPOC.POC(igood == 1),Total_aDOC.POC(igood == 1));
Tmp1= Tmp1(2,1) %Value = 0.7149 (24/03/22)
Tmp2= Tmp2(2,1) %Value = 5.8360e-52 (24/03/22) Reported as P<.001

clearvars 'Tmp1' 'Tmp2'

[Tmp1,Tmp2] = corrcoef(Total_uPOC.POC(iP_good),Total_aDOC.POC(iP_good),'Rows','pairwise');
Tmp1= Tmp1(2,1) %Value = 0.7149 (24/03/22)
Tmp2= Tmp2(2,1) %Value = 5.8360e-52 (24/03/22) Reported as P<.001

clearvars 'Tmp1' 'Tmp2'

[Tmp1,Tmp2] = corrcoef(Total_uPOC.POC(iM_good),Total_aDOC.POC(iM_good),'Rows','pairwise');
Tmp1= Tmp1(2,1) %Value = 0.7149 (24/03/22)
Tmp2= Tmp2(2,1) %Value = 5.8360e-52 (24/03/22) Reported as P<.001

clearvars 'Tmp1' 'Tmp2'

[Tmp1,Tmp2] = corrcoef(Total_uPOC.POC(Paul_Results.Province == 1 & igood == 1),Total_aDOC.POC(Paul_Results.Province == 1 & igood == 1));
Tmp1= Tmp1(2,1) %Value = 0.8582 (24/03/22)
Tmp2= Tmp2(2,1) %Value = 2.6569e-09 (24/03/22) Reported as P<.001

clearvars 'Tmp1' 'Tmp2'

[Tmp1,Tmp2] = corrcoef(Total_uPOC.TrueCarbonMass(Paul_Results.Province == 6 & igood == 1),Total_aDOC.TrueCarbonMass(Paul_Results.Province == 6 & igood == 1));
Tmp1= Tmp1(2,1) %Value = 0.6807 (24/03/22)
Tmp2= Tmp2(2,1) %Value = 5.1303e-12 (24/03/22) Reported as P<.001
clearvars 'Tmp1' 'Tmp2'

[Tmp1,Tmp2] = corrcoef(Total_uPOC.POC(Paul_Results.Province == 6 & igood == 1),Total_aDOC.POC(Paul_Results.Province == 6 & igood == 1));
Tmp1= Tmp1(2,1) %Value = 0.8866 (24/03/22)
Tmp2= Tmp2(2,1) %Value = 1.6598e-27 (24/03/22) Reported as P<.001
clearvars 'Tmp1' 'Tmp2'

[Tmp1,Tmp2] = corrcoef(Total_aDOC.Volumen,Total_aDOC.TrueCarbonMass);
Tmp1= Tmp1(2,1) %Value = -0.2616 (24/03/22)
Tmp2= Tmp2(2,1) %Value = 1.4843e-07 (24/03/22) Reported as P<.001
clearvars 'Tmp1' 'Tmp2'

% Total_aDOC.MaDOC_MuPOC_percent = (Total_aDOC.TrueCarbonMass./ Total_uPOC.TrueCarbonMass)*100;
% Total_aDOC.MaDOC_MuPOC_median_productive_percent = median(Total_aDOC.MaDOC_MuPOC_percent(iP_good)); %Value = 9.1180 (24/03/22)
% Total_aDOC.MaDOC_MuPOC_std_productive_percent = mad(Total_aDOC.MaDOC_MuPOC_percent(iP_good),1); %Value = 1.9924 (24/03/22)
% Total_aDOC.MaDOC_MuPOC_median_mesopelagic_percent = median(Total_aDOC.MaDOC_MuPOC_percent(iM_good)); %Value = 11.6025 (24/03/22)
% Total_aDOC.MaDOC_MuPOC_std_mesopelagic_percent = mad(Total_aDOC.MaDOC_MuPOC_percent(iM_good),1); %Value = 2.9885 (24/03/22)

Total_aDOC.MaDOC_MuPOC_percent = (Total_aDOC.TrueCarbonMass./ Total_uPOC.TrueCarbonMass)*100;
Total_aDOC.MaDOC_MuPOC_median_productive_percent = median(Total_aDOC.MaDOC_MuPOC_percent(iP_good)); %Value = 9.1180 (24/03/22)
Total_aDOC.MaDOC_MuPOC_std_productive_percent = prcrng(Total_aDOC.MaDOC_MuPOC_percent(iP_good)); %Value = 3.1078 (24/03/22)
Total_aDOC.MaDOC_MuPOC_median_mesopelagic_percent = median(Total_aDOC.MaDOC_MuPOC_percent(iM_good)); %Value = 11.6025 (24/03/22)
Total_aDOC.MaDOC_MuPOC_std_mesopelagic_percent = prcrng(Total_aDOC.MaDOC_MuPOC_percent(iM_good)); %Value = 4.9897 (24/03/22)


%% Calculation of POC concentration for section 3.7 New manuscript (25/03/22)
% Total_aDOC.aDOC_concentration_corrected = POC_concentration(Total_aDOC.TrueCarbonMass,Total_aDOC.Volumen); %POC concentration of aDOC filters using the true mass
% Total_aDOC.aDOC_concentration_corrected_median_good_productive = median(Total_aDOC.POC(iP_good)); %Value = 2.1785 (25/03/22)
% Total_aDOC.aDOC_concentration_corrected_std_good_productive = mad(Total_aDOC.POC(iP_good),1); %Value = 0.9559 (25/03/22)
% Total_aDOC.aDOC_concentration_corrected_median_good_mesopelagic = median(Total_aDOC.POC(iM_good)); %Value = 0.9922 (25/03/22)
% Total_aDOC.aDOC_concentration_corrected_std_good_mesopelagic = mad(Total_aDOC.POC(iM_good),1); %Value = 0.2311 (25/03/22)
% 
% Total_aDOC.aDOC_concentration = POC_concentration(Total_aDOC.CarbonMass,Total_aDOC.Volumen); %POC concentration of aDOC filters using the carbon mass
% Total_aDOC.aDOC_concentration_median_good = median(Total_aDOC.aDOC_concentration(igood)); %Value = 1.9012 (25/03/22)
% Total_aDOC.aDOC_concentration_std_good = mad(Total_aDOC.aDOC_concentration (igood)); %Value = 1.2535 (25/03/22)
% Total_aDOC.aDOC_concentration_median_good_productive = median(Total_aDOC.aDOC_concentration(iP_good)); %Value = 2.6765 (25/03/22)
% Total_aDOC.aDOC_concentration_std_good_productive = mad(Total_aDOC.aDOC_concentration(iP_good),1); %Value = 1.1028 (25/03/22)
% Total_aDOC.aDOC_concentration_median_good_mesopelagic = median(Total_aDOC.aDOC_concentration(iM_good)); %Value = 1.3963 (25/03/22)
% Total_aDOC.aDOC_concentration_std_good_mesopelagic = mad(Total_aDOC.aDOC_concentration(iM_good),1); %Value = 0.2913 (25/03/22)
% 
% Total_uPOC.uPOC_concentration_corrected = POC_concentration(Total_uPOC.TrueCarbonMass,Total_uPOC.Volumen); %POC concentration of uPOC filters using the true mass
% Total_uPOC.uPOC_concentration_corrected_median_good_productive = median(Total_uPOC.uPOC_concentration_corrected(iP_good)); %Value = 20.7175 (25/03/22)
% Total_uPOC.uPOC_concentration_corrected_std_good_productive = mad(Total_uPOC.uPOC_concentration_corrected(iP_good),1); %Value = 9.4717 (25/03/22)
% Total_uPOC.uPOC_concentration_corrected_median_good_mesopelagic = median(Total_uPOC.uPOC_concentration_corrected(iM_good)); %Value = 8.4073 (25/03/22)
% Total_uPOC.uPOC_concentration_corrected_std_good_mesopelagic = mad(Total_uPOC.uPOC_concentration_corrected(iM_good),1); %Value = 1.7381 (25/03/22)
% 
% Total_uPOC.uPOC_concentration_corrected_POC = ((Total_uPOC.uPOC_concentration_corrected - Total_uPOC.POC)./Total_uPOC.POC)*100; % Difference uPOC conc and POC %
% Total_uPOC.uPOC_concentration_corrected_POC_median_good_productive = median(Total_uPOC.uPOC_concentration_corrected_POC(iP_good)); %Value = 10.0329(25/03/22)
% Total_uPOC.uPOC_concentration_corrected_POC_std_good_productive  = mad(Total_uPOC.uPOC_concentration_corrected_POC(iP_good),1); %Value = 2.4111 (25/03/22)
% Total_uPOC.uPOC_concentration_corrected_POC_median_good_mesopelagic = median(Total_uPOC.uPOC_concentration_corrected_POC(iM_good)); %Value = 12.9557 (25/03/22)
% Total_uPOC.uPOC_concentration_corrected_POC_std_good_mesopelagic = mad(Total_uPOC.uPOC_concentration_corrected_POC(iM_good),1); %Value = 3.8492 (25/03/22)
% 
% Total_uPOC.uPOC_concentration = POC_concentration(Total_uPOC.CarbonMass,Total_uPOC.Volumen); %POC concentration of uPOC filters using the carbon mass
% Total_uPOC.uPOC_concentration_POC = ((Total_uPOC.uPOC_concentration - Total_uPOC.POC)./Total_uPOC.POC)*100; % Difference uPOC conc and POC %
% Total_uPOC.uPOC_concentration_POC_median_good_productive = median(Total_uPOC.uPOC_concentration_POC(iP_good)); %Value = 12.5629(25/03/22)
% Total_uPOC.uPOC_concentration_POC_std_good_productive  = mad(Total_uPOC.uPOC_concentration_POC(iP_good),1); %Value = 2.9557(25/03/22)
% Total_uPOC.uPOC_concentration_POC_median_good_mesopelagic = median(Total_uPOC.uPOC_concentration_POC(iM_good)); %Value = 18.1276(25/03/22)
% Total_uPOC.uPOC_concentration_POC_std_good_mesopelagic = mad(Total_uPOC.uPOC_concentration_POC(iM_good),1); %Value = 4.9878 (25/03/22)
% 
% Total_uPOC.MassPOC_singleblank = (Total_uPOC.TrueCarbonMass - Total_aDOC.MedianTrueCarbonMassM); % Mass of POC corrected using a single blank median of values >200 m
% Total_uPOC.POC_singleblank = POC_concentration(Total_uPOC.MassPOC_singleblank,Total_uPOC.Volumen); %POC concentration using a single blank
% Total_uPOC.POC_singleblank_POC = ((Total_uPOC.POC_singleblank - Total_uPOC.POC)./Total_uPOC.POC)*100; % Difference POC_singleBlank conc and POC in %
% Total_uPOC.POC_singleblank_POC_median_good_productive = median(Total_uPOC.POC_singleblank_POC(iP_good)); %Value = 3.3125 (25/03/22) in %
% Total_uPOC.POC_singleblank_POC_std_good_productive  = mad(Total_uPOC.POC_singleblank_POC(iP_good),1); %Value = 2.0504 (25/03/22) in %
% Total_uPOC.POC_singleblank_POC_median_good_mesopelagic = median(Total_uPOC.POC_singleblank_POC(iM_good)); %Value = -0.0515 (25/03/22) in %
% Total_uPOC.POC_singleblank_POC_std_good_mesopelagic = mad(Total_uPOC.POC_singleblank_POC(iM_good),1); %Value = 3.0547 (25/03/22) in %

Total_aDOC.aDOC_concentration_corrected_median_good_productive = median(Total_aDOC.POC(iP_good)); %Value = 2.1785 (25/03/22)
Total_aDOC.aDOC_concentration_corrected_std_good_productive = prcrng(Total_aDOC.POC(iP_good)); %Value = 1.3616 (25/03/22)
Total_aDOC.aDOC_concentration_corrected_median_good_mesopelagic = median(Total_aDOC.POC(iM_good)); %Value = 0.9922 (25/03/22)
Total_aDOC.aDOC_concentration_corrected_std_good_mesopelagic = prcrng(Total_aDOC.POC(iM_good)); %Value = 0.4119 (25/03/22)

Total_aDOC.aDOC_concentration = POC_concentration(Total_aDOC.CarbonMass,Total_aDOC.Volumen); %POC concentration of aDOC filters using the carbon mass
Total_aDOC.aDOC_concentration_median_good = median(Total_aDOC.aDOC_concentration(igood)); %Value = 1.9012 (25/03/22)
Total_aDOC.aDOC_concentration_std_good = prcrng(Total_aDOC.aDOC_concentration (igood)); %Value = 1.3076 (25/03/22)
Total_aDOC.aDOC_concentration_median_good_productive = median(Total_aDOC.aDOC_concentration(iP_good)); %Value = 2.6765 (25/03/22)
Total_aDOC.aDOC_concentration_std_good_productive = prcrng(Total_aDOC.aDOC_concentration(iP_good)); %Value = 1.6072 (25/03/22)
Total_aDOC.aDOC_concentration_median_good_mesopelagic = median(Total_aDOC.aDOC_concentration(iM_good)); %Value = 1.3963 (25/03/22)
Total_aDOC.aDOC_concentration_std_good_mesopelagic = prcrng(Total_aDOC.aDOC_concentration(iM_good)); %Value = 0.5206 (25/03/22)

Total_uPOC.uPOC_concentration_corrected = POC_concentration(Total_uPOC.TrueCarbonMass,Total_uPOC.Volumen); %POC concentration of uPOC filters using the true mass
Total_uPOC.uPOC_concentration_corrected_median_good_productive = median(Total_uPOC.uPOC_concentration_corrected(iP_good)); %Value = 20.7175 (25/03/22)
Total_uPOC.uPOC_concentration_corrected_std_good_productive = prcrng(Total_uPOC.uPOC_concentration_corrected(iP_good)); %Value = 20.4734 (25/03/22)
Total_uPOC.uPOC_concentration_corrected_median_good_mesopelagic = median(Total_uPOC.uPOC_concentration_corrected(iM_good)); %Value = 8.4073 (25/03/22)
Total_uPOC.uPOC_concentration_corrected_std_good_mesopelagic = prcrng(Total_uPOC.uPOC_concentration_corrected(iM_good)); %Value = 2.7810 (25/03/22)

Total_uPOC.uPOC_concentration_corrected_POC = ((Total_uPOC.uPOC_concentration_corrected - Total_uPOC.POC)./Total_uPOC.POC)*100; % Difference uPOC conc and POC %
Total_uPOC.uPOC_concentration_corrected_POC_median_good_productive = median(Total_uPOC.uPOC_concentration_corrected_POC(iP_good)); %Value = 10.0329(25/03/22)
Total_uPOC.uPOC_concentration_corrected_POC_std_good_productive  = prcrng(Total_uPOC.uPOC_concentration_corrected_POC(iP_good)); %Value = 3.8108 (25/03/22)
Total_uPOC.uPOC_concentration_corrected_POC_median_good_mesopelagic = median(Total_uPOC.uPOC_concentration_corrected_POC(iM_good)); %Value = 12.9557 (25/03/22)
Total_uPOC.uPOC_concentration_corrected_POC_std_good_mesopelagic = prcrng(Total_uPOC.uPOC_concentration_corrected_POC(iM_good)); %Value = 6.5670 (25/03/22)

Total_uPOC.uPOC_concentration = POC_concentration(Total_uPOC.CarbonMass,Total_uPOC.Volumen); %POC concentration of uPOC filters using the carbon mass
Total_uPOC.uPOC_concentration_POC = ((Total_uPOC.uPOC_concentration - Total_uPOC.POC)./Total_uPOC.POC)*100; % Difference uPOC conc and POC %
Total_uPOC.uPOC_concentration_POC_median_good_productive = median(Total_uPOC.uPOC_concentration_POC(iP_good)); %Value = 12.5629(25/03/22)
Total_uPOC.uPOC_concentration_POC_std_good_productive  = prcrng(Total_uPOC.uPOC_concentration_POC(iP_good)); %Value = 5.3215 (25/03/22)
Total_uPOC.uPOC_concentration_POC_median_good_mesopelagic = median(Total_uPOC.uPOC_concentration_POC(iM_good)); %Value = 18.1276(25/03/22)
Total_uPOC.uPOC_concentration_POC_std_good_mesopelagic = prcrng(Total_uPOC.uPOC_concentration_POC(iM_good)); %Value = 8.5492 (25/03/22)

Total_uPOC.MassPOC_singleblank = (Total_uPOC.TrueCarbonMass - Total_aDOC.MedianTrueCarbonMassM); % Mass of POC corrected using a single blank median of values >200 m
Total_uPOC.POC_singleblank = POC_concentration(Total_uPOC.MassPOC_singleblank,Total_uPOC.Volumen); %POC concentration using a single blank
Total_uPOC.POC_singleblank_POC = ((Total_uPOC.POC_singleblank - Total_uPOC.POC)./Total_uPOC.POC)*100; % Difference POC_singleBlank conc and POC in %
Total_uPOC.POC_singleblank_POC_median_good_productive = median(Total_uPOC.POC_singleblank_POC(iP_good)); %Value = 3.3125 (25/03/22) in %
Total_uPOC.POC_singleblank_POC_std_good_productive  = prcrng(Total_uPOC.POC_singleblank_POC(iP_good)); %Value = 3.3325 (25/03/22) in %
Total_uPOC.POC_singleblank_POC_median_good_mesopelagic = median(Total_uPOC.POC_singleblank_POC(iM_good)); %Value = -0.0515 (25/03/22) in %
Total_uPOC.POC_singleblank_POC_std_good_mesopelagic = prcrng(Total_uPOC.POC_singleblank_POC(iM_good)); %Value = 4.9965 (25/03/22) in %

a = cell2mat(Total_uPOC.ID);% converting cell array into an ordinary array
flds = string(unique(a(:,9:end), 'rows')); % getting a string with the number of the stations ".XX"
for ides = 1:length(flds)
    strfld = flds(ides); % getting station by station 
    istation = contains (Total_aDOC.ID, strfld);% logical index of station containing strfld 
    [V, imaxdepth] = max(contains (Total_aDOC.ID, strfld)); %getting the index of the deepest aDOC blank measured at the corresponding station
    aDOC_TrueCarbonMass_maxdepth = Total_aDOC.TrueCarbonMass(imaxdepth); % getting the value for the previous index
    uPOC_POCMass_maxdepth(istation,1) = Total_uPOC.TrueCarbonMass(istation) - aDOC_TrueCarbonMass_maxdepth;
    Total_uPOC.POCMass_maxdepth(istation,1) = Total_uPOC.TrueCarbonMass(istation) - aDOC_TrueCarbonMass_maxdepth;
    
end

% Total_uPOC.POC_maxdepth = POC_concentration(Total_uPOC.POCMass_maxdepth,Total_uPOC.Volumen);%POC concentration using a aDOC at deepest depth
% Total_uPOC.POC_maxdepth_POC = ((Total_uPOC.POC_maxdepth - Total_uPOC.POC)./Total_uPOC.POC)*100; % Difference POC_singleBlank conc and POC %
% Total_uPOC.POC_maxdepth_POC_median_good_productive = median(Total_uPOC.POC_maxdepth_POC(iP_good)); %Value = 2.4826 (25/03/22) in %
% Total_uPOC.POC_maxdepth_POC_std_good_productive  = mad(Total_uPOC.POC_maxdepth_POC(iP_good),1); %Value = 2.3550 (25/03/22) in %
% Total_uPOC.POC_maxdepth_POC_median_good_mesopelagic = median(Total_uPOC.POC_maxdepth_POC(iM_good)); %Value = 0 (25/03/22) in %
% Total_uPOC.POC_maxdepth_POC_std_good_mesopelagic = mad(Total_uPOC.POC_maxdepth_POC(iM_good)); %Value = 3.7014 (25/03/22) in %

Total_uPOC.POC_maxdepth = POC_concentration(Total_uPOC.POCMass_maxdepth,Total_uPOC.Volumen);%POC concentration using a aDOC at deepest depth
Total_uPOC.POC_maxdepth_POC = ((Total_uPOC.POC_maxdepth - Total_uPOC.POC)./Total_uPOC.POC)*100; % Difference POC_singleBlank conc and POC %
Total_uPOC.POC_maxdepth_POC_median_good_productive = median(Total_uPOC.POC_maxdepth_POC(iP_good)); %Value = 2.4826 (25/03/22) in %
Total_uPOC.POC_maxdepth_POC_std_good_productive  = prcrng(Total_uPOC.POC_maxdepth_POC(iP_good)); %Value = 3.5407 (25/03/22) in %
Total_uPOC.POC_maxdepth_POC_median_good_mesopelagic = median(Total_uPOC.POC_maxdepth_POC(iM_good)); %Value = 0 (25/03/22) in %
Total_uPOC.POC_maxdepth_POC_std_good_mesopelagic = prcrng(Total_uPOC.POC_maxdepth_POC(iM_good)); %Value = 3.5141 (25/03/22) in %

% plot(Total_aDOC.TrueCarbonMass(istation),Total_aDOC.Depth(istation), '--b');
% set(gca,'Ydir','reverse');
% xlabel ('Carbon mass aDOC ($\mu g$)', 'FontSize', 16, 'FontAngle', 'italic', 'Interpreter', 'latex');
% ylabel ('Depth (m)', 'FontSize', 16, 'FontAngle', 'italic');
% ylim([0 600]);
% hold on
%     plot(Total_aDOC.TrueCarbonMass(imaxdepth),Total_aDOC.Depth(imaxdepth), 'or');
%     plot(Total_aDOC.TrueCarbonMass(istation),Total_aDOC.Depth(istation), '-g');
%     plot(Total_aDOC.TrueCarbonMass(imaxdepth),Total_aDOC.Depth(imaxdepth), 'xk');
%     plot(Total_aDOC.TrueCarbonMass(istation),Total_aDOC.Depth(istation), '.-y');
%     plot(Total_aDOC.TrueCarbonMass(imaxdepth),Total_aDOC.Depth(imaxdepth), 'dc')
% legend('Station 68', 'Max depth station 68', 'Station 7', 'Max depth station 7','Station 1', 'Max depth station 1' );

%% Estimation of loss of particles from uPOC to aDOC
% Total_aDOC.Diff_TrueCarbonMass_allaDOC_median200deep = (Total_aDOC.TrueCarbonMass - Total_aDOC.MedianTrueCarbonMassM);
% Total_uPOC.MassPOC_lossofparticles = (Total_uPOC.MassPOC + Total_aDOC.Diff_TrueCarbonMass_allaDOC_median200deep);
% Total_uPOC.POC_lossofparticles = POC_concentration(Total_uPOC.MassPOC_lossofparticles,Total_uPOC.Volumen);
% Total_uPOC.POC_lossofparticles_POC = ((Total_uPOC.POC_lossofparticles - Total_uPOC.POC)./Total_uPOC.POC)*100; % Difference POC_singleBlank conc and POC %
% Total_uPOC.POC_lossofparticles_POC_median_good_productive = median(Total_uPOC.POC_lossofparticles_POC(iP_good)); %Value = 3.3135 (25/03/22) in %
% Total_uPOC.POC_lossofparticles_POC_std_good_productive  = mad(Total_uPOC.POC_lossofparticles_POC(iP_good)); %Value = 4.1709 (25/03/22) in %
% Total_uPOC.POC_lossofparticles_POC_median_good_mesopelagic = median(Total_uPOC.POC_lossofparticles_POC(iM_good)); %Value = -0.0515 (25/03/22) in %
% Total_uPOC.POC_lossofparticles_POC_std_good_mesopelagic = mad(Total_uPOC.POC_lossofparticles_POC(iM_good),1); %Value = 3.0547 (25/03/22) in %

Total_aDOC.Diff_TrueCarbonMass_allaDOC_median200deep = (Total_aDOC.TrueCarbonMass - Total_aDOC.MedianTrueCarbonMassM);
Total_uPOC.MassPOC_lossofparticles = (Total_uPOC.MassPOC + Total_aDOC.Diff_TrueCarbonMass_allaDOC_median200deep);
Total_uPOC.POC_lossofparticles = POC_concentration(Total_uPOC.MassPOC_lossofparticles,Total_uPOC.Volumen);
Total_uPOC.POC_lossofparticles_POC = ((Total_uPOC.POC_lossofparticles - Total_uPOC.POC)./Total_uPOC.POC)*100; % Difference POC_singleBlank conc and POC %
Total_uPOC.POC_lossofparticles_POC_median_good_productive = median(Total_uPOC.POC_lossofparticles_POC(iP_good)); %Value = 3.3135 (25/03/22) in %
Total_uPOC.POC_lossofparticles_POC_std_good_productive  = prcrng(Total_uPOC.POC_lossofparticles_POC(iP_good)); %Value = 3.3325 (25/03/22) in %
Total_uPOC.POC_lossofparticles_POC_median_good_mesopelagic = median(Total_uPOC.POC_lossofparticles_POC(iM_good)); %Value = -0.0515 (25/03/22) in %
Total_uPOC.POC_lossofparticles_POC_std_good_mesopelagic = prcrng(Total_uPOC.POC_lossofparticles_POC(iM_good)); %Value = 4.9965 (25/03/22) in %


%% Calculating positive bias by neglecting to use aDOC blanks
load(strcat(TotalDir, fname1)); %Loading data from duplicate analysis. 

% PositiveBias = Total_aDOC.CarbonMass - FilterAcidifiedTotal.MedianCarbonMass;
% PositiveBiasAll = median(PositiveBias, 'omitnan'); %Value = 6.8478 (25/03/22)
% PositiveBiasProductive = median(PositiveBias(iP_good),'omitnan'); %Value = 7.9599 (25/03/22) 
% PositiveBiasProductive_STD = mad(PositiveBias(iP_good),1); %Value = 2.7477 (25/03/22) 
% PositiveBiasMesopelagic = median(PositiveBias(iM_good),'omitnan'); %Value = 4.9942 (25/03/22) 
% PositiveBiasMesopelagic_STD = mad(PositiveBias(iM_good),1); %Value = 1.7183 (25/03/22) 
% 
% Total_uPOC.uPOC_Acidified = Total_uPOC.CarbonMass - FilterAcidifiedTotal.MedianCarbonMass;
% Total_uPOC.uPOC_Acidified_concentration = POC_concentration(Total_uPOC.uPOC_Acidified,Total_uPOC.Volumen); %POC concentration of uPOC filters using the carbon mass minus Acidified
% 
% Total_uPOC.uPOC_Acidified_concentration_POC = ((Total_uPOC.uPOC_Acidified_concentration - Total_uPOC.POC)./Total_uPOC.POC)*100; % Difference uPOC conc and POC %
% Total_uPOC.uPOC_Acidified_concentration_POC_median_Prod = median(Total_uPOC.uPOC_Acidified_concentration_POC(iP_good)); %Value = 7.7869 (25/03/22) in %
% Total_uPOC.uPOC_Acidified_concentration_POC_std_Prod  = mad(Total_uPOC.uPOC_Acidified_concentration_POC(iP_good),1); %Value = 1.9764 (25/03/22) in %
% Total_uPOC.uPOC_Acidified_concentration_POC_median_meso = median(Total_uPOC.uPOC_Acidified_concentration_POC(iM_good)); %Value = 9.8470 (25/03/22) in %
% Total_uPOC.uPOC_Acidified_concentration_POC_std_meso = mad(Total_uPOC.uPOC_Acidified_concentration_POC(iM_good),1); %Value = 4.1563 (25/03/22) in %
% 
% Total_uPOC.uPOC_NonAcidified = Total_uPOC.CarbonMass - FilterNonAcidifiedTotal.MedianCarbonMass;
% Total_uPOC.uPOC_NonAcidified_concentration = POC_concentration(Total_uPOC.uPOC_NonAcidified,Total_uPOC.Volumen); %POC concentration of uPOC filters using the carbon mass minus NonAcidified
% 
% Total_uPOC.uPOC_NonAcidified_concentration_POC = ((Total_uPOC.uPOC_NonAcidified_concentration - Total_uPOC.POC)./Total_uPOC.POC)*100; % Difference uPOC conc and POC %
% Total_uPOC.uPOC_NonAcidified_concentration_POC_median_Prod = median(Total_uPOC.uPOC_NonAcidified_concentration_POC(iP_good)); %Value = 8.8421 (25/03/22) in %
% Total_uPOC.uPOC_NonAcidified_concentration_POC_std_Prod  = mad(Total_uPOC.uPOC_NonAcidified_concentration_POC(iP_good),1); %Value = 2.0721 (25/03/22) in %
% Total_uPOC.uPOC_NonAcidified_concentration_POC_median_meso = median(Total_uPOC.uPOC_NonAcidified_concentration_POC(iM_good)); %Value = 11.7023 (25/03/22) in %
% Total_uPOC.uPOC_NonAcidified_concentration_POC_std_meso = mad(Total_uPOC.uPOC_NonAcidified_concentration_POC(iM_good),1); %Value = 4.3127 (25/03/22) in %

PositiveBias = Total_aDOC.CarbonMass - FilterAcidifiedTotal.MedianCarbonMass;
PositiveBiasAll = median(PositiveBias, 'omitnan'); %Value = 6.8478 (25/03/22)
PositiveBiasProductive = median(PositiveBias(iP_good),'omitnan'); %Value = 7.9599 (25/03/22) 
PositiveBiasProductive_STD = prcrng(PositiveBias(iP_good)); %Value = 4.3762 (25/03/22) 
PositiveBiasMesopelagic = median(PositiveBias(iM_good),'omitnan'); %Value = 4.9942 (25/03/22) 
PositiveBiasMesopelagic_STD = prcrng(PositiveBias(iM_good)); %Value = 2.9559 (25/03/22) 

Total_uPOC.uPOC_Acidified = Total_uPOC.CarbonMass - FilterAcidifiedTotal.MedianCarbonMass;
Total_uPOC.uPOC_Acidified_concentration = POC_concentration(Total_uPOC.uPOC_Acidified,Total_uPOC.Volumen); %POC concentration of uPOC filters using the carbon mass minus Acidified

Total_uPOC.uPOC_Acidified_concentration_POC = ((Total_uPOC.uPOC_Acidified_concentration - Total_uPOC.POC)./Total_uPOC.POC)*100; % Difference uPOC conc and POC %
Total_uPOC.uPOC_Acidified_concentration_POC_median_Prod = median(Total_uPOC.uPOC_Acidified_concentration_POC(iP_good)); %Value = 7.7869 (25/03/22) in %
Total_uPOC.uPOC_Acidified_concentration_POC_std_Prod  = prcrng(Total_uPOC.uPOC_Acidified_concentration_POC(iP_good)); %Value = 3.5771 (25/03/22) in %
Total_uPOC.uPOC_Acidified_concentration_POC_median_meso = median(Total_uPOC.uPOC_Acidified_concentration_POC(iM_good)); %Value = 9.8470 (25/03/22) in %
Total_uPOC.uPOC_Acidified_concentration_POC_std_meso = prcrng(Total_uPOC.uPOC_Acidified_concentration_POC(iM_good)); %Value = 7.4427 (25/03/22) in %

Total_uPOC.uPOC_NonAcidified = Total_uPOC.CarbonMass - FilterNonAcidifiedTotal.MedianCarbonMass;
Total_uPOC.uPOC_NonAcidified_concentration = POC_concentration(Total_uPOC.uPOC_NonAcidified,Total_uPOC.Volumen); %POC concentration of uPOC filters using the carbon mass minus NonAcidified

Total_uPOC.uPOC_NonAcidified_concentration_POC = ((Total_uPOC.uPOC_NonAcidified_concentration - Total_uPOC.POC)./Total_uPOC.POC)*100; % Difference uPOC conc and POC %
Total_uPOC.uPOC_NonAcidified_concentration_POC_median_Prod = median(Total_uPOC.uPOC_NonAcidified_concentration_POC(iP_good)); %Value = 8.8421 (25/03/22) in %
Total_uPOC.uPOC_NonAcidified_concentration_POC_std_Prod  = prcrng(Total_uPOC.uPOC_NonAcidified_concentration_POC(iP_good)); %Value = 3.6963 (25/03/22) in %
Total_uPOC.uPOC_NonAcidified_concentration_POC_median_meso = median(Total_uPOC.uPOC_NonAcidified_concentration_POC(iM_good)); %Value = 11.7023 (25/03/22) in %
Total_uPOC.uPOC_NonAcidified_concentration_POC_std_meso = prcrng(Total_uPOC.uPOC_NonAcidified_concentration_POC(iM_good)); %Value = 7.7242 (25/03/22) in %

%% Plot aDOC vs POC, with points coloured by depth and with one symbol for surface and one symbol for mesopelagic
% compute the correlation between aDOC and POC for the surface only and for the mesopelagic only

tmp = corrcoef(Total_POC.POC(iP_good), Total_aDOC.POC(iP_good));
tmp = tmp(2,1) %Value = 0.6627 (25/03/22)
clearvars tmp


mdl_aDOC_vs_POC_productive = fitlm (Total_POC.POC(iP_good), Total_aDOC.POC(iP_good), 'linear', 'RobustOpts', 'on');
IpValue = mdl_aDOC_vs_POC_productive.Coefficients{1,4};
if IpValue > 0.05
    mdl_aDOC_vs_POC_productive = fitlm (Total_POC.POC(iP_good), Total_aDOC.POC(iP_good), 'linear','Intercept',false, 'RobustOpts', 'on');
end

tmp = corrcoef(Total_POC.POC(iM_good), Total_aDOC.POC(iM_good));
tmp = tmp(2,1) %Value = 0.2535 (25/03/22)
clearvars tmp

mdl_aDOC_vs_POC_mesopelagic = fitlm (Total_POC.POC(iM_good), Total_aDOC.POC(iM_good), 'linear', 'RobustOpts', 'on');
IpValue = mdl_aDOC_vs_POC_mesopelagic.Coefficients{1,4};
if IpValue > 0.05
    mdl_aDOC_vs_POC_mesopelagic = fitlm (Total_POC.POC(iM_good), Total_aDOC.POC(iM_good), 'linear','Intercept',false, 'RobustOpts', 'on');
end


POC_mock = linspace(0.8,100,1000);
POC_mock = POC_mock';

aDOC_productive_mock = mdl_aDOC_vs_POC_productive.Coefficients{1,1} + (POC_mock * mdl_aDOC_vs_POC_productive.Coefficients{2,1});   
aDOC_mesopelagic_mock = mdl_aDOC_vs_POC_mesopelagic.Coefficients{1,1} + (POC_mock * mdl_aDOC_vs_POC_mesopelagic.Coefficients{2,1});  


Label = ' Depth (m)';
Range = [0 500];

f10 = figure(10);
clf
h1= scatter(Total_POC.POC(iP_good), Total_aDOC.POC(iP_good), 45, Total_POC.Depth(iP_good), 'o', 'filled', 'MarkerEdgeColor','k', 'LineWidth', 0.25, 'MarkerEdgeAlpha', 0.5);
set(gca, 'TickLength',[0 0], 'FontSize', 16, 'box' ,'on', 'LineWidth', 1, 'Xcolor' ,'k', 'Ycolor', 'k');
xlabel('$POC$ $(mg/m^{3})$', 'FontSize', 18, 'FontAngle', 'italic', 'Interpreter', 'latex');
ylabel('$aDOC$ $(mg/m^{3})$', 'FontSize', 15, 'FontAngle', 'italic', 'Interpreter', 'latex');
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlim([2 100])
ylim([0.2 15])
yticks([0.2, 0.5, 1, 2, 5, 10, 15])
xticks([2, 5, 10, 25, 50, 100])
c = colorbar;
c.Label.String = Label;
caxis(Range)
c.Label.FontSize = 16;
c.Label.FontAngle = 'italic';
colormap(flipud(othercolor('Spectral11')));  
x0=1;
y0=1;
width=20;
height=20;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
    hold on
    
%   h2 = plot(POC_mock, aDOC_productive_mock, 'r--', 'LineWidth', 1,'DisplayName','PZ');
    h3 = scatter(Total_POC.POC(iM_good), Total_aDOC.POC(iM_good), 45, Total_POC.Depth(iM_good), 's', 'filled', 'MarkerEdgeColor','none', 'LineWidth', 0.25, 'MarkerEdgeAlpha', 0.5);
%     xlim([2 100])
%     ylim([0 15])
%   h4 = plot(POC_mock, aDOC_mesopelagic_mock, 'b-.', 'LineWidth', 1, 'DisplayName','MZ');
%    legend([h2,h4],'box', 'off', 'location', 'southeast');

f11 = figure(11);
clf
scatter(Total_POC.POC(iP_good), Total_aDOC.POC(iP_good), 45, Total_POC.Depth(iP_good), 'o', 'filled', 'MarkerEdgeColor','k', 'LineWidth', 0.25);
set(gca, 'TickLength',[0 0], 'FontSize', 16, 'box' ,'on', 'LineWidth', 1, 'Xcolor' ,'k', 'Ycolor', 'k');
%set(gca, 'XDir','reverse','YDir','reverse');
%plot (Total_POC.POC(iP), Total_aDOC.POC(iP),'s','MarkerSize',6); %'MarkerFaceColor','k'
xlabel('$POC$ $(mg/m^{3})$', 'FontSize', 20, 'FontAngle', 'italic', 'Interpreter', 'latex');
%set(gca, 'FontSize', 16, 'box' ,'on', 'LineWidth', 1, 'Xcolor' ,'k', 'Ycolor', 'k');
ylabel('$aDOC$ $(mg/m^{3})$', 'FontSize', 20, 'FontAngle', 'italic', 'Interpreter', 'latex');
%text(5,15,'(a)','FontSize', 20, 'Color', 'k');
xlim([0 100])
ylim([0 15])
%yline(0,'-.k');
c = colorbar;
c.Label.String = Label;
caxis(Range)
c.Label.FontSize = 16;
c.Label.FontAngle = 'italic';
colormap(flipud(othercolor('Spectral11')));  
    hold on    
    plot(POC_mock, aDOC_productive_mock, 'r--', 'LineWidth', 1);
    
    axes('Position',[.63 .65 .21 .25])
    box on 
    scatter(Total_POC.POC(iM_good), Total_aDOC.POC(iM_good), 45, Total_POC.Depth(iM_good), 's', 'filled', 'MarkerEdgeColor','k', 'LineWidth', 0.25);
    xlim([0 40])
    ylim([0 8])
    xlabel('$POC$ $(mg/m^{3})$', 'FontSize', 15, 'FontAngle', 'italic', 'Interpreter', 'latex');
    ylabel('$aDOC$ $(mg/m^{3})$', 'FontSize', 15, 'FontAngle', 'italic', 'Interpreter', 'latex');
        hold on
        plot(POC_mock, aDOC_mesopelagic_mock, 'b-.', 'LineWidth', 1);
%     legend ('PZ', 'MZ', 'box', 'off'); 

% x0=10;
% y0=10;
% width=12;
% height=12
% set(gcf,'units','points','position',[x0,y0,width,height])
%% Saving figure

f10name = 'POC_vs_aDOC_correlation_nofittedlines';
f11name = 'POC_vs_aDOC_correlation_inset';

print (f10, '-dpng', strcat(figDir,f10name), '-r350');
print (f10, '-depsc', strcat(figDir,f10name));
print (f11, '-dpng', strcat(figDir,f11name), '-r350');
print (f11, '-depsc', strcat(figDir,f11name));

%%  
%Saving the data
newdataDir = "D:\Measuring POC_Paper\poc\POC\Data\Total\New Data\";
filename = 'Paul_New_Results.mat';
fnmat = strcat (newdataDir, filename);
save(fnmat,'Paul_Results', 'Total_aDOC', 'Total_uPOC', 'Total_POC');
