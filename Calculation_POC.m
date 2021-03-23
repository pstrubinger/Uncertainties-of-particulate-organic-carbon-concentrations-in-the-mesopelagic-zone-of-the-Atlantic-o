clear all
%% To do list %%
%% load the .mat file
%% calculate the mean value in the capsules, and acidified and non-acidified filters
%% calculate the true mass of carbon in Samples, Blanks, R_samples and R_blanks taking into account acidified, non-acidified, capsules
%% Calculate the POC concentration using the volumen of filtered seawater 

%% Load data from mat file.
%% select current directory
% directories
dataDir = "D:\Measuring POC_Paper\poc\POC\Data\";
fn = dir(strcat(dataDir, "*mat"));
%% 
for ifn = 1:length(fn)
    filename = fn(ifn).name;
    load(strcat(dataDir, filename));

%% calculating the mean value in the capsules, and acidified and non-acidified filters
Capsules.M_capsules = mean(Capsules.CarbonMass); % average carbon mass in the tin capsules.
Capsules.meanPI = mean(Capsules.PI); % average PI in the tin capsules.

Filter_non_acidified.M_Filter_non_acidified = mean(Filter_non_acidified.CarbonMass); % Average carbon mass in the non-acidified filters
Filter_non_acidified.mean_PI = mean(Filter_non_acidified.PI); % Average PI in the non-acidified filters
Filter_non_acidified.std_error = std(Filter_non_acidified.CarbonMass /sqrt(length(Filter_non_acidified.CarbonMass)));

Volume_bottle = 2.2; % in litres 
Uncer_volume = 0.01; % uncertainty of volumetric measurements of each bottle in litres
Samples.bottles = floor(Samples.Volumen / Volume_bottle); % number of bottles used
Samples.Comb_Unc_Volume = sqrt(Samples.bottles * Uncer_volume^2); % combined uncertainty of volume
Blanks.bottles = floor(Samples.Volumen / Volume_bottle); % number of bottles used
Blanks.Comb_Unc_Volume = sqrt(Blanks.bottles * Uncer_volume^2); % combined uncertainty of volume
%% calculating the true mass of carbon POC uncertainty(uPOC, aDOC) in Samples, Blanks, R_samples and R_blanks 
a = cell2mat(Samples.ID);
flds = string(unique(a(:,3:5), 'rows')); % Obtaining the name of the desiccator ex: a01, a02.... 
s = 'mean_PI_';
w = 'std_error_';
for ides = 1:length(flds)
    strfld = flds(ides);
    strfld1 = strcat(s, strfld);  
    strfld2 = strcat(w, strfld); 
    ifilter = contains (Filter_acidified.ID, strfld);
    isample = contains(Samples.ID, strfld);
    iblanks = contains(Blanks.ID, strfld);
    Filter_acidified.(strfld) = mean(Filter_acidified.CarbonMass(ifilter));
    Filter_acidified.(strfld1) = mean(Filter_acidified.PI(ifilter));
    Filter_acidified.(strfld2) = std(Filter_acidified.CarbonMass(ifilter))/sqrt(length(Filter_acidified.CarbonMass(ifilter)));
   
    %out = removeCoffset(sample,capsule,non_acidified,acidified);
    Samples.(strfld) = removeCoffset(Samples.CarbonMass(isample), Capsules.M_capsules, Filter_acidified.(strfld), Filter_non_acidified.M_Filter_non_acidified);

    %out = total_uncertainty_u_POC(Volume, PI_uPOC, PI_cap, PI_acidified,PI_nonacidified, M_uPOC, M_cap,M_acidified, M_nonacidified sigma_volume); 
    uPOCuncertainty.(strfld) = total_uncertainty_u_POC(Samples.Volumen(isample), Samples.PI(isample), Capsules.meanPI, Filter_acidified.(strfld1), Filter_non_acidified.mean_PI, Samples.CarbonMass(isample), Capsules.M_capsules, Filter_acidified.(strfld), Filter_non_acidified.M_Filter_non_acidified, Samples.Comb_Unc_Volume(isample));

    %out = uncertainty_handling(std_mean_ac, std_mean_non, volume)
    uncertainty_hand.(strfld) = uncertainty_handling(Filter_acidified.(strfld2), Samples.Volumen(isample));
    
    %out = removeCoffset(sample,capsule,non_acidified,acidified);
    Blanks.(strfld) = removeCoffset(Blanks.CarbonMass(iblanks), Capsules.M_capsules, Filter_acidified.(strfld), Filter_non_acidified.M_Filter_non_acidified);
    
    %out = total_uncertainty_a_DOC(Volume, PI_aDOC, PI_cap, PI_acidified, PI_nonacidified, M_aDOC, M_cap,M_acidified, M_nonacidified , sigma_volume)
    aDOCuncertainty.(strfld) = total_uncertainty_a_DOC(Blanks.Volumen(iblanks), Blanks.PI(iblanks), Capsules.meanPI, Filter_acidified.(strfld1), Filter_non_acidified.mean_PI, Blanks.CarbonMass(iblanks), Capsules.M_capsules, Filter_acidified.(strfld), Filter_non_acidified.M_Filter_non_acidified, Blanks.Comb_Unc_Volume(iblanks));
   
    if ides==1
        Samples.TrueCarbonMass = Samples.(strfld);
        Blanks.TrueCarbonMass = Blanks.(strfld);
        Samples.uPOCuncertainty = uPOCuncertainty.(strfld);
        Blanks.aDOCuncertainty = aDOCuncertainty.(strfld);
        Samples.uncertainty_handling = uncertainty_hand.(strfld);
        Blanks.uncertainty_handling = uncertainty_hand.(strfld);
        
    else
        Samples.TrueCarbonMass = vertcat(Samples.TrueCarbonMass,Samples.(strfld));
        Blanks.TrueCarbonMass = vertcat(Blanks.TrueCarbonMass,Blanks.(strfld));
        Samples.uPOCuncertainty = vertcat(Samples.uPOCuncertainty,uPOCuncertainty.(strfld));
        Blanks.aDOCuncertainty = vertcat(Blanks.aDOCuncertainty,aDOCuncertainty.(strfld));
        Samples.uncertainty_handling = vertcat(Samples.uncertainty_handling, uncertainty_hand.(strfld));
        Blanks.uncertainty_handling = vertcat(Blanks.uncertainty_handling, uncertainty_hand.(strfld));
         
    end
end
 
Samples = rmfield(Samples,flds);
Blanks = rmfield(Blanks,flds);

if ifn == 16
    R_samples.TrueCarbonMass = removeCoffset(R_samples.CarbonMass, Capsules.M_capsules, Filter_acidified.a08, Filter_non_acidified.M_Filter_non_acidified);
%   R_samples.TrueCarbonMass = R_samples.CarbonMass - Capsules.M_capsules - Filter_acidified.a08 - Filter_non_acidified.M_Filter_non_acidified;
    
    R_blanks.TrueCarbonMass = removeCoffset(R_blanks.CarbonMass, Capsules.M_capsules, Filter_acidified.a08, Filter_non_acidified.M_Filter_non_acidified);
%   R_blanks.TrueCarbonMass = R_blanks.CarbonMass - Capsules.M_capsules - Filter_acidified.a08 - Filter_non_acidified.M_Filter_non_acidified;
else
    R_samples.TrueCarbonMass = removeCoffset(R_samples.CarbonMass, Capsules.M_capsules, Filter_acidified.a04, Filter_non_acidified.M_Filter_non_acidified);   
%   R_samples.TrueCarbonMass = R_samples.CarbonMass - Capsules.M_capsules - Filter_acidified.a04 - Filter_non_acidified.M_Filter_non_acidified;
    
    R_blanks.TrueCarbonMass = removeCoffset(R_blanks.CarbonMass, Capsules.M_capsules, Filter_acidified.a04, Filter_non_acidified.M_Filter_non_acidified);
%   R_blanks.TrueCarbonMass = R_blanks.CarbonMass - Capsules.M_capsules - Filter_acidified.a04 - Filter_non_acidified.M_Filter_non_acidified; 
end
    
%% Calculating the mass of particulate organic carbon

if Blanks.TrueCarbonMass > 0 
    Samples.MassPOC = Samples.TrueCarbonMass - Blanks.TrueCarbonMass;
else
    Samples.MassPOC = Samples.TrueCarbonMass;
end

if R_blanks.TrueCarbonMass > 0 
    R_samples.MassPOC = R_samples.TrueCarbonMass - R_blanks.TrueCarbonMass;
else
    R_samples.MassPOC = R_samples.TrueCarbonMass;
end

%% Calculating the POC concentration (mg m^-3) in samples.
for ipoc = 1:length(Samples.ID)
    Samples.POC(ipoc,1) = POC_concentration(Samples.MassPOC(ipoc),Samples.Volumen(ipoc));%Sample.MassPOC is in ug; Sample.Volumen is in litres
end

Blanks.MassPOC = Blanks.TrueCarbonMass;
for ibpoc = 1:length(Blanks.ID)
    Blanks.POC (ibpoc,1) = POC_concentration(Blanks.MassPOC(ibpoc),Blanks.Volumen(ibpoc));%Sample.MassPOC is in ug; Sample.Volumen is in litres 
end

for irpoc = 1:length(R_samples.ID)
    R_samples.POC (irpoc,1) = POC_concentration(R_samples.MassPOC(irpoc),R_samples.Volumen(irpoc));%Sample.MassPOC is in ug; Sample.Volumen is in litres  
end

R_blanks.MassPOC = R_blanks.TrueCarbonMass;
for irbpoc = 1:length(R_blanks.ID)
    R_blanks.POC (irbpoc,1) = POC_concentration(R_samples.MassPOC(irbpoc),R_samples.Volumen(irbpoc)); %Sample.MassPOC is in ug; Sample.Volumen is in litres
end

for irpoc = 1:length(R_samples.ID)
    R_samples.POC (irpoc,1) = POC_concentration(R_samples.MassPOC(irpoc),R_samples.Volumen(irpoc));%Sample.MassPOC is in ug; Sample.Volumen is in litres
end
%% Saving the data
fnmat = strcat (dataDir, filename);
save(fnmat, 'Standards', 'Samples', 'Blanks', 'R_samples', 'R_blanks', 'Capsules', 'Filter_acidified', 'Filter_non_acidified', 'mdl')
end

