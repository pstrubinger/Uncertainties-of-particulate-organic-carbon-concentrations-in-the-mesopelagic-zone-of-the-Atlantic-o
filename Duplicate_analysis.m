clear all;
dataDir = "D:\Measuring POC_Paper\poc\POC\Data\";
figDir = "D:\Measuring POC_Paper\poc\POC\Figures\Duplicates\";
fn = dir(strcat(dataDir, "*mat"));

flds = {'RUN', 'ID', 'Depth', 'Volumen','CarbonMass','TrueCarbonMass','MassPOC','PI', 'POC'};%, 'Duplicate'};
% for iflds = 1:length(flds)
%     Duplicate_Samples.(flds{iflds}) = []; %Creating a empty structure to store data
%     Duplicate_Blanks.(flds{iflds}) = [];
% end
iD = 1; % index for the duplicates of samples
for ifn = 1:length (fn) % calling each file .mat and loading it.
    filename = fn(ifn).name;
    load(strcat(dataDir, filename));
    
    for iR = 1:length(R_samples.ID)% index of the number of replicate in each run       
       
        % new variable containing just the string (last 6 numbers xxx.xx) for comparison the R_sample and samples 
        Sfld = R_samples.ID{iR}(end-5:end);
        Bfld = R_blanks.ID{iR}(end-5:end);

            % for loop to find in the Sample.ID which sample contains replicate
            for iID = 1: length(Samples.ID)
                if ~isempty(strfind(Samples.ID{iID}, Sfld))
                    break
                end %endif
            end %endfor third for
            
            for iIDB = 1: length(Blanks.ID)
                if ~isempty(strfind(Samples.ID{iIDB}, Bfld))
                    break
                end %endif         
            
            end %endfor fourth for
         
        
            for ifld = 1:length(flds)
                Duplicate_Samples.(flds{ifld})(iD,:) = Samples.(flds{ifld})(iID); 
                Duplicate_Blanks.(flds{ifld})(iD,:) = Blanks.(flds{ifld})(iIDB); 
            end %endfor fifth for  
            
           Duplicate_Samples.POC_R(iD,:) = R_samples.POC(iR);
           Duplicate_Samples.CarbonMass_R(iD,:) = R_samples.CarbonMass(iR);
           Duplicate_Samples.TrueCarbonMass_R(iD,:) = R_samples.TrueCarbonMass(iR);
           Duplicate_Samples.MassPOC_R(iD,:) = R_samples.MassPOC(iR);
           
           Duplicate_Blanks.POC_R(iD,:) = R_blanks.POC(iR);
           Duplicate_Blanks.CarbonMass_R(iD,:) = R_blanks.CarbonMass(iR);
           Duplicate_Blanks.TrueCarbonMass_R(iD,:) = R_blanks.TrueCarbonMass(iR);
           Duplicate_Blanks.MassPOC_R(iD,:) = R_blanks.MassPOC(iR);
           iD = iD + 1;
           
    end %endfor second for
    
end %endfor first for
clearvars 'Standards' 'Samples' 'Blanks' 'R_samples' 'R_blanks' 'Capsules' 'Filter_acidified' 'Filter_non_acidified' 'mdl' 'ifn' 'ifld' 'iR' 'Sfld' 'iD' 'iR' 'iID' 'iIDB' 'Bfld' ;
%% Dividing the samples into zone, Productive and Mesopelagic.

for idepth = 1:length(Duplicate_Samples.ID)
       if Duplicate_Samples.Depth(idepth) <= 200 %Zones : 1 for Productive zone; 2 for mesopelagic zone
          Duplicate_Samples.Zone(idepth) = 1;
       else
          Duplicate_Samples.Zone(idepth)= 2;
       end
end
Duplicate_Samples.Zone = Duplicate_Samples.Zone';
iP = find (Duplicate_Samples.Zone == 1);
iM = find (Duplicate_Samples.Zone == 2);
%% Calculate concentration of uPOC and aDOC in mg/m3
Duplicate_Samples.uPOC_conc = POC_concentration(Duplicate_Samples.CarbonMass,Duplicate_Samples.Volumen);
Duplicate_Samples.uPOC_conc_R = POC_concentration(Duplicate_Samples.CarbonMass_R,Duplicate_Samples.Volumen);
Duplicate_Blanks.aDOC_conc = POC_concentration(Duplicate_Blanks.CarbonMass,Duplicate_Blanks.Volumen);
Duplicate_Blanks.aDOC_conc_R = POC_concentration(Duplicate_Blanks.CarbonMass_R,Duplicate_Blanks.Volumen);

%% Calculating the difference and relative difference between duplicates

Duplicate_Samples.differenceabs = abs((Duplicate_Samples.POC - Duplicate_Samples.POC_R)/sqrt(2)); %absolute values of Scaled arithmetic difference
Duplicate_Samples.differenceabs_median = median (Duplicate_Samples.differenceabs);
Duplicate_Samples.differenceabs_std = mad(Duplicate_Samples.differenceabs);
Duplicate_Samples.difference_noabs = (Duplicate_Samples.POC - Duplicate_Samples.POC_R)/sqrt(2); %Scaled arithmetic difference
Duplicate_Samples.mean = ((Duplicate_Samples.POC + Duplicate_Samples.POC_R)/2);
Duplicate_Samples.median = median(Duplicate_Samples.mean);
Duplicate_Samples.relative_diffeabs = Duplicate_Samples.differenceabs(:,1)./ Duplicate_Samples.mean(:,1); %absolute values of Scaled relative difference
Duplicate_Samples.relative_diffen_noabs = Duplicate_Samples.difference_noabs(:,1)./ Duplicate_Samples.mean(:,1); % Scaled relative difference

Duplicate_Blanks.differenceabs = abs((Duplicate_Blanks.POC - Duplicate_Blanks.POC_R)/sqrt(2));
Duplicate_Blanks.differenceabs_median = median(Duplicate_Blanks.differenceabs);
Duplicate_Blanks.differenceabs_std = mad(Duplicate_Blanks.differenceabs);
Duplicate_Blanks.difference_noabs = (Duplicate_Blanks.POC - Duplicate_Blanks.POC_R)/sqrt(2);
Duplicate_Blanks.mean = ((Duplicate_Blanks.POC + Duplicate_Blanks.POC_R)/2);
Duplicate_Blanks.relative_diffeabs = Duplicate_Blanks.differenceabs(:,1)./ Duplicate_Blanks.mean(:,1);
Duplicate_Blanks.relative_diffen_noabs = Duplicate_Blanks.difference_noabs(:,1)./ Duplicate_Blanks.mean(:,1);

%% Calculating the relative uncertainty of duplicate estimates using percentil precision

% RUD = relative_uncertainty_duplicates(Std_RD, mean_D, D1, D2)
% RUD = Std_RD * mean_D * sqrt(2) / sqrt(D1.^2 + D2.^2);

%Productive Zone
for iRUD_P = 1:length(Duplicate_Samples.ID(iP))

    % relative uncertainty of each duplicate pair in the productive zone. 
    RUD_P_PP(iRUD_P, 1) = relative_uncertainty_duplicates(PP_P, Duplicate_Samples.mean(iRUD_P),Duplicate_Samples.POC(iRUD_P), Duplicate_Samples.POC_R(iRUD_P));
    
end

%Mesopelagic Zone
for iRUD_M = 1:length(Duplicate_Samples.ID(iM))

    % relative uncertainty of each duplicate pair in the mesopelagic zone. 
    RUD_M_PP(iRUD_M, 1) = relative_uncertainty_duplicates(PP_M, Duplicate_Samples.mean(iRUD_M),Duplicate_Samples.POC(iRUD_M), Duplicate_Samples.POC_R(iRUD_M));
    
end

%% Calculating the total uncertainty of the POC concentrations using Percentile precision

RU_P_PP = prctile(RUD_P_PP,50); %50th percentile in the relative uncertainty of each duplicate pair in the productive zone distribution 
RU_M_PP = prctile(RUD_M_PP,50); %50th percentile in the relative uncertainty of each duplicate pair in the mesopelagic zone distribution 

Duplicate_Samples.POCerr_PP = nan(size(Duplicate_Samples.POC)); 
Duplicate_Samples.POCerr_PP(iP) = RU_P_PP*Duplicate_Samples.POC(iP);
Duplicate_Samples.POCerr_PP(iM) = RU_M_PP*Duplicate_Samples.POC(iM);
