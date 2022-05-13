%Script for extract samples and blanks and their replicates into two new structures called.....
%...Duplicate_samples and Duplicate_Blanks.

%To do:
%Create a new structure with 7 fields RUN,ID,Depth,Volumen,PI,POC,Duplicate
%Create two new fields called difference and relative_diffe

%% 
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
clearvars 'Standards' 'Samples' 'Blanks' 'R_samples' 'R_blanks' 'Capsules' 'Filter_acidified' 'Filter_non_acidified' 'mdl' 'ifn' 'ifld' 'iR' 'Sfld' 'iD' 'iR' 'iID' 'iIDB' 'Bfld' 'filename' 'fn' 'flds' ;
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
clearvars 'idepth';

%% Calculate concentration of uPOC and aDOC in mg/m3
cd 'D:\Measuring POC_Paper\poc\POC\Scripts'

Duplicate_Samples.uPOC_conc = POC_concentration(Duplicate_Samples.CarbonMass,Duplicate_Samples.Volumen);
Duplicate_Samples.uPOC_conc_R = POC_concentration(Duplicate_Samples.CarbonMass_R,Duplicate_Samples.Volumen);

Duplicate_Blanks.aDOC_conc = POC_concentration(Duplicate_Blanks.CarbonMass,Duplicate_Blanks.Volumen);
Duplicate_Blanks.aDOC_conc_R = POC_concentration(Duplicate_Blanks.CarbonMass_R,Duplicate_Blanks.Volumen);

cd('D:\Measuring POC_Paper\poc\POC\Data');
cd

%% including time in the structure % added after create Paul_Results, so Paul_Results must be added
filename = 'Paul_Results.mat';
cd('D:\Measuring POC_Paper\poc\POC\Data\Total');
load(filename);

itime = [];
for iT = 1:length(Duplicate_Samples.ID)
    timefld = Duplicate_Samples.ID{iT}(end-11:end);
for i = 1:length(Paul_Results.ID)
    if contains (Paul_Results.ID(i),timefld)
        itime = [itime, i];
    end
end
end
Duplicate_Samples.Time = Paul_Results.Time(itime);

clearvars 'itime' 'iT' 'timefld' 'i' 'Paul_Results' Total_aDOC Total_POC Total_uPOC;

cd('D:\Measuring POC_Paper\poc\POC\Data');

%% creating an index for pre-dawn and noon casts

Duplicate_Samples.TimeNumber = hour(Duplicate_Samples.Time);
for itime = 1:length(Duplicate_Samples.ID)
       if Duplicate_Samples.TimeNumber(itime)<= 7 %{'07:00'} %Zones : 1 for Pre-Dawn; 2 for Noon
          Duplicate_Samples.Cast(itime) = 1;
       else
          Duplicate_Samples.Cast(itime)= 2;
       end
end
Duplicate_Samples.Cast = Duplicate_Samples.Cast';
idawn = find (Duplicate_Samples.Cast == 1);
inoon = find (Duplicate_Samples.Cast == 2);

%clearvars itime inoon idawn
clearvars itime 

%% Excluding bad runs
Bad_runs = ["Results_11"; "Results_12"];
ibad = contains(Duplicate_Samples.RUN, Bad_runs);
igood = ~ibad;

%% %% Calculating the difference and relative difference between duplicates

% Values of mean, median and std for all runs
Duplicate_Samples.mean = ((Duplicate_Samples.POC + Duplicate_Samples.POC_R)/2);
Duplicate_Samples.median = median(Duplicate_Samples.mean);
Duplicate_Samples.medianDawn = median(Duplicate_Samples.mean(idawn));
Duplicate_Samples.medianNoon = median(Duplicate_Samples.mean(inoon));
Duplicate_Samples.medianDawn_std = mad(Duplicate_Samples.mean(idawn),1);
Duplicate_Samples.medianNoon_std = mad(Duplicate_Samples.mean(inoon),1);

%absolute values of Scaled arithmetic difference for all runs
Duplicate_Samples.differenceabs = abs((Duplicate_Samples.POC - Duplicate_Samples.POC_R)/sqrt(2)); %absolute values of Scaled arithmetic difference
Duplicate_Samples.differenceabs_median = median(Duplicate_Samples.differenceabs);
Duplicate_Samples.differenceabs_std = mad(Duplicate_Samples.differenceabs,1);

%absolute values of Scaled arithmetic difference for good runs
Duplicate_Samples.differenceabs_good = nan(size(Duplicate_Samples.POC)); 
Duplicate_Samples.differenceabs_good(igood) = abs((Duplicate_Samples.POC(igood) - Duplicate_Samples.POC_R(igood))/sqrt(2));
Duplicate_Samples.differenceabs_median_good = median (Duplicate_Samples.differenceabs_good,'omitnan');
Duplicate_Samples.differenceabs_std_good = mad(Duplicate_Samples.differenceabs_good,1);

%Scaled arithmetic difference for all runs
Duplicate_Samples.difference_noabs = (Duplicate_Samples.POC - Duplicate_Samples.POC_R)/sqrt(2); %Scaled arithmetic difference

%Scaled arithmetic difference for good runs
Duplicate_Samples.difference_noabs_good = nan(size(Duplicate_Samples.POC)); 
Duplicate_Samples.difference_noabs_good(igood) = (Duplicate_Samples.POC(igood) - Duplicate_Samples.POC_R(igood))/sqrt(2);
Duplicate_Samples.difference_noabs_median_good = median (Duplicate_Samples.difference_noabs_good,'omitnan');
Duplicate_Samples.difference_noabs_std_good = mad(Duplicate_Samples.difference_noabs_good,1);


% Values of mean, median and std for good runs
Duplicate_Samples.mean_good = nan(size(Duplicate_Samples.POC)); 
Duplicate_Samples.mean_good(igood) = ((Duplicate_Samples.POC(igood) + Duplicate_Samples.POC_R(igood))/2);

%absolute and non-absolute values of Scaled relative difference for all runs
Duplicate_Samples.relative_diffeabs = Duplicate_Samples.differenceabs(:,1)./ Duplicate_Samples.mean(:,1); %absolute values of Scaled relative difference
Duplicate_Samples.relative_diffeabs_median = median (Duplicate_Samples.relative_diffeabs);
Duplicate_Samples.relative_diffeabs_median_std = mad(Duplicate_Samples.relative_diffeabs,1);

Duplicate_Samples.relative_diffen_noabs = Duplicate_Samples.difference_noabs(:,1)./ Duplicate_Samples.mean(:,1); % Scaled relative difference
Duplicate_Samples.relative_diffen_noabs_median = median (Duplicate_Samples.relative_diffen_noabs);
Duplicate_Samples.relative_diffen_noabs_median_std = mad(Duplicate_Samples.relative_diffen_noabs,1);

%absolute and non-absolute values of Scaled relative difference for good runs
Duplicate_Samples.relative_diffeabs_good = nan(size(Duplicate_Samples.POC)); 
Duplicate_Samples.relative_diffeabs_good(igood) = Duplicate_Samples.differenceabs(igood,1)./ Duplicate_Samples.mean_good(igood,1);
Duplicate_Samples.relative_diffeabs_good_median = median (Duplicate_Samples.relative_diffeabs_good,'omitnan');
Duplicate_Samples.relative_diffeabs_good_median_std = mad(Duplicate_Samples.relative_diffeabs_good,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Eq.8 in the new manuscript(16/03/2022) used to calculate Sigma r (Eq.9 (robust
%standard deviation of the relative duplicates in each zone)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Duplicate_Samples.relative_diffen_noabs_good = nan(size(Duplicate_Samples.POC));
Duplicate_Samples.relative_diffen_noabs_good(igood) = Duplicate_Samples.difference_noabs_good(igood,1)./ Duplicate_Samples.mean_good(igood,1);
Duplicate_Samples.relative_diffen_noabs_good_median = median (Duplicate_Samples.relative_diffen_noabs_good,'omitnan');
Duplicate_Samples.relative_diffen_noabs_good_median_std = mad(Duplicate_Samples.relative_diffen_noabs_good,1);

%values of Scaled arithmetic difference for blanks
Duplicate_Blanks.differenceabs = abs((Duplicate_Blanks.POC - Duplicate_Blanks.POC_R)/sqrt(2));
Duplicate_Blanks.differenceabs_median = median(Duplicate_Blanks.differenceabs);
Duplicate_Blanks.differenceabs_std = mad(Duplicate_Blanks.differenceabs,1);
Duplicate_Blanks.difference_noabs = (Duplicate_Blanks.POC - Duplicate_Blanks.POC_R)/sqrt(2);
Duplicate_Blanks.mean = ((Duplicate_Blanks.POC + Duplicate_Blanks.POC_R)/2);
Duplicate_Blanks.relative_diffeabs = Duplicate_Blanks.differenceabs(:,1)./ Duplicate_Blanks.mean(:,1);
Duplicate_Blanks.relative_diffen_noabs = Duplicate_Blanks.difference_noabs(:,1)./ Duplicate_Blanks.mean(:,1);


% [rho_p, pval_p] = corr(Duplicate_Samples.mean(iP), Duplicate_Samples.differenceabs(iP));
% [rho_m, pval_m] = corr(Duplicate_Samples.mean(iM), Duplicate_Samples.differenceabs(iM));


%% plotting the difference and relative difference
iP = find (Duplicate_Samples.Zone == 1);
iM = find (Duplicate_Samples.Zone == 2);

pd_iP = fitdist(Duplicate_Samples.relative_diffeabs(iP),'Normal'); %Create a normal distribution object by fitting it to the data.
pd_iM = fitdist(Duplicate_Samples.relative_diffeabs(iM),'Normal'); %Create a normal distribution object by fitting it to the data.

x_values = 0:0.001:1.5; %Define the input vector x_values to contain the values at which to calculate the pdf

y_iP = pdf(pd_iP,x_values);%Compute the pdf (Probability density function) values for the standard normal distribution at the values in x_values.
y_iM = pdf(pd_iM,x_values);

% mdl_productive_duplicate = fitlm (Duplicate_Samples.mean(iP), Duplicate_Samples.differenceabs(iP), 'linear', 'RobustOpts','andrews');
% IpValue = mdl_productive_duplicate.Coefficients{1,4}
% if IpValue > 0.05
%     mdl_productive_duplicate = fitlm (Duplicate_Samples.mean(iP), Duplicate_Samples.differenceabs(iP), 'linear','Intercept',false, 'RobustOpts','andrews');
% end
% 
% mdl_mesopelagic_duplicate = fitlm (Duplicate_Samples.mean(iM), Duplicate_Samples.differenceabs(iM), 'linear', 'RobustOpts','andrews');
% IpValue = mdl_mesopelagic_duplicate.Coefficients{1,4}
% if IpValue > 0.05
%     mdl_mesopelagic_duplicate = fitlm (Duplicate_Samples.mean(iM), Duplicate_Samples.differenceabs(iM), 'linear','Intercept',false, 'RobustOpts','andrews');
% end
% opt.RobustWgtFun = 'bisquare';
% opt.Tune = 6; % optional

mdl_productive_duplicate = fitlm (Duplicate_Samples.mean(iP), Duplicate_Samples.differenceabs(iP), 'linear', 'RobustOpts', 'on');
IpValue = mdl_productive_duplicate.Coefficients{1,4}
if IpValue > 0.05
    mdl_productive_duplicate = fitlm (Duplicate_Samples.mean(iP), Duplicate_Samples.differenceabs(iP), 'linear','Intercept',false, 'RobustOpts', 'on');
end

% mdl_mesopelagic_duplicate = fitlm (Duplicate_Samples.mean(iM), Duplicate_Samples.differenceabs(iM), 'linear', 'RobustOpts',opt);
% if IpValue > 0.05
%     mdl_mesopelagic_duplicate = fitlm (Duplicate_Samples.mean(iM), Duplicate_Samples.differenceabs(iM), 'linear','Intercept',false, 'RobustOpts', opt);
% end

mean_duplicate_prod = Duplicate_Samples.mean(iP);
mean_duplicate_meso = Duplicate_Samples.mean(iM);
Diff_Abs_meso = Duplicate_Samples.differenceabs(iM); % 1, 15, 29 (1, 31, 16)

mdl_mesopelagic_duplicate = fitlm (Duplicate_Samples.mean(iM), Duplicate_Samples.differenceabs(iM), 'linear', 'RobustOpts', 'on');
IpValue = mdl_mesopelagic_duplicate.Coefficients{1,4}
if IpValue > 0.05
    mdl_mesopelagic_duplicate = fitlm (Duplicate_Samples.mean(iM), Duplicate_Samples.differenceabs(iM) , 'linear','Intercept',false, 'RobustOpts', 'on');
end

% mdl_mesopelagic_duplicate = fitlm (mean_duplicate_meso, Diff_Abs_meso, 'linear', 'RobustOpts', 'on');
% IpValue = mdl_mesopelagic_duplicate.Coefficients{1,4}
% if IpValue > 0.05
%     mdl_mesopelagic_duplicate = fitlm (mean_duplicate_meso, Diff_Abs_meso , 'linear','Intercept',false, 'RobustOpts', 'on');
% end

POC_mock = linspace(0.8,100,1000);
POC_mock = POC_mock';

Duplicate_differenceabs_prod_mock = POC_mock * mdl_productive_duplicate.Coefficients{1,1};  
% Duplicate_differenceabs_prod_mock = POC_mock * mdl_productive_duplicate.Coefficients{1,1} + (POC_mock * mdl_productive_duplicate.Coefficients{2,1});  
Duplicate_differenceabs_meso_mock = POC_mock * mdl_mesopelagic_duplicate.Coefficients{1,1};  
% Duplicate_differenceabs_meso_mock = mdl_mesopelagic_duplicate.Coefficients{1,1} + (POC_mock * mdl_mesopelagic_duplicate.Coefficients{2,1});  


f2 = figure(2);
clf
subplot(1,2,1);
plot (Duplicate_Samples.mean(iP), Duplicate_Samples.differenceabs(iP), 'ks', 'MarkerFaceColor','k', 'MarkerSize',6);
% title('Difference versus Average POC Concentration');
xlabel('$\bar{D}$ $(mg/m^{3})$', 'FontSize', 20, 'FontAngle', 'italic', 'Interpreter', 'latex');
set(gca, 'FontSize', 16, 'box' ,'on', 'LineWidth', 1, 'Xcolor' ,'k', 'Ycolor', 'k');
ylabel('$\frac{|D_1 - D_2|}{\sqrt 2}$ $(mg/m^{3})$', 'FontSize', 20, 'FontAngle', 'italic', 'Interpreter', 'latex');
text(5,15,'(a)','FontSize', 20, 'Color', 'k');
xlim([0 60])
ylim([-0.5 16])
yline(0,'-.k');
    hold on
    plot (Duplicate_Samples.mean(iM), Duplicate_Samples.differenceabs(iM), 'ks', 'MarkerSize',6);
    plot (Duplicate_Samples.mean(ibad), Duplicate_Samples.differenceabs(ibad), 'ks', 'MarkerEdgeColor','m','LineWidth',0.8, 'MarkerSize',6);
% set (gcf, 'paperposition', [-7.57104234666667 5.70864534333333 36.1420833333333 18.2827083333333]);
    plot(POC_mock, Duplicate_differenceabs_prod_mock, 'r--', 'LineWidth', 1);
    plot(POC_mock, Duplicate_differenceabs_meso_mock, 'b-.', 'LineWidth', 1);
%     plot(mdl_productive_duplicate,'Marker','none', 'MarkerFaceColor', 'none');
%     plot(mdl_mesopelagic_duplicate, 'Marker','none', 'MarkerFaceColor', 'none');
%     xlabel('POC (mg/m^{3})', 'FontSize', 16, 'FontAngle', 'italic','Interpreter', 'tex');
%     set(gca, 'FontSize', 16, 'box' ,'on', 'LineWidth', 1, 'Xcolor' ,'k', 'Ycolor', 'k');
%     ylabel('$\frac{|D_1 - D_2|}{\sqrt 2}$ $(mg/m^{3})$', 'FontSize', 20, 'FontAngle', 'italic', 'Interpreter', 'latex');
%     legend('off')
%     title('')
    
subplot(1,2,2);
plot (Duplicate_Samples.mean(iP), Duplicate_Samples.relative_diffeabs(iP), 'ks', 'MarkerFaceColor','k', 'MarkerSize',6);
% xlabel('$\bar{D}$ (mg/m^{3})', 'FontSize', 16, 'FontAngle', 'italic', 'Interpreter', 'latex');
xlabel('$\bar{D}$ $(mg/m^{3})$', 'FontSize', 20, 'FontAngle', 'italic', 'Interpreter', 'latex');
set(gca,'YTick',(0:0.2:1.5),'FontSize', 16, 'box' ,'on', 'LineWidth', 1, 'Xcolor' ,'k', 'Ycolor', 'k');
ylabel('$\frac{|D_1 - D_2|}{\bar{D}\sqrt 2}$', 'FontSize', 20, 'FontAngle', 'italic', 'Interpreter', 'latex');
text(5,1.2,'(b)','FontSize', 20, 'Color', 'k');
xlim([0 60])
ylim([-0.1 1.3])
yline(0,'-.k');
    hold on
    plot (Duplicate_Samples.mean(iM), Duplicate_Samples.relative_diffeabs(iM), 'ks','MarkerSize',6);
    plot (Duplicate_Samples.mean(ibad), Duplicate_Samples.relative_diffeabs(ibad), 'ks','MarkerEdgeColor','m', 'LineWidth',0.8, 'MarkerSize',6);
    axes('Position',[.7 .7 .19 .19])
    box on 
    plot(y_iP, x_values, 'r--');
    set(gca,'YTick',(0:0.2:1.3), 'FontSize', 12, 'box' ,'on', 'LineWidth', 1, 'Xcolor' ,'k', 'Ycolor', 'k');
    ylabel('$\frac{|D_1 - D_2|}{\bar{D}\sqrt 2}$', 'FontSize', 16, 'FontAngle', 'italic', 'Interpreter', 'latex');
    hold on 
    plot(y_iM, x_values, 'b-.');
    legend ('PZ', 'MZ', 'box', 'off'); 
set (gcf, 'paperposition', [-7.57104234666667 5.70864534333333 36.1420833333333 18.2827083333333]);

% print (f1,'-dpng', strcat(figDir, 'difference_vs_POC_FittedLines_distribution'), '-r350');
% print (f1, '-depsc', strcat(figDir, 'difference_vs_POC_FittedLines_distribution'));
% print (f2,'-dpng', strcat(figDir, 'difference_vs_POC_FittedLines_nodistribution'), '-r350');
% print (f2, '-depsc', strcat(figDir, 'difference_vs_POC_FittedLines_nodistribution'));

%% Calculating the difference and relative difference between duplicates in aDOC and uPOC

Duplicate_Samples.differenceabs_uPOC_good = nan(size(Duplicate_Samples.POC)); 
Duplicate_Samples.differenceabs_uPOC_good(igood) = abs((Duplicate_Samples.uPOC_conc(igood) - Duplicate_Samples.uPOC_conc_R(igood))/sqrt(2)); %absolute values of Scaled arithmetic difference
Duplicate_Samples.differenceabs_median_uPOC_good = median (Duplicate_Samples.differenceabs_uPOC_good,'omitnan');
Duplicate_Samples.differenceabs_std_uPOC_good = mad(Duplicate_Samples.differenceabs_uPOC_good,1);

Duplicate_Samples.differenceabs_median_uPOC_good_productive = median(Duplicate_Samples.differenceabs_uPOC_good(iP),'omitnan');
Duplicate_Samples.differenceabs_std_uPOC_good_productive = mad(Duplicate_Samples.differenceabs_uPOC_good(iP),1);
Duplicate_Samples.differenceabs_median_uPOC_good_mesopelagic = median(Duplicate_Samples.differenceabs_uPOC_good(iM),'omitnan');
Duplicate_Samples.differenceabs_std_uPOC_good_mesopelagic = mad(Duplicate_Samples.differenceabs_uPOC_good(iM),1);

Duplicate_Samples.difference_noabs_uPOC_good = nan(size(Duplicate_Samples.POC)); 
Duplicate_Samples.difference_noabs_uPOC_good(igood) = (Duplicate_Samples.uPOC_conc(igood) - Duplicate_Samples.uPOC_conc_R(igood))/sqrt(2); %Scaled arithmetic difference
Duplicate_Samples.difference_noabs_median_uPOC_good = median (Duplicate_Samples.difference_noabs_uPOC_good,'omitnan');
Duplicate_Samples.difference_noabs_std_uPOC_good = mad(Duplicate_Samples.difference_noabs_uPOC_good,1);

Duplicate_Samples.difference_noabs_median_uPOC_good_productive = median(Duplicate_Samples.difference_noabs_uPOC_good(iP),'omitnan');
Duplicate_Samples.difference_noabs_std_uPOC_good_productive = mad(Duplicate_Samples.difference_noabs_uPOC_good(iP),1);
Duplicate_Samples.difference_noabs_median_uPOC_good_mesopelagic = median(Duplicate_Samples.difference_noabs_uPOC_good(iM),'omitnan');
Duplicate_Samples.difference_noabs_std_uPOC_good_mesopelagic = mad(Duplicate_Samples.difference_noabs_uPOC_good(iM),1);

tmp = corrcoef(Duplicate_Samples.difference_noabs_uPOC_good,Duplicate_Samples.Volumen,'Rows','pairwise');
tmp= tmp(2,1);% value = 0.2779
clearvars 'tmp'
tmp = corrcoef(Duplicate_Samples.differenceabs_uPOC_good,Duplicate_Samples.Volumen,'Rows','pairwise');
tmp= tmp(2,1); % value = -0.3002
clearvars 'tmp'

Duplicate_Blanks.differenceabs_aDOC_good = nan(size(Duplicate_Samples.POC)); 
Duplicate_Blanks.differenceabs_aDOC_good(igood) = abs((Duplicate_Blanks.aDOC_conc(igood) - Duplicate_Blanks.aDOC_conc_R(igood))/sqrt(2));
Duplicate_Blanks.differenceabs_median_aDOC_good = median(Duplicate_Blanks.differenceabs_aDOC_good,'omitnan');
Duplicate_Blanks.differenceabs_std_aDOC_good = mad(Duplicate_Blanks.differenceabs_aDOC_good,1);

Duplicate_Blanks.differenceabs_median_aDOC_good_productive = median(Duplicate_Blanks.differenceabs_aDOC_good(iP),'omitnan');
Duplicate_Blanks.differenceabs_std_aDOC_good_productive = mad(Duplicate_Blanks.differenceabs_aDOC_good(iP),1);
Duplicate_Blanks.differenceabs_median_aDOC_good_mesopelagic = median(Duplicate_Blanks.differenceabs_aDOC_good(iM),'omitnan');
Duplicate_Blanks.differenceabs_std_aDOC_mesopelagic = mad(Duplicate_Blanks.differenceabs_aDOC_good(iM),1);

 %% Calculating the relative uncertainty in POC (Sigma r) using fuction prcng in the relative difference of duplicate estimates
cd 'D:\Measuring POC_Paper\poc\POC\Scripts'

 %function out = prcrng(x)
Sr_P =  prcrng(Duplicate_Samples.relative_diffen_noabs_good(iP)); %robust standard deviation of relative duplicate differences for productive zone
Sr_M =  prcrng(Duplicate_Samples.relative_diffen_noabs_good(iM)); %robust standard deviation of relative duplicate differences for mesopelagic zone

cd('D:\Measuring POC_Paper\poc\POC\Data');
cd

%% Saving Data
dataDir = 'D:\Measuring POC_Paper\poc\POC\Data\Duplicates\';
fn = 'New_Duplicates';
fnmat = strcat (dataDir, fn);
save(fnmat, 'Duplicate_Samples','Duplicate_Blanks','Sr_P','Sr_M');
