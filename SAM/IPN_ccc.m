function patients ()
   
PatientFolder = '/home/kh/ShareWindows/data/patients/patients_SAM';
ControlsFolder = '/home/kh/ShareWindows/data/controls/controls_SAM';

for_all (PatientFolder, 'patients')
for_all (ControlsFolder, 'controls')

end


function for_all (Folder, group)


    DIR = dir (Folder)
    isub = [DIR(:).isdir]; %  returns logical vector
    nameFolds = {DIR(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];
    PearsonR=zeros(size(nameFolds),5);
    CCC=zeros(size(nameFolds),5);

    for i= 1:size(nameFolds)
   

    TimeInt = [.32, .6];
      
   [PearsonR, CCC]=kh_corr(strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, i, TimeInt, group, 'Verbgeneration', PearsonR, CCC, 'Broca_left_dil', 1)  
   [PearsonR, CCC]=kh_corr(strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, i, TimeInt, group, 'Verbgeneration', PearsonR, CCC, 'Broca_right_dil', 2)  
   [PearsonR, CCC]=kh_corr(strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, i, TimeInt, group, 'Verbgeneration', PearsonR, CCC, 'Wernicke_left_dil', 3)  
   [PearsonR, CCC]=kh_corr(strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, i, TimeInt, group, 'Verbgeneration', PearsonR, CCC, 'Wernicke_right_dil', 4)  
   [PearsonR, CCC]=kh_corr(strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, i, TimeInt, group, 'Verbgeneration', PearsonR, CCC, 'ROIs_complete_MEG', 5)  
   
    end

    PathPearsonR=strcat('/home/kh/ShareWindows/Corr', filesep, 'PersonR_', group,'_ROIs_MEG_abs.mat')
    save(PathPearsonR, 'PearsonR')
    PathCCC=strcat('/home/kh/ShareWindows/Corr', filesep, 'CCC_', group,'_ROIs_MEG_abs.mat')
    save(PathCCC, 'CCC')
    
end


function [PearsonR, CCC]=kh_corr(SubjectPath, SubjectName, i, TimeInt, group, Task, PearsonR, CCC, ROI, A)

if 1==strcmp(SubjectName, 'Pat_02_13008rh_1') || 1==strcmp(SubjectName, 'Pat_02_13008rh_2') || 1==strcmp(SubjectName,'Pat_03_13014bg_1') || 1==strcmp(SubjectName,'Pat_03_13014bg_2')
    return
end

PathMask = strcat ('/home/kh/ShareWindows/data/patients/', ROI, '+tlrc');

if 1==strcmp(group, 'patients')
    Path2Subj=strcat(SubjectPath, filesep, 'TimeIntervalls');
else
    Path2Subj=strcat(SubjectPath, filesep, 'SAM');
end

cd(Path2Subj)

% if ~exist ('fMRI_comparison', 'dir')
%     mkdir('fMRI_comparison')
% end

if 1==strcmp (SubjectName,'Pat_02_13008rh') ||  1==strcmp (SubjectName,'Pat_03_13014bg')
    PathERF = strcat('BothRuns_br_z_transf_brain01ERF_noise_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's', '+tlrc');
elseif 1==strcmp(group, 'Patients')
    PathERF = strcat('br_z_transf_brain01ERF_noise_abs_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's_', SubjectName, '+tlrc');
else
    PathERF = strcat('br_z_transf_brain01ERF_noise_abs_', num2str(TimeInt(1)), '-', num2str(TimeInt(2)), 's_', SubjectName, '+tlrc');
end

[V_ERF, Info_ERF] = BrikLoad (PathERF);

[V_Mask, Info_Mask] = BrikLoad (PathMask); 

Act_MEG = V_Mask.*V_ERF;

Path_fMRI=strcat('/media/truecrypt7/kirsten_thesis/data/', group, filesep, SubjectName,'/fMRI/statistics/', Task);

[V_fMRI, Info_fMRI] = BrikLoad (strcat(Path_fMRI, filesep, 'wspmT_0001_downsampled_z_transf+tlrc'));

Act_fMRI=V_Mask.*V_fMRI;

Act_MEG=squeeze(Act_MEG(:));
Act_fMRI=squeeze(Act_fMRI(:));
[ind]=find(Act_MEG==0);
Act_MEG(ind)=[];
ind=[];
[ind]=find(Act_fMRI==0);
Act_fMRI(ind)=[];


%%


% [ind_Act_MEG_p05]=find(Act_MEG>1.644853)
% Act_MEG_p05=Act_MEG(ind_Act_MEG_p05)
% 
% [ind_Act_fMRI_p05]=find(Act_fMRI>1.644853)
% Act_fMRI_p05=Act_fMRI(ind_Act_fMRI_p05)
% 
% Pcorr=corr(Act_fMRI_p05, Act_MEG_p05); % normale Person Correlation
% 
% ccc=[2*corr(Act_fMRI,Act_MEG).*std(Act_fMRI).*std(Act_MEG)]./[std(Act_fMRI)^2+std(Act_MEG)^2+(mean(Act_fMRI)-mean(Act_MEG))^2 ];


%%


Pcorr=corr(Act_fMRI, Act_MEG); % normale Person Correlation

ccc=[2*corr(Act_fMRI,Act_MEG).*std(Act_fMRI).*std(Act_MEG)]./[std(Act_fMRI)^2+std(Act_MEG)^2+(mean(Act_fMRI)-mean(Act_MEG))^2 ];

PearsonR(i,A)=Pcorr;
CCC(i,A)=ccc;

% X=LeftAct_fMRI;
% Y=[Act_fMRI Act_MEG]
% ssd = IPN_ssd(X)
% rc = kh_IPN_ccc(Y, ssd)



end



%% Computes the concordance correlation coefficient for evaluating reproducibility.
function rc = kh_IPN_ccc(Y,ssd)
% INPUT:
%   Y - a N*R data matrix
%
% REFERENCE:
%   Lin, L.I. 1989. A Corcordance Correlation Coefficient to Evaluate
%   Reproducibility. Biometrics 45, 255-268.
%
% XINIAN ZUO 2008
% zuoxinian@gmail.com

Ybar = mean(Y);S = cov(Y,1);R = size(Y,2);
tmp = triu(S,1);
rc = 2*sum(tmp(:))/((R-1)*trace(S)+IPN_ssd(Ybar));

end


%% Computes the sum square distance for evaluating reproducibility.
function ssd = IPN_ssd(X)
% INPUT:
%   X - a 1*R data vector
%
% REF:
%   Lin, L.I. 1989. A Corcordance Correlation Coefficient to Evaluate
%   Reproducibility. Biometrics 45, 255-268.
%
% XINIAN ZUO 2008
% zuoxinian@gmail.com

R=length(X);ssd=0;
for k=1:R-1
    ssd=ssd+sum((X(k+1:R)-X(k)).*(X(k+1:R)-X(k)));
end


end