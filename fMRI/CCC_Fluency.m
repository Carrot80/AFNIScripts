function patients ()
   
PatientFolder = '/media/truecrypt7/kirsten_thesis/data/patients';
ControlsFolder = '/media/truecrypt7/kirsten_thesis/data/controls';

for_all (PatientFolder, 'patients')
% for_all (ControlsFolder, 'controls')

end


function for_all (Folder, group)


    DIR = dir (Folder)
    isub = [DIR(:).isdir]; %  returns logical vector
    nameFolds = {DIR(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];
    PearsonR=zeros(size(nameFolds),5);
    CCC=zeros(size(nameFolds),5);

    for i= 1:size(nameFolds)
      
   [PearsonR, CCC]=kh_corr(strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, i, group, 'Verbgeneration', PearsonR, CCC, 'Broca_left_dil', 1)  
   [PearsonR, CCC]=kh_corr(strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, i, group, 'Verbgeneration', PearsonR, CCC, 'Broca_right_dil', 2)  
   [PearsonR, CCC]=kh_corr(strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, i, group, 'Verbgeneration', PearsonR, CCC, 'Wernicke_left_dil', 3)  
   [PearsonR, CCC]=kh_corr(strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, i, group, 'Verbgeneration', PearsonR, CCC, 'Wernicke_right_dil', 4)  
   [PearsonR, CCC]=kh_corr(strcat(Folder, filesep, nameFolds{i,1}), nameFolds{i}, i, group, 'Verbgeneration', PearsonR, CCC, 'ROIs_complete', 5)  
   
    end

    PathPearsonR=strcat('/home/kh/ShareWindows/Corr', filesep, 'PersonR_', group,'_ROIs_fMRI_Fluency_VG.mat')
    save(PathPearsonR, 'PearsonR')
    PathCCC=strcat('/home/kh/ShareWindows/Corr', filesep, 'CCC_', group,'_ROIs_fMRI_Fluency_VG.mat')
    save(PathCCC, 'CCC')
    
end


function [PearsonR, CCC]=kh_corr(SubjectPath, SubjectName, i, group, Task, PearsonR, CCC, ROI, A)


PathMask = strcat ('/home/kh/ShareWindows/Masks_fMRI/', ROI, '_fMRI+tlrc');

% Maske auf T-Verteilung anpassen:
%     PathMask = strcat ('/home/kh/ShareWindows/data/patients/ROIs_complete_MEG+tlrc');
%     eval(['!3dresample -master wspmT_0001+tlrc', ' -prefix ', ' Broca_left_dil_fMRI', ' -inset ', 'Broca_left_dil+tlrc' ]) 

%%
% PathFluency = strcat(SubjectPath, filesep, 'fMRI/statistics/Fluency/');
% cd(PathFluency)
% 
% delete 'wspmT_0001+tlrc.BRIK'
% delete 'wspmT_0001+tlrc.HEAD'
% ! 3dcopy wspmT_0001.hdr wspmT_0001
% 
% PathVG = strcat(SubjectPath, filesep, 'fMRI/statistics/Verbgeneration/');
% cd(PathVG)
% % 
% delete 'wspmT_0001+tlrc.BRIK'
% delete 'wspmT_0001+tlrc.HEAD'
% ! 3dcopy wspmT_0001.hdr wspmT_0001
%  eval(['!3dresample -master ', PathMask, ' -prefix ', ' wspmT_0001_res', ' -inset ', 'wspmT_0001+tlrc' ]) 
% delete 'wspmT_0001+tlrc.BRIK'
% delete 'wspmT_0001+tlrc.HEAD'


%%
PathFluency = strcat(SubjectPath, filesep, 'fMRI/statistics/Fluency/wspmT_0001+tlrc');

[V_Fluency, Info_Fluency] = BrikLoad (PathFluency);

[V_Mask, Info_Mask] = BrikLoad (PathMask); 

Act_Fluency = V_Mask.*V_Fluency;

PathVG=strcat(SubjectPath, filesep, 'fMRI/statistics/Verbgeneration/wspmT_0001+tlrc')

[V_VG, Info_VG] = BrikLoad (PathVG);

Act_VG=V_Mask.*V_VG;

Act_Fluency=squeeze(Act_Fluency(:));
Act_VG=squeeze(Act_VG(:));
[ind_Fluency]=find(Act_Fluency==0);
[ind_VG]=find(Act_VG==0);

if 1==length(ind_Fluency)~=length(ind_VG)
    Act_Fluency(ind_Fluency)=NaN;
    Act_VG(ind_VG)=NaN;
    Pcorr=corr(Act_VG, Act_Fluency, 'rows', 'pairwise'); 
    Act_Fluency(ind_Fluency)=[];
    Act_VG(ind_VG)=[];
    ccc=[2*Pcorr.*std(Act_VG).*std(Act_Fluency)]./[std(Act_VG)^2+std(Act_Fluency)^2+(mean(Act_VG)-mean(Act_Fluency))^2 ];
%stimmt noch etwas nicht!
else
    Act_Fluency(ind_Fluency)=[];
    Act_VG(ind_VG)=[];
    Pcorr=corr(Act_VG, Act_Fluency); % normale Person Correlation
    ccc=[2*corr(Act_VG,Act_Fluency).*std(Act_VG).*std(Act_Fluency)]./[std(Act_VG)^2+std(Act_Fluency)^2+(mean(Act_VG)-mean(Act_Fluency))^2 ];

end
 

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