% this function classifies patients according to Wada test results;
% Patients with bilateral LIs are classified as True negativ
% input = List of LIs of all 24 patients (one column)

function [NBroca,NWernicke,Class_Broca,Class_Wernicke]=Threshold(LIs_rawdata)
  
[LIs_rawdata]=createDataMatrix()
Thr=0.3
[corClass, pat_wada, Sensitivity, Specificity, PPV, NPV]=classification_incl_sensitivity_atypical(LIs_rawdata.Broca, Thr)
Class_Broca=[Sensitivity Specificity PPV NPV] % Parameters according to atpyical lateralization (right, bilateral)
NBroca=[corClass.Total.cases corClass.left.cases corClass.right.cases] % classification according to right lateralization (bilaterals are not considered)
[corClass, pat_wada, Sensitivity, Specificity, PPV, NPV]=classification_incl_sensitivity_atypical(LIs_rawdata.Wernicke, Thr)
Class_Wernicke=[Sensitivity Specificity PPV NPV]
NWernicke=[corClass.Total.cases corClass.left.cases corClass.right.cases]
end


function [LIs_rawdata]=createDataMatrix()
    PathResults='D:\Arbeit\LinuxExchange\Results';
    TargetColumn = 1;
    
    % LI_max sqrt:
    ColumnSize = 4;
    load (strcat(PathResults, '\LI_maxAct\LI_all_noise.mat'))
    LIs_rawdata.Broca(:,TargetColumn:ColumnSize)=LI_All_noise(:,1:4);
    LIs_rawdata.Wernicke(:,TargetColumn:ColumnSize)=LI_All_noise(:,5:8);
    for i=1:4
        LIs_rawdata.Name{i,1}=GetLineName(i)
    end
    clear LI_All_noise
    TargetColumn = TargetColumn + ColumnSize;
    
    %LI single Trials 320-600ms:
    ColumnSize = 1;
    load(strcat(PathResults, '\LI_MEG_singleTrials\patients_TTest_Int320To600ms.mat'))
    LI.Int320To600ms.Broca([2:3,5:6],:)=[];
    LI.Int320To600ms.Wernicke([2:3,5:6],:)=[];
    LIs_rawdata.Broca(:,TargetColumn)=LI.Int320To600ms.Broca(:,1);
    LIs_rawdata.Wernicke(:,TargetColumn)=LI.Int320To600ms.Wernicke(:,1);
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    clear LI
    TargetColumn = TargetColumn + ColumnSize;
    
    % LI single Trials 320-470ms:
    ColumnSize = 1;
    load(strcat(PathResults, '\LI_MEG_singleTrials\patients_TTest_Int320To470ms.mat'))
    LI.Int320To470ms.Broca([2:3,5:6],:)=[];
    LI.Int320To470ms.Wernicke([2:3,5:6],:)=[];
    LIs_rawdata.Broca(:,TargetColumn)=LI.Int320To470ms.Broca(:,1);
    LIs_rawdata.Wernicke(:,TargetColumn)=LI.Int320To470ms.Wernicke(:,1);
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI
    
    %LI_avg10 trials 320-600ms:
    ColumnSize = 1;
    load( strcat(PathResults, '\LI_UTest\patients\SumOfSignVoxels_dil_0.32_0.6_s.mat'))
    LIs_rawdata.Broca(:,TargetColumn)=LI_patients.Broca(:,1);
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_patients.Wernicke(:,1);
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_patients
    
    %LI_avg10 trials 320-470ms:
    ColumnSize = 1;
    load( strcat(PathResults, '\LI_UTest\patients\SumOfSignVoxels_dil_0.32_0.47_s.mat'))
    LIs_rawdata.Broca(:,TargetColumn)=LI_patients.Broca(:,1);
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_patients.Wernicke(:,1);
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_patients
    
    %LI RelAct_0.32_0.6s BL100ms no sqrt
    ColumnSize = 1;
    load (strcat(PathResults, '\LI_RelativeActSubstractionBL100ms\LI_All_RelAct_0.32_0.6s_Patients_BL100ms.mat'))
    LIs_rawdata.Broca(:,TargetColumn)=LI_All_RelAct(:,1)
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_All_RelAct(:,3)
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_All_RelAct
    
    %LI RelAct_0.32_0.47s BL100ms no sqrt
    ColumnSize = 1;
    load (strcat(PathResults, '\LI_RelativeActSubstractionBL100ms\LI_All_RelAct_0.32_0.47s_Patients_BL100ms.mat'))
    LIs_rawdata.Broca(:,TargetColumn)=LI_All_RelAct(:,1)
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_All_RelAct(:,3)
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_All_RelAct
    
     %LI RelAct_0.32_0.6s BL100ms sqrt
    ColumnSize = 1;
    load (strcat(PathResults, '\LI_RelativeActSubstractionBL100ms\LI_All_RelAct_0.32_0.6s_Patients_BL100ms.mat'))
    LIs_rawdata.Broca(:,TargetColumn)=LI_All_RelAct(:,2)
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_All_RelAct(:,4)
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_All_RelAct
    
    %LI RelAct_0.32_0.47s BL100ms sqrt
    ColumnSize = 1;
    load (strcat(PathResults, '\LI_RelativeActSubstractionBL100ms\LI_All_RelAct_0.32_0.47s_Patients_BL100ms.mat'))
    LIs_rawdata.Broca(:,TargetColumn)=LI_All_RelAct(:,2)
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_All_RelAct(:,4)
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_All_RelAct
    
    %%
    %LI RelAct_0.32_0.6s BL350ms no sqrt
    ColumnSize = 1;
    load (strcat(PathResults, '\LI_RelativeAct\LI_All_RelAct_0.32_0.6s_Patients.mat'))
    LI_All_RelAct([2,3,6],:)=[];
    LIs_rawdata.Broca(:,TargetColumn)=LI_All_RelAct(:,1)
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_All_RelAct(:,3)
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_All_RelAct
    
    %LI RelAct_0.32_0.47s BL350ms no sqrt
    ColumnSize = 1;
    load (strcat(PathResults, '\LI_RelativeAct\LI_All_RelAct_0.32_0.47s_Patients.mat'))
    LI_All_RelAct([2,3,5],:)=[];
    LIs_rawdata.Broca(:,TargetColumn)=LI_All_RelAct(:,1)
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_All_RelAct(:,3)
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_All_RelAct
    %%
    %LI RelAct Substraction 0.32_0.6s BL350ms no sqrt
    ColumnSize = 1;
    load (strcat(PathResults, '\LI_RelativeActSubstraction\LI_All_RelAct_0.32_0.6s_Patients.mat'))
    LIs_rawdata.Broca(:,TargetColumn)=LI_All_RelAct(:,1)
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_All_RelAct(:,3)
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_All_RelAct
    
    %LI RelAct Substraction 0.32_0.47s BL350ms no sqrt
    ColumnSize = 1;
    load (strcat(PathResults, '\LI_RelativeActSubstraction\LI_All_RelAct_0.32_0.47s_Patients.mat'))
    LIs_rawdata.Broca(:,TargetColumn)=LI_All_RelAct(:,1)
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_All_RelAct(:,3)
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_All_RelAct
    
   %% 
    %LI RelAct_0.32_0.6s BL350ms sqrt
    ColumnSize = 1;
    load (strcat(PathResults, '\LI_RelativeAct\LI_All_RelAct_0.32_0.6s_Patients.mat'))
    LI_All_RelAct([2,3,6],:)=[];
    LIs_rawdata.Broca(:,TargetColumn)=LI_All_RelAct(:,2)
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_All_RelAct(:,4)
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_All_RelAct
    
    %LI RelAct_0.32_0.47s BL350ms sqrt
    ColumnSize = 1;
    load (strcat(PathResults, '\LI_RelativeAct\LI_All_RelAct_0.32_0.47s_Patients.mat'))
    LI_All_RelAct([2,3,6],:)=[];
    LIs_rawdata.Broca(:,TargetColumn)=LI_All_RelAct(:,2)
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_All_RelAct(:,4)
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_All_RelAct
    
    %%
       %LI RelAct Substraction 0.32_0.6s BL350ms sqrt
    ColumnSize = 1;
    load (strcat(PathResults, '\LI_RelativeActSubstraction\LI_All_RelAct_0.32_0.6s_Patients.mat'))
    LIs_rawdata.Broca(:,TargetColumn)=LI_All_RelAct(:,2)
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_All_RelAct(:,4)
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_All_RelAct
    
    %LI RelAct Substraction 0.32_0.47s BL350ms sqrt
    ColumnSize = 1;
    load (strcat(PathResults, '\LI_RelativeActSubstraction\LI_All_RelAct_0.32_0.47s_Patients.mat'))
    LIs_rawdata.Broca(:,TargetColumn)=LI_All_RelAct(:,2)
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_All_RelAct(:,4)
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_All_RelAct
    
    %%
    %LI RelAct Z-scores voxelvalues 0.32_0.6s BL100ms 
    ColumnSize = 3;
    load (strcat(PathResults, '\LI_RelativeActSubstractionBL100ms\zscores_LI_All_RelAct_0.32_0.6s_Patients_BL100ms.mat'))
    LI_ALL_zscores([2 4 5 6],:)=[];
    LIs_rawdata.Broca(:,TargetColumn:(TargetColumn+2))=LI_ALL_zscores(:,[1 9 17])
    LIs_rawdata.Wernicke(:,TargetColumn:(TargetColumn+2))=LI_ALL_zscores(:,[5 13 21])
    for i=TargetColumn:(TargetColumn+2)
        LIs_rawdata.Name{i,1}=GetLineName(i)
    end
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_ALL_zscores
    
    %LI RelAct Z-scores voxelvalues 0.32_0.47s BL100ms 
    ColumnSize = 3;
    load (strcat(PathResults, '\LI_RelativeActSubstractionBL100ms\zscores_LI_All_RelAct_0.32_0.47s_Patients_BL100ms.mat'))
    LI_ALL_zscores([2 4 5 6],:)=[];
    LIs_rawdata.Broca(:,TargetColumn:(TargetColumn+2))=LI_ALL_zscores(:,[1 9 17])
    LIs_rawdata.Wernicke(:,TargetColumn:(TargetColumn+2))=LI_ALL_zscores(:,[5 13 21])
    for i=TargetColumn:(TargetColumn+2)
        LIs_rawdata.Name{i,1}=GetLineName(i)
    end
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_ALL_zscores
    %%
    %LI RelAct Z-scores voxelcount 0.32_0.6s BL100ms 
    ColumnSize = 3;
    load (strcat(PathResults, '\LI_RelativeActSubstractionBL100ms\zscores_LI_All_RelAct_0.32_0.6s_Patients_BL100ms.mat'))
    LI_ALL_zscores([2 4 5 6],:)=[];
    LIs_rawdata.Broca(:,TargetColumn:(TargetColumn+2))=LI_ALL_zscores(:,[2 10 18])
    LIs_rawdata.Wernicke(:,TargetColumn:(TargetColumn+2))=LI_ALL_zscores(:,[6 14 22])
    for i=TargetColumn:(TargetColumn+2)
        LIs_rawdata.Name{i,1}=GetLineName(i)
    end
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_ALL_zscores
    
    %LI RelAct Z-scores voxelcount 0.32_0.47s BL100ms 
    ColumnSize = 3;
    load (strcat(PathResults, '\LI_RelativeActSubstractionBL100ms\zscores_LI_All_RelAct_0.32_0.47s_Patients_BL100ms.mat'))
    LI_ALL_zscores([2 4 5 6],:)=[];
    LIs_rawdata.Broca(:,TargetColumn:(TargetColumn+2))=LI_ALL_zscores(:,[2 10 18])
    LIs_rawdata.Wernicke(:,TargetColumn:(TargetColumn+2))=LI_ALL_zscores(:,[6 14 22])
    for i=TargetColumn:(TargetColumn+2)
        LIs_rawdata.Name{i,1}=GetLineName(i)
    end
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_ALL_zscores
end


function [LineName] = GetLineName( LineIndex )
switch LineIndex
    case 1
        LineName = 'Sqrt of Max Amp. 320-600ms'
    case 2
        LineName = 'Sqrt of Sum of Max Amp. 320-600ms'
    case 3
        LineName = 'Sqrt of Sum of Max Amp. 320-470ms'
    case 4
        LineName = 'Sqrt of Sum of Max Amp. 400-600ms'
    case 5
        LineName = 'T-test Single Trials 320-600ms' 
    case 6
        LineName = 'T-test Single Trials 320-470ms'
    case 7
        LineName = 'Wilkoxon Test Avg 10 Trials 320-600ms'
    case 8
        LineName = 'Wilkoxon Test Avg 10 Trials 320-470ms'
    case 9
        LineName = 'BL 100 no sqrtMax of RelAct.NoiseEst + [Sum(abs(VG)-sum(abs(BL))]  320-600ms'
    case 10
        LineName = 'BL 100 no sqrtMax of RelAct. NoiseEst + [Sum(abs(VG)-sum(abs(BL))]  320-470 ms'
    case 11
        LineName = 'BL 100 sqrt Max of RelAct. NoiseEst + [Sum(abs(VG)-sum(abs(BL))]  320-600ms'
    case 12
        LineName = 'BL 100 sqrt Max of RelAct.  NoiseEst + [Sum(abs(VG)-sum(abs(BL))]  320-470 ms'
    case 13
        LineName = 'BL 350 no sqrt Max of RelAct. [mean(abs(VG)./mean(abs(BL))]  320-600ms'
    case 14
        LineName = 'BL 350 no sqrt Max of RelAct.  [mean(abs(VG)./mean(abs(BL))]   320-470ms'
    case 15
        LineName = 'BL 350 no sqrt Max of RelAct.NoiseEst + [Sum(abs(VG)-sum(abs(BL))]  320-600ms'
    case 16
        LineName = 'BL 350 no sqrt Max of RelAct. NoiseEst + [Sum(abs(VG)-sum(abs(BL))]  320-470 ms'
    case 17
        LineName = 'BL 350 sqrt Max of RelAct. [mean(abs(VG)./mean(abs(BL))]  320-600ms'
    case 18
        LineName = 'BL 350 sqrt Max of RelAct.  [mean(abs(VG)./mean(abs(BL))]   320-470ms'
    case 19
        LineName = 'BL 350 sqrt Max of RelAct. NoiseEst + [Sum(abs(VG)-sum(abs(BL))]  320-600ms'
    case 20
        LineName = 'BL 350 sqrt Max of RelAct.  NoiseEst + [Sum(abs(VG)-sum(abs(BL))]  320-470 ms'
    case 21
        LineName = 'z-distribution voxelvalues p<.05 320-600ms'
    case 22
        LineName = 'z-distribution voxelvalues p<.01 320-600ms'
    case 23
        LineName = 'z-distribution voxelvalues p<.001  320-600ms'
    case 24
        LineName = 'z-distribution voxelvalues p<.05 320-470ms'
    case 25
        LineName = 'z-distribution voxelvalues p<.01 320-470ms'
    case 26
        LineName = 'z-distribution voxelvalues p<.001  320-470ms'
    case 27
        LineName = 'z-distribution voxelcount p<.05 320-600ms'
    case 28
        LineName = 'z-distribution voxelcount p<.01 320-600ms'
    case 29
        LineName = 'z-distribution voxelcount p<.001  320-600ms'
    case 30
        LineName = 'z-distribution voxelcount p<.05 320-470ms'
    case 31
        LineName = 'z-distribution voxelcount p<.01 320-470ms'
    case 32
        LineName = 'z-distribution voxelcount p<.001  320-470ms'
    otherwise
        LineName = strcat('FEHLER: GetLineName mit unbekanntem Index aufgerufen: ', num2str(LineIndex))
    end
end


function [corClass, pat_wada, Sensitivity, Specificity, PPV, NPV]=classification_incl_sensitivity_atypical (values, CutOff)

    CutOff_left=CutOff;
    CutOff_right=CutOff*(-1);
    
    if size(values,1)>24
        disp ('Warning: length of input data is larger than 24 - aborting')
        return
    end
    
    % classification wada test Patient 1-24:

    Wada= {'right' 'left' 'left' 'left' 'kein Wada' 'right' 'kein Wada' 'right' 'kein Wada' 'right' 'kein Wada' ...
        'left' 'left' 'right' 'kein Wada' 'kein Wada' 'left' 'kein Wada' 'kein Wada' 'left' 'right' 'kein Wada' 'left' 'kein Wada'}
    wada_no=[5 7 9 11 15 16 18 19 22 24];
    wada_yes=[1:4, 6, 8, 10, 12:14, 17, 20, 21, 23];
    result_wada=Wada(wada_yes)';
    pat_wada.data=values(wada_yes, :);
    
    % categorisation of LI in left, right and bilateral according to values:
    
    for i=1:size(pat_wada.data,1) % length of row
        for j=1:size(pat_wada.data,2) % size of columns of data
            if pat_wada.data(i,j)<=CutOff_right
                pat_wada.nom{i,j}='right';
            elseif pat_wada.data(i,j)>=CutOff_left
                pat_wada.nom{i,j}='left';
            else
                pat_wada.nom{i,j}='bilateral';
            end
        end
    end
   % classification according to correctness:
   for i=1:size(pat_wada.data,1)
       for j=1:size(pat_wada.data,2)
           if 1==strcmp(result_wada{i}, pat_wada.nom{i,j})
               pat_wada.classification{i,j}='correct';
           else
               pat_wada.classification{i,j}='incorrect';
           end
       end
   end
   
  % Percentage correct classified 
 
  corClass   = [];
  corClass.Total.cases(:,1)=sum(ismember(pat_wada.classification, 'correct'))
  corClass.Total.percent(:,1)=(sum(ismember(pat_wada.classification, 'correct')))/length(pat_wada.classification)
  
  % Percentage correct classified of Left dominant patients (according to
  % Wada): 
 [corClass] = classification_dominance (result_wada, pat_wada, 'left', corClass) 
 [corClass] = classification_dominance (result_wada, pat_wada, 'right', corClass) 
 [Sensitivity, Specificity, PPV, NPV]=sensitivity (corClass, pat_wada, result_wada, CutOff)
    
end

function [corClass] = classification_dominance (result_wada, pat_wada, hem, corClass )
    [ind]= ismember(result_wada, hem)
    row=find(ind==1)
    for i=1:size(pat_wada.classification,2)
        Wada_hem=pat_wada.classification(row,i)
        corClass.(hem).cases(i,1)=sum(ismember(Wada_hem, 'correct'))
        corClass.(hem).percent(i,1)=(sum(ismember(Wada_hem, 'correct'))/length(Wada_hem))
    end
end



function [Sensitivity, Specificity, PPV, NPV]=sensitivity (corClass,  pat_wada, result_wada, CutOff)

    % corrClass: column 1: correct classified total, column 2: correct
    % classified left wada; column 3: correct classified atypical patients
    CuffOff_left=CutOff;
    CutOff_right=CutOff*(-1);
    
    for i=1:length(corClass.left.cases)
        ntotal=14;
        nleft=8;
        natypical=6;
        
        [ind_leftwada]=ismember(result_wada, 'left')
        LIs_Wada_left=pat_wada.data(find(ind_leftwada==1), i);
        [ind_rightwada]=ismember(result_wada, 'right')
        LIs_Wada_right=pat_wada.data(find(ind_rightwada==1), i);
        
        TP = sum(LIs_Wada_right<CuffOff_left); % True positive,  all patients detected as atypical by fMRI that have atypical Wada test
        FN = sum(LIs_Wada_right>=CuffOff_left);    %natypical-TP;   % false negative
        TN = corClass.left.cases(i,1); % True negative, all patients with left fMRI that also have left wada test
        FP = sum(LIs_Wada_left<CuffOff_left)       %       nleft-TN;       % False positive
        Sensitivity(i,1) = (TP/(TP+FN));
        Specificity(i,1) = (TN/(FP+TN));
        PPV(i,1) = TP/(TP+FP);     % positive predictive value
        NPV(i,1) = TN/(TN+FN);     % negative predictive value
    end
end