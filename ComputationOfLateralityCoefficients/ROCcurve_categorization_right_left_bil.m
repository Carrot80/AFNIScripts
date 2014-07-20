% this function classifies patients according to Wada test results;
% patients with bilateral LIs are classified as True negativ
% input = List of LIs of all 24 patients (one column)

function [NBroca,NWernicke,Class_Broca,Class_Wernicke]=Threshold(LIs_rawdata)
    group='patients'; % only works for patients!      
    [LIs_rawdata]=createDataMatrix(group)
    corClass=[];
    Sensitivity=[];
    Specificity=[];
    PPV=[];
    NPV=[];
    CutOffsAll= 0.0:0.05:1% 0.00:0.05:1 %change
    sizeCutOffsAll=length(CutOffsAll)
    IntCutOffs=1
    for CutOff=CutOffsAll
        [corClass, Sensitivity, Specificity, PPV, NPV]=classification_incl_sensitivity_atypical(LIs_rawdata.Broca, CutOff, sizeCutOffsAll, IntCutOffs, corClass, Sensitivity, Specificity, PPV, NPV)
        IntCutOffs=IntCutOffs+1
    end
    Broca.corClass=corClass;
    Broca.Sensitivity=Sensitivity;
    Broca.Specificity=Specificity;
    Broca.PPV=PPV;
    Broca.NPV=NPV;
        
    % Das gleiche f�r Wernicke:
    corClass=[];
    Sensitivity=[];
    Specificity=[];
    PPV=[];
    NPV=[];
    IntCutOffs=1    
    for CutOff=CutOffsAll
        [corClass, Sensitivity, Specificity, PPV, NPV]=classification_incl_sensitivity_atypical(LIs_rawdata.Wernicke, CutOff, sizeCutOffsAll, IntCutOffs, corClass, Sensitivity, Specificity, PPV, NPV)
    IntCutOffs=IntCutOffs+1
    end
    Wernicke.corClass=corClass;
    Wernicke.Sensitivity=Sensitivity;
    Wernicke.Specificity=Specificity;
    Wernicke.PPV=PPV;
    Wernicke.NPV=NPV;
      
    % finds N with no activation (out of N=24):
    [No_Act]=findNoAct(LIs_rawdata) 
    plotROC (Broca, Wernicke, LIs_rawdata.Name, CutOffsAll)
end


function plotROC (Broca, Wernicke, Method, CutOffsAll)

% Broca:

Mean_Broca_Sens=mean(Broca.Sensitivity)
Mean_Broca_Spec=mean(Broca.Specificity)
Mean_Wernicke_Sens=mean(Wernicke.Sensitivity)
Mean_Wernicke_Spec=mean(Wernicke.Specificity)

Max_Broca_Sens=max(Broca.Sensitivity)
Max_Broca_Spec=max(Broca.Specificity)
Max_Wernicke_Sens=max(Wernicke.Sensitivity)
Max_Wernicke_Spec=max(Wernicke.Specificity)

[row, column]=find(Broca.Sensitivity==1)
Broca.Specificity(row, column)

% Plot four best methods (Rel Act):
for i=1:4
    ind_maxBroca=[12 24 25 31]
    subplot(2,2,i)
    plot(1-Broca.Specificity(:,ind_maxBroca(i)), Broca.Sensitivity(:,ind_maxBroca(i)))
    xlim([0 1])
    ylim([0 1])
    hold on
    plot(1-Wernicke.Specificity(:,ind_maxBroca(i)), Wernicke.Sensitivity(:,ind_maxBroca(i)), 'r')
    hold on
    % bisecting angle:
    bistect=[0:0.1:1]
    plot(bistect,bistect, 'k')
    xlabel('1 - Specificity');
    ylabel('Sensitvity');   
%     title(Method{12})
end

Youden_Broca=Broca.Sensitivity+Broca.Specificity-1
[MaxYouden_Broca MaxYouden_BrocaInd]=max(Youden_Broca)

Youden_Wernicke=Broca.Sensitivity+Wernicke.Specificity-1;
[MaxYouden_Wernicke MaxYouden_WernickeInd]=max(Youden_Wernicke);

% find corresponding cut off value for max Youden index:
for i=1:length(MaxYouden_BrocaInd)
    CuttOffValue_Broca(i)=CutOffsAll(MaxYouden_BrocaInd(i))
end

for i=1:length(MaxYouden_WernickeInd)
    CuttOffValue_Wernicke(i)=CutOffsAll(MaxYouden_WernickeInd(i));
end

%%
% Werte mit h�chsten Sensitivit�t Broca:
[Max_Sensitivity_Broca, indMax_Sensitivity_Broca]=max(Broca.Sensitivity);
for i=1:32   
Rel_specificity(i)=Broca.Specificity(indMax_Sensitivity_Broca(i),i)
end
J_Broca=Max_Sensitivity_Broca+Rel_specificity-1 % Youden's Index

for i=1:length(indMax_Sensitivity_Broca)
    CuttOffValue_Sens_maximized_Broca(i)=CutOffsAll(indMax_Sensitivity_Broca(i))
end

% Werte mit h�chsten Sensitivit�t Wernicke:
[Max_Sensitivity_Wernicke, indMax_Sensitivity_Wernicke]=max(Wernicke.Sensitivity);
for i=1:32   
Rel_specificity_Wernicke(i)=Wernicke.Specificity(indMax_Sensitivity_Wernicke(i),i)
end
J_Wernicke=Max_Sensitivity_Wernicke+Rel_specificity_Wernicke-1 % Youden's Index
for i=1:length(indMax_Sensitivity_Wernicke)
    CuttOffValue_Sens_maximized_Wernicke(i)=CutOffsAll(indMax_Sensitivity_Wernicke(i))
end

% compare Time intervalls:
[p,h,stats]=signrank(MaxYouden_Broca(1, [2,5,7,9,11,21,22,23,27,28,29]),MaxYouden_Broca(1, [3,6,8,10,12,24,25,26,30,31,32])) 


figure
for i=1:size(Broca.Sensitivity,2)
    Sensitivity=Broca.Sensitivity(:,i)
    Specificity=Broca.Specificity(:,i)
    
    subplot(4,6,i)
    plot(1-Specificity, Sensitivity)
    xlim([0 1])
    ylim([0 1])
    hold on
    plot(1-Wernicke.Specificity(:,i), Wernicke.Sensitivity(:,i), 'r')
    hold on
    % bisecting angle:
    bistect=[0:0.1:1]
    plot(bistect,bistect, 'k')
    xlabel('1 - Specificity');
    ylabel('Sensitvity');   
end
end


function [No_Act]=findNoAct(LIs_rawdata)

No_Act.Broca(:,1)=sum(isnan(LIs_rawdata.Broca));
No_Act.Wernicke(:,1)=sum(isnan(LIs_rawdata.Wernicke));
% [Row, Col]=find(isnan(LIs_rawdata.Broca)==1)
end


function [LIs_rawdata]=createDataMatrix(group)
    PathResults='D:\Arbeit\LinuxExchange\Results';
    TargetColumn = 1;
    
    % LI_max sqrt:
    ColumnSize = 4;
    load (strcat(PathResults, '\LI_maxAct\LI_all_noise_', group, '.mat'))
    LIs_rawdata.Broca(:,TargetColumn:ColumnSize)=LI_All_noise(:,1:4);
    LIs_rawdata.Wernicke(:,TargetColumn:ColumnSize)=LI_All_noise(:,5:8);
    for i=1:4
        LIs_rawdata.Name{i,1}=GetLineName(i)
    end
    clear LI_All_noise
    TargetColumn = TargetColumn + ColumnSize;
    
    %LI single Trials 320-600ms:
    ColumnSize = 1;
    load(strcat(PathResults, '\LI_MEG_singleTrials\', group, '_TTest_Int320To600ms.mat'))
    if 1==strcmp(group, 'patients');
        LI.Int320To600ms.Broca([2:3,5:6],:)=[];
        LI.Int320To600ms.Wernicke([2:3,5:6],:)=[];
    end
    LIs_rawdata.Broca(:,TargetColumn)=LI.Int320To600ms.Broca(:,1);
    LIs_rawdata.Wernicke(:,TargetColumn)=LI.Int320To600ms.Wernicke(:,1);
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    clear LI
    TargetColumn = TargetColumn + ColumnSize;
    
    % LI single Trials 320-470ms:
    ColumnSize = 1;
    load(strcat(PathResults, '\LI_MEG_singleTrials\', group, '_TTest_Int320To470ms.mat'))
    if 1==strcmp(group, 'patients');
        LI.Int320To470ms.Broca([2:3,5:6],:)=[];
        LI.Int320To470ms.Wernicke([2:3,5:6],:)=[];
    end
    LIs_rawdata.Broca(:,TargetColumn)=LI.Int320To470ms.Broca(:,1);
    LIs_rawdata.Wernicke(:,TargetColumn)=LI.Int320To470ms.Wernicke(:,1);
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI
    
    %LI_avg10 trials 320-600ms:
    ColumnSize = 1;
    load( strcat(PathResults, '\LI_UTest\', group, '\SumOfSignVoxels_dil_0.32_0.6_s_', group,'.mat'))
    LIs_rawdata.Broca(:,TargetColumn)=LI.Broca(:,1);
    LIs_rawdata.Wernicke(:,TargetColumn)=LI.Wernicke(:,1);
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI
    
    %LI_avg10 trials 320-470ms:
    ColumnSize = 1;
    load( strcat(PathResults, '\LI_UTest\', group, '\SumOfSignVoxels_dil_0.32_0.47_s_', group, '.mat'))
    LIs_rawdata.Broca(:,TargetColumn)=LI.Broca(:,1);
    LIs_rawdata.Wernicke(:,TargetColumn)=LI.Wernicke(:,1);
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI
    
    %LI RelAct_0.32_0.6s BL100ms no sqrt
    ColumnSize = 1;
    load (strcat(PathResults, '\LI_RelativeActSubstractionBL100ms\LI_All_RelAct_0.32_0.6s_', group, '_BL100ms.mat'))
    LIs_rawdata.Broca(:,TargetColumn)=LI_All_RelAct(:,1)
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_All_RelAct(:,3)
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_All_RelAct
    
    %LI RelAct_0.32_0.47s BL100ms no sqrt
    ColumnSize = 1;
    load (strcat(PathResults, '\LI_RelativeActSubstractionBL100ms\LI_All_RelAct_0.32_0.47s_', group, '_BL100ms.mat'))
    LIs_rawdata.Broca(:,TargetColumn)=LI_All_RelAct(:,1)
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_All_RelAct(:,3)
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_All_RelAct
    
     %LI RelAct_0.32_0.6s BL100ms sqrt
    ColumnSize = 1;
    load (strcat(PathResults, '\LI_RelativeActSubstractionBL100ms\LI_All_RelAct_0.32_0.6s_', group, '_BL100ms.mat'))
    LIs_rawdata.Broca(:,TargetColumn)=LI_All_RelAct(:,2)
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_All_RelAct(:,4)
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_All_RelAct
    
    %LI RelAct_0.32_0.47s BL100ms sqrt
    ColumnSize = 1;
    load (strcat(PathResults, '\LI_RelativeActSubstractionBL100ms\LI_All_RelAct_0.32_0.47s_', group, '_BL100ms.mat'))
    LIs_rawdata.Broca(:,TargetColumn)=LI_All_RelAct(:,2)
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_All_RelAct(:,4)
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_All_RelAct
    
    %%
    %LI RelAct_0.32_0.6s BL350ms no sqrt
    ColumnSize = 1;
    load (strcat(PathResults, '\LI_RelativeAct\LI_All_RelAct_0.32_0.6s_', group, '.mat'))
    if 1==strcmp(group, 'patients');
        LI_All_RelAct([2,3,6],:)=[];
    end
    LIs_rawdata.Broca(:,TargetColumn)=LI_All_RelAct(:,1)
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_All_RelAct(:,3)
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_All_RelAct
    
    %LI RelAct_0.32_0.47s BL350ms no sqrt
    ColumnSize = 1;
    load (strcat(PathResults, '\LI_RelativeAct\LI_All_RelAct_0.32_0.47s_', group, '.mat'))
    if 1==strcmp(group, 'patients');
        LI_All_RelAct([2,3,5],:)=[];
    end
    LIs_rawdata.Broca(:,TargetColumn)=LI_All_RelAct(:,1)
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_All_RelAct(:,3)
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_All_RelAct
    %%
    %LI RelAct Substraction 0.32_0.6s BL350ms no sqrt
    ColumnSize = 1;
    load (strcat(PathResults, '\LI_RelativeActSubstraction\LI_All_RelAct_0.32_0.6s_', group, '.mat'))
    LIs_rawdata.Broca(:,TargetColumn)=LI_All_RelAct(:,1)
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_All_RelAct(:,3)
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_All_RelAct
    
    %LI RelAct Substraction 0.32_0.47s BL350ms no sqrt
    ColumnSize = 1;
    load (strcat(PathResults, '\LI_RelativeActSubstraction\LI_All_RelAct_0.32_0.47s_', group, '.mat'))
    LIs_rawdata.Broca(:,TargetColumn)=LI_All_RelAct(:,1)
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_All_RelAct(:,3)
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_All_RelAct
    
   %% 
    %LI RelAct_0.32_0.6s BL350ms sqrt
    ColumnSize = 1;
    load (strcat(PathResults, '\LI_RelativeAct\LI_All_RelAct_0.32_0.6s_', group, '.mat'))
    if 1==strcmp(group, 'patients');
        LI_All_RelAct([2,3,6],:)=[];
    end
    LIs_rawdata.Broca(:,TargetColumn)=LI_All_RelAct(:,2)
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_All_RelAct(:,4)
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_All_RelAct
    
    %LI RelAct_0.32_0.47s BL350ms sqrt
    ColumnSize = 1;
    load (strcat(PathResults, '\LI_RelativeAct\LI_All_RelAct_0.32_0.47s_', group, '.mat'))
    if 1==strcmp(group, 'patients');
        LI_All_RelAct([2,3,6],:)=[];
    end
    LIs_rawdata.Broca(:,TargetColumn)=LI_All_RelAct(:,2)
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_All_RelAct(:,4)
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_All_RelAct
    
    %%
       %LI RelAct Substraction 0.32_0.6s BL350ms sqrt
    ColumnSize = 1;
    load (strcat(PathResults, '\LI_RelativeActSubstraction\LI_All_RelAct_0.32_0.6s_', group,'.mat'))
    LIs_rawdata.Broca(:,TargetColumn)=LI_All_RelAct(:,2)
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_All_RelAct(:,4)
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_All_RelAct
    
    %LI RelAct Substraction 0.32_0.47s BL350ms sqrt
    ColumnSize = 1;
    load (strcat(PathResults, '\LI_RelativeActSubstraction\LI_All_RelAct_0.32_0.47s_', group, '.mat'))
    LIs_rawdata.Broca(:,TargetColumn)=LI_All_RelAct(:,2)
    LIs_rawdata.Wernicke(:,TargetColumn)=LI_All_RelAct(:,4)
    LIs_rawdata.Name{TargetColumn,1}=GetLineName(TargetColumn)
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_All_RelAct
    
    %%
    %LI RelAct Z-scores voxelvalues 0.32_0.6s BL100ms 
    ColumnSize = 3;
    load (strcat(PathResults, '\LI_RelativeActSubstractionBL100ms\zscores_LI_All_RelAct_0.32_0.6s_', group, '_BL100ms.mat'))
    if 1==strcmp(group, 'patients');
        LI_ALL_zscores([2 4 5 6],:)=[];
    end
    LIs_rawdata.Broca(:,TargetColumn:(TargetColumn+2))=LI_ALL_zscores(:,[1 9 17])
    LIs_rawdata.Wernicke(:,TargetColumn:(TargetColumn+2))=LI_ALL_zscores(:,[5 13 21])
    for i=TargetColumn:(TargetColumn+2)
        LIs_rawdata.Name{i,1}=GetLineName(i)
    end
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_ALL_zscores
    
    %LI RelAct Z-scores voxelvalues 0.32_0.47s BL100ms 
    ColumnSize = 3;
    load (strcat(PathResults, '\LI_RelativeActSubstractionBL100ms\zscores_LI_All_RelAct_0.32_0.47s_', group, '_BL100ms.mat'))
    if 1==strcmp(group, 'patients');
        LI_ALL_zscores([2 4 5 6],:)=[];
    end
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
    load (strcat(PathResults, '\LI_RelativeActSubstractionBL100ms\zscores_LI_All_RelAct_0.32_0.6s_', group, '_BL100ms.mat'))
    if 1==strcmp(group, 'patients');
        LI_ALL_zscores([2 4 5 6],:)=[];
    end
    LIs_rawdata.Broca(:,TargetColumn:(TargetColumn+2))=LI_ALL_zscores(:,[2 10 18])
    LIs_rawdata.Wernicke(:,TargetColumn:(TargetColumn+2))=LI_ALL_zscores(:,[6 14 22])
    for i=TargetColumn:(TargetColumn+2)
        LIs_rawdata.Name{i,1}=GetLineName(i)
    end
    TargetColumn = TargetColumn + ColumnSize;
    clear LI_ALL_zscores
    
    %LI RelAct Z-scores voxelcount 0.32_0.47s BL100ms 
    ColumnSize = 3;
    load (strcat(PathResults, '\LI_RelativeActSubstractionBL100ms\zscores_LI_All_RelAct_0.32_0.47s_', group, '_BL100ms.mat'))
    if 1==strcmp(group, 'patients');
        LI_ALL_zscores([2 4 5 6],:)=[];
    end
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


function [corClass, Sensitivity, Specificity, PPV, NPV]=classification_incl_sensitivity_atypical (values, CutOff, sizeCutOffsAll, IntCutOffs, corClass, Sensitivity, Specificity, PPV, NPV)

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
 
%   Thr=strcat('CutOff_', num2str(CutOff*100))

    corClass.Total.cases(IntCutOffs,:,1)=sum(ismember(pat_wada.classification, 'correct'))
    corClass.Total.percent(IntCutOffs, :,1)=(sum(ismember(pat_wada.classification, 'correct')))/length(pat_wada.classification)
    corClass.NoAct(:,1)=sum(isnan(pat_wada.data));

  % Percentage correct classified of Left dominant patients (according to
  % Wada): 
 [corClass] = classification_dominance (result_wada, pat_wada, 'left', corClass, CutOff, IntCutOffs) 
 [corClass] = classification_dominance (result_wada, pat_wada, 'right', corClass, CutOff, IntCutOffs) 
 [Sensitivity, Specificity, PPV, NPV]=sensitivity (corClass, pat_wada, result_wada, CutOff, sizeCutOffsAll, IntCutOffs, Sensitivity, Specificity, PPV, NPV)

end

function [corClass] = classification_dominance (result_wada, pat_wada, hem, corClass, CutOff, IntCutOffs )
%     Thr=strcat('CutOff_', num2str(CutOff*100))    
    [ind]= ismember(result_wada, hem)
    row=find(ind==1)
    for i=1:size(pat_wada.classification,2)
        Wada_hem=pat_wada.classification(row,i)
        corClass.(hem).cases(IntCutOffs,i,1)=sum(ismember(Wada_hem, 'correct'))
        corClass.(hem).percent(IntCutOffs,i,1)=(sum(ismember(Wada_hem, 'correct'))/length(Wada_hem))
    end
end



function [Sensitivity, Specificity, PPV, NPV]=sensitivity (corClass,  pat_wada, result_wada, CutOff, sizeCutOffsAll, IntCutOffs, Sensitivity, Specificity, PPV, NPV)

    % corrClass: column 1: correct classified total, column 2: correct
    % classified left wada; column 3: correct classified atypical patients
    CuffOff_left=CutOff;
    CutOff_right=CutOff*(-1);
%     Thr=strcat('CutOff_', num2str(CutOff*100))
    
    for i=1:length(corClass.left.cases)
        ntotal=14;
        nleft=8;
        natypical=6;
        
        [ind_leftwada]=ismember(result_wada, 'left')
        LIs_Wada_left=pat_wada.data(find(ind_leftwada==1), i);
        [ind_rightwada]=ismember(result_wada, 'right')
        LIs_Wada_right=pat_wada.data(find(ind_rightwada==1), i);
        
        TP = sum(LIs_Wada_right<=CutOff_right); % True positive,  all patients detected as right lateralized by MEG that have atypical Wada test
        FN = sum(LIs_Wada_right>CutOff_right);    %False negative, rechtslaterale, deren LI gr��er ist als CutOff rechts
        TN = corClass.left.cases(IntCutOffs,i,1); % True negative, all patients with left fMRI that also have left wada test
        FP = sum(LIs_Wada_left<=CutOff_right)      % False positive
        Sensitivity(IntCutOffs,i,1) = (TP/(TP+FN));
        Specificity(IntCutOffs,i,1) = (TN/(FP+TN));
        PPV(IntCutOffs, i,1) = TP/(TP+FP);     % positive predictive value
        NPV(IntCutOffs, i,1) = TN/(TN+FN);     % negative predictive value
    end
end