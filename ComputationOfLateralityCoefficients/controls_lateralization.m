
function [NBroca,NWernicke,Class_Broca,Class_Wernicke]=Threshold(LIs_rawdata)
    group='controls'; % only works for controls
    Thr=0.25           % change!
         
    [LIs_rawdata]=createDataMatrix(group)

    [Lat_Broca]=leftlat(LIs_rawdata.Broca, Thr)
    
    % Das gleiche für Wernicke:
    [Lat_Wernicke]=leftlat(LIs_rawdata.Wernicke, Thr)
   
     % finds N with no activation (out of N=24):
    [No_Act]=findNoAct(LIs_rawdata) 
end

function [No_Act]=findNoAct(LIs_rawdata)

No_Act.Broca(:,1)=sum(isnan(LIs_rawdata.Broca));
No_Act.Wernicke(:,1)=sum(isnan(LIs_rawdata.Wernicke));
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


function [corClass]=leftlat(data, Thr)


for i=1:size(data,2)
right(i)=sum(data(:,i)<=-Thr)
left(i)=sum(data(:,i)>=Thr)
bil(i)=sum(data(:,i)<Thr & data(:,i)>-Thr)

corClass(i,1:3)=[left(i) right(i) bil(i)]
end




end