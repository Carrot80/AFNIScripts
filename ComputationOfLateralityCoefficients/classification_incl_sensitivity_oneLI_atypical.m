%Classifies patients according to :
% 1 ROI bilateral + 1 ROI right => right lateralization
% 1 ROI bilateral + 1 ROI left => left lateralization
% 2 ROIs bialteral = bilateral
% Rest = dissociated language
% Input: LIs

function forALL (LI_Broca, LI_Wernicke)

% load 'LI_All_RelAct_0.32_0.6s_Patients.mat'
% LI_All_RelAct([2 3 6],:)=[];
% LI_Broca=LI_All_RelAct(:,2)
% LI_Wernicke=LI_All_RelAct(:,4)
[Classification_2, corClass_2, ~]=classification_incl_sensitivity_oneLI (LI_Broca, LI_Wernicke, .2)
[Classification_15, corClass_15, ~]=classification_incl_sensitivity_oneLI (LI_Broca, LI_Wernicke, .15)
[Classification_1, corClass_1, ~]=classification_incl_sensitivity_oneLI (LI_Broca, LI_Wernicke, .1)
end

function [Classification, corClass, Sensitivity, Specificity, PPV, NPV]=classification_incl_sensitivity_oneLI(LI_Broca, LI_Wernicke, CuttOff)

    % corClass= percentage of correct classified cases
    %Total
    %left=patients with left wada test result
    
    CuttOff_left=CuttOff;
    CuttOff_right=CuttOff*-1
    if length(LI_Broca)>24
        disp ('Warning: length of input data is larger than 24 - aborting')
        return
    end
   
    if length(LI_Broca)<24
        disp ('Warning: length of input data is smaller than 24 - aborting')
        return
    end
    for i=1:length(LI_Broca)
        if LI_Broca(i)<CuttOff_left && LI_Broca(i)>CuttOff_right && LI_Wernicke(i)<=CuttOff_right
            Laterality{i,1}='atypical'
        elseif LI_Broca(i)<=CuttOff_right && LI_Wernicke(i)<CuttOff_left && LI_Wernicke(i)>CuttOff_right
            Laterality{i,1}='atypical'
        elseif LI_Broca(i)<=CuttOff_right && LI_Wernicke(i)<=CuttOff_right
            Laterality{i,1}='atypical'
        elseif LI_Broca(i)<CuttOff_left && LI_Broca(i)>CuttOff_right && LI_Wernicke(i)>=CuttOff_left
            Laterality{i,1}='dissociate'
        elseif LI_Broca(i)>=CuttOff_left && LI_Wernicke(i)<CuttOff_left && LI_Wernicke(i)>CuttOff_right
            Laterality{i,1}='left'
        elseif LI_Broca(i)>=CuttOff_left && LI_Wernicke(i)>=CuttOff_left
            Laterality{i,1}='left'
        elseif LI_Broca(i)<CuttOff_left && LI_Wernicke(i)<CuttOff_left && LI_Wernicke(i)>CuttOff_right
            Laterality{i,1}='atypical'        
        elseif LI_Broca(i)<=CuttOff_right && LI_Wernicke(i)<=CuttOff_right
        else Laterality{i,1}='dissociate'
        end
    end
    % classification wada test Patient 1-24:

    Wada= {'atypical' 'left' 'left' 'left' 'kein Wada' 'atypical' 'kein Wada' 'atypical' 'kein Wada' 'atypical' 'kein Wada' ...
        'left' 'left' 'atypical' 'kein Wada' 'kein Wada' 'left' 'kein Wada' 'kein Wada' 'left' 'atypical' 'kein Wada' 'left' 'kein Wada'}
    wada_no=[5 7 9 11 15 16 18 19 22 24];
    wada_yes=[1:4, 6, 8, 10, 12:14, 17, 20, 21, 23];
    result_wada=Wada(wada_yes)';
    fMRI_result=Laterality(wada_yes');
    fMRI_noWada=Laterality(wada_no);
    
   % classification according to correctness:
   for i=1:length(fMRI_result)
        if 1==strcmp(fMRI_result{i}, result_wada{i})
               classification{i}='correct';
           else
              classification{i}='incorrect';
        end
   end
        
  % Percentage correct classified 
  corClass   = [];
  corClass.Total.cases=sum(ismember(classification, 'correct'));
  corClass.Total.percent=(sum(ismember(classification, 'correct')))/length(classification);
  
  % Percentage correct classified of Left dominant patients (according to
  % Wada): 
 [corClass] = classification_dominance (classification, result_wada, fMRI_result, 'left', corClass) 
 [corClass] = classification_dominance (classification, result_wada, fMRI_result, 'atypical', corClass) 
 [Sensitivity, Specificity, PPV, NPV]=sensitivity (corClass, fMRI_result, result_wada, CuttOff)
 Classification=[Sensitivity, Specificity, PPV, NPV];
end

function [corClass] = classification_dominance (classification, result_wada, fMRI_result, hem, corClass )
    [ind]= ismember(result_wada, hem)
    row=find(ind==1)
    Wada_hem=classification(row)
    corClass.(hem).cases=sum(ismember(Wada_hem, 'correct'))
    corClass.(hem).percent=(sum(ismember(Wada_hem, 'correct'))/length(Wada_hem))
end

function [Sensitivity, Specificity, PPV, NPV]=sensitivity (corClass,  fMRI_result, result_wada, CuttOff)

    % corrClass: column 1: correct classified total, column 2: correct
    % classified left wada; column 3: correct classified atypical patients

    CuttOff_left=CuttOff;
    CuttOff_right=CuttOff*-1;
    for i=1:length(corClass.left.cases)
        ntotal=14;
        nleft=8;
        natypical=6;
        
        [ind_leftwada]=ismember(result_wada, 'left')
        Lateratliy_fMRI_left=fMRI_result(find(ind_leftwada==1));
        [ind_atypicalwada]=ismember(result_wada, 'atypical')
        Lateratliy_fMRI_atypical=fMRI_result(find(ind_atypicalwada==1));
        
        TP = sum(ismember(Lateratliy_fMRI_atypical, 'atypical')); % True positive
        FN = natypical-sum(ismember(Lateratliy_fMRI_atypical, 'atypical'));   % false negative
        TN = sum(ismember(Lateratliy_fMRI_left, 'left')); % True negative
        FP = nleft-sum(ismember(Lateratliy_fMRI_left, 'left'));   %       nleft-TN;       % False positive
        Sensitivity(i,1) = (TP/(TP+FN));
        Specificity(i,1) = (TN/(FP+TN));
        PPV(i,1) = TP/(TP+FP);     % positive predictive value
        NPV(i,1) = TN/(TN+FN);     % negative predictive value
    end
end