% this function classifies patients according to Wada test results;
% Patients with bilateral LIs are classified as True negativ
% input = List of LIs of all 24 patients (one column)
function [Classification, corClass, pat_wada, Sensitivity, Specificity, PPV, NPV]=classification_incl_sensitivity_atypical (values, CutOff)

    CutOff_left=CutOff;
    CutOff_right=CutOff*(-1);
    
    if length(values)>24
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
    for i=1:length(pat_wada.data)
        for j=1:size(pat_wada.data,2)
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
   for i=1:length(pat_wada.data)
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
  corClass.Total.cases=sum(ismember(pat_wada.classification, 'correct'))
  corClass.Total.percent=(sum(ismember(pat_wada.classification, 'correct')))/length(pat_wada.classification)
  
  % Percentage correct classified of Left dominant patients (according to
  % Wada): 
 [corClass] = classification_dominance (result_wada, pat_wada, 'left', corClass) 
 [corClass] = classification_dominance (result_wada, pat_wada, 'right', corClass) 
 [Sensitivity, Specificity, PPV, NPV]=sensitivity (corClass, pat_wada, result_wada, CutOff)
 if length(Sensitivity)==4
     Classification=[Sensitivity(2), Specificity(2), PPV(2), NPV(2),Sensitivity(4), Specificity(4), PPV(4), NPV(4)];
 else Classification = [Sensitivity(1), Specificity(1), PPV(1), NPV(1),Sensitivity(2), Specificity(2), PPV(2), NPV(2)];
 end
end

function [corClass] = classification_dominance (result_wada, pat_wada, hem, corClass )
    [ind]= ismember(result_wada, hem)
    row=find(ind==1)
    for i=1:size(pat_wada.classification,2)
        Wada_hem=pat_wada.classification(row,i)
        corClass.(hem).cases(1,i)=sum(ismember(Wada_hem, 'correct'))
        corClass.(hem).percent(1,i)=(sum(ismember(Wada_hem, 'correct'))/length(Wada_hem))
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
        TN = corClass.left.cases(1,i); % True negative, all patients with left fMRI that also have left wada test
        FP = sum(LIs_Wada_left<CuffOff_left)       %       nleft-TN;       % False positive
        Sensitivity(i,1) = (TP/(TP+FN));
        Specificity(i,1) = (TN/(FP+TN));
        PPV(i,1) = TP/(TP+FP);     % positive predictive value
        NPV(i,1) = TN/(TN+FN);     % negative predictive value
    end
end