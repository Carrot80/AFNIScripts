% calculates Sensitivity and specificity for detection of typical language :

function [Sensitivity, Specificity, PPV, NPV]=sensitivity_typical_Language (corrClass)

    % corrClass: column 1: correct classified total, column 2: correct
    % classified left wada; column 3: correct classified atypical patients

    for i=1:length(corrClass)
        ntotal=14;
        nleft=8;
        natypical=6;

        TP = corrClass(i,2); % True positive
        FN = nleft-TP;   % false negative
        TN = corrClass(i,3); % True negative
        FP = natypical-TN;       % False positive
        Sensitivity(i,1) = (TP/(TP+FN));
        Specificity(i,1) = (TN/(FP+TN)) ;
        PPV(i,1) = TP/(TP+FP);     % positive predictive value
        NPV(i,1) = TN/(TN+FN);     % negative predictive value
    end
end