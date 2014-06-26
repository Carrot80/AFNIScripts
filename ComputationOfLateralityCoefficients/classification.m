function [corClass, pat_wada]=classification (values)

    % corClass= percentage of correct classified cases
    %Total
    %left=patients with left wada test result
    
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
            if pat_wada.data(i,j)<=-.20
                pat_wada.nom{i,j}='right';
            elseif pat_wada.data(i,j)>=.20
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