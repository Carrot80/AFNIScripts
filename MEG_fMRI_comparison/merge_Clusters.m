
function for_all ()
    
    mergeTasks('Patients')

end

function mergeTasks (group)

%% Broca

load(strcat('D:\Arbeit\LinuxExchange\Results\Comparison_fMRI_MEG\Clusters_and_Locations_Controls_fMRI_MEG_comparison.mat'));

Table_Both = zeros(size(Table.MEG.Broca.Coord,1) ,16);

j=1

for i=1:size(Table.MEG.Broca.Coord,1)
   
    Table_Both(j,:)=Table.MEG.Broca.Coord(i, :);
    j=j+2;
end


 j=2

for i=1:size(Table.MEG.Broca.Coord,1)
   
    Table_Both(j,:)=Table.fMRI.Broca.Coord(i, :);
    j=j+2;
end

xlswrite('Table_clust_merged_Broca_controls', Table_Both)

%% Location 


j=1

for i=1:size(Table.MEG.Broca.Coord,1)
   
    Location{j,:}=Table.MEG.Broca.Location{i, 1} ;
    j=j+2;
end


 j=2

for i=1:size(Table.MEG.Broca.Coord,1)
   
    Location{j,:}=Table.fMRI.Broca.Location{i, 1} ;
    j=j+2;
end




%% Wernicke

Table_Both = zeros(size(Table.MEG.Wernicke.Coord,1) ,16);

j=1

for i=1:size(Table.MEG.Wernicke.Coord,1)
   
    Table_Both(j,:)=Table.MEG.Wernicke.Coord(i, :);
    j=j+2;
end


 j=2

for i=1:size(Table.MEG.Wernicke.Coord,1)
   
    Table_Both(j,:)=Table.fMRI.Wernicke.Coord(i, :);
    j=j+2;
end

xlswrite('Table_clust_merged_Wernicke_controls', Table_Both)

%% Location 


j=1

for i=1:size(Table.MEG.Wernicke.Coord,1)
   
    Location{j,:}=Table.MEG.Wernicke.Location{i, 1} ;
    j=j+2;
end


 j=2

for i=1:size(Table.MEG.Wernicke.Coord,1)
   
    Location{j,:}=Table.fMRI.Wernicke.Location{i, 1} ;
    j=j+2;
end


for k=1:length(Location)
    if isempty(strfind(Location{k,1}, 'Left'))
        Hem{k,1}='R';
    else
       Hem{k,1}='L';
    end
    
end


%%
















end