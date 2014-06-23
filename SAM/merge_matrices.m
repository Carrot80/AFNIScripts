
function for_all ()
    
    mergeTasks('Patients')

end

function mergeTasks (group)

%% Broca

Table_Fluency=load(strcat('D:\Arbeit\LinuxExchange\data\Cluster_and_Locations_fMRI_t-values_', group,'_Fluency_ROIs.mat'));
Table_VG=load(strcat('D:\Arbeit\LinuxExchange\data\Cluster_and_Locations_fMRI_t-values_', group,'_Verbgeneration_ROIs.mat'));

Table_Both = zeros(size(Table_Fluency.Table.Broca.Coord,1) ,16);

j=1

for i=1:size(Table_Fluency.Table.Broca.Coord,1)
   
    Table_Both(j,:)=Table_Fluency.Table.Broca.Coord(i, :);
    j=j+2;
end


 j=2

for i=1:size(Table_Fluency.Table.Broca.Coord,1)
   
    Table_Both(j,:)=Table_VG.Table.Broca.Coord(i, :);
    j=j+2;
end

xlswrite(strcat('Table_clust_merged_Broca_', group), Table_Both)

%% Location 


j=1

for i=1:size(Table_Fluency.Table.Broca.Coord,1)
   
    Location{j,:}=Table_Fluency.Table.Broca.Location{i, 1} ;
    j=j+2;
end


 j=2

for i=1:size(Table_Fluency.Table.Broca.Coord,1)
   
    Location{j,:}=Table_VG.Table.Broca.Location{i, 1} ;
    j=j+2;
end




%% Wernicke

Table_Both = zeros(size(Table_Fluency.Table.Broca.Coord,1) ,16);

j=1

for i=1:size(Table_Fluency.Table.Broca.Coord,1)
   
    Table_Both(j,:)=Table_Fluency.Table.Wernicke.Coord(i, :);
    j=j+2;
end


 j=2

for i=1:size(Table_Fluency.Table.Broca.Coord,1)
   
    Table_Both(j,:)=Table_VG.Table.Wernicke.Coord(i, :);
    j=j+2;
end

xlswrite(strcat('Table_clust_merged_Wernicke_', group), Table_Both)

%% Location 


j=1

for i=1:size(Table_Fluency.Table.Broca.Coord,1)
   
    Location{j,:}=Table_Fluency.Table.Wernicke.Location{i, 1} ;
    j=j+2;
end


 j=2

for i=1:size(Table_Fluency.Table.Broca.Coord,1)
   
    Location{j,:}=Table_VG.Table.Wernicke.Location{i, 1} ;
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