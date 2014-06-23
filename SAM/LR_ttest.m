
function for_all ()

% Script written by Yuval:

ControlsFolder = '/home/kh/ShareWindows/data/controls/controls_SAM';

DIR = dir (ControlsFolder)
isub = [DIR(:).isdir]; %  returns logical vector
nameFolds = {DIR(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];


TimeIntervall (nameFolds, ControlsFolder, 0.32, 0.47 )
% TimeIntervall (nameFolds, ControlsFolder, 0.25, 0.55 )
% TimeIntervall (nameFolds, ControlsFolder, 0.25, 0.65 )


end


function TimeIntervall (nameFolds, ControlsFolder, TimeBeg, TimeEnd)

VlrAll= [];
Vall=[];


for i= 1:size(nameFolds)
    SubjectPath = strcat(ControlsFolder, filesep, nameFolds{i,1});
    SubjectName = nameFolds{i};

    
    [V, Info, VlrAll, Vall, TimeBeg, TimeEnd] = get_V (SubjectPath, SubjectName, TimeBeg, TimeEnd, VlrAll, Vall, i);
    
    
end

     UtestLR (SubjectPath, SubjectName, Info, VlrAll, Vall, TimeBeg, TimeEnd, V);
     
     
     

end

%%




function [V, Info, VlrAll, Vall, TimeBeg, TimeEnd] = get_V (SubjectPath, SubjectName, TimeBeg, TimeEnd, VlrAll, Vall, i)

    Path = strcat(SubjectPath, filesep, 'SAM');
    cd (Path)

    FileName = strcat('br_z_transf_brain01ERF_', num2str(TimeBeg), '-', num2str(TimeEnd), 's', '_', SubjectName, '+tlrc');

    [V, Info] = BrikLoad (FileName);
    Vlr=flipdim(V,1);
    if i==1
        Vall=V;
        VlrAll=Vlr;
    else
        Vall(:,:,:,i)=V;
        VlrAll(:,:,:,i)=Vlr;
    end
end

%%

function UtestLR (SubjectPath, SubjectName, Info, VlrAll, Vall, TimeBeg, TimeEnd, V);
 
    U=zeros(size(V));
    for i=1:36
        for j=1:44
            for k=1:36
                if Vall(i,j,k)>0
                    p = ranksum(squeeze(Vall(i,j,k,:)),squeeze(VlrAll(i,j,k,:)));
                    u=1-p;
                    dif=mean(squeeze(Vall(i,j,k,:)))-mean(squeeze(VlrAll(i,j,k,:)));
                    if dif>0
                        R=1;% R is 1 when current is larger than other side
                    else
                        R=-1;
                    end
                    U(i,j,k)=R*u;
                end
            end
        end
        i
    end

    PathAllControls = '/home/kh/ShareWindows/data/controls/AVGcontrols';
    
    cd (PathAllControls)
    
    OptTSOut.Scale = 1;
    OptTSOut.Prefix = strcat('Utest', '_', 'LR', '_', num2str(TimeBeg), '_', num2str(TimeEnd), 's');
    OptTSOut.verbose = 1;
    OptTSOut.View = '+tlrc'
    %Vsymm=double(Vlr+V>0);
    WriteBrik (U, Info, OptTSOut);
    
end


