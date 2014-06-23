
% Read trigger:
fileName=('hb_lf_c,rfhp0.1Hz') ;

trig=readTrig_BIU(fileName);
figure;plot(trig)
unique(trig)
trig4224=bitand(uint16(trig),4224);
hold on
plot(trig4224,'r')

trig=double(bitand(uint16(trig),4224));
trigSh=zeros(size(trig));
trigSh(2:end)=trig(1:end-1);
trigOnset=find((trig-trigSh)>0);
offset=round(0.5*1017.25);
trl=trigOnset'-offset;
trl(:,2)=trl+1017;
trl(:,3)=-offset;

save trl trl

  
hdr                     = ft_read_header(fileName) ;
     cfg = [] ;
     cfg.dataset             = fileName ;
     cfg.channel             = 'MEG' ;

    cfg.trl     = trl;
    cfg.dataset = fileName;
    cfg.demean  ='yes';
    cfg.baselinewindow = [-0.3 0];
    cfg.channel     ='MEG';
    cfg.hpfilter    = 'yes';
    cfg.hpfreq      = 60;
    cfg.padding     = 1;
    dataHp60=ft_preprocessing(cfg);
    
    cfg=[];
    cfg.method='abs';
    cfg.criterion='sd';
    cfg.critval=3;
    [good,bad]=badTrials(cfg,dataHp60,1)
    save BadTrials bad
    save GoodTrials good