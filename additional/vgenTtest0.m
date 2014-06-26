function vgenTtest0(gaData,xlim,mirror,pThr,gaRaw)
if ischar(gaData)
    load (gaData)
    eval(['gaData=',gaData,';']);
end
if isfield(gaData,'time')
    xSamp=nearest(gaData.time,xlim)
    [~,p] = ttest(gaData.individual(:,:,xSamp));
else
    xSamp=nearest(gaData.freq,xlim)
    [~,p] = ttest(gaData.powspctrm(:,:,xSamp));
end

cfg=[];
if strcmp(gaData.label{1,1},'Fp1')
    cfg.layout='WG32.lay';
else
    cfg.layout='4D248.lay';
end
cfg.xlim=[xlim xlim];
% cfg.marker  =  'labels';
cfg.highlight          = 'marker';
cfg.highlightchannel   =  gaData.label(p<pThr);
cfg.interactive='yes';
if ~mirror
    if isfield(gaData,'individual')
        avg=mean(gaData.individual(:,:,xSamp));
        gaData.individual(:,avg<0,xSamp)=0;
    else
        avg=mean(gaData.powspctrm(:,:,xSamp));
        gaData.powspctrm(:,avg<0,xSamp)=0;
    end
    cfg.zlim=[-max(avg) max(avg)];
end
figure;
ft_topoplotER(cfg,gaData)
if exist('gaRaw','var')
    figure;
    ft_topoplotER(cfg,gaRaw)
end
    