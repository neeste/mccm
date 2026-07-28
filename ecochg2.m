% ECochG protocol #2 - simulate intra-cochlear ECochG tone-burst response
function [pr,S]=ecochg2(var,val)
pr.fr=0.5;      % tone frequency (kHz)
pr.lv=100;      % tone level (dB SPL)
pr.td=50;       % stimulus duration (ms)
pr.tr=2;        % ramp duration (ms)
pr.ts=0.5;      % stimulus start (ms)
pr.tp=60;       % plotted time (ms)
pr.xe=[24 23];  % electrode positions (mm)
pr.wvc=1;       % volume-conduction width (mm)
pr.hca=0.999;   % high-pass-filter coefficient
pr.hcc=750;     % hair-cell currrent
pr.nrc=0;       % neural-rate currrent
pr.sumlim=[-40 40]; % summed range
pr.diflim=[-40 40]; % difference range
if (nargin)
    if (var==1), pr.gam=val;  end % OHC motility
    if (var==2), pr.xtap=val; end % apical taper
    if (var==3), pr.chl=val;  end % conductive loss
    pr.xe=16:26;                  % electrode positions (mm)
    if (nargout)                  % summarize across frequency and position
        pr.fr=[0.25*(1:4) 2]; % tone frequency (kHz)
    end
else
    var=0;
end
S=tdm25('ecochg',pr);
if (nargout==0)
    if (var==0)
        if (length(pr.xe)==2)
            print('-dpng',1,'','ecochg2.png');
            print('-dpng',2,'','ecochg2a.png');
            print('-dpng',3,'','ecochg2b.png');
        else
            print('-dpng',3,'','ecochg2c.png');
        end
    end
    clear pr
end
end
