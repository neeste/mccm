% ECochG protocol #5 - simulate ECochG across days post op
function ecochg5
fr=0.5;
ne=4;
xes=[18 20 22 24];
chl=[1 8 4 2];
fig=4;
nday=4;
am=zeros(nday,1);
gd=zeros(nday,1);
ph=zeros(nday,1);
gn=zeros(nday,1);
figure(fig);clf
for ke=1:ne
    xe=xes(ke);
    lab=sprintf('%d',xe);
    for kd=1:nday
        [am(kd),gd(kd),ph(kd),gn(kd)]=ecochg_dpo(fr,xe,chl(kd));
    end
    xx=1:nday;
    ph=unwrap(pi-ph)/(2*pi);
    plot_dpo(xx,am,gd,ph,gn,ke,lab,fig);
end
print('-dpng',fig,'','ecochg5.png');
end

%==========================================================

% simulate intra-cochlear ECochG on day of surgery
function [am,gd,ph,gn]=ecochg_dpo(fr,xe,chl)
pr.fr=0.5;      % tone frequency (kHz)
pr.lv=100;      % tone level (dB SPL)
pr.td=13;       % stimulus duration (ms)
pr.tr=1;        % ramp duration (ms)
pr.ts=0.5;      % stimulus start (ms)
pr.tp=20;       % plotted time (ms)
pr.xe=[24 23];  % electrode positions (mm)
pr.wd=0;        % write data to TXT file ???
pr.wvc=1;       % volume-conduction width (mm)
pr.hca=0.999;   % high-pass-filter coefficient
pr.hcc=750*8;   % hair-cell currrent
pr.nrc=0;       % neural-rate currrent
pr.xe=16:26;    % electrode positions (mm)
pr.diflim=[-40 40]*8; % difference-plot limits
pr.fr=fr;       % tone frequency (kHz)
pr.xe=xe;       % electrode positions (mm)
pr.chl=chl;     % CHL damping factor
gn=me_gain(pr); % middle-ear gain 
S=tdm25('ecochg',pr);
sd=S.td/2;      % stimulus delay is half its duration
am=S.am;        % amplitude
gd=S.gd-sd;     % group delay
ph=S.ph;        % phase
end

function gn=me_gain(pr)
S1=tdm25(0,1);        % run click to obtain default parameters
pa=S1.pa;             % default parameters
pa.rst=pa.rst*pr.chl; % modify stapes damping parameter (rst)
S2=tdm25('param',pa); % rerun click with modified rst
% compute average ME gain across 1-4 kHz for modified rst
f=S2.f;                    % frequency
g=abs(S2.gn);              % magnitude of ME gain
gn=mean(g((f>=1)&(f<=4))); % average ME gain
end

function plot_dpo(xx,am,gd,ph,gn,ke,lab,fig)
ne=4;
nr=4;
figure(fig)
subplot(ne,nr,ke)
plot(xx,am)
ylim([0 40])
if (ke==1), ylabel('amplitude'); end
title(lab)
subplot(ne,nr,ke+ne)
plot(xx,gd)
ylim([0 5])
if (ke==1), ylabel('delay (ms)'); end
subplot(ne,nr,ke+ne*2)
plot(xx,ph)
ylim([-0.1 1.1])
if (ke==1), ylabel('phase (cyc)'); end
subplot(ne,nr,ke+ne*3)
plot(xx,gn)
ylim([0 5])
if (ke==1), ylabel('ME gain'); end
xlabel('days post op')
drawnow
end
