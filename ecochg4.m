% ECochG protocol #4 - simulate ECochG on day of surgery
function ecochg4
fig=4;
%-----------------------------
% compute WNR, spectrogram, VC, & single-frequency ECochG waveforms
ecochg_dos;
%-----------------------------
% compute multiple-frequency ECochG waveforms
S=ecochg_dos;
% prepare to plot 
nfrq=size(S.am,1);
fr=S.fr;
sd=S.td/2; % stimulus delay is half its duration
am1=S.am;
gd1=S.gd-sd;
for k=1:nfrq
    ii=am1(k,:)<0.1; % don't plot small amplitudes
    am1(k,ii)=nan;
    gd1(k,ii)=nan;
end
lt1='-'; % line type
% specify plotting limits
ff=[  0.2        2     ];
aa=[  0   max(max(am1))];
dd=[  0   max(max(gd1))]*1.05;
% plot ECochG amplitude & delay
figure(fig);clf
subplot(2,1,1)
reset_color_index;
semilogx(fr,am1,lt1)
axis([ff aa])
ylabel('amplitude (\muv)')
title('ECochG difference')
subplot(2,1,2)
reset_color_index;
semilogx(fr,gd1,lt1)
axis([ff dd])
ylabel('delay (ms)')
xlabel('stimulus frequency (kHz)')
drawnow
print('-dpng',fig,'','ecochg4c.png');
end

%==========================================================

% simulate intra-cochlear ECochG on day of surgery
function S=ecochg_dos
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
pr.dsp1=4;      % upper panel [1=Vr 2=Pe 3=Ps 4=WNR]
pr.dsp2=4;      % lower panel [1=HB 2=NR 3=LD 4=HC]
if (nargout)
    pr.fr=[0.25*(1:4) 2]; % set of tone frequencies (kHz)
end
S=tdm25('ecochg',pr);
if (isscalar(pr.fr))
    print('-dpng',1,'','ecochg4a.png');
    print('-dpng',3,'','ecochg4b.png');
end
end

%==========================================================

function reset_color_index
ax=gca();
set(ax,'ColorOrderIndex',1);
end % return
