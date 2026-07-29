% ECochG protocol #3 - various influences on ECochG
function ecochg3(var)
if (nargin<1), var=0; end
if (var==3), val=10; else, val=0; end % CHL
ecochg_comp(var,val,4);
end

function ecochg_comp(var,val,fig)
%-----------------------------
if (var==0), ecochg2(0); end % capture electrode-location plot
% Compute tuning curves
[S1,S2]=cmptc(var,val); % without and with influence
% plot tuning curves without and with influence
figure(fig);clf
plot_tc(S1.pa,S1,S1.ts,S1.cp,'-')
plot_tc(S2.pa,S2,S2.ts,S2.cp,'--')
print('-dpng',fig,'',sprintf('tc_%d.png',var));
%-----------------------------
% Compute ECochG
[~,S3]=ecochg2(0);           % initial condition
if (var)
    [~,S4]=ecochg2(var,val); % influence condition
else
    S4=S3;                   % no influence
end
gd_stim=S3.td/2;             % stimulus center is half its duration
% plot ECochG without and with influence 
nfrq=size(S3.am,1);
xe=S3.xe;
fr=S3.fr;
lg=S3.lg;
am1=S3.am;
am0=S4.am;
gd1=S3.gd-gd_stim;
gd0=S4.gd-gd_stim;
for k=1:nfrq
    ii=am1(k,:)<0.6;           % don't plot small amplitudes
    am1(k,ii)=nan;
    gd1(k,ii)=nan;
    am0(k,ii)=nan;
    gd0(k,ii)=nan;
end
xx=[min(xe) max(xe)];
ff=[  0.2      2   ];
aa=[   0          max(max(am1))];
dd=[min(min(gd1)) max(max(gd1))];
lt0='--';
lt1='-';
figure(fig);clf
subplot(2,2,1)
reset_color_index;
plot(xe',am1',lt1); hold on
reset_color_index;
plot(xe',am0',lt0); hold off
axis([xx aa])
ylabel('amplitude (\muv)')
legend(lg,'location','northwest')
subplot(2,2,3)
reset_color_index;
plot(xe',gd1',lt1); hold on
reset_color_index;
plot(xe',gd0',lt0); hold off
axis([xx dd])
ylabel('delay (ms)')
xlabel('electrode location (mm)')
subplot(2,2,2)
reset_color_index;
semilogx(fr,am1,lt1); hold on
reset_color_index;
semilogx(fr,am0,lt0); hold off
axis([ff aa])
ylabel('amplitude (\muv)')
subplot(2,2,4)
reset_color_index;
semilogx(fr,gd1,lt1); hold on
reset_color_index;
semilogx(fr,gd0,lt0); hold off
axis([ff dd])
ylabel('delay (ms)')
xlabel('stimulus frequency (kHz)')
drawnow
print('-dpng',fig,'',sprintf('echchg_%d.png',var));
end

%==========================================================

% Compare tuning curves between gam=1 and gam=0
function [S1,S2]=cmptc(var,val)
S1=tdm25(0,1);
pa=S1.pa;
if (var==1), pa.gam=val;  end
if (var==2), pa.xtap=val; end
if (var==3), pa.rst=pa.rst*val; end
S2=tdm25('param',pa);
end

% plot displacement tuning curves
function plot_tc(pa,sav,ts,cp,lt)
f=sav.f;
s=2i*pi*f;
nf=length(f);
nsv=length(pa.isv);
spl=zeros(nf,nsv);
phv=zeros(nf,nsv);
gdv=zeros(nf,nsv);
dsp_ref = 140;    % dB re 1 nm displacement
spl_ref = 0.0002; % SPL reference pressure (rms Pa)
j1=ceil(nf*(0.1/max(f))); % first freqeuncy above 0.1 kHz
for i=1:nsv
    [d1,d2,pe]=fetch_sav(i,sav);
    dh = hbmix(cp.hb(i,:), d1, d2); % HB displacement (2-col; see hbmix)
    ii=isnan(dh)|(abs(dh)<1e-16);
    dh(ii)=1e-16;
    vh = dh.*s;  % HB velocity (cm/s)
    pe = pe / spl_ref;
    spl(:,i)=20*log10(abs(pe./vh))+pa.hbt-dsp_ref;
    phv(:,i)=unwrap(angle(vh./pe))/(2*pi);
    gdv(:,i)=delay(vh./pe,f);
    for j=j1:nf
        if (gdv(j,i)<1e-9)
            phv(j:end,i)=nan;
            gdv(j:end,i)=nan;
            break;
        end
    end
end
spl(f>2,1)=nan;
spl(f>4,2)=nan;
% plot HB threshold tuning curves
subplot(2,1,1)
reset_color_index;
semilogx(f,spl,lt); hold on
axis([0.1 20 0 100])
ylabel('Pe (dB SPL)')
title('HB velocity tuning curves')
text(10,15,ts.lab)
subplot(2,1,2)
reset_color_index;
semilogx(f,phv,lt); hold on
axis([0.1 20 -9 1])
ylabel('phase (cyc)')
text(10,-7,pa.parlab)
xlabel('frequency (kHz)')
drawnow
end % return

function [d1,d2,pe,ps,vr,vs,ve,nr]=fetch_sav(i,sav)
d1 = fft(sav.d1(:,i));
d2 = fft(sav.d2(:,i));
nr = fft(sav.wnr);
vs = fft(sav.vst);
ps = fft(sav.pst);
pe = fft(sav.ped);
ve = fft(sav.ved);
vr = fft(sav.vep);
nf=length(sav.f);
ii=1:nf;
d1 = d1(ii);
d2 = d2(ii);
nr = nr(ii);
vs = vs(ii);
ps = ps(ii);
pe = pe(ii);
ve = ve(ii);
vr = vr(ii);
end % return;

% group delay
function gd=delay(R,f)
[n,m] = size(R);
ph = unwrap(angle(R))/(2*pi);
gd = zeros([n,m]);
for k=1:m
   gd(:,k) = -cdif(ph(:,k))./cdif(f(:));
end
end % return

% centered difference
function dx=cdif(x)
n=length(x);
dx=zeros(size(x));
dx(1)=x(2)-x(1);
dx(2:(n-1))=(x(3:n)-x(1:(n-2)))/2;
dx(n)=x(n)-x(n-1);
end % return

%----------------------------------------------------------

function reset_color_index
ax=gca();
set(ax,'ColorOrderIndex',1);
end % return
