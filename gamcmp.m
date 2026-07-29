%==========================================================

% Compare tuning curves between gam=1 and gam=0
function gamcmp
S=tdm25(0,1);
fig=4;
figure(fig);clf
plot_tc(S.pa,S,S.ts,S.cp,'-')
S.pa.gam=0;
S=tdm25('param',S.pa);
figure(fig)
plot_tc(S.pa,S,S.ts,S.cp,'--')
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
