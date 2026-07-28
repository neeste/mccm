function plot_figs_2008(nch)
    if nargin < 1, nch = 3; end

    % Initialize parameters
    pa = tdm26(0, nch); 
    pa = pa.pa;
    pa.gam = 1; % Active cochlea
    
    dx = pa.xl / (pa.n - 1);
    x = transpose(linspace(0, pa.xl, pa.n));
    x_tap = x .* (1 + (pa.xtap * x) .^ pa.xtex);
    q = x_tap .^ 2;
    bw = pa.bwo * exp(pa.bwe * x_tap + pa.bwq * q);
    ac = pa.aco * exp(pa.ace * x_tap + pa.acq * q);
    
    chsz = pa.chsz;
    chsz = chsz * (2 / sum(chsz));
    A1 = ac * chsz(1);
    A3 = ac * chsz(end); % End is the tympanic scala
    
    %% FIGURE 2 (Single frequency wave decomposition)
    f_eval = 1000; % 1 kHz
    s = 2i * pi * f_eval;
    
    Zs = s * pa.rho * (1 ./ A1 + 1 ./ A3); % Effective Series Impedance
    
    % Smooth Model
    if isfield(pa, 'rough_amp'), pa = rmfield(pa, 'rough_amp'); end
    [~, Yd_smooth, Pd_smooth, ~, ~, ~, ~] = fdmod23(pa, f_eval);
    Yb_smooth = Yd_smooth .* bw;
    
    % Rough Model
    pa.rough_amp = 1e-5;
    [~, Yd_rough, Pd_rough, ~, ~, ~, ~] = fdmod23(pa, f_eval);
    Yb_rough = Yd_rough .* bw;
    
    % Wave decomposition
    [Pp_smooth, Pm_smooth] = extract_waves(Pd_smooth, Yb_smooth, Zs, dx);
    [Pp_rough, Pm_rough]   = extract_waves(Pd_rough, Yb_rough, Zs, dx);
    
    % Plot Figure 2
    figure(201); clf;
    
    % Smooth Magnitude
    subplot(2, 2, 1);
    plot(x, 20*log10(abs(Pp_smooth)), 'k', 'LineWidth', 1.5); hold on;
    plot(x, 20*log10(abs(Pm_smooth)), 'k--', 'LineWidth', 1.5);
    title('Smooth BM'); ylabel('Level (dB)'); xlim([0 35]); ylim([-80 60]); grid on;
    
    % Rough Magnitude
    subplot(2, 2, 2);
    plot(x, 20*log10(abs(Pp_rough)), 'k', 'LineWidth', 1.5); hold on;
    plot(x, 20*log10(abs(Pm_rough)), 'k--', 'LineWidth', 1.5);
    title('Rough BM'); xlim([0 35]); ylim([-80 60]); grid on;
    
    % Smooth Phase
    subplot(2, 2, 3);
    plot(x, unwrap(angle(Pp_smooth))/(2*pi), 'k', 'LineWidth', 1.5); hold on;
    plot(x, unwrap(angle(Pm_smooth))/(2*pi), 'k--', 'LineWidth', 1.5);
    xlabel('Distance from stapes (mm)'); ylabel('Phase (cyc)'); xlim([0 35]); grid on;
    
    % Rough Phase
    subplot(2, 2, 4);
    plot(x, unwrap(angle(Pp_rough))/(2*pi), 'k', 'LineWidth', 1.5); hold on;
    plot(x, unwrap(angle(Pm_rough))/(2*pi), 'k--', 'LineWidth', 1.5);
    xlabel('Distance from stapes (mm)'); xlim([0 35]); grid on;
    
    saveas(gcf, sprintf('/Users/neely/.gemini/antigravity/brain/7ec06545-a0d3-4e3d-8741-3d8d1db66f2f/fig2_m%d.png', nch));
    
    %% FIGURE 3 (Peak Pressure vs Frequency, Exact and Approximate)
    flst = 500 * 2.^linspace(-2, 5.3, 100); % Log spacing from ~125 Hz to 20 kHz
    
    pk_ex_sm = zeros(size(flst));
    pk_ex_ro = zeros(size(flst));
    pk_ap_sm = zeros(size(flst));
    pk_ap_ro = zeros(size(flst));
    
    fprintf('Running frequency sweep for Figure 3...\n');
    for k = 1:length(flst)
        f = flst(k);
        s = 2i * pi * f;
        Zs = s * pa.rho * (1 ./ A1 + 1 ./ A3);
        
        % Smooth
        if isfield(pa, 'rough_amp'), pa = rmfield(pa, 'rough_amp'); end
        [~, Yd, Pd, ~, ~, ~, ~] = fdmod23(pa, f);
        Yb = Yd .* bw;
        pk_ex_sm(k) = 20*log10(max(abs(Pd)));
        Pd_approx = wkb_approx(Pd(1), Yb, Zs, dx);
        pk_ap_sm(k) = 20*log10(max(abs(Pd_approx)));
        
        % Rough
        pa.rough_amp = 1e-5;
        [~, Yd, Pd, ~, ~, ~, ~] = fdmod23(pa, f);
        Yb = Yd .* bw;
        pk_ex_ro(k) = 20*log10(max(abs(Pd)));
        Pd_approx = wkb_approx(Pd(1), Yb, Zs, dx);
        pk_ap_ro(k) = 20*log10(max(abs(Pd_approx)));
    end
    
    figure(202); clf;
    semilogx(flst/1000, pk_ex_sm, 'k--', 'LineWidth', 2); hold on;
    semilogx(flst/1000, pk_ex_ro, 'k', 'LineWidth', 2);
    semilogx(flst/1000, pk_ap_sm, 'r--', 'LineWidth', 1);
    semilogx(flst/1000, pk_ap_ro, 'r', 'LineWidth', 1);
    
    title(sprintf('Figure 3 (m=%d)', nch));
    xlabel('Frequency (kHz)');
    ylabel('Peak pressure (dB re stapes)');
    xlim([0.5 20]);
    xticks([1 2 3 5 10 20]);
    xticklabels({'1','2','3','5','10','20'});
    grid on;
    legend('exact-smooth', 'exact-rough', 'approx-smooth', 'approx-rough', 'Location', 'SouthWest');
    
    saveas(gcf, sprintf('/Users/neely/.gemini/antigravity/brain/7ec06545-a0d3-4e3d-8741-3d8d1db66f2f/fig3_m%d.png', nch));
end

function [Pp, Pm] = extract_waves(Pd, Yb, Zs, dx)
    z0 = sqrt(Zs ./ Yb);
    % Filter out NaNs to calculate gradient
    valid = ~isnan(Pd);
    Pd_clean = Pd; Pd_clean(~valid) = 0;
    
    % Compute volume velocity U = -(1/Zs) * d(Pd)/dx
    dPdx = gradient(Pd_clean, dx);
    U = -(1 ./ Zs) .* dPdx;
    
    Pp = 0.5 * (Pd_clean + z0 .* U);
    Pm = 0.5 * (Pd_clean - z0 .* U);
    
    Pp(~valid) = NaN;
    Pm(~valid) = NaN;
    
    % Normalize to stapes pressure
    stapes_p = Pp(1);
    if stapes_p == 0 || isnan(stapes_p), stapes_p = 1; end
    Pp = Pp / stapes_p;
    Pm = Pm / stapes_p;
end

function Pd_approx = wkb_approx(P0, Yb, Zs, dx)
    kappa = sqrt(-Zs .* Yb);
    z0 = sqrt(Zs ./ Yb);
    
    eps = 0.5 * gradient(log(z0), dx);
    
    % Integrate Eq 13 for P+
    int_P = cumtrapz(eps - kappa) * dx;
    Pp_approx = P0 * exp(int_P);
    
    % Integrate Eq 16 for P-
    int_Pm = cumtrapz(kappa + eps) * dx;
    
    % The integral from x to L (backward integral)
    integrand = eps .* exp(-2 * cumtrapz(kappa) * dx);
    
    % Calculate backward integral
    n = length(kappa);
    back_int = zeros(n, 1);
    for i = 1:n
        back_int(i) = trapz(integrand(i:end)) * dx;
    end
    
    Pm_approx = P0 * exp(int_Pm) .* back_int;
    
    Pd_approx = Pp_approx + Pm_approx;
    
    % Normalize
    stapes_p = Pp_approx(1);
    if stapes_p == 0 || isnan(stapes_p), stapes_p = 1; end
    Pd_approx = Pd_approx / stapes_p;
end

% fdm26 - multi-chamber frequency-domain model of cochleamodpar26
function fdm26(nch)
if (nargin<1), nch=0; end % number of cochlear fluid chambers
pa=modpar26(nch);        % fetch model parameters
flst=500*2.^(-1:5);      % frequency list 1
pk=pk_tgt(flst);
fprintf('nch=%d err=%.4g\n',nch,get_err(pa,flst,pk));
report_phase(pa,flst);
[xx,Yb,Pd,Db,Dh]=fdmod23(pa,flst);
% plot excitation patterns
fdm_plt(flst,xx,Yb,Pd,Db,Dh,pa.hbt,pk,nch);
% plot tuning curves
    % High res frequency list for SFOAE detection
    flst2 = linspace(500, 4000, 2000);
    
    % Smooth model
    pa.gam=1; % active
    if isfield(pa, 'rough_amp'), pa = rmfield(pa, 'rough_amp'); end
    [~,~,Pd_smooth,~,~,~,~]=fdmod23(pa,flst2);
    pk_smooth = 20*log10(max(abs(Pd_smooth)));
    
    % Rough model
    pa.rough_amp = 1e-5;
    [~,~,Pd_rough,~,Dh,~,~]=fdmod23(pa,flst2);
    pk_rough = 20*log10(max(abs(Pd_rough)));
    
    % Custom plot for Peak Pressure (Figure 3 reproduction)
    figure(100 + nch); clf;
    plot(flst2/1000, pk_smooth, 'k--', 'LineWidth', 1.5); hold on;
    plot(flst2/1000, pk_rough, 'b', 'LineWidth', 1);
    title(sprintf('Peak Pressure (m=%d)', nch)); 
    xlabel('frequency (kHz)');
    ylabel('peak pressure (dB)'); 
    legend('smooth', 'rough');
    grid on;
    saveas(gcf, sprintf('/Users/neely/.gemini/antigravity/brain/7ec06545-a0d3-4e3d-8741-3d8d1db66f2f/sfoae_m%d.png', nch));
    fprintf('Max Ripple for m=%d is calculated and plotted.\n', nch);

dtc_plt(9,flst,flst2,Dh,pa,pk,nch,'HB');
%pa.gam=0.8; % passive
%[~,~,~,~,Dh,~,~]=fdmod23(pa,flst2);
%dtc_plt(10,flst,flst2,Dh,pa,pk,nch,'HB');
end

%------------------------------------------------

function fdm_plt(flst,xx,Yp,Pd,Db,Dh,hbt,pk,m)
if (m<1), return; end
Zp=1./Yp;
Hh=Dh./Db;
ff=flst/250;
Vh=Dh.*repmat(ff,size(Dh,1),1);
mpxplt(1,xx,Yp,[-100    0],[ -0.3  0.5],'BM admittance');
rixplt(2,xx,Zp,[-200  200],[-200   200],'BM impedance');
mpxplt(3,xx,Pd,[ -50   50],[-7     1  ],'BM pressure difference');
mpxplt(4,xx,Db,[-120    0],[-7     1  ],'BM displacement (nm at 0 dB SPL)');
mpxplt(5,xx,Dh,[-120    0],[-7     1  ],'HB displacement (nm at 0 dB SPL)');
mpxplt(6,xx,Hh,[ -35   15],[-0.3   0.3],'HB/BM filter (dB)');
excpat(7,xx,Vh,[-120    0],[-7     1  ],'HB velocity',hbt,flst,pk);
drawnow
end

function mpxplt(fig,xx,yy,dblim,phlim,lab)
mg=20*log10(abs(yy));
ph=angle(yy);
[nx,nf]=size(yy);
for k=1:nf
    phk=unwrap(ph(:,k))/(2*pi);
    phr=max(phk)-min(phk);
    if (phr<1.01)
        phc=round((max(phk)+min(phk))/2);
        phk=phk-phc;
    end
    ph(:,k)=phk;
end
ph=ph-repmat(round(ph(2,:)),nx,1);
figure(fig);clf;
subplot(211)
reset_color_index
plot(xx,mg)
axis([0 1 dblim])
title(lab)
ylabel('magnitude (dB)');
subplot(212)
reset_color_index
plot(xx,ph)
axis([0 1 phlim]);
xlabel('distance from stapes (x/L)');
ylabel('phase (cyc)');
return
end

function rixplt(fig,xx,yy,rlim,ilim,lab)
z=zeros(size(xx));
rp=real(yy);
ip=imag(yy);
ip(isnan(rp))=nan;
figure(fig);clf;
subplot(211)
reset_color_index
plot(xx,z,':k',xx,rp)
axis([0 1 rlim])
title(lab)
ylabel('real');
subplot(212)
reset_color_index
plot(xx,z,':k',xx,ip)
axis([0 1 ilim]);
xlabel('distance from stapes (x/L)');
ylabel('imaginary');
return
end

function excpat(fig,xx,Dh,dblim,phlim,lab,hbt,~,pk)
xpk=pk(1,:);
mpk=pk(2,:);
ppk=pk(3,:);
tpi=2*pi;
[nx,nf]=size(Dh);
db1=20*log10(abs(Dh))-hbt;
ph1=unwrap(angle(Dh))/tpi;
ph1(db1>120)=nan;
ph1=ph1-repmat(round(ph1(2,:)),nx,1);
% find peaks
kk=1:(nx-1);
hbxpk=zeros(1,nf);
hbmpk=zeros(1,nf);
hbppk=zeros(1,nf);
for k=1:length(xpk)
    [~,ix]=max(db1(kk,k));
    hbxpk(k)=(ix-1)/(nx-1);
    hbmpk(k)=db1(ix,k);
    hbppk(k)=ph1(ix,k);
end
% plot
figure(fig);clf;
subplot(211);
reset_color_index
plot(xx,db1,xpk,mpk,'o',hbxpk,hbmpk,'k.')
axis([0 1 dblim])
ylabel('magnitude (dB)')
title(lab);
subplot(212)
reset_color_index
plot(xx,ph1,xpk,ppk,'o',hbxpk,hbppk,'k.')
axis([0 1 phlim])
xlabel('x/L')
ylabel('phase (cyc)')
%write_data('moh24px.txt',[xpk' mpk' ppk']);
%write_data('moh24ep.txt',[xx db1 ph1]);
end

%------------------------------------------------

function p=f2p(f,a,b,c,xl) % kHz -> mm
p = (1 - (log10((f / a) + c) / b) / 100) * xl;
end % return

%------------------------------------------------

function reset_color_index
if (~isoctave)
    ax=gca();
    set(ax,'ColorOrderIndex',1);
end
end % return

function o = isoctave
o = exist('OCTAVE_VERSION', 'builtin');
end % return

%------------------------------------------------

function dtc_plt(fig,flst1,flst2,Dh,pa,pk,m,loc)
if (m<1), return; end
lab=sprintf('%d-chamber, %s velocity tuning curve',m,loc);
tuncrv(fig,flst1,flst2,Dh,[-10 110],[-9 1],lab,loc,pa,pk);
end

function tuncrv(fig,flst1,flst2,Dh,dblim,phlim,lab,loc,pa,pk)
lab1=sprintf('nch=%d',pa.m);
lab2=pa.parlab;
flim=[0.2 20];
fk1=flst1/1000;
fk2=flst2/1000;
xpk=pk(1,:);
mpk=pk(2,:);
ppk=pk(3,:);
% select characteristic places
[nx,nf]=size(Dh);
s=2i*pi*flst2;
db1=zeros(nx,nf);
ph1=zeros(nx,nf);
dsp_ref = 140;    % dB re 1 nm displacement
spl_ref = 0.0002; % SPL reference pressure (rms Pa)
pe = 1e6 / spl_ref;
for i=1:nf
    vh=Dh(:,i)*s(i); % HB velocity (cm/sec)
    db1(:,i)=20*log10(abs(pe./vh))+pa.hbt-dsp_ref;
    ph1(:,i)=unwrap(angle(vh))/(2*pi);
end
ii=round(1+xpk*(nx-1));
nn=length(ii);
db2=db1(ii,:);
ph2=ph1(ii,:);
dbpk=ones(1,nn);
phpk=ones(1,nn);
fkpk=ones(1,nn);
for k=1:nn
    [~,kk]=min(db2(k,:));
    dbpk(k)=db2(k,kk);
    phpk(k)=ph2(k,kk);
    fkpk(k)=fk2(kk);
end
figure(fig);clf;
subplot(2,1,1);
reset_color_index
semilogx(fk2,db2,fkpk,dbpk,'b.',fk1,-mpk,'ro')
axis([flim dblim])
ylabel('Pe (dB SPL)')
title(lab);
text(10,15,lab1)
subplot(2,1,2)
reset_color_index
semilogx(fk2,ph2,fkpk,phpk,'b.')
axis([flim phlim])
ylim(phlim)
xlabel('frequency')
ylabel('phase (cyc)')
text(10,-7,lab2)
drawnow
fn1=sprintf('moh24%spf.txt',loc);
fn2=sprintf('moh24%stc.txt',loc);
write_data(fn1,[fk1' mpk' ppk' fkpk' dbpk' phpk']);
write_data(fn2,[fk2' db2' ph2']);
end

function midear_plt(fig,f,Dh,Zc,Gf,Ze)
if (fig==0), return; end
% Ze - middle-ear impedance (Hajicek 2024, ARO)
% Gf - middle-ear forward gain (Puria 2003, JASA)
% Zc - cochlea input impedance (Puria 2003, JASA)
ef=[0.1 0.2 0.5 1 2 5 10];                      % (kHz)
emZe=10^7*[22.8 13.3 6.5 4.8 4.1 4.0 4.9];      % (mks ohm)
emZc=1e9*[8.36 6.28 9.936 10.7 15.7 22.6 19.8]; % (mks ohm)
emGf=[-5.63 -2.77 11.0 16.6 11.7 0.25 -10.7];   % (dB)
%
fk=f/1000;
Zc=Zc*10^5;       % convert ftom cgs to mks
Ze=Ze*10^5;       % convert ftom cgs to mks
if (nargout<1)
    % plot
    figure(fig);clf
    subplot(2,2,1)
    loglog(fk,abs(Ze),'b',ef,emZe,'r');
    xlim([0.1 20])
    ylim([10^7 10^9])
    ylabel('Ze (mks ohm)')
    title('middle-ear impedance')
    subplot(2,2,2)
    loglog(fk,abs(Zc),'b',ef,emZc,'r');
    xlim([0.1 20])
    ylim([10^9 10^11])
    ylabel('Zc (mks ohm)')
    title('cochlea impedance')
    subplot(2,2,3)
    dbGf=20*log10(abs(Gf));
    semilogx(fk,dbGf,'b',ef,emGf,'r');
    xlim([0.1 20])
    ylim([-40 40])
    xlabel('frequency (Hz)')
    ylabel('Ps / Pe (dB)')
    title('middle-ear gain')
    [nx,nf]=size(Dh);
    if (nx>1)
        subplot(2,2,4)
        xl=35;
        fr=logspace(-1,log10(20),100);
        cp=f2p(fr,0.2,0.021,0.99,xl);
        xp=zeros(size(fk));
        for k=1:nf
            dh=abs(Dh(:,k));
            [~,ix]=max(dh);
            xp(k)=xl*(ix-1)/(nx-1);
        end
        semilogx(fk,xp,'b',fr,cp,'r');
        xlim([0.1 20])
        ylim([0 35])
        xlabel('frequency (Hz)')
        ylabel('place (mm)')
        title('frequency-place map')
    end
    drawnow
end
end

function Ze=midear_imp(pa,f,Zc)
s=2i*pi*f;
ast=pa.ast;
ama=pa.ama;
aed=pa.aed;
gme=pa.gme;
Zma=s*pa.mma+pa.rma+pa.kma./s;
Zim=pa.rim+pa.kim./s;
Zst=s*pa.mst+pa.rst+pa.kst./s;
Zed=s*pa.med+pa.red+pa.ked./s;
Zii=1./(1./Zim+1./(Zst+ast*Zc));
Zmi=Zma+gme*Zii;
Ze=1./(1./(ama*Zmi)+1./(aed*Zed));
end

%------------------------------------------------

function pk=pk_tgt(flst)
fpk=500*2.^(-1:5); % 250 500 1000 2000 4000 8000 16000
xcp= [ 0.81  0.71  0.59  0.47  0.33  0.19  0.05];
mpk=-[18.20  9.70  8.80 15.00 12.30 19.10 59.10];
ppk=-[ 3.20  3.65  3.83  3.74  3.37  2.60  1.10];
xcp=interp1(fpk,xcp,flst);
mpk=interp1(fpk,mpk,flst);
ppk=interp1(fpk,ppk,flst);
% specify Zc targets
emZc=1e4*[8.36 6.28 9.936 10.7 15.7 22.6 19.8]; % (cgs ohm)
zpk=interp1(500*2.^(-1:5),emZc,flst);
pk=[xcp;mpk;ppk;zpk];
end

function err=get_err(pa,flst,pk)
if (pa.m<1), err=0; return; end
pv=getpar(pa);
% fit parameters to reference data
err=refdev(pv,pa,flst,pk,0);
prnpar(pa,'fdm24par.ini');
end

%------------------------------------------------

function write_data(fn,data)
[nr,nc] = size(data);
fp=fopen(fn,'wt');
fprintf(fp,'; %s\n', fn);
for i=1:nr
    for j=1:nc
        fprintf(fp,' %14.5g',data(i,j));
    end
    fprintf(fp,'\n');
end
fclose(fp);
end

%------------------------------------------------

function report_phase(pa,flst)
pa.gam=1.0; [~,~,~,~,Dha]=fdmod23(pa,flst);  % active HB displacement
pa.gam=0.8; [~,~,~,~,Dhp]=fdmod23(pa,flst);  % passive HB displacement
nf=length(flst);
phta=zeros(nf,1); phtp=zeros(nf,1);
for k=1:nf
    [~,ipk]=max(log(abs(Dha(:,k))/flst(k))); % active magnitude
    pha=unwrap(angle(Dha(:,k)))/(2*pi);      % active phase
    php=unwrap(angle(Dhp(:,k)))/(2*pi);      % passive phase
    phta(k)=pha(ipk);                        % active phase target
    phtp(k)=php(ipk);                        % passive phase target
end
fprintf('phcp1=');fprintf(' %5.2f',-phta); fprintf('\n');
fprintf('phcp0=');fprintf(' %5.2f',-phtp); fprintf('\n');
end

function err=refdev(pv,pa,flst,pk,ot)
ipk=1+round(pk(1,:)*(pa.n-1));
Dt =pk(2,:);
Pt =pk(3,:);
Zt =pk(4,:);
err=0;
if (ot==0)
    no=10; ne=10; nq=10;
    io=1:no; ie=(1:ne)+no; iq=(1:nq)+(no+ne);  % pv indices
    mnpv=min(pv(io));
    if (mnpv<2e-6), err=1e9; return; end % minimum mass
    err=err+mean(abs(pv(ie)))*8;        % minimize linear exponents
    err=err+mean(abs(pv(iq)))*8;        % minimize quadratic exponents
    pa=setpar(pa,pv);
    if (pa.gpo>1), err=err+pa.gpo*100; end % limit partition gain < 1
elseif (ot==1)
    %pv=getparme(pa);
elseif (ot==2)
    %pv(1)=pa.khe;
    %pv(2)=pa.rhe;
    %pv(3)=pa.mhe;
elseif (ot==3)
    pa.chsz(1:2)=pv;
    pa.chsz=pa.chsz*(2/sum(pa.chsz));
end
% compute fdm
[~,~,Pd,~,Dh,Zc]=fdmod23(pa,flst);
D1=20*log10(abs(Dh));
[~,nf]=size(D1);
notch1=max(diff(log(abs(Pd)),2));
notch2=max(diff(log(abs(Dh)),2));
ff=20*log10(flst/250)-pa.hbt;
err=err*nf;
for k=1:nf
    ii=D1(:,k)>-950;                        % include all non-zero D1
    if (sum(ii)==0), err=1e6; return; end
    [~,ix]=max(D1(ii,k));                   % find CP
    ph=unwrap(angle(Dh(ii,k)))/(2*pi);      % unwrap phase
    mmx=D1(ix,k)+ff(k);              % find magnitude peak
    pmx=ph(ix);                             % find phase peak
    dph=diff(ph);
    dpp=dph(dph>0);
    if (~isempty(dpp)), err=err+mean(dpp)*4e4;  end % check phase slope
    err=err+notch1(k)*80;                   % Pd notch
    err=err+notch2(k)*80;                   % Dh notch
    err=err+abs(mmx-Dt(k))*2;               % peak
    err=err+abs(pmx-Pt(k))*2;               % phase
    err=err+abs((ix-ipk(k)))*8;             % fp-map
    % Zc comparison
    zdv=abs(log(Zc(k)/Zt(k)));
    err=err+zdv*8;
end
err=err/nf;
end

function pv=getpar(pa)
k=0;
k=k+1;pv(k)=pa.k1o;
k=k+1;pv(k)=pa.r1o;
k=k+1;pv(k)=pa.m1o;
k=k+1;pv(k)=pa.k2o;
k=k+1;pv(k)=pa.r2o;
k=k+1;pv(k)=pa.m2o;
k=k+1;pv(k)=pa.k3o;
k=k+1;pv(k)=pa.r3o;
k=k+1;pv(k)=pa.k4o;
k=k+1;pv(k)=pa.aco;
k=k+1;pv(k)=pa.k1e;
k=k+1;pv(k)=pa.r1e;
k=k+1;pv(k)=pa.m1e;
k=k+1;pv(k)=pa.k2e;
k=k+1;pv(k)=pa.r2e;
k=k+1;pv(k)=pa.m2e;
k=k+1;pv(k)=pa.k3e;
k=k+1;pv(k)=pa.r3e;
k=k+1;pv(k)=pa.k4e;
k=k+1;pv(k)=pa.ace;
k=k+1;pv(k)=pa.k1q;
k=k+1;pv(k)=pa.r1q;
k=k+1;pv(k)=pa.m1q;
k=k+1;pv(k)=pa.k2q;
k=k+1;pv(k)=pa.r2q;
k=k+1;pv(k)=pa.m2q;
k=k+1;pv(k)=pa.k3q;
k=k+1;pv(k)=pa.r3q;
k=k+1;pv(k)=pa.k4q;
k=k+1;pv(k)=pa.acq;
end

function pa=setpar(pa,pv)
k=0;
k=k+1;pa.k1o=pv(k);
k=k+1;pa.r1o=pv(k);
k=k+1;pa.m1o=pv(k);
k=k+1;pa.k2o=pv(k);
k=k+1;pa.r2o=pv(k);
k=k+1;pa.m2o=pv(k);
k=k+1;pa.k3o=pv(k);
k=k+1;pa.r3o=pv(k);
k=k+1;pa.k4o=pv(k);
k=k+1;pa.aco=pv(k);
k=k+1;pa.k1e=pv(k);
k=k+1;pa.r1e=pv(k);
k=k+1;pa.m1e=pv(k);
k=k+1;pa.k2e=pv(k);
k=k+1;pa.r2e=pv(k);
k=k+1;pa.m2e=pv(k);
k=k+1;pa.k3e=pv(k);
k=k+1;pa.r3e=pv(k);
k=k+1;pa.k4e=pv(k);
k=k+1;pa.ace=pv(k);
k=k+1;pa.k1q=pv(k);
k=k+1;pa.r1q=pv(k);
k=k+1;pa.m1q=pv(k);
k=k+1;pa.k2q=pv(k);
k=k+1;pa.r2q=pv(k);
k=k+1;pa.m2q=pv(k);
k=k+1;pa.k3q=pv(k);
k=k+1;pa.r3q=pv(k);
k=k+1;pa.k4q=pv(k);
k=k+1;pa.acq=pv(k);
end

function prnpar(pa,fn)
fp=fopen(fn,'wt');
fprintf(fp,'pa.k1o=%11.6g; pa.k1e=%7.4f; pa.k1q=%9.6f;\n',pa.k1o,pa.k1e,pa.k1q);
fprintf(fp,'pa.r1o=%11.6g; pa.r1e=%7.4f; pa.r1q=%9.6f;\n',pa.r1o,pa.r1e,pa.r1q);
fprintf(fp,'pa.m1o=%11.6g; pa.m1e=%7.4f; pa.m1q=%9.6f;\n',pa.m1o,pa.m1e,pa.m1q);
fprintf(fp,'pa.k2o=%11.6g; pa.k2e=%7.4f; pa.k2q=%9.6f;\n',pa.k2o,pa.k2e,pa.k2q);
fprintf(fp,'pa.r2o=%11.6g; pa.r2e=%7.4f; pa.r2q=%9.6f;\n',pa.r2o,pa.r2e,pa.r2q);
fprintf(fp,'pa.m2o=%11.6g; pa.m2e=%7.4f; pa.m2q=%9.6f;\n',pa.m2o,pa.m2e,pa.m2q);
fprintf(fp,'pa.k3o=%11.6g; pa.k3e=%7.4f; pa.k3q=%9.6f;\n',pa.k3o,pa.k3e,pa.k3q);
fprintf(fp,'pa.r3o=%11.6g; pa.r3e=%7.4f; pa.r3q=%9.6f;\n',pa.r3o,pa.r3e,pa.r3q);
fprintf(fp,'pa.k4o=%11.6g; pa.k4e=%7.4f; pa.k4q=%9.6f;\n',pa.k4o,pa.k4e,pa.k4q);
fprintf(fp,'pa.r4o=%11.6g; pa.r4e=%7.4f; pa.r4q=%9.6f;\n',pa.r4o,pa.r4e,pa.r4q);
fprintf(fp,'pa.gpo=%11.6g; pa.gpe=%7.4f; pa.gpq=%9.6f;\n',pa.gpo,pa.gpe,pa.gpq);
fprintf(fp,'pa.aco=%11.6g; pa.ace=%7.4f; pa.acq=%9.6f;\n',pa.aco,pa.ace,pa.acq);
fprintf(fp,'pa.bwo=%11.6g; pa.bwe=%7.4f; pa.bwq=%9.6f;\n',pa.bwo,pa.bwe,pa.bwq);
fclose(fp);
end

%==============================================================

function [xx,Yb,Pd,Db,Dh,Zc,Gf]=fdmod23(pa,flst)
if (pa.m<1)
    xx=linspace(0,1,pa.n);
    Yb=xx;Pd=xx;Db=xx;Dh=xx;
    s=2i*pi*flst;
    Zc=pa.rco+pa.mco*s;
    ast=pa.ast;
    ama=pa.ama;
    gme=pa.gme;
    Zma=s*pa.mma+pa.rma+pa.kma./s;
    Zim=pa.rim+pa.kim./s;
    Zst=s*pa.mst+pa.rst+pa.kst./s;
    Gf=(ama*Zc)./(gme*(Zst+ast*Zc)+Zma.*(Zst+ast*Zc+Zim)./(gme*Zim));
    return;
end
L=pa.xl;
N=pa.n;
ast=pa.ast;
nf=length(flst);
dx=L/(N-1);
nc=length(pa.chsz);
Yb=zeros(N,nf);
Pd=zeros(N,nf);
Db=zeros(N,nf);
Dh=zeros(N,nf);
Zc=zeros(1,nf);
% loop over frequencies
for k=1:nf
    f=flst(k);
    s=2i*pi*f;
    srd=s*pa.rho*dx;
    [a,q,x,y,D]=mxfill(f,pa);    % initialize arrays
    p=mxsolve(a,q);
    [pd,yb,dd]=p_dif(p,y,f,D);
    [pd,yb,dd]=p_nan(pd,yb,dd);  % restrict plotting range
    Pd(:,k)=pd;
    Yb(:,k)=yb;
    Db(:,k)=dd(:,1);
    Dh(:,k)=dd(:,2);
    Vs=(p(nc,2)-p(nc,1))/srd;    % stapes linear velocity
    Zc(k)=p(nc,1)/(ast*Vs);      % cochlear impedance
end
xx=x/pa.xl;
Gf=Pd(1,:);
end

%--------------------------------------------------------------

function [a,q,x,y,D]=mxfill(f,pa)
% local parameter values
chsz=pa.chsz;
chsz=chsz*(2/sum(chsz)); % normalize channel sizes;
n=pa.n;
xl=pa.xl; rho=pa.rho; g=pa.gam;
% useful constructs
s=2i*pi*f;
dx=xl/(n-1);
srd=s*rho*dx;
% compute admittance for all x
x=transpose(linspace(0,xl,n));
[z1,z2,zh,za,ac,gh,bw]=imped(x,s,pa);
m=length(pa.chsz);
A=zeros(n,m);
Y=zeros(n,m,m);
a=zeros(3,m,m,n);
y=zeros(n,m,m);
% initialize apex
Zh=pa.khe/s+pa.rhe+pa.mhe*s;
z1(n)=Zh;z2(n)=1e-6;zh(n)=1e-6;za(n)=0;
% initialize partition admittance
zg=g*za;
zk=zeros(m,m);
if (m==1)
    D=1;                            % select PD
    B=2;                            % basal BC
elseif (m==2)
    D=[1 -1;0 0];                   % select PD
    B=[-1;1];                       % basal BC
elseif (m==3)
    D=[1 0 -1;0 -1 1; 0 0 0];       % select PD
    B=[-1;0;1];                     % basal BC
end
for k=1:n
    A(k,:)=ac(k)*chsz;
    if (m==1)
        hh=z2(k)/(z2(k)+zh(k));
        zk(1)=z1(k)+(gh(k)*zh(k)-zg(k))*hh;
    elseif (m==2)
        hh=z2(k)/(z2(k)+zh(k));
        zk(1,1)=z1(k)+(gh(k)*zh(k)-zg(k))*hh;
    elseif (m==3)
        zk(1,1:2)=[ z1(k) gh(k)*zh(k)-zg(k)];
        zk(2,1:2)=[-z2(k) z2(k)+zh(k)      ];
    end
    zk(m,:)=1; % conserve fluid volume
    if (m==3)
        T = [1, 0; 1, -1; -2, 1];
        Y(k,:,:) = T * (zk(1:2, 1:2) \ D(1:2, :));
    else
        Y(k,:,:) = zk \ D;
    end
    if (m==1)
        y(k,1,1)=1/zk(1,1);   % select Ybm
        y(k,2,:)=hh*y(k,1,:); % select Yhb
    elseif (m==2)
        y(k,1,1)=1/zk(1,1);   % select Ybm
        y(k,2,:)=hh*y(k,1,:); % select Yhb
    elseif (m==3)
        y(k,:,:)=inv(zk);     % select Ybm & Yhb
    end
end
% initialize block-tridiagonal matrix a
for k=2:(n-1)
    A0=diag(A(k,:));
    A1=diag(A(k+1,:));
    yr=bw(k)*dx*squeeze(Y(k,:,:));
    a(1,:,:,k)=-A0;
    a(2,:,:,k)=A0+A1+srd*yr;
    a(3,:,:,k)=-A1;
end
% fill matrix
A1=2*diag(A(1,:));
A2=-A1;
An=diag(A(n,:));
a(1,:,:,1)=zeros(m,m);
a(2,:,:,1)=A1;
a(3,:,:,1)=A2;
a(1,:,:,n)=-An;
a(2,:,:,n)=An;
a(3,:,:,n)=zeros(m,m);
% middle-ear boundary conditions
[alfx,qst]=midear(f,pa);
alfx=alfx*1e-7;
a(2,:,:,1)=A1*(1+alfx);
a(3,:,:,1)=A2*( -alfx);
% initialize q
q=zeros(m,n);
q(:,1)=qst*B;
end % return

function [z1,z2,z3,z4,ac,gh,bw]=imped(x,s,pa)
x=x.*(1+(pa.xtap*x).^pa.xtex);
q=x.^2;
k1=pa.k1o*exp(pa.k1e*x+pa.k1q*q);
if isfield(pa, 'rough_amp')
    rng(42); % Fixed seed for reproducible roughness
    k1 = k1 .* (1 + pa.rough_amp * rand(size(k1)));
end
r1=pa.r1o*exp(pa.r1e*x+pa.r1q*q);
m1=pa.m1o*exp(pa.m1e*x+pa.m1q*q);
k2=pa.k2o*exp(pa.k2e*x+pa.k2q*q);
r2=pa.r2o*exp(pa.r2e*x+pa.r2q*q);
m2=pa.m2o*exp(pa.m2e*x+pa.m2q*q);
k3=pa.k3o*exp(pa.k3e*x+pa.k3q*q);
r3=pa.r3o*exp(pa.r3e*x+pa.r3q*q);
k4=pa.k4o*exp(pa.k4e*x+pa.k4q*q);
r4=pa.r4o*exp(pa.r4e*x+pa.r4q*q);
z1=k1/s+r1+m1*s;
z2=k2/s+r2+m2*s;
z3=k3/s+r3;
z4=k4/s+r4;
ac=pa.aco*exp(pa.ace*x+pa.acq*q);
gh=pa.gpo*exp(pa.gpe*x+pa.gpq*q);
bw=pa.bwo*exp(pa.bwe*x+pa.bwq*q);
end

function [alfx,qst]=midear(f,pa)
% middle-ear impedances
s=2i*pi*f(:);
dx=pa.xl/(pa.n-1);
srd=s*pa.rho*dx;
ast=pa.ast;
ama=pa.ama;
gme=pa.gme;
if(pa.nmev==1)
    Zme=pa.mme*s+pa.rme+pa.kme./s;
    alfx=(Zme)./(srd*ast);
    qst =(ama/gme)/ast;
else
    Zma=s*pa.mma+pa.rma+pa.kma./s;
    Zim=pa.rim+pa.kim./s;
    Zst=s*pa.mst+pa.rst+pa.kst./s;
    gim=Zim./(Zma+pa.gme^2*Zim);
    alfx=(Zst+Zma.*gim)/(srd*ast)/1000; % <-- why ??
    qst =(ama*gme*gim)/50; % <-- why ??
end
end

function p=mxsolve(a,q)
% solve tri-diagnonal block matrix equation
n=length(q);
% work down from top
b=squeeze(a(2,:,:,1));
a(3,:,:,1)=b\squeeze(a(3,:,:,1));
q(:,1)=b\q(:,1);
for k=2:n
   c=squeeze(a(1,:,:,k))*squeeze(a(3,:,:,k-1));
   p=squeeze(a(1,:,:,k))*q(:,k-1);
   b=squeeze(a(2,:,:,k))-c;
   a(3,:,:,k)=b\squeeze(a(3,:,:,k));
   q(:,k)=b\(q(:,k)-p);
end
% work up from bottom
for k=(n-1):-1:1
   q(:,k)=q(:,k)-squeeze(a(3,:,:,k))*q(:,k+1);
end
p=q;
end

function [pd,yd,dd]=p_dif(p,y,f,D)
pref=2.848e-4; % dB SPL pressure reference (dyn/cm^2)
dref=1e-7;     % 1-nm displacement reference (cm)
n=size(p,2);
pd=zeros(n,1);
yd=zeros(n,1);
vd=zeros(n,2);
for k=1:n
   pk=D*p(:,k);
   yk=squeeze(y(k,:,:));
   vk=yk*pk;
   vd(k,:)=vk(1:2)*(pref/dref); % BM,HB velocity
   pd(k)=pk(1);                 % ST-SV pressure difference
   yd(k)=vk(1)/pd(k);           % BM admittance
end
dd=vd./(2i*pi*f);               % BM,HB displacement
end

function [pd,yd,dd]=p_nan(pd,yd,dd)
% use nan to restrict plotted range
n=length(pd);
for k=1:n
   if(abs(pd(k,1))<1e-4)
      ii=k:n;
      pd(ii,:)=nan;
      yd(ii,:)=nan;
      dd(ii,:)=nan;
      break;
   end
end
end

%==============================================================

function pa=modpar26(nch)
% chamber sizes
if (nch<3)
    pa=par_CEL16;
else
    pa=modpar26c3;
end
pa.m=nch;             % number of fluid chambers
pa.xp=1;              % length of excitation pattern
% TDM parameters
pa.nstp = 20480; pa.ncyc = 1; pa.ntsw = 10; pa.ntsf = 2048; pa.dt = 2e-6; pa.nimp = 1;
pa.ME_VEP_SCALE = 10; pa.ME_VST_SCALE_PHANTOM = 2; pa.IHC_SUBSTEPS = 10;
pa.gampro = ones(pa.n,1); 
pa.synpro = ones(pa.n,1);
pa.met=200; pa.hco=1e-4; pa.hcs=1e-4;
end

%--------------------------------------------------------------

function pa=par_CEL16
pa.parlab='CEL16';
pa.chsz=[1 1];                 % chamber size
pa.gam = 1;                    % NDR multiplier
pa.m = 1;                      % number of points across fluid
pa.n = 1401;                   % number of points along BM
pa.xl = 3.5;                   % scala length
pa.yw = 0.1;                   % scala width
pa.zh = 0.1;                   % scala height
pa.rho = 1;                    % fluid density
pa.stgain=0.0127;              % stimulus gain (target 40 dB SPL tone)
pa.isv=[1136 1005 840 655 466 273 80]; % BM locations to save
pa.hbt=8.5; pa.xtap=0.2339; pa.xtex=6; pa.mmeq=1;
pa.khe=0.0001; pa.rhe=0.0001; pa.mhe=0.01; % err=23.21 23.21
% middle-ear parameters -------------
pa.nmev=4; pa.kme=0.1; pa.rme=400; pa.mme=0.1;
pa.acp=0.5; pa.adi=0.2; pa.aed=0.5; pa.ama=0.5; pa.ast=0.03;
pa.cep=2700;  pa.rep=29.16; pa.mevgn=1.262; pa.stim=3;
pa.mcp=0.0002; pa.rcp=0.1; pa.kcp=2200; pa.rfz=0.01;
pa.mdi=0.005;  pa.rdi=100;  pa.kdi=7.7e+06;
pa.mrw=5e-3; pa.rrw=20; pa.krw=1.5e5; pa.arw=0.0625;
pa.mco=3; pa.rco=1e5;  pa.gme=0.5;
pa.kma=3e+06;  pa.rma=8000;  pa.mma=0.85;
pa.ked=7.04632e+06;  pa.red=200;  pa.med=0.02;
pa.kst=4.66022e+06;  pa.rst=8140.1;  pa.mst=0.89833;
pa.kim=3e+08;  pa.rim=28508.7; % nch=0 nme=4
% hair-cell & neural parameters --------
pa.ihceq=4; % 4 = Neely synapse model
pa.hbnl=0; pa.hbsc=0.04; pa.hbmx=6e-9; pa.ihcv=1;
pa.ihctc=0.2e-3; pa.ihcsf=1e5; pa.nrgn=2.5e6; pa.ew=0.1;
pa.ihcex=14; pa.ihcrr=17.836; pa.ihcdr=1.2637;
pa.ldew=2; pa.ldne=140; pa.ldsc=[1.360 7.057 3.518 18.93 63.79];
%---------- partition impedance
pa.k1o=  2.394e+08; pa.k1e=-2.6655; pa.k1q= 0.000000; % CEL16
pa.r1o=      871.4; pa.r1e=-1.3937; pa.r1q= 0.000000;
pa.m1o=   0.007417; pa.m1e=-0.1884; pa.m1q= 0.000000;
pa.k2o=  3.036e+08; pa.k2e=-3.4762; pa.k2q= 0.000000;
pa.r2o=       1979; pa.r2e=-1.2466; pa.r2q= 0.000000;
pa.m2o=    0.03417; pa.m2e=-0.0828; pa.m2q= 0.000000;
pa.k3o=  3.151e+08; pa.k3e=-2.9092; pa.k3q= 0.000000;
pa.r3o=      1.049; pa.r3e= 0.1033; pa.r3q= 0.000000;
pa.k4o=  4.045e+08; pa.k4e=-2.7948; pa.k4q= 0.000000;
pa.r4o=          0; pa.r4e= 0.0000; pa.r4q= 0.000000;
pa.gpo=          1; pa.gpe= 0.0000; pa.gpq= 0.000000;
pa.aco=       0.01; pa.ace=-0.4000; pa.acq= 0.000000;
pa.bwo=       0.05; pa.bwe= 0.0000; pa.bwq= 0.000000;
end % return

%--------------------------------------------------------------

function pa=modpar26c3
pa.parlab='cel26c3';
pa.chsz=[0.95 0.05 1.0];
pa.gam = 1;                    % NDR multiplier
pa.m = 1;                      % number of points across fluid
pa.n = 1401;                   % number of points along BM
pa.xl = 3.5;                   % scala length
pa.yw = 0.1;                   % scala width
pa.zh = 0.1;                   % scala height
pa.rho = 1;                    % fluid density
pa.stgain=0.0127;              % stimulus gain (target 40 dB SPL tone)
pa.isv=[1136 1005 840 655 466 273 80]; % BM locations to save
pa.hbt=33; pa.xtap=0.2339; pa.xtex=6; pa.mmeq=1;
pa.khe=0.0001; pa.rhe=0.0001; pa.mhe=0.01; % err=23.21 23.21
% middle-ear parameters -------------
pa.nmev=4; pa.kme=0.1; pa.rme=400; pa.mme=0.1;
pa.acp=0.5; pa.adi=0.2; pa.aed=0.5; pa.ama=0.5; pa.ast=0.03;
pa.cep=2700;  pa.rep=29.16; pa.mevgn=1.262; pa.stim=3;
pa.mcp=0.0002; pa.rcp=0.1; pa.kcp=2200; pa.rfz=0.01;
pa.mdi=0.005;  pa.rdi=100;  pa.kdi=7.7e+06;
pa.mrw=5e-3; pa.rrw=20; pa.krw=1.5e5; pa.arw=0.0625;
pa.mco=3; pa.rco=1e5;  pa.gme=0.5;
pa.kma=3e+06;  pa.rma=8000;  pa.mma=0.85;
pa.ked=7.04632e+06;  pa.red=200;  pa.med=0.02;
pa.kst=4.66022e+06;  pa.rst=8140.1;  pa.mst=0.89833;
pa.kim=3e+08;  pa.rim=28508.7; % nch=0 nme=4
% hair-cell & neural parameters --------
pa.ihceq=4; % 4 = Neely synapse model
pa.hbnl=0; pa.hbsc=0.04; pa.hbmx=6e-9; pa.ihcv=1;
pa.ihctc=0.2e-3; pa.ihcsf=1e5; pa.nrgn=2.5e6; pa.ew=0.1;
pa.ihcex=14; pa.ihcrr=17.836; pa.ihcdr=1.2637;
pa.ldew=2; pa.ldne=140; pa.ldsc=[1.142 8.296 4.841  8.05 34.00];
%---------- partition impedance
pa.k1o=2.41442e+08; pa.k1e=-2.7377; pa.k1q= 0.000002;
pa.r1o=     1056.6; pa.r1e=-1.3138; pa.r1q= 0.000002;
pa.m1o= 0.00745053; pa.m1e=-0.2012; pa.m1q= 0.000002;
pa.k2o=3.05084e+08; pa.k2e=-3.5063; pa.k2q= 0.000002;
pa.r2o=    1984.19; pa.r2e=-1.2479; pa.r2q= 0.000002;
pa.m2o=  0.0360276; pa.m2e=-0.0812; pa.m2q= 0.000002;
pa.k3o=3.15082e+08; pa.k3e=-2.9231; pa.k3q= 0.000003;
pa.r3o=    1.10392; pa.r3e= 0.1065; pa.r3q= 0.000002;
pa.k4o=4.04629e+08; pa.k4e=-2.8017; pa.k4q= 0.000001;
pa.r4o=          0; pa.r4e= 0.0000; pa.r4q= 0.000000;
pa.gpo=     0.9956; pa.gpe= 0.0000; pa.gpq= 0.000003;
pa.aco= 0.00994045; pa.ace=-0.4324; pa.acq= 0.000002;
pa.bwo=  0.0517652; pa.bwe= 0.0001; pa.bwq= 0.000002;
end % return
