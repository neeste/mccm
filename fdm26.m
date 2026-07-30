% fdm26 - multi-chamber frequency-domain model of cochleamodpar26
function R=fdm26(nch)
if (nargin<1), nch=0; end % number of cochlear fluid chambers
% pa.tmrigid IS NOT IMPLEMENTED HERE YET. tdm26 clamps d2 kinematically
% (a(i2) = 0 in fold_p); the frequency-domain equivalent is to eliminate V2 from
% zk, which is a real restructuring of every m branch, not a flag.
% REFUSE rather than silently score an UNCLAMPED model while tdm26 clamps. That
% exact fdm/tdm divergence cost a full session on 2026-07-29: fdm26 had no third
% DOF at m=3 while tdm26 had one, and g4_maperr_m3b reported PASS on the
% disagreement for months.
if (isstruct(nch) && isfield(nch,'pa') && isstruct(nch.pa) ...
        && isfield(nch.pa,'tmrigid') && nch.pa.tmrigid)
    error('fdm26:tmrigid', ...
        ['pa.tmrigid (TM clamped) is not implemented in fdm26. tdm26 clamps d2 ' ...
         'but fdm26 would score it FREE, so the two solvers would disagree ' ...
         'silently.\nEliminate V2 from zk in the relevant m branch first, or ' ...
         'measure tmrigid in the time domain only.']);
end
if (isstruct(nch))                                        % struct arg => analysis mode
    if (isfield(nch,'xprofile')),    R=xprofile(nch);
    elseif (isfield(nch,'xreflect')),R=xreflect(nch);
    elseif (isfield(nch,'fdsolve')), R=fd_solve(nch);
    elseif (isfield(nch,'mreflect')),R=modal_reflectance(nch);
    elseif (isfield(nch,'modal')),   R=modal_wave_analysis(nch);
    elseif (isfield(nch,'macroD')),  R=macroD_export(nch);
    elseif (isfield(nch,'cfmap')),   R=cfmap_analysis(nch);
    elseif (isfield(nch,'reflect')), R=reflect_analysis(nch);
    else,                            R=grpdelay_analysis(nch); end
    return;
end
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
    if (isfield(pa,'chsznorm') && pa.chsznorm), pa.chsz=pa.chsz*(2/sum(pa.chsz)); end
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

% Frequency-domain cochlear group delay (forward-latency analog), Month-2 pivot.
% tau(f_cf) = -d(phase)/d(omega) of the HB response at each place's CF; slope d
% (tau ~ f^-d) is the ABR-latency slope to steer to 0.413, with NO stability wall.
function R=grpdelay_analysis(pr)
if (isfield(pr,'pa')), pa=pr.pa;
elseif (isfield(pr,'nch')), pa=modpar26(pr.nch);
else, pa=modpar26(1); end
if (isfield(pr,'flst')), flst=pr.flst; else, flst=logspace(log10(200),log10(20000),700); end
[~,~,~,~,Dh]=fdmod23(pa,flst);
[N,nf]=size(Dh); w=2*pi*flst(:);
places=unique(round(linspace(round(0.04*N),round(0.96*N),140)));
f_cf=nan(1,numel(places)); tau=nan(1,numel(places)); gain=nan(1,numel(places));
for idx=1:numel(places)
    p=places(idx); mag=abs(Dh(p,:)); [pkmag,kk]=max(mag);
    if (kk<3||kk>nf-2), continue; end
    ph=unwrap(angle(Dh(p,:)));                       % radians
    dphdw=(ph(kk+1)-ph(kk-1))/(w(kk+1)-w(kk-1));
    tau(idx)=-dphdw*1000;                            % ms
    f_cf(idx)=flst(kk)/1000;                         % kHz
    gain(idx)=20*log10(pkmag/(abs(Dh(p,1))+eps));    % tip re basal-tail (dB)
end
sel=(f_cf>=0.5)&(f_cf<=4)&isfinite(tau)&(tau>0);
if (sum(sel)>=2)                                     % need >=2 points for a slope fit
    pf=polyfit(log10(f_cf(sel)),log10(tau(sel)),1); R.d=-pf(1);
else
    R.d=NaN;                                         % fd slope undefined (e.g. 3-chamber tip-tail collapses the window)
end
R.f=f_cf; R.tau=tau; R.gain=gain; R.sel=sel;
R.tipgain=max(gain(isfinite(gain)));
flmap=500*2.^(-1:5); R.maperr=tuning_err(pa,flmap,pk_tgt(flmap));  % MAP/tuning-fit constraint
end

% MOH2008 stapes reflectance (the SFOAE / retrograde wave), from the model's
% wave variables:  Zs = s*rho/A (A = ac*chsz1),  Yv = bw*Yb,  kappa = sqrt(Zs*Yv),
% z0 = sqrt(Zs/Yv),  eps = 0.5 d/dx ln z0,  R(f) = INT eps*exp(-2 INT kappa) dx.
% Emissions arise where z0 varies rapidly (eps large): imposed BM roughness in the
% 1-chamber model, or -- the MOH2027 hypothesis -- the abrupt tip-tail transition
% of the 3-chamber model, with no roughness.
function R=reflect_analysis(pr)
if (isfield(pr,'pa')), pa=pr.pa; elseif (isfield(pr,'nch')), pa=modpar26(pr.nch); else, pa=modpar26(1); end
if (isfield(pr,'flst')), flst=pr.flst; else, flst=linspace(500,4000,1500); end
if (pa.m>=3), error('reflect_analysis: 3-chamber wave decomposition not yet implemented (needs the 2-DOF extension)'); end
n=pa.n; xl=pa.xl; dx=xl/(n-1); x=linspace(0,xl,n)';
chsz=pa.chsz; if (isfield(pa,'chsznorm') && pa.chsznorm), chsz=chsz*(2/sum(chsz)); end
rho=pa.rho; gam=pa.gam;
ii=2:(n-1);                                     % drop the boundary points
nf=numel(flst); refl=zeros(1,nf);
[~,kref]=min(abs(flst-median(flst)));           % reference freq for the kappa-profile diagnostic
for kf=1:nf
    s=2i*pi*flst(kf);
    [z1,z2,zh,za,ac,gh,bw]=imped(x,s,pa);       % unmasked; imped applies the x-warp
    hh=z2./(z2+zh);
    Yb=1./(z1 + (gh.*zh - gam.*za).*hh);        % 1-/2-chamber BM admittance
    A=ac*chsz(1);                               % scala area
    Zs=s*rho./A; Yv=bw.*Yb;                     % scala impedance, shunt admittance (per unit length)
    kap=sqrt(Zs.*Yv);                           % wavenumber
    fl=imag(kap)<0; kap(fl)=-kap(fl);           % forward wave propagates apically: Im(kappa)>=0
    z0=sqrt(Zs./Yv);                            % characteristic impedance
    eps_=0.5*gradient(log(z0))/dx;              % eps = 0.5 d/dx ln z0
    cumk=cumsum(kap)*dx;                        % INT_0^x kappa
    refl(kf)=sum(eps_(ii).*exp(-2*cumk(ii)))*dx;% Eq 18
    if (kf==kref), R.kap=kap; R.eps=eps_; R.x=x; R.fref=flst(kf); end
end
w=2*pi*flst(:)';
R.f=flst; R.refl=refl; R.pa=pa;
R.mag=20*log10(abs(refl)+eps); R.ripple=R.mag-movmean(R.mag,50);
R.gd=-gradient(unwrap(angle(refl)))./gradient(w)*1000;   % round-trip group delay (ms)
end

% ---- Multi-chamber wave decomposition (the 3-chamber SFOAE extension) ----------
% The FD coupling  -A0 p(k-1) + (A0+A1+srd*bw*dx*Y) p(k) - A1 p(k+1) = 0  is the
% discretization of the COUPLED transmission line  d/dx[A dp/dx] = s*rho*bw*Y*p,
% i.e. the m-vector first-order system
%       dP/dx = -Zs*U ,  dU/dx = -Yv*P ,   Zs = s*rho*A^-1 (diag) ,  Yv = bw*Y  (m x m).
% The modal wavenumbers are kappa_i = sqrt(eig(Zs*Yv)) (branch Im>=0). For m=1 this
% reduces to the scalar kappa; for m>=3 there are m coupled modes and the ACOUSTIC
% (slow, tonotopic) mode is the SFOAE carrier. modal_wave_analysis returns the modal
% dispersion, continuity-tracked, with the acoustic mode identified.
function R=cfmap_analysis(pr)
% CF MAP: peak PLACE per frequency -- the model's frequency-place map.
%
% Exists because fdmod23/tuning_err/imped/mxfill are LOCAL subfunctions and are
% unreachable from an external script, so a cross-check against tdm26 cannot be
% assembled from outside this file. Four external harnesses failed on that.
% This computes the map exactly as tuning_err does: max over PLACE of |Dh| (and
% |Db|) at each frequency in flst.
%
%   pr.pa   or pr.nch   model (default 3-chamber)
%   pr.flst frequency list, Hz (default 500*2.^(-1:5) = the maperr grid)
% Returns peak place as a normalized position x/L in [0,1] (larger = apical),
% directly comparable to tdm26's isv/pa.n convention.
if (isfield(pr,'pa')), pa=pr.pa; elseif (isfield(pr,'nch')), pa=modpar26(pr.nch); else, pa=modpar26(3); end
if (isfield(pr,'flst')), flst=pr.flst; else, flst=500*2.^(-1:5); end
[xx,~,~,Db,Dh,~]=fdmod23(pa,flst);
nf=numel(flst);
R.f=flst; R.xpk_bm=nan(1,nf); R.xpk_hb=nan(1,nf); R.amp_bm=nan(1,nf);
for k=1:nf
    a=abs(Db(:,k)); a(~isfinite(a))=0; [pk,ip]=max(a);
    if (pk>0), R.xpk_bm(k)=ip/pa.n; R.amp_bm(k)=20*log10(pk); end
    b=abs(Dh(:,k)); b(~isfinite(b))=0; [pk2,ip2]=max(b);
    if (pk2>0), R.xpk_hb(k)=ip2/pa.n; end
end
R.pa=pa; R.xx=xx; R.n=pa.n;
end

function R=modal_wave_analysis(pr)
if (isfield(pr,'pa')), pa=pr.pa; elseif (isfield(pr,'nch')), pa=modpar26(pr.nch); else, pa=modpar26(3); end
if (isfield(pr,'f')), f=pr.f; else, f=2000; end
n=pa.n; rho=pa.rho; s=2i*pi*f;
[~,~,x,~,~,Y,A,bw]=mxfill(f,pa);                 % coupling admittance Y (n x m x m), areas A (n x m)
m=size(A,2);                                     % true coupling dimension (=length(chsz), not pa.m)
kall=zeros(n,m); Vall=zeros(n,m,m);
for k=1:n
    Zs=s*rho*diag(1./A(k,:));                    % m x m series impedance (diagonal)
    Yv=bw(k)*squeeze(Y(k,:,:));                  % m x m shunt admittance
    M=Zs*Yv;
    [V,Dg]=eig(M);
    kap=sqrt(diag(Dg)); fl=imag(kap)<0; kap(fl)=-kap(fl);   % forward branch Im(kappa)>=0
    for j=1:m, V(:,j)=V(:,j)/norm(V(:,j)); end   % unit-norm mode shapes
    kall(k,:)=kap.'; Vall(k,:,:)=V;
end
% continuity-track the modes apically by eigenvector overlap
for k=2:n
    pV=squeeze(Vall(k-1,:,:)); cV=squeeze(Vall(k,:,:));
    ov=abs(pV'*cV); perm=zeros(1,m); used=false(1,m);
    for i=1:m                                    % greedy best-overlap assignment
        [~,ord]=sort(ov(i,:),'descend');
        for c=ord, if (~used(c)), perm(i)=c; used(c)=true; break; end, end
    end
    kall(k,:)=kall(k,perm); Vall(k,:,:)=cV(:,perm);
    % keep eigenvector phase continuous (align sign to previous)
    for j=1:m, if (real(sum(conj(pV(:,j)).*squeeze(Vall(k,:,j)')))<0), Vall(k,:,j)=-Vall(k,:,j); end, end
end
% acoustic mode = smallest NON-TRIVIAL |kappa| at a basal reference. The symmetric
% scalae leave a trivial common/compression mode (kappa~0, p1=p2, no BM drive) that
% must be excluded; among the propagating modes the acoustic wave is the longest-
% wavelength (smallest |kappa|) one at the base -- it carries the tonotopic map.
kref=3; kabs=abs(kall(kref,:));
cand=find(kabs > 1e-3*max(kabs));               % drop the trivial kappa~0 mode(s)
[~,j]=min(kabs(cand)); ia=cand(j);
ka=kall(:,ia);
% CF place for this f: where the acoustic wavenumber peaks (turning point)
[~,icf]=max(real(ka));
R.f=f; R.x=x; R.kall=kall; R.ka=ka; R.iac=ia; R.icf=icf; R.pa=pa; R.m=m;
R.Vall=Vall;   % n x m x m mode SHAPES in chamber-pressure space. Needed to say
               % WHICH chambers each mode occupies: the primary acoustic mode is
               % an ST-vs-SV differential, and a secondary SS<->CL wave (RL as
               % the partition, the structural analog of BM between SV and ST)
               % would appear as an SS-vs-CL differential.
end

% EXACT stapes reflectance for the coupled model, following fdm12a's press_ex:
% solve the actual pressure field, then decompose the acoustic (differential)
% pressure into forward + backward waves and read r = P-(base)/P+(base). The DEFAULT
% is a ROBUST windowed least-squares fit over a basal window: the field Pd(x) is fit
% to A*exp(-kappa_a*(x-xc)) + B*exp(+kappa_a*(x-xc)) and r=B/A. Fitting over many
% nodes averages out the single-node finite-difference noise AND the second-DOF mode
% admixture (which does not match the acoustic basis) -- both of which corrupt the
% 3-chamber single-node estimate (pr.win=0 selects the old single-node method).
function R=xreflect(pr)
if (isfield(pr,'pa')), pa=pr.pa; elseif (isfield(pr,'nch')), pa=modpar26(pr.nch); else, pa=modpar26(3); end
if (isfield(pr,'flst')), flst=pr.flst; else, flst=linspace(500,4000,600); end
rho=pa.rho; n=pa.n; nf=numel(flst);
robust = ~(isfield(pr,'win') && pr.win==0);
kb=max(6,round(0.02*n)); if (isfield(pr,'kb')), kb=pr.kb; end   % basal reference node
refl=zeros(1,nf); Ppv=zeros(1,nf); Pmv=zeros(1,nf);
kapv=zeros(1,nf); havek=isfield(pr,'kap_fixed');
for kf=1:nf
    f=flst(kf); s=2i*pi*f;
    [a,q,x,y,D,Y,A,bw]=mxfill(f,pa); p=mxsolve(a,q);         % raw m x n scala pressures
    Pd=(D(1,:)*p).';                                         % differential (acoustic) pressure, unmasked
    dx=x(2)-x(1);
    kref=acoustic_kappa(A,bw,Y,kb,s,rho);
    if (robust)
        khi=min(round(0.22*n), kb+max(20,round(pi/max(imag(kref),1)/dx)));  % ~1/2 wavelength, below CF
        wij=(kb:khi).'; ns=min(numel(wij),12); ss=round(linspace(kb,khi,ns));
        kas=arrayfun(@(kk) acoustic_kappa(A,bw,Y,kk,s,rho), ss);            % kappa_a(x), subsampled
        kaw=interp1(x(ss),kas,x(wij),'linear','extrap');
        phi=cumsum(kaw)*dx; phi=phi-phi(round(numel(wij)/2));              % WKB phase, ref at window center
        Mmat=[exp(-phi), exp(phi)]; AB=Mmat\Pd(wij);                       % least-squares fwd/bwd fit
        Pp=AB(1); Pm=AB(2); kap=kaw(round(numel(wij)/2));
    else
        if (havek), kref=pr.kap_fixed(kf); end
        kap=kref;
        dPd=(Pd(kb+1)-Pd(kb-1))/(2*dx);
        Pp=(Pd(kb)-dPd/kap)/2; Pm=(Pd(kb)+dPd/kap)/2;       % single-node decomposition
    end
    refl(kf)=Pm/Pp; Ppv(kf)=Pp; Pmv(kf)=Pm; kapv(kf)=kap;
end
w=2*pi*flst(:)';
R.Pp=Ppv; R.Pm=Pmv; R.kap=kapv;
R.f=flst; R.refl=refl; R.pa=pa; R.kb=kb;
R.mag=20*log10(abs(refl)+eps); R.ripple=R.mag-movmean(R.mag,50);
R.gd=-gradient(unwrap(angle(refl)))./gradient(w)*1000;       % round-trip group delay (ms)
R.gdcyc=-gradient(unwrap(angle(refl)))./gradient(log(flst(:)'))/(2*pi);  % delay in cycles (fdm12a units)
end

% Acoustic-mode wavenumber at node k: smallest NON-trivial sqrt(eig(Zs*Yv)), branch Im>=0.
function kap=acoustic_kappa(A,bw,Y,k,s,rho)
Zs=s*rho*diag(1./A(k,:)); Yv=bw(k)*squeeze(Y(k,:,:));
[~,Dg]=eig(Zs*Yv); ka=sqrt(diag(Dg)); fl=imag(ka)<0; ka(fl)=-ka(fl);
kab=abs(ka); cand=find(kab>1e-3*max(kab)); [~,j]=min(kab(cand)); kap=ka(cand(j));
end

% Full-array forward/backward decomposition at one frequency, to LOCALIZE where the
% backward wave originates (basal boundary artifact vs distributed intracochlear
% reflection). Pd+-(x) = (Pd -+ Pd'/kappa_a)/2 with kappa_a(x) the tracked acoustic mode.
function R=xprofile(pr)
if (isfield(pr,'pa')), pa=pr.pa; elseif (isfield(pr,'nch')), pa=modpar26(pr.nch); else, pa=modpar26(1); end
if (isfield(pr,'f')), f=pr.f; else, f=2000; end
rho=pa.rho; n=pa.n; s=2i*pi*f;
[a,q,x,y,D,Y,A,bw]=mxfill(f,pa); p=mxsolve(a,q); Pd=(D(1,:)*p).'; dx=x(2)-x(1);
m=size(A,2); kap=zeros(n,m); Wall=zeros(n,m,m);
for k=1:n
    Zs=s*rho*diag(1./A(k,:)); Yv=bw(k)*squeeze(Y(k,:,:));
    [Wk,Dg]=eig(Zs*Yv); ka=sqrt(diag(Dg)); fl=imag(ka)<0; ka(fl)=-ka(fl);
    kap(k,:)=ka.'; Wall(k,:,:)=Wk;
end
for k=2:n                                        % continuity-track
    pV=squeeze(Wall(k-1,:,:)); cV=squeeze(Wall(k,:,:)); ov=abs(pV'*cV);
    perm=zeros(1,m); used=false(1,m);
    for i=1:m, [~,ord]=sort(ov(i,:),'descend');
        for c=ord, if(~used(c)), perm(i)=c; used(c)=true; break; end, end, end
    kap(k,:)=kap(k,perm); Wall(k,:,:)=cV(:,perm);
end
kb=abs(kap(3,:)); cand=find(kb>1e-3*max(kb)); [~,j]=min(kb(cand)); ia=cand(j);
ka=kap(:,ia);
dPd=gradient(Pd)/dx;
Pp=(Pd - dPd./ka)/2; Pm=(Pd + dPd./ka)/2;
R.x=x; R.Pd=Pd; R.Pp=Pp; R.Pm=Pm; R.ka=ka; R.f=f; R.pa=pa; R.ratio=Pm./Pp;
end

% Expose the model's actual frequency-domain field (fdmod23) for wave diagnostics.
function R=fd_solve(pr)
if (isfield(pr,'pa')), pa=pr.pa; elseif (isfield(pr,'nch')), pa=modpar26(pr.nch); else, pa=modpar26(3); end
if (isfield(pr,'flst')), flst=pr.flst; else, flst=2000; end
[xx,Yb,Pd,Db,Dh]=fdmod23(pa,flst);
R.x=xx(:); R.Yb=Yb; R.Pd=Pd; R.Db=Db; R.Dh=Dh; R.pa=pa; R.f=flst;
end

% ---- Multi-chamber modal reflectance (the 3-chamber SFOAE, Stage 2) -------------
% Generalizes the scalar R = INT eps*exp(-2 INT kappa) dx to the coupled system by
% projecting onto the ACOUSTIC mode. Using the right/left eigenvectors w_a,u_a of
% M=Zs*Yv (u_a^T w_a = 1), the acoustic characteristic impedance is the SCALING-
% INVARIANT scalar  z0_a = kappa_a / (u_a^T*Yv*w_a)  (reduces to sqrt(Zs/Yv) for m=1),
% and eps_a = 0.5 d/dx ln z0_a is the local reflection coefficient. The tip-tail
% transition of the 3-chamber model makes z0_a vary rapidly near the CF place, so
% eps_a is large THERE (CF-place reflection, ms round-trip) -- the MOH2027 emission.
function R=modal_reflectance(pr)
if (isfield(pr,'pa')), pa=pr.pa; elseif (isfield(pr,'nch')), pa=modpar26(pr.nch); else, pa=modpar26(3); end
if (isfield(pr,'flst')), flst=pr.flst; else, flst=linspace(500,4000,600); end
rho=pa.rho; n=pa.n; nf=numel(flst); refl=zeros(1,nf);
[~,kref]=min(abs(flst-median(flst)));
for kf=1:nf
    f=flst(kf); s=2i*pi*f;
    [~,~,x,~,~,Y,A,bw]=mxfill(f,pa); m=size(A,2); dx=x(2)-x(1); ii=2:(n-1);
    kap=zeros(n,m); z0=zeros(n,m); Wall=zeros(n,m,m);
    for k=1:n
        Zs=s*rho*diag(1./A(k,:)); Yv=bw(k)*squeeze(Y(k,:,:));
        Mk=Zs*Yv; [Wk,Dg]=eig(Mk); lam=diag(Dg);
        ka=sqrt(lam); fl=imag(ka)<0; ka(fl)=-ka(fl);     % forward branch Im>=0
        Uk=Wk\eye(m);                                    % left eigenvectors (rows), u^T w = 1
        yv=sum((Uk*Yv).*Wk.',2);                         % u_j^T Yv w_j  (scaling-invariant)
        z0(k,:)=(ka./yv).';  kap(k,:)=ka.'; Wall(k,:,:)=Wk;
    end
    % continuity-track modes apically by eigenvector overlap
    for k=2:n
        pV=squeeze(Wall(k-1,:,:)); cV=squeeze(Wall(k,:,:));
        ov=abs(pV'*cV); perm=zeros(1,m); used=false(1,m);
        for i=1:m, [~,ord]=sort(ov(i,:),'descend');
            for c=ord, if (~used(c)), perm(i)=c; used(c)=true; break; end, end
        end
        kap(k,:)=kap(k,perm); z0(k,:)=z0(k,perm); Wall(k,:,:)=cV(:,perm);
    end
    kb=abs(kap(3,:)); cand=find(kb>1e-3*max(kb)); [~,j]=min(kb(cand)); ia=cand(j);
    ka=kap(:,ia); za=z0(:,ia);
    eps_=0.5*gradient(log(za))/dx; cumk=cumsum(ka)*dx;
    refl(kf)=sum(eps_(ii).*exp(-2*cumk(ii)))*dx;
    if (kf==kref), R.kap=ka; R.z0=za; R.eps=eps_; R.x=x; R.fref=f; R.cumk=cumk; end
end
w=2*pi*flst(:)';
R.f=flst; R.refl=refl; R.pa=pa;
R.mag=20*log10(abs(refl)+eps); R.ripple=R.mag-movmean(R.mag,50);
R.gd=-gradient(unwrap(angle(refl)))./gradient(w)*1000;
end

% Data-fit-only tuning error (MAP peak, phase, CF-place map, notches, Zc), the
% data terms of refdev() WITHOUT the exponent regularization -- so it measures
% "does the tuning still fit the data", independent of how the params got there.
function err=tuning_err(pa,flst,pk)
ipk=1+round(pk(1,:)*(pa.n-1)); Dt=pk(2,:); Pt=pk(3,:); Zt=pk(4,:);
[~,~,Pd,~,Dh,Zc]=fdmod23(pa,flst);
D1=20*log10(abs(Dh)); [~,nf]=size(D1);
notch1=max(diff(log(abs(Pd)),2)); notch2=max(diff(log(abs(Dh)),2));
ff=20*log10(flst/250)-pa.hbt; err=0;
for k=1:nf
    ii=D1(:,k)>-950;
    if (sum(ii)==0), err=1e6; return; end
    [~,ix]=max(D1(ii,k));
    ph=unwrap(angle(Dh(ii,k)))/(2*pi);
    mmx=D1(ix,k)+ff(k); pmx=ph(ix);
    dph=diff(ph); dpp=dph(dph>0);
    if (~isempty(dpp)), err=err+mean(dpp)*4e4; end   % phase-slope (causality) penalty
    err=err+notch1(k)*80+notch2(k)*80;               % tuning-shape notches
    err=err+abs(mmx-Dt(k))*2+abs(pmx-Pt(k))*2+abs(ix-ipk(k))*8;  % MAP peak, phase, fp-map
    err=err+abs(log(Zc(k)/Zt(k)))*8;                 % cochlear impedance
end
err=err/nf;
end

%--------------------------------------------------------------

function [a,q,x,y,D,Y,A,bw]=mxfill(f,pa)
% local parameter values
chsz=pa.chsz;
if (isfield(pa,'chsznorm') && pa.chsznorm), chsz=chsz*(2/sum(chsz)); end % see tdm26 note
n=pa.n;
xl=pa.xl; rho=pa.rho; g=pa.gam;
% useful constructs
s=2i*pi*f;
dx=xl/(n-1);
srd=s*rho*dx;
% compute admittance for all x
x=transpose(linspace(0,xl,n));
[z1,z2,zh,za,ac,gh,bw,z5,zv,z5c]=imped(x,s,pa); % z5 = OC-height, zv = CL vent shunt
m=length(pa.chsz);
A=zeros(n,m);
Y=zeros(n,m,m);
a=zeros(3,m,m,n);
y=zeros(n,m,m);
% initialize apex
Zh=pa.khe/s+pa.rhe+pa.mhe*s;
z1(n)=Zh;z2(n)=1e-6;zh(n)=1e-6;za(n)=0;
if (~isempty(z5)), z5(n)=1e-6; end   % apical terminal for the OC-height DOF,
                                     % matching the z2/zh treatment (m=4 only)
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
elseif (m==4)
    % 4-chamber [ST SS SV CL], 3 DOFs, matching tdm26's topology:
    %   d1 (BM)        ST <-> SV
    %   d2 (shear)     SV <-> SS
    %   d3 (OC height) CL <-> SS      <- the secondary SS/CL partition (RL)
    % Last row is the volume-conservation slot, zeroed as in m==3.
    D=[1 0 -1 0; 0 -1 1 0; 0 -1 0 1; 0 0 0 0];   % select PD
    % NESTED topology: d1 moves from ST<->SV to ST<->CL, so CL sits BETWEEN the
    % BM and the rest rather than beside them. Only row 1 changes; d2 and d3 are
    % untouched. This must match tdm26's fold_p pickup (P_ST - P_CL) or the two
    % solvers silently describe different models and maperr is meaningless.
    if (isfield(pa,'nested') && pa.nested)
        D(1,:)=[1 0 0 -1];          % d1: ST <-> CL
    end
    % Basal BC: the stapes drives ST and SV only; SS and CL have no direct
    % stapes path (same treatment as tdm26.m:424-425).
    B=[-1;0;1;0];                   % basal BC
end
for k=1:n
    A(k,:)=ac(k)*chsz;
    if (m==1)
        hh=z2(k)/(z2(k)+zh(k));
        zk(1)=z1(k)+(gh(k)*zh(k)-zg(k))*hh;
    elseif (m==2)
        hh=z2(k)/(z2(k)+zh(k));
        zk(1,1)=z1(k)+(gh(k)*zh(k)-zg(k))*hh;
        % THIRD DOF at m<3, mirroring micro26 (SN: "when m<3 the micro-mechanics
        % should mimic m=4 via interaction between SS and CL without longitudinal
        % coupling"). Without this fdm26 has no z5 below m=3, so the distillation
        % objective was structurally BLIND to k5/r5/m5 -- measured as exactly
        % 0.000e+00 sensitivity, which would have let a fit run happily while the
        % third DOF sat inert.
        %
        % NOTE THIS IS THE m==2 BRANCH, and it serves pa.m=1 as well: mxfill sets
        % m=length(pa.chsz) (fdm26.m:895, and fdm24.m:552 before it), and
        % modpar26(1) ships chsz=[1 1]. The m==1 branch above requires
        % numel(chsz)==1 and is DEAD for every standard parameter set. That is
        % also why m=1 == m=2 holds exactly -- it is the same model twice.
        %
        % zk3 is the m==4 3-DOF form. It is built separately rather than into zk
        % because zk(m,:)=1 below overwrites row m for volume conservation, which
        % at m=2 is row 2 -- one of the rows the condensation needs.
        use3 = isfield(pa,'dof') && pa.dof>=3 && ~isempty(z5);
        if (use3)
            zact3 = gh(k)*zh(k)-zg(k);
            zk3 = [ z1(k)+z5(k)  -zact3        -z5(k)
                   -z2(k)         z2(k)+zh(k)   0
                   -z5(k)         zact3         z5(k)];
            % ypart, NOT v3. This is the 3-element PARTITION RESPONSE VECTOR,
            % not the velocity of DOF 3 -- which is what "v3" reads as, and is
            % the same naming collision that made the algebraic shear (now dh in
            % micro26) confusable with the dynamical third DOF. Renamed
            % 2026-07-29. ypart(1)=BM, ypart(2)=hair bundle, ypart(3)=DOF 3.
            ypart = zk3 \ [1;0;0];   % partition response to unit BM-row pressure
        end
    elseif (m==3)
        % pa.m3form MIRRORS tdm26's switch of the same name (tdm26 m==3 branch).
        % Without this mirror the frequency domain SILENTLY IGNORED m3form, so a
        % model run with the m=2-style force law was still scored for maperr with
        % the d2-only law -- which is why nesting_test configs B and E returned
        % identical maperr (499.3) and the force law could not be tested at all.
        %   0 (default): FDM form, active force driven by d2 ALONE
        %   1          : m<3 form, active force driven by the RELATIVE
        %                displacement d3 = d1 - d2
        % Derived from tdm26's m<3 branch with zact = gh*zh - zg:
        %   s1 = -(z1*V1 + zact*(V1-V2))  -> row1 = [z1 + zact, -zact]
        %   s2 = -(z2*V2 - zh*(V1-V2))    -> row2 = [-zh,       z2+zh]
        % versus the d2-only form row1 = [z1, +zact], row2 = [-z2, z2+zh].
        % The recovered zact*V1 term is exactly what tdm26's m==3 comment calls
        % the "erroneous k_act*d1 stiffening" -- it legitimately belongs to the
        % m<3 law, and only its presence in the d2-only form was erroneous.
        % CONTINUOUS alpha in [0,1] (mirrors tdm26). Operator interpolation, so
        % both endpoints are reproduced EXACTLY:
        %   alpha=0 -> [z1, zact ; -z2, z2+zh]     (FDM d2-only, default)
        %   alpha=1 -> [z1+zact, -zact ; -zh, z2+zh] (m=2 relative-displacement)
        al = 0; if (isfield(pa,'m3form')), al = pa.m3form; end
        zact = gh(k)*zh(k)-zg(k);
        zmix = (1-al)*z2(k) + al*zh(k);
        % THIRD DOF AT m=3b, CONDENSED (2026-07-29). tdm26's m<4 third DOF is
        % INTERNAL -- no fluid compartment -- so it adds no pressure unknown and
        % can be eliminated analytically instead of adding a row. That matters
        % here because zk is m x m and row 3 is already the volume-conservation
        % row, leaving only two partition rows.
        %
        % Partition rows, with z5c = k5/s + r5 the coupling (spring and damper,
        % NO mass) and z5 = z5c + m5*s the full DOF-3 impedance:
        %   row1: (z1 + al*zact + z5c)*V1 + (1-2al)*zact*V2 - z5c*V3 = F1
        %   row3: -z5c*V1                                   + z5*V3  = 0
        % Row 3 has no active term: at m<4 the OHC force reaches the BM only
        % (micro26, and cochlea_proc.docx eq. 24). Solving row 3 for V3 gives
        % V3 = (z5c/z5)*V1, and substituting leaves the BM impedance raised by
        %   dz1 = z5c - z5c^2/z5 = z5c*(z5 - z5c)/z5 = z5c*(m5*s)/z5
        % the SERIES combination of the spring-damper with d3's mass, which is
        % what an attached mass presents as a load. Limits, all checked:
        %   m5 -> 0    dz1 -> 0      a massless d3 loads nothing
        %   m5 -> inf  dz1 -> z5c    d3 pinned in the lab, k5 a spring to ground
        %                            (exactly what breaks arch_gate's g2 freeze)
        %   z5 -> 0    dz1 has a POLE at d3's resonance: the abrupt change in
        %                            effective BM stiffness
        % Without this fdm26 had NO third DOF at m=3 at all, so g4_maperr_m3b
        % scored a 2-DOF model while tdm26 ran a 3-DOF one and reported PASS.
        dz1 = 0;
        if (isfield(pa,'d3int') && pa.d3int && ~isempty(z5) && ~isempty(z5c))
            dz1 = z5c(k) .* (z5(k) - z5c(k)) ./ z5(k);
        end
        zk(1,1:2)=[ z1(k)+al*zact+dz1   (1-2*al)*zact];
        zk(2,1:2)=[-zmix                 z2(k)+zh(k) ];
    elseif (m==4)
        % Derived from tdm26's m>=4 force_cp. NOTE THE SIGN FLIP vs m==3: the
        % m==3 BM equation carries +k_act*d2, whereas m>=4 carries -act (the
        % OHC pair's REACTION on the BM), with +act on the OC-height row.
        %   tdm26  s1 = -(k1 d1 + r1 v1 - act)          -> BM   : -zact*V2
        %          s2 = -(-k2 d1 - r2 v1 + (k2+k3) d2 + (r2+r3) v2)  (as m==3)
        %          s3 = -(k5 dc + r5 vc + act)          -> OCht : +zact*V2
        % with zact = gh*zh - zg the same active impedance as m==3.
        % LIMITATION: fdm26 has no ohcgain/ohcsgn, so this is the tdm26 default
        % (sgn=+1, fsp=1). A non-default pair will NOT be reproduced here.
        zact = gh(k)*zh(k)-zg(k);
        % m=4 IS UNCHANGED FROM HEAD, DELIBERATELY (2026-07-29, SN: "split the
        % change: keep the m<4 half, revert m=4"; "revert m=4 in fdm26 too").
        %
        % The relative-attachment form -- zk(1,1:3)=[z1+z5 -zact -z5],
        % zk(3,1:3)=[-z5 zact z5] -- was applied here to mirror micro26, and the
        % matching tdm26 change made m=4 DIVERGE in the time domain. Reverting
        % tdm26 alone left fdm26 carrying it, so g4_maperr_m4 stayed at 601.5
        % instead of returning to 525.2 -- the two solvers silently disagreeing,
        % which is the failure mode arch_gate's new stability gate exists to stop.
        %
        % NOTE what the old form costs, so it is not mistaken for correct: the
        % third column is [0;0;z5], so rows 1-2 close WITHOUT V3 and d3 is a
        % driven observer with no back-action. That is the frequency-domain twin
        % of the tdm26 inertness, and it is why k5/r5/m5 measure EXACTLY
        % 0.000e+00 in any objective built on this response. m=4 is therefore
        % reverted to a KNOWN-LIMITED form, not to a correct one; fixing it needs
        % the fluid-coupled d3 question settled first (see micro26's m>=4 note).
        zk(1,1:3)=[ z1(k)  -zact        0    ];
        zk(2,1:3)=[-z2(k)   z2(k)+zh(k) 0    ];
        zk(3,1:3)=[ 0       zact        z5(k)];
    end
    if (~exist('use3','var')), use3 = false; end   % only the m==2 branch sets it
    zk(m,:)=1; % conserve fluid volume
    if (m==3)
        T = [1, 0; 1, -1; -2, 1];
        Y(k,:,:) = T * (zk(1:2, 1:2) \ D(1:2, :));
    elseif (m==4)
        % T maps DOF velocities to chamber volume flows; columns MUST sum to
        % zero (volume conservation), as they do for m==3. The ST/SS/SV rows for
        % V1,V2 are unchanged from m==3 (same partitions, same convention); CL
        % couples only to V3, and V3 exchanges between SS and CL, which fixes
        % the third column up to ONE SIGN. pa.T4sgn selects it; the default is
        % validated against tdm26 (see mxfill4_check.m). Do not change without
        % re-running that check.
        % T4sgn = -1 VALIDATED (mxfill4_final.m, 2026-07-23): against tdm26's
        % click CF map, -1 matches to 0.010 in x/L at 4 kHz where +1 errs by
        % 0.112. Calibrated on m=2, where fdm26 and tdm26 agree to 0.000-0.006
        % across 250 Hz-8 kHz. Do not change without re-running that check.
        s4=-1; if (isfield(pa,'T4sgn')), s4=pa.T4sgn; end
        T = [ 1,  0,   0
              1, -1,  s4
             -2,  1,   0
              0,  0, -s4];
        Y(k,:,:) = T * (zk(1:3, 1:3) \ D(1:3, :));
    else
        Y(k,:,:) = zk \ D;
    end
    if (m==1)
        y(k,1,1)=1/zk(1,1);   % select Ybm
        y(k,2,:)=hh*y(k,1,:); % select Yhb
    elseif (m==2)
        if (use3)
            % 3-DOF condensation: BM and hair-bundle admittances come from the
            % full partition rather than the 2-DOF hh shortcut. Reduces to the
            % 2-DOF form when the third DOF is absent, which arch_gate checks.
            y(k,1,1)=ypart(1);   % select Ybm
            y(k,2,1)=ypart(2);   % select Yhb
        else
            y(k,1,1)=1/zk(1,1);   % select Ybm
            y(k,2,:)=hh*y(k,1,:); % select Yhb
        end
    elseif (m==3)
        y(k,:,:)=inv(zk);     % select Ybm & Yhb
    elseif (m==4)
        y(k,:,:)=inv(zk);     % select Ybm & Yhb & Yoc (same as m==3)
    end
    % CL VENT as a direct chamber-to-chamber shunt. It is NOT a partition, so it
    % does not belong in D; it is an admittance Yv=1/zv between CL (chamber 4)
    % and pa.clvtgt, stamped onto Y with the standard four-entry conductance
    % pattern -- the exact frequency-domain counterpart of the a2 stamp tdm26
    % writes at entries (t,t),(4,4),(t,4),(4,t).
    % REQUIRES pa.nested, matching tdm26, whose vent stamp lives INSIDE its
    % nested branch and is therefore ignored when nested=0. The port originally
    % applied the vent unconditionally, and fdm_port_gate caught it: the
    % appended+vent row moved maperr 1020.3 -> 1065.2 while amp and maxRe stayed
    % bit-identical to appended-without-vent, because tdm26 was not applying a
    % vent at all. Two solvers describing different models is exactly what makes
    % maperr meaningless, which is the fault this port existed to remove.
    if (m==4 && ~isempty(zv) && isfield(pa,'clvent') && pa.clvent>0 && ...
        isfield(pa,'nested') && pa.nested)
        vt=3; if (isfield(pa,'clvtgt')), vt=pa.clvtgt; end
        Yv=1/zv(k);
        Yk=squeeze(Y(k,:,:));
        Yk(4,4)=Yk(4,4)+Yv;   Yk(vt,vt)=Yk(vt,vt)+Yv;
        Yk(4,vt)=Yk(4,vt)-Yv; Yk(vt,4)=Yk(vt,4)-Yv;
        Y(k,:,:)=Yk;
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
if (isfield(pa,'alfx_scale')), alfx=alfx*pa.alfx_scale; else, alfx=alfx*1e-7; end
a(2,:,:,1)=A1*(1+alfx);
a(3,:,:,1)=A2*( -alfx);
% initialize q
q=zeros(m,n);
q(:,1)=qst*B;
end % return

function [z1,z2,z3,z4,ac,gh,bw,z5,zv,z5c]=imped(x,s,pa)
% z5 (8th output, 4-chamber only) is the OC-HEIGHT / cortilymph-pump impedance,
% the frequency-domain counterpart of tdm26's cp.k5/cp.r5/cp.m5. Appended LAST so
% the two existing 7-output callers (fdm26.m:617, :864) are unaffected.
x=x.*(1+(pa.xtap*x).^pa.xtex);
q=x.^2;
k1=pa.k1o*exp(pa.k1e*x+pa.k1q*q);
if isfield(pa, 'rough_amp')
    sd=42; if (isfield(pa,'rough_seed')), sd=pa.rough_seed; end
    rng(sd); % reproducible roughness (seed selectable for realization statistics)
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
% OC-height (cortilymph pump) impedance -- 4-chamber only. Mirrors tdm26.m:475-477
% cp.k5/cp.r5/cp.m5, and carries the MASS term (m5*s) as z1 does, because the
% OC-height DOF has its own inertia (tdm26.m:630 a(i3)=a(i1)+s3./cp.m5).
z5=[]; z5c=[];
if (isfield(pa,'k5o'))
    k5=pa.k5o*exp(pa.k5e*x+pa.k5q*q);
    r5=pa.r5o*exp(pa.r5e*x+pa.r5q*q);
    m5=pa.m5o*exp(pa.m5e*x+pa.m5q*q);
    z5=k5/s+r5+m5*s;
    % z5c is the COUPLING impedance: spring and damper only, NO mass. The
    % attachment transmits k5/r5 between BM and d3, while m5 belongs to d3
    % alone. Needed to condense the internal third DOF at m<4 (see m==3).
    z5c=k5/s+r5;
end
% VENT impedance -- the CL shunt, mirroring tdm26's resonant vent. The channel
% inertance is mv = m5/clvent, which is NOT a free choice: tdm26's a2 stamp
% G = clvent*mu3 = clvent*m1/m5 is exactly the coefficient m1/mv of a DOF of
% mass m5/clvent, so the pure-inertance vent already had that mass implicitly.
% clvk (cortilymph compliance) and clvr (channel resistance) are the added
% elements that make the CL resonance placeable independently of the coupling
% strength. clvk=clvr=0 leaves zv=mv*s, the pure inertance, which reproduces the
% previous frequency-flat behaviour.
zv=[];
if (isfield(pa,'clvent') && pa.clvent>0 && ~isempty(z5))
    mv=m5/pa.clvent;
    kv=0; if (isfield(pa,'clvk')), kv=pa.clvk; end
    rv=0; if (isfield(pa,'clvr')), rv=pa.clvr; end
    % pa.clvoct places the CL resonance a fixed number of octaves BELOW the
    % local BM resonance sqrt(k1/m1), and takes precedence over a raw clvk. It
    % has to be resolved HERE rather than in modpar26 because the required
    % stiffness is place-dependent (it tracks k1/m1/m5), so a scalar set in the
    % parameter file would be correct at exactly one place. Resolving it from
    % the solver's own k1/m1 also keeps it right if those are ever refit.
    if (isfield(pa,'clvoct') && isfinite(pa.clvoct))
        kv=(m5/pa.clvent).*(k1./m1)*4^(-pa.clvoct);
    end
    kv=kv(:).*ones(size(m5)); rv=rv(:).*ones(size(m5));
    zv=kv/s+rv+mv*s;
end
ac=pa.aco*exp(pa.ace*x+pa.acq*q);
gh=pa.gpo*exp(pa.gpe*x+pa.gpq*q);
bw=pa.bwo*exp(pa.bwe*x+pa.bwq*q);
if (isfield(pa,'z0unif') && pa.z0unif)
    % Uniform-z0 area compensation (fdm12a-style): retaper the scala area ac and
    % BM width bw by the SAME factor (=> bw/ac unchanged => kappa=sqrt(Zs*Yv), hence
    % CF map / tuning / latency, preserved to leading order) so that z0 becomes
    % uniform in the stiffness tail, killing the spurious basal reflection.
    % Co-scaling ac,bw by exp(sz*x) scales z0 by exp(-sz*x); sz is the z0 tail
    % log-slope. pa.z0slope = MEASURED acoustic-mode z0_a slope (correct for the
    % 3-chamber's 2-DOF coupling); default falls back to the k1-scalar prediction
    % sz=(k1e-ace-bwe)/... i.e. force ac*bw ~ k1 (exact only when Yb ~ s/k1).
    if (isfield(pa,'z0slope'))
        sz=pa.z0slope;  szq=0; if (isfield(pa,'z0slopeq')), szq=pa.z0slopeq; end
        ace=pa.ace+sz; bwe=pa.bwe+sz; acq=pa.acq+szq; bwq=pa.bwq+szq;
    else
        d =pa.bwe-pa.ace; ace=(pa.k1e-d)/2; bwe=(pa.k1e+d)/2;
        dq=pa.bwq-pa.acq; acq=(pa.k1q-dq)/2; bwq=(pa.k1q+dq)/2;
    end
    ac=pa.aco*exp(ace*x+acq*q);
    bw=pa.bwo*exp(bwe*x+bwq*q);
end
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

%% modpar26, par_CEL16, modpar26c3 moved to shared files modpar26.m / modpar26c3.m

function R = macroD_export(pr)
%MACROD_EXPORT  Expose this solver's D and B for the fdm/tdm parity check.
%
% The capstone design note requires "fdm/tdm parity checked at the partition
% boundary". That check needs fdm26's coupling matrix, which is built inside
% mxfill and is therefore unreachable from outside. This dispatch rebuilds it
% from the same definitions and returns it, touching no solve.
%
% IT EXPORTS THE TOPOLOGY, which is the invariant that must never diverge:
% which chamber pair each DOF spans, including the nested variant. It does NOT
% claim the matrices should be numerically EQUAL to tdm26's. fdm26 poses
% Y = T*(zk\D) in an impedance formulation; tdm26 poses s = s_int - Df*p with
% the mass ratios applied separately. So the scalings legitimately differ --
% m=2 by a factor of two, the d3 row by clcouple, B by sign convention.
% Conflating "same topology" with "same matrix" would make the parity test
% either vacuous or permanently red, and neither is worth having.
pa = pr.pa; m = pa.m;
if (m==1)
    D=1;                      B=2;
elseif (m==2)
    D=[1 -1;0 0];             B=[-1;1];
elseif (m==3)
    D=[1 0 -1;0 -1 1; 0 0 0]; B=[-1;0;1];
else
    D=[1 0 -1 0; 0 -1 1 0; 0 -1 0 1; 0 0 0 0];
    if (isfield(pa,'nested') && pa.nested), D(1,:)=[1 0 0 -1]; end
    B=[-1;0;1;0];
end
R.D = D; R.B = B; R.m = m; R.topo = sign(D);
end
