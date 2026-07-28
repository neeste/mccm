function sone_fit
pa.ldew=2; pa.ldne=70; pa.ldsc=[6.913 3.823 0.454 4.499]; % 8.2%
slv=[0 20 40 60 80 100];
ldn=zeros(size(slv));
lcu=zeros(size(slv));
lph=zeros(size(slv));
ew=pa.ldew;
ne=pa.ldne;
sc=pa.ldsc;
mbw=0; % model bandwidth
nlv=length(slv);
nrs=cell(1,nlv);
for k=1:nlv
    fn = sprintf('bw%03d/ldngrw%03d.mat',mbw*100,round(slv(k)));
    load(fn,'nr','fr','bw');
    nt=size(nr,1);
    i1=1+round((80/120)*(nt-1)); % ramp=70, steady=50
    nr=nr(i1:nt,:)';
    nrs{k}=nr;
end
if (1)
    op=optimset('MaxFunEvals',1000);
    sc=fminsearch(@(sc) ld_err(sc,ew,ne,nrs),sc,op);
    fprintf('sc=[%5.3f %5.3f %5.3f %5.3f];\n',sc) 
end
ewn=ew/35; % normalize by xl=35 mm
fprintf('  fr   lv   bw   ldph\n');
for k=1:length(slv)
    lv=slv(k);
    nr=nrs{k};
    pld = loudness(nr,ewn,ne,sc);
    ldso = mean(pld); % sone
    ldcu = sone2cu(ldso);
    ldph = sone2phon(ldso);
    ldn(k) = ldso;
    lcu(k) = ldcu;
    lph(k) = ldph;
    fprintf('%5.2f %3.0f %5.2f %5.2f\n',fr,lv,bw,ldph);
end
% CLS data - Rasetshwane et al. 2015
lvd=[0.17 20.14 40.10 60.07 80.03 100.00];
ldd=[0.96  6.64 11.69 18.00 25.19  47.15];
ldp=cu2phon(ldd);
% Fletcher & Munson
[son,bet]=fm33; 
% tdm24 comparison
figure(1); clf
subplot(1,2,1)
semilogy(slv,ldn,bet,son,'k:')
xlabel('level (dB SPL)')
ylabel('loudness (sone)')
legend('tdm24','Fletcher','Location','SouthEast')
axis([-5 105 0.001 100])
subplot(1,2,2)
plot(slv,lph,lvd,ldp,'k:')
xlabel('level (dB SPL)')
ylabel('loudness (phon)')
axis([-5 105 -5 105])
legend('tdm24','Rasetshwane','Location','NorthWest')
drawnow
end

function err=ld_err(sc,ew,ne,nrs)
lv=[0 20 40 60 80 100];
son=fm33;
ew=ew/35; % normalize by xl=35 mm
err=0;
for k=1:length(lv)
    kk=lv(k)+11;
    fld=son(kk);
    nr=nrs{k};
    pld=loudness(nr,ew,ne,sc);
    ld=mean(pld);
    err=err+abs(log10(ld/fld));
end
fprintf('sc=[%5.3f %5.3f %5.3f %5.3f]; %% err=%5.3f\n',sc,err) 
end

function [pld,enr]=loudness(nr,ewn,ne,c)
nx=size(nr,1);
kw=round(nx*ewn)-1;
enr=zeros(ne,1);
p=1/3; % log-loudness slope
for k=1:ne
    k1=1+round((k-1)*(nx-1)/(ne-1));
    k2=min(k1+kw,nx);
    kk=k1:k2;
    enr(k)=mean(mean(abs(nr(kk,:)).^p));
end
N=enr.^(1/p);
pld=(N.^c(1)./(1+N.^c(2)/c(3)))*c(4);
end

function phon=cu2phon(CU)
% reference phon for NH1k listeners from Rasetshwane et al. (2015)
NH1k=[0 18.309 33.668 48.953 63.718 75.942 82.27 87.093 91.878 96.643 102];
CU(CU<00)=00;
CU(CU>50)=50;
phon=zeros(size(CU));
for k=1:length(CU)
    phon(k)=interp1(0:5:50,NH1k,CU(k));
end
end

% sone-to-CU conversion (Jesteadt et al. 2017)
function cu=sone2cu(sone)
cu=2.6253*log10(sone+0.0887).^3+0.7799*log10(sone+0.0887).^2+8.0856*log10(sone+0.0887)+13.4493;
end

% sone-to-phon conversion
function phon=sone2phon(sone)
[son,bet]=fm33;
phon=zeros(size(sone));
for k=1:length(sone)
    phon(k)=interp1(son,bet,sone(k));
end
end % return

% Fletcher & Muson 1933 (sone)
function [son,bet]=fm33
bet=-10:129;
flu=[0.015     0.025      0.04      0.06      0.09      0.14      0.22      0.32 ...
      0.45       0.7         1       1.4       1.9      2.51       3.4      4.43 ...
       5.7      7.08         9      11.2      13.9      17.2      21.4      26.6 ...
      32.6      39.3      47.5      57.5      69.5      82.5      97.5       113 ...
       131       151       173       197       222       252       287       324 ...
       360       405       455       505       555       615       675       740 ...
       810       890       975      1060      1155      1250      1360      1500 ...
      1640      1780      1920      2070      2200      2350      2510      2680 ...
      2880      3080      3310      3560      3820      4070      4350      4640 ...
      4950      5250      5560      5870      6240      6620      7020      7440 ...
      7950      8510      9130      9850  1.06e+04  1.14e+04  1.24e+04  1.35e+04 ...
  1.46e+04  1.58e+04  1.71e+04  1.84e+04  1.98e+04  2.14e+04  2.31e+04   2.5e+04 ...
  2.72e+04  2.96e+04  3.22e+04   3.5e+04   3.8e+04  4.15e+04   4.5e+04   4.9e+04 ...
   5.3e+04   5.7e+04   6.2e+04  6.75e+04   7.4e+04   8.1e+04   8.8e+04   9.7e+04 ...
  1.06e+05  1.16e+05  1.26e+05  1.38e+05   1.5e+05  1.64e+05   1.8e+05  1.97e+05 ...
  2.15e+05  2.35e+05   2.6e+05  2.88e+05  3.16e+05  3.46e+05   3.8e+05  4.18e+05 ...
   4.6e+05  5.06e+05  5.56e+05  6.09e+05  6.68e+05  7.32e+05     8e+05  8.75e+05 ...
  9.56e+05 1.047e+06  1.15e+06 1.266e+06];
son=flu/flu(51); % convert FLU to sone
end % return
