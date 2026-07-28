% Which transduction drive should the WNR latency be measured from?
%   default : cur.hb = v1 - a*v2   (mixed; free-standing IHC bundle in sub-tectorial flow)
%   bm      : cur.hb = v1          (BM velocity only -- what the re-baseline used)
%   dof2    : cur.hb = -v2         (TM-RL SHEAR only -- matches the OHC wiring d3=d2,
%                                   and the SS-area/shear account of bundle deflection)
% All three leave the MECHANICS identical (cur.hb feeds only the neural path), so this
% is purely a question of which quantity the Wave-V analogue should be read from.
L=load('parfit26_recip.mat'); pa0=L.R.pa;
modes={'default','bm','dof2'};
lbl={'v1-a*v2 (default)','v1 (BM only)','-v2 (shear only)'};
fprintf('\n%-20s %7s %8s %9s %7s %9s %9s\n','drive','slope','level_c','bandRMS','Delta','worstStep','anchor');
fprintf('%s\n',repmat('-',1,76));
G=cell(1,3);
for k=1:3
    pa=pa0;
    if (k==2), pa.hbmode='bm'; elseif (k==3), pa.hbmode='dof2'; end
    m=abr_metric(pa,false);
    if (~m.ok), fprintf('%-20s FAILED: %s\n',lbl{k},m.msg); continue; end
    f=m.f(:); slv=m.slv(:); [F,I]=ndgrid(f,slv/100);
    t88=12.90*(5.00.^(-I)).*(F.^(-0.413)); t13=12.63*(5.34.^(-I)).*(F.^(-0.390));
    tlo=min(t88,t13); thi=max(t88,t13);
    Lm=m.lat; ok=isfinite(Lm)&Lm>0;
    ds=linspace(0,1.5,301); best=inf; D=0;
    for j=1:numel(ds)
        V=max(0,tlo-(Lm+ds(j)))+max(0,(Lm+ds(j))-thi);
        r=sqrt(mean(V(ok).^2)); if(r<best),best=r; D=ds(j); end
    end
    st=0;
    for a=1:numel(f), for b=1:numel(slv)-1
        if (ok(a,b)&&ok(a,b+1)), st=max(st,abs(100*(Lm(a,b+1)-Lm(a,b))/Lm(a,b))); end
    end, end
    fi=find(abs(f-1)<0.01,1); li=find(abs(slv-20)<0.1,1); anc=NaN;
    if (~isempty(fi)&&~isempty(li)), anc=Lm(fi,li); end
    fprintf('%-20s %7.3f %8.2f %9.3f %7.2f %8.1f%% %9.2f\n', lbl{k}, m.slope, m.level_c, best, D, st, anc);
    G{k}=Lm;
end
fprintf('\ntargets: slope 0.39-0.41 | level_c 5.0-5.34 | bandRMS ->0 | steps ~28%%/20dB | anchor 9.03-9.35 ms\n');
for k=1:3
    if (isempty(G{k})), continue; end
    fprintf('\nlatency grid -- %s  (rows 0.5/1/2/4 kHz, cols 20/40/60/80 dB)\n',lbl{k});
    disp(round(G{k},2));
end
save('dof2_compare.mat','G','modes');
disp('DOF2_COMPARE_DONE');
