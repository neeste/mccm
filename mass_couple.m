% MASS-COUPLING LEVER: does lightening the BM load decompress the m=4 CF map?
%
% Hypothesis: the extra shear (d2) and OC-height (d3) DOFs load the BM through
% the mass couplings mu2 = cp.m1./cp.m2 and mu3 = cp.m1./cp.m5 (tdm26.m:406-407),
% compressing the tonotopic range to ~3.4 oct (vs 5.9 for m=2). Raising m2o
% (lowers mu2) and/or m5o (lowers mu3) lightens the load.
%
% MECHANISM PRIOR: the 2-chamber ALSO has m2/shear and spans a full 5.9 oct, so
% mu2 alone cannot be the compressor -- the m=4-specific element is the OC-height
% DOF, i.e. mu3/m5. Both are tested; mu3 is expected to dominate.
%
% PRE-FLIGHT (standing rule): confirm the lever moves the OUTPUT (map RANGE in
% octaves), not merely the intermediate mu2/mu3 -- an inert or cancelled lever
% shows as an unchanged range, as chsz(4)=0 did for m5.
%
% Primary metric: BF-map RANGE from the linear click (fast). Stability (coupeig,
% ~4 min/eval for m=4) is checked only on the best config(s) afterwards.

ISV=[1136 1005 840 655 466 273 80];
b2=modpar26(4).m2o; b5=modpar26(4).m5o;
fprintf('native m2o=%.4g  m5o=%.4g   (mult scales these)\n', b2, b5);

rangeof = @(pa) local_range(pa, ISV);

% ---- PRE-FLIGHT: does each lever move the RANGE at all? ----
fprintf('\n=== PRE-FLIGHT: range at mult 1 vs 8 ===\n');
p2a=modpar26(4); p2b=modpar26(4); p2b.m2o=b2*8;
p5a=modpar26(4); p5b=modpar26(4); p5b.m5o=b5*8;
[r2a,~,~]=rangeof(p2a); [r2b,~,~]=rangeof(p2b);
[r5a,~,~]=rangeof(p5a); [r5b,~,~]=rangeof(p5b);
fprintf('  m2o x1 -> x8 : range %.2f -> %.2f oct  (delta %.2f)\n', r2a, r2b, abs(r2a-r2b));
fprintf('  m5o x1 -> x8 : range %.2f -> %.2f oct  (delta %.2f)\n', r5a, r5b, abs(r5a-r5b));
if (abs(r2a-r2b)<0.01 && abs(r5a-r5b)<0.01)
    fprintf('  *** ABORT: neither lever moves the range. ***\n'); disp('MASS_COUPLE_ABORTED'); return
end
fprintf('  -> at least one lever moves the range. Proceeding.\n');

mult=[1 2 4 8 16 32];
fprintf('\n=== BF-map RANGE (oct) vs mass multiplier  [target m=2 = 5.90] ===\n');
fprintf('  mult |  m2o only        |  m5o only        |  BOTH\n');
fprintf('       | range  mono  fin | range  mono  fin | range  mono  fin\n');
for mm=mult
    pa2=modpar26(4); pa2.m2o=b2*mm;
    pa5=modpar26(4); pa5.m5o=b5*mm;
    pab=modpar26(4); pab.m2o=b2*mm; pab.m5o=b5*mm;
    [r2,m2s,f2]=rangeof(pa2); [r5,m5s,f5]=rangeof(pa5); [rb,mbs,fb]=rangeof(pab);
    fprintf('  %4g | %5.2f  %-4s %-3s | %5.2f  %-4s %-3s | %5.2f  %-4s %-3s\n', ...
        mm, r2,m2s,f2, r5,m5s,f5, rb,mbs,fb);
end
fprintf(['\nDECOMPRESSION if range rises toward 5.90 as mass rises. If mu3/m5\n' ...
         'dominates (expected), the m5o and BOTH columns move while m2o stays flat.\n' ...
         'If NEITHER rises, the compression is not BM mass-loading by these DOFs.\n']);
disp('MASS_COUPLE_DONE');

function [rng,mono,finl]=local_range(pa, ISV)
rng=NaN; mono='?'; finl='yes'; pa.isv=ISV;
try, evalc('S=tdm26(0,pa,0,0);'); catch, finl='THRW'; return; end
if (any(~isfinite(S.d1(:)))), finl='DIVG'; return; end
nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf); bf=nan(1,numel(ISV));
for i=1:numel(ISV)
    D=fft(S.d1(:,i)); H=D(1:nf)./max(abs(P),eps); a=abs(H).*f; a(~(f>0.15&f<18))=0;
    [pk,ip]=max(a); if(pk>0), bf(i)=f(ip); end
end
v=bf(isfinite(bf)&bf>0);
if (numel(v)>=2), rng=log2(max(v)/min(v)); end
bfx=bf(isfinite(bf));
if (numel(bfx)>=3), mono=tern(all(diff(bfx)>0),'ok','FOLD'); end
end
function s=tern(c,a,b), if c, s=a; else, s=b; end, end
