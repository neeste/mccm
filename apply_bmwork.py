# Add the BM-WORK diagnostic to tdm26. Apply ONLY when no job is using tdm26.
# ohcp measures work on the RELATIVE (BM-to-RL) coordinate. What decides whether
# the OHC force AMPLIFIES THE BM is the work it does against v1, which is a
# different quantity and is not currently measured.
p='tdm26.m'; s=open(p).read()
E=[
 ("cur.stp = 0; cur.cyc = 0; cur.wnr = 0; cur.ohcp = 0;",
  "cur.stp = 0; cur.cyc = 0; cur.wnr = 0; cur.ohcp = 0; cur.ohcbm = 0;"),
 ("sav.ohcp = zeros(nt,1);",
  "sav.ohcp = zeros(nt,1);\nsav.ohcbm = zeros(nt,1);   % OHC work against BM velocity"),
 ("[ss,ii,gam,ohcp]=force_cp(pa,cp,st);",
  "[ss,ii,gam,ohcp,ohcbm]=force_cp(pa,cp,st);"),
 ("cur.ohcp = ohcp;",
  "cur.ohcp = ohcp; cur.ohcbm = ohcbm;"),
 ("   sav.ohcp(i) = cur.ohcp;",
  "   sav.ohcp(i) = cur.ohcp; sav.ohcbm(i) = cur.ohcbm;"),
 ("function [ss,ii,gam,ohcp]=force_cp(pa,cp,st)",
  "function [ss,ii,gam,ohcp,ohcbm]=force_cp(pa,cp,st)"),
 ("n = pa.n; i1 = 1:n; i2 = (n+1):(2*n); gam = cp.gm; ohcp = 0;",
  "n = pa.n; i1 = 1:n; i2 = (n+1):(2*n); gam = cp.gm; ohcp = 0; ohcbm = 0;"),
 ("    ohcp = -sum(wpow(isfinite(wpow)));",
  """    ohcp = -sum(wpow(isfinite(wpow)));
    % BM WORK. ohcp above is the work on the RELATIVE (BM-to-RL) coordinate,
    % which is the right measure for the internal force pair. It does NOT say
    % whether the BM is being amplified. The BM receives +act (the reaction in
    % s1 below), so the power delivered to BM motion is act.*v1. ohcbm>0 means
    % the OHC force opposes BM damping and drives the travelling wave; ohcbm~0
    % with ohcp>0 means the force is injecting energy into the OC-height
    % coordinate only and never reaching the BM, which would explain an
    % amplifier that stays near 2.5 dB regardless of gain or damping.
    bpow  = act .* v1;
    ohcbm = sum(bpow(isfinite(bpow)));"""),
 ("dgn.ohcP=mean(sav.ohcp(isfinite(sav.ohcp)));                 % >0 = energy injected",
  """dgn.ohcP=mean(sav.ohcp(isfinite(sav.ohcp)));                 % >0 = energy injected
dgn.ohcBM=mean(sav.ohcbm(isfinite(sav.ohcbm)));              % >0 = BM amplified
dgn.ohcBMW=sum(sav.ohcbm(isfinite(sav.ohcbm)))*dtms/1000;    % integrated BM work"""),
]
for old,new in E:
    assert old in s, 'anchor missing: %.60s' % old
    s=s.replace(old,new,1)
open(p,'w').write(s)
print('tdm26: added ohcbm (BM-work) diagnostic, %d edits' % len(E))
