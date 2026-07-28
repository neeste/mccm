# Make score26's place indices n-relative. Apply ONLY when no job is using score26.
import re
p='score26.m'; s=open(p).read()
old="ISV = 1391:-10:11; XLO = 0.05; XHI = 0.85; FP = [1 2 4 8];"
new=("% Place grid as FRACTIONS of pa.n, not fixed indices. The original\n"
     "% 1391:-10:11 was valid only for n=1401; at any other grid size those\n"
     "% indices fall outside the array. Fractions keep the SAME physical places\n"
     "% under a change of n, which is what makes a convergence test meaningful.\n"
     "ISVFRAC = (1391:-10:11)/1401;\n"
     "ISV = unique(max(1, min(pa.n, round(ISVFRAC*pa.n))), 'stable');\n"
     "XLO = 0.05; XHI = 0.85; FP = [1 2 4 8];")
assert old in s, 'anchor not found (already patched?)'
s=s.replace(old,new)
open(p,'w').write(s)
print('score26: place indices are now n-relative')
