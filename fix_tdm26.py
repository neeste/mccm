import re

with open("tdm26.m", "r") as f:
    content = f.read()

# Revert xpnd_q
content = content.replace(
    "qx(j2) =  s1 + s2 .* mu;\n    qx(j3) = -2*s1 - s2 .* mu - qst;",
    "qx(j2) = -s2 .* mu;\n    qx(j3) = -s1 + s2 .* mu - qst;"
)

# Revert cp_imped a2 modification
content = content.replace(
    "a2(k,4) = 1; a2(k,5) = L2_p + L2_c - mu(k); a2(k,6) = mu(k) - 1;\n        % Chamber 3 (SV): Symmetric Coupling\n        a1(k,9) = -L3_p; a3(k,9) = -L3_c;\n        a2(k,7) = -2; a2(k,8) = mu(k); a2(k,9) = L3_p + L3_c + 2 - mu(k);",
    "a2(k,4) = 0; a2(k,5) = L2_p + L2_c + mu(k); a2(k,6) = -mu(k);\n        % Chamber 3 (SV): Symmetric Coupling\n        a1(k,9) = -L3_p; a3(k,9) = -L3_c;\n        a2(k,7) = -1; a2(k,8) = -mu(k); a2(k,9) = L3_p + L3_c + 1 + mu(k);"
)

# Revert cp_imped boundary conditions
content = content.replace(
    "a2(1,4) = 1;                    a2(1,5) = L2_c - mu(1); a2(1,6) = mu(1) - 1;\n    a2(1,7) = -2;                   a2(1,8) = mu(1);  a2(1,9) = L3_c + 2 - mu(1) + 2*cp.alfx;",
    "a2(1,4) = 0;                    a2(1,5) = L2_c + mu(1); a2(1,6) = -mu(1);\n    a2(1,7) = -1;                   a2(1,8) = -mu(1); a2(1,9) = L3_c + 1 + mu(1) + 2*cp.alfx;"
)

with open("tdm26.m", "w") as f:
    f.write(content)
