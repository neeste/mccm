% THE SS-VENT LIMIT, WITH THE CONNECTIVITY FIXED TOO.
%
% ss_limit FAILED and the reason was my error, not the model's. I merged the
% AREAS (chsz [.95 .10 .95]) but left m=3b's CONNECTIVITY, which puts SV in the
% middle of the chain. maxRe came out +20.0 against the target's +3.1.
%
% THE CHAIN. Verified from the tdm26 stamp: chambers are 1=ST, 2=SS, 3=SV, d1
% couples 1-3, d2 couples 2-3. So the chain is
%     ST(1) --d1-- SV(3) --d2-- SS(2)
% and chamber 3 is the MIDDLE. Merging CL into SS under the nested chain gives
%     ST --d1-- (SS+CL) --d2-- SV
% with the small merged pool in the middle instead. With three chambers and two
% partitions there is only one structural degree of freedom -- which chamber is
% in the middle -- and it is chosen entirely by which SLOT holds which area.
% So the nested 3-chamber topology needs no new code: put the merged pool in
% slot 3 and the real SV in slot 2.
%
%     chsz = [0.95  0.95  0.10]
%             ST    SV    SS+CL      <- slot 2 holds SV, slot 3 holds the merge
%
% WHY THIS IS WORTH A SECOND ATTEMPT. If it matches, the nested construction is
% a permuted 3-chamber model, fdm26 scores its map TODAY, and the clean-map half
% of SN's bar stops being blocked on porting nested/clvent into fdm26. Every
% nested row so far has reported maperr 1020.3, which is the appended model's
% number for a different topology and means nothing here.
%
% TARGET nested vent 100 -> SS: d1 +67.27  d2 +80.31  d3 +72.00, maxRe +3.1.
% Stock m=3b for contrast: +81.15 / +84.17 / +86.33, maxRe +19.3, maperr 499.3.
%
% The sweep around 0.10 is there because the effective merged area need not be
% exactly the arithmetic sum: the vent is finite (100, not infinite) and d1 also
% couples CL to ST, so some of CL's area may still act separately.

TGT = [67.27 80.31 72.00 3.1];

cfg = { 'permuted [.95 .95 .10]', [0.95 0.95 0.10]
        'permuted [.95 .95 .05]', [0.95 0.95 0.05]
        'permuted [.95 .95 .15]', [0.95 0.95 0.15]
        'permuted [.95 1.0 .10]', [0.95 1.00 0.10]
        'stock    [.95 .05 1.0]', [0.95 0.05 1.00] };

fprintf('\n  TARGET d1 %+.2f  d2 %+.2f  d3 %+.2f | maxRe %+.1f\n', TGT(1), TGT(2), TGT(3), TGT(4));
fprintf('  slot order is [ST SS SV]; d1 couples 1-3, d2 couples 2-3, so slot 3 is the MIDDLE\n\n');
fprintf('  config                 | amp d1  amp d2  amp d3 | maxRe     | maperr | range mono\n');
fprintf('%s\n', repmat('-',1,88));
for i = 1:size(cfg,1)
    pa = modpar26(3); pa.chsz = cfg{i,2};
    try
        S = score26(pa, 'fast', false);
    catch e
        fprintf('  %-22s | FAILED: %s\n', cfg{i,1}, e.message); continue
    end
    mk = '';
    da = abs(S.amp_gain-TGT(1)); dr = abs(S.maxRe-TGT(4));
    if (da < 3 && dr < 3), mk = ' <== MATCHES THE VENT LIMIT'; end
    fprintf('  %-22s | %+6.2f %+7.2f %+7.2f | %+9.1f | %6.1f | %5.2f %-4s%s\n', ...
        cfg{i,1}, S.amp_gain, S.amp_d2, S.amp_d3, S.maxRe, S.maperr, ...
        S.bf_range, S.bf_mono, mk);
end
fprintf(['\n  IF A ROW MATCHES: its maperr is the FIRST live measurement of the map\n' ...
         '  half of SN''s bar for this session''s construction. Read it against\n' ...
         '  m=3b 499.3 and appended 1020.3. Below 499.3 would mean the nested\n' ...
         '  construction improves the map as well as the amplifier.\n' ...
         '  IF NO ROW MATCHES: d3 keeps a fluid role that no 3-chamber model can\n' ...
         '  represent, and the fdm26 port is genuinely required.\n']);
disp('SS_LIMIT2_DONE');
