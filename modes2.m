% IS THERE A SECONDARY PROPAGATION MODE BETWEEN SS AND CL?
%
% SN: traditional cochlear models propagate the traveling wave by cross-section
% AREA EXCHANGE between SV and ST (BM as the partition). It would make sense if
% a SECONDARY propagation mode existed between SS and CL.
%
% The structural analog is already present: RL separates SS from CL just as BM
% separates SV from ST, the OC-height coordinate d3 is that partition, and
% chsz=[0.95 0.05 1.0 0.05] makes SS/CL a MATCHED pair (0.05/0.05) exactly as
% ST/SV are matched (0.95/1.0). Open question: does a pair at 20x smaller
% cross-section actually PROPAGATE, or is it evanescent (energy accumulating
% locally rather than travelling)? The latter would also explain the STATIC,
% f=0, BASAL instability found in the OC-height DOF.
%
% modal_wave_analysis (fdm26.m:645) gives modal wavenumbers kappa_i =
% sqrt(eig(Zs*Yv)) with forward branch Im>=0, plus (now) the mode SHAPES Vall in
% chamber-pressure space [ST SS SV CL].
%   PROPAGATING : |Im(kappa)| >> |Re(kappa)|  (phase accumulates, little decay)
%   EVANESCENT  : |Re(kappa)| dominant        (decays without travelling)
% Mode identity is read off the eigenvector: an ST-vs-SV differential is the
% primary acoustic wave; an SS-vs-CL differential is the secondary SN proposes;
% an all-same-sign shape is the trivial compression mode (kappa ~ 0).

lab3={'ST','SS','SV'}; lab4={'ST','SS','SV','CL'};
for nch=[3 4]
    lab = lab3; if (nch==4), lab=lab4; end
    pa=modpar26(nch);
    for f=[1000 2000 4000]
        try
            R=fdm26(struct('modal',1,'pa',pa,'f',f));
        catch e
            fprintf('\nnch=%d f=%g FAILED: %s\n', nch, f, e.message); continue
        end
        m=R.m; x=R.x; n=numel(x);
        kb=round(0.15*n); kc=R.icf; if (isempty(kc)||kc<1||kc>n), kc=round(0.5*n); end
        fprintf('\n===== nch=%d  f=%g Hz   (acoustic mode index %d, CF place x/L=%.3f) =====\n', ...
                nch, f, R.iac, x(kc)/max(x));
        fprintf('  place    mode  |kappa|   Re(k)     Im(k)   Im/|Re|  char        shape [%s]\n', strjoin(lab,' '));
        for tag=1:2
            if (tag==1), k=kb; nm='basal '; else, k=kc; nm='CF    '; end
            for j=1:m
                kp=R.kall(k,j); V=squeeze(R.Vall(k,:,j)); V=V/max(abs(V));
                rr=abs(imag(kp))/max(abs(real(kp)),eps);
                if (abs(kp) < 1e-3*max(abs(R.kall(k,:)))), ch='trivial';
                elseif (rr>3), ch='PROPAGATE';
                elseif (rr<0.33), ch='evanesc';
                else, ch='mixed'; end
                sv=sprintf('%+5.2f ', real(V));
                fprintf('  %s %3d  %9.3g %9.3g %9.3g %8.2f  %-10s %s\n', ...
                        nm, j, abs(kp), real(kp), imag(kp), rr, ch, sv);
            end
        end
    end
end
fprintf(['\nREAD: a secondary SS<->CL wave shows as a mode whose shape is an\n' ...
         'SS-vs-CL DIFFERENTIAL (opposite signs on those two, small ST/SV) AND\n' ...
         'whose character is PROPAGATE. If the SS/CL-dominant mode is evanescent,\n' ...
         'the secondary path exists structurally but does not carry a wave --\n' ...
         'consistent with the static basal accumulation in the OC-height DOF.\n']);
disp('MODES2_DONE');
