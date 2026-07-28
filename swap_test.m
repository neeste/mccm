% PARAMETER-SWAP TEST. Is the 3-chamber's collapsing high-CF tip caused by the
% THIRD CHAMBER, or inherited from the cel26c3 parameter set?
%   native2 : m=2, par_CEL16   (reference, d=0.445, tip sustained)
%   native3 : m=3, cel26c3     (reference, d=0.774, tip collapses)
%   swap2   : m=2, cel26c3     <-- decisive: c3 params, only two chambers
%   swap3   : m=3, par_CEL16   <-- complement: CEL16 params, three chambers
p2=modpar26(2); p3=modpar26(3);
s2=p3; s2.m=2;                        % cel26c3 params, 2 chambers (chsz unused for m<3)
s3=p2; s3.m=3; s3.chsz=[1 1 1];       % CEL16 params, 3 chambers (pad chsz)
cfg={{'native2 (m=2,CEL16)',p2},{'native3 (m=3,c3)',p3}, ...
     {'swap2   (m=2,c3)',s2},  {'swap3   (m=3,CEL16)',s3}};
for c=1:numel(cfg)
    nm=cfg{c}{1}; pa=cfg{c}{2};
    fprintf('\n===== %s =====\n',nm);
    try
        evalc('S=tdm26(0,pa,0,0);');
        nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf);
        BF=[]; NN=[];
        fprintf('  place    BF    contrast   Q3     N=tau*BF\n');
        for i=1:size(S.d1,2)
            D=fft(S.d1(:,i)-S.d2(:,i)); D=D(1:nf);
            H=D./max(abs(P),eps); A=abs(H).*f(:);
            ok=f>0.15 & f<18; A(~ok)=0;
            [pk,ip]=max(A); if (pk<=0), continue; end
            bf=f(ip); if (bf<0.6 || bf>16), continue; end
            [~,it]=min(abs(f-bf/4)); con=20*log10(pk/max(A(it),eps));
            q3=bf/max(deal_bw(A,f,ip,pk/10^(3/20)),eps);
            ph=unwrap(angle(H)); w=2*pi*f*1000; gd=-gradient(ph,w);
            tau=median(gd(max(1,ip-3):min(nf,ip+3)))*1000;
            fprintf('  %5d %7.3f %8.1f %6.2f %9.2f\n', pa.isv(i), bf, con, q3, tau*bf);
            BF(end+1)=bf; NN(end+1)=tau*bf; %#ok<AGROW>
        end
        if (numel(BF)>=3)
            b=polyfit(log(BF),log(NN),1);
            fprintf('  --> N ~ f^%.3f   =>  d = %.3f   (target 0.39-0.41)\n', b(1), 1-b(1));
        else
            fprintf('  --> too few valid places for a slope\n');
        end
    catch e
        fprintf('  FAILED: %s\n', e.message);
    end
end
disp('SWAP_TEST_DONE');
