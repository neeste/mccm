% ldngrw_st - loudness summation @ 60 dB SPL for each stimulus type
ddn={'qolr_tone','qolr_flat','qolr_lono','qolr_amod'};
nrc=[1 4 4 1];                % repetition count
for k=1:4                     % make subfolders ???
    if (~exist(ddn{k},'dir')), mkdir(ddn{k}); end
end
if (~exist('vep','dir')), mkdir('vep'); end
tdm24(0,1,0,0)                % Ped calibration
tdm24('ldngrw')               % single tone only
for k=1:4
    st=k-1; % stimulus type: 0=five_tone 1=flat_noise 2=lono_noise 3=AM_tone
    for rc=1:nrc(k)               % repetition count
        tdm24('ldngrw',1,st,rc)   % loudness summation: level & bandwidth
        dd=ddn{k};                % destination folder
        tfn='ldnsumbw.txt';       % text file
        dfn=sprintf('ldnsumbw%d.txt',rc);
        dpn=fullfile(ddn{k},dfn); % destination path name
        sc=copyfile(tfn,dpn);     % copy txt file
        fprintf('%s->%s %d\n',tfn,dpn,sc);
        wfn=fullfile('vep','vep*.mat');
        sc=movefile(wfn,dd);      % move vep file
        fprintf('%s->%s\n',wfn,dpn);
    end
end
