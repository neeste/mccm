% ldngrw_av - loudness summation average across stimulus types
function ldngrw_av
ddn={'qolr_tone','qolr_flat','qolr_lono','qolr_amod'};
nrc=[1 4 4 1];               % repetition count
for k=1:4
    m=zeros(6,7,nrc(k));
    for rc=1:nrc(k)               % repetition count
        dfn=sprintf('ldnsumbw%d.txt',rc);
        tpn=fullfile(ddn{k},dfn); % text path name
        m(:,:,rc)=read_data(tpn);
    end
    dfn=fullfile(ddn{k},'ldnsumbw.txt'); % average
    mav=mean(m,3);
    msd=std(m,[],3);
    write_data(dfn,[mav msd])
end
return

function m=read_data(fn)
fp=fopen(fn,'r');
if (fp == -1)
    fprintf('unable to open file %s\n',fn)
    return;
end
%-----------------------
m=[];
% get first line of text
txt  = fgetl(fp); 
% skip over lines that begin with ;
while(txt(1) == ';')
    txt  = fgetl(fp);
    if (txt == -1), break; end
end
% add data as columns
while(txt(1) ~= ';')
    p = m;            % previous m
    m = [p sscanf(txt,'%f')];
    txt  = fgetl(fp); % get another line of text
    if (txt == -1), break; end
end
if (~isempty(m)), m = transpose(m); end
%-----------------------
fclose(fp);
return

function write_data(fn,data)
[nr,nc] = size(data);
fp=fopen(fn,'wt');
fprintf(fp,'; %s\n', fn);
for i=1:nr
    for j=1:nc
        fprintf(fp,' %14.5g',data(i,j));
    end
    fprintf(fp,'\n');
end
fclose(fp);
return

