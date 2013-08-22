function save_log_series(curlogfile,i,measure);
% функция сохранения каждой конкретной меры (строки в соотв. лог-файле)
n=length(measure);
fp=fopen(curlogfile,'a');
if (i==1)
    fprintf(fp,'t\t');
    for j=1:n-1
        fprintf(fp,'%d\t',j);
    end
    fprintf(fp,'%d\n',n);
end
fprintf(fp,'%d\t',i);
for j=1:n-1
    fprintf(fp,'%d\t',measure(j));
end
fprintf(fp,'%d\n',measure(n));
fclose(fp);
