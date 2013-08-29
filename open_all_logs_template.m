function open_all_logs_template(logfile,docalcentr);
if (nargin<2)
    docalcentr=1;
end

 baselogfile=logfile;
%logfile1=strrep(baselogfile,'.txt','_degrees.txt');
%open_log_series(logfile1,'degree',docalcentr);
logfile1=strrep(baselogfile,'.txt','_GrSpec.txt');
open_log_series(logfile1,'Graph spectrum',docalcentr);
%logfile1=strrep(baselogfile,'.txt','_FiedVectn.txt');
%open_log_series(logfile1,'FiedVectn',docalcentr);
%logfile1=strrep(baselogfile,'.txt','_Eigence.txt');
%open_log_series(logfile1,'Eigence',docalcentr);

logfile1=strrep(baselogfile,'.txt','_GrSpec_adjmatr.txt');
open_log_series(logfile1,'Graph spectrum adjacency matrix',docalcentr);
