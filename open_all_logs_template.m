function open_all_logs_template(logfile);

 baselogfile=logfile;
%logfile1=strrep(baselogfile,'.txt','_degrees.txt');
%open_log_series(logfile1,'degree');
logfile1=strrep(baselogfile,'.txt','_GrSpec.txt');
open_log_series(logfile1,'Graph spectrum');
%logfile1=strrep(baselogfile,'.txt','_FiedVectn.txt');
%open_log_series(logfile1,'FiedVectn');
%logfile1=strrep(baselogfile,'.txt','_Eigence.txt');
%open_log_series(logfile1,'Eigence');
