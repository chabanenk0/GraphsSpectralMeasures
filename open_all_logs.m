function open_all_logs(logfile);

 baselogfile=logfile;
logfile1=strrep(baselogfile,'.txt','_degrees.txt');
open_log_series(logfile1,'degree');
logfile1=strrep(baselogfile,'.txt','_nb.txt');
open_log_series(logfile1,'node betweenness');
logfile1=strrep(baselogfile,'.txt','_eb.txt');
open_log_series(logfile1,'Edge betweenness');
logfile1=strrep(baselogfile,'.txt','_c3.txt');
open_log_series(logfile1,'Local clustering coef');
logfile1=strrep(baselogfile,'.txt','_wcc.txt');
open_log_series(logfile1,'weighted cluster coefficient');
logfile1=strrep(baselogfile,'.txt','_vert_ecc.txt');
open_log_series(logfile1,'vertex eccentrality');

