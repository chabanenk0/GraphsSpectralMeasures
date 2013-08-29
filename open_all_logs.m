function open_all_logs(logfile,docalcentr);
if (nargin<2)
    docalcentr=1;
end

 baselogfile=logfile;
logfile1=strrep(baselogfile,'.txt','_degrees.txt');
open_log_series(logfile1,'degree',docalcentr);
logfile1=strrep(baselogfile,'.txt','_nb.txt');
open_log_series(logfile1,'node betweenness',docalcentr);
logfile1=strrep(baselogfile,'.txt','_eb.txt');
open_log_series(logfile1,'Edge betweenness',docalcentr);
logfile1=strrep(baselogfile,'.txt','_c3.txt');
open_log_series(logfile1,'Local clustering coef',docalcentr);
logfile1=strrep(baselogfile,'.txt','_wcc.txt');
open_log_series(logfile1,'weighted cluster coefficient',docalcentr);
logfile1=strrep(baselogfile,'.txt','_vert_ecc.txt');
open_log_series(logfile1,'vertex eccentrality',docalcentr);

