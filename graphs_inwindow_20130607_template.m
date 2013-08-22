filename='dj_19_10_1987.txt';
y=dlmread(filename);
n=length(y);
wind=500;
tstep=1;
thr_fr=0.5;
crp_emb_dimm=1;
crp_delay=1;
crp_epsilon=0.4; %0.1 - Marwan default
logfile='tcp_vg_500_1_spec.txt';
graph_type=2;%1 - crp, 2 - visibility
%graph_type=2;%visibility
% 20130601 Graph Entropy
Q=5;
types=round((1:wind+1)/(wind/(Q-1))+1);% вектор классов узлов (каждый из узлов должен быть отнесен к какому-то классу)
Nrunmax=3;% кол-во запусков процедуры
Nd=10;% кількість градацій відстані для ентропії відстаней
   % l=logspace(log(d_min)/log(10),log(Dmax)/log(10),Nd);
   %l=linspace(-14,55,30);
 % l=linspace(0.2,9,Nd);
skip_complex_measures=1;
% end 20130601 Graph Entropy

% Создание пустых массивов для динам. мер
% В случае глобальных мер результат - одно число. (например, мера "диаметр")
% mas_diam=[];
% В случае локальных мер на выходе - вектор, характеризующий меру для
% каждого узла или ребра. Например, степень вершин (для каждой вершины -
% своя степень)
% Массив степеней (матрица)
mas_degr=[];
% Статистические характеристики степеней (макс, мин, среднее)
mas_maxdegr=[];
mas_mindegr=[];
mas_meandegr=[];

mas_graphSpectrum=[];
mas_graphSpectrum_delta=[];
mas_min_graphSpectrum=[];
mas_max_graphSpectrum=[];
mas_mean_graphSpectrum=[];

mas_algebraic_connectivity = [];

mas_Alg_Connect=[];
mas_fiedler_vector = [];
mas_min_fiedler_vector=[];
mas_max_fiedler_vector=[];
mas_mean_fiedler_vector=[];

mas_eigencentrality = [];
mas_min_eigencentrality=[];
mas_max_eigencentrality=[];
mas_mean_eigencentrality=[];
mas_eigence_maxlambda=[];
mas_max_lambda=[];
mas_graph_energy = [];

% Создание лог-файла. Запись заголовков.
% При добавлении новых мер, нужно добавить соотв. заголовок. \t - символ
% табуляции, разделяющий ячейки будущей таблицы
fp=fopen(logfile,'a');
fprintf(fp,'filename\tt\tmax node degree\tmin node degree\tmean node degree\tmaxGrSpectr\tminGrSpec\tmeanGrSpec\tAlgConnect\tmaxFiedVectn\tminFiedVectn\tmeanFiedVectn\tEigence\tGraph_ener\tgraphsSpectralDelta\n');
fclose(fp);
for i=1:tstep:n-wind % главный цикл
    i
    n-wind
    y_fragm=y(i:i+wind); %вырезание фрагмента (окна) исходного ряда
    % вычисление матрицы смежности (или срп, или граф видимости)
    if(graph_type==1)
        Adj=double(crp(y_fragm,crp_emb_dimm,crp_delay,crp_epsilon));
    else
        Adj=ts2visgraph(y_fragm);
    end
    
    % Вычисление меры диаметра (также можно добавить свою функцию, которая 
    % возвращает одно число для глобальных или вектор для локальных мер)
    % [d,dij]=diameter(Adj);
    % d=0;
    % mas_diam=[mas_diam;d]; % добавление посчитанного значения в массив

    % Вычисление меры степеней вершин (на выход - вектор)
    degr=degrees(Adj);
    mas_degr=[mas_degr;degr];% Сохранение всего вектора как строки матрицы
    maxdegr=max(degr);% вычисление максимума степеней вершин (получается глобальная мера на основе локальной)
    mas_maxdegr=[mas_maxdegr;maxdegr]; % сохранение максимума в массиве
    mindegr=min(degr); % аналогично минимуму
    mas_mindegr=[mas_mindegr;mindegr]; % Сохранение (аналогично как выше, в отдельные массивы)
    meandegr=mean(degr); % среднее значение
    mas_meandegr=[mas_meandegr;meandegr]; % сохранение среднего значения...
    % В конце каждой итерации результаты вычисления добавляются в лог-файл
    % При добавлении своих мер, нужно соотв. заголовкам (см. fprintf выше
    % строка 41 этой программы, добавить выведение новых мер (еще один %d и
    % имя переменной, откуда вывести результат в конец списка параметров
    % fprintf
    gr_sp=graph_spectrum_unsort(Adj);
    mas_graphSpectrum=[mas_graphSpectrum;gr_sp];
    min_gr_sp=min(gr_sp);
    mas_min_graphSpectrum=[mas_min_graphSpectrum,min_gr_sp];
    max_gr_sp=max(gr_sp);
    mas_max_graphSpectrum=[mas_max_graphSpectrum,max_gr_sp];
    mean_gr_sp=mean(gr_sp);
    mas_mean_graphSpectrum=[mas_mean_graphSpectrum,mean_gr_sp];
    
    % 20130717 eigenv(max)-eigenv(max2)  второе по величине после макс
    gr_sp_sort=-sort(-gr_sp);
    gr_sp_delta=gr_sp_sort(1)-gr_sp_sort(2);
    mas_graphSpectrum_delta=[mas_graphSpectrum_delta,gr_sp_delta];
    % end 20130717
    Alg_Connect=algebraic_connectivity(Adj);
    mas_Alg_Connect=[mas_Alg_Connect;Alg_Connect];
    
    FiedVect=fiedler_vector(Adj);
    mas_fiedler_vector=[mas_fiedler_vector;FiedVect];
    min_FiedVect=min(FiedVect);
    mas_min_fiedler_vector=[mas_min_fiedler_vector,min_FiedVect];
    max_FiedVect=max(FiedVect);
    mas_max_fiedler_vector=[mas_max_fiedler_vector,max_FiedVect];
    mean_FiedVect=mean(FiedVect);
    mas_mean_fiedler_vector=[mas_mean_fiedler_vector,mean_FiedVect];
    
    Eigence=eigencentrality(Adj);
    mas_eigencentrality=[mas_eigencentrality;Eigence];
    min_Eigence=min(Eigence); 
    mas_min_eigencentrality=[mas_min_eigencentrality,min_Eigence];
    max_Eigence=max(Eigence);
    mas_max_eigencentrality=[mas_max_eigencentrality,max_Eigence];
    mean_Eigence=mean(Eigence);
    mas_mean_eigencentrality=[mas_mean_eigencentrality,mean_Eigence];
    Lambda=maxlambda(Adj);
    max_Lambda=max(Lambda);
    mas_max_lambda=[mas_max_lambda,max_Lambda];
    
    
    Graph_ener=graph_energy(Adj);
    mas_graph_energy =[mas_graph_energy ;Graph_ener];
  
    fp=fopen(logfile,'a');
    fprintf(fp,'%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',filename,i,maxdegr,mindegr,meandegr,min_gr_sp,max_gr_sp,mean_gr_sp,Alg_Connect,min_FiedVect,max_FiedVect,mean_FiedVect,min_Eigence,max_Eigence,mean_Eigence,max_Lambda,Graph_ener,gr_sp_delta );
    fclose(fp);
    % Сохранение локальных мер выполняется в функции save_logs_template. При
    % добавлении новой локальной меры, нужно: добавить в параметр этой
    % функции (следующая строка) новую меру, а также подкорректировать саму
    % функцию (открыть для редактирования save_logs_template.m  например командой edit save_logs_template
    % Дальнейшие инструкции по исправлению см. в save_logs_template.m
    save_logs_template (logfile,i,degr,gr_sp,FiedVect,Eigence);
end

% Вывод на экран глобальной меры. Для своей меры нужно скопировать такой же
% фрагмент, но исправить имя переменной-массива и подписи.
% figure;
% plot(mas_diam); % сюда вместо mas_diam нужно имя своего массива
% ('diameter'); %сюда - имя заголовка графика
% outfile=strrep(filename,'.txt','_diameter.txt');% тут можно подправить файл, куда сохранится массив как ряд 
% dlmwrite(outfile,mas_diam,'\n');% тут тоже подкорректировать имя переменной



% Для глобальных мер - аналогично, но выводится 3 графика при
% необходимости: максимум, минимум и среднее значение локальной меры
figure;
plot(mas_degr);
title('node degree dynamics');

figure;
plot(mas_maxdegr);
title('max node degree dynamics');
outfile=strrep(filename,'.txt','_maxdegr.txt');
dlmwrite(outfile,mas_maxdegr,'\n');
figure;
plot(mas_mindegr);
title('min node degree dynamics');
outfile=strrep(filename,'.txt','_mindegr.txt');
dlmwrite(outfile,mas_mindegr,'\n');
figure;
plot(mas_meandegr);
title('mean node degree dynamics');
outfile=strrep(filename,'.txt','_meandegr.txt');
dlmwrite(outfile,mas_meandegr,'\n');

figure;
plot(mas_graphSpectrum);
title('graphSpectrum');

figure;
plot(mas_max_graphSpectrum);
title('max graphSpectrum');
outfile=strrep(filename,'.txt','_max_graphSpectrum.txt');
dlmwrite(outfile,mas_max_graphSpectrum,'\n');

figure;
plot(mas_min_graphSpectrum);
title('min graphSpectrum');
outfile=strrep(filename,'.txt','_min_graphSpectrum.txt');
dlmwrite(outfile,mas_min_graphSpectrum,'\n');

figure;
plot(mas_mean_graphSpectrum);
title('mean graphSpectrum');
outfile=strrep(filename,'.txt','_mean_graphSpectrum.txt');
dlmwrite(outfile,mas_mean_graphSpectrum,'\n');

figure;
plot(mas_graphSpectrum_delta);
title('GraphSpectrum delta');
outfile=strrep(filename,'.txt','_graphSpectrum_delta.txt');
dlmwrite(outfile,mas_graphSpectrum_delta,'\n');


figure;
plot(mas_Alg_Connect);
title('algebraic connectivity');
outfile=strrep(filename,'.txt','_algebraic connectivity.txt');
dlmwrite(outfile,mas_Alg_Connect,'\n');

figure;
plot(mas_fiedler_vector);
title('fiedler vector');

figure;
plot(mas_max_fiedler_vector);
title('max fiedler vector');
outfile=strrep(filename,'.txt','_max_fiedler_vector.txt');
dlmwrite(outfile,mas_max_fiedler_vector,'\n');

figure;
plot(mas_min_fiedler_vector);
title('min fiedler vector');
outfile=strrep(filename,'.txt','_min_fiedler_vector.txt');
dlmwrite(outfile,mas_min_fiedler_vector,'\n');

figure;
plot(mas_mean_fiedler_vector);
title('mean fiedler vector');
outfile=strrep(filename,'.txt','_mean_fiedler_vector.txt');
dlmwrite(outfile,mas_mean_fiedler_vector,'\n');

figure;
plot(mas_eigencentrality);
title('eigencentrality');

figure;
plot(mas_max_eigencentrality);
title('max eigencentrality');
outfile=strrep(filename,'.txt','_max_eigencentrality.txt');
dlmwrite(outfile,mas_max_eigencentrality,'\n');

figure;
plot(mas_min_eigencentrality);
title('min eigencentrality');
outfile=strrep(filename,'.txt','_min_eigencentrality.txt');
dlmwrite(outfile,mas_min_eigencentrality,'\n');

figure;
plot(mas_mean_eigencentrality);
title('mean eigencentrality');
outfile=strrep(filename,'.txt','_mean_eigencentrality.txt');
dlmwrite(outfile,mas_mean_eigencentrality,'\n');

figure;
plot(mas_max_lambda);
title('max lambda');
outfile=strrep(filename,'.txt','_max_lambda.txt');
dlmwrite(outfile,mas_max_lambda,'\n');

figure;
plot(mas_graph_energy);
title('graph energy');
outfile=strrep(filename,'.txt','_graph energy.txt');
dlmwrite(outfile,mas_graph_energy,'\n');




% Пример вывода максимума, минимума и среднего локальной меры на одном
% графике

% figure;
% plot(mas_max_closeness,'k');
% hold on;
% plot(mas_min_closeness,'k');
% plot(mas_mean_closeness,'r');
% title('closeness');
% outfile=strrep(filename,'.txt','_MaxCloseness.txt');
% dlmwrite(outfile,mas_max_closeness,'\n');
% outfile=strrep(filename,'.txt','_MinCloseness.txt');
% dlmwrite(outfile,mas_min_closeness,'\n');
% outfile=strrep(filename,'.txt','_MeanCloseness.txt');
% dlmwrite(outfile,mas_mean_closeness,'\n');

% Далее вызывается функция open_all_logs, которая открывает сохраненные
% матрицы локальных мер и выводит графики. При добавлении новых локальных
% мер, необходимо добавить для новой меры последовательность операторов,
% аналогичные добавленным.
open_all_logs_template(logfile);
