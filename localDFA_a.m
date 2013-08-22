function [a]=localDFA_a(ARG,stepBegin,stepEnd,step, ratio);

if (size(ARG,1)>size(ARG,2))
    ARG=transpose(ARG);
end

% ¬ычисление нормализированных прибыльностей
n=length(ARG);
ARG=(ARG(2:n)-ARG(1:n-1))./(0.5*(ARG(2:n)+ARG(1:n-1)));
sigm=std(ARG);
ARG=ARG/sigm;

 stepBegin=8;
 stepEnd=floor(length(ARG)/2);
 step=5;
 ratio=1.05;
    [x,y]=DFA(ARG,stepBegin,stepEnd,step,ratio);
    [a,b]=MinSquare(log(x),log(y));
    
end




