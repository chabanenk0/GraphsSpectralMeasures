function [e]=permen(y,m,L)
%sd = std(y);
%y = (y-mean(y))/sd;
ly = length(y);
    
permlist = perms(1:m);
c(1:length(permlist))=0;
    
for i=1:ly-(m-1)*L
    y1=[];
    for k=1:m
        y1=[y1 y(i+(k-1)*L)];
    end;
    [a,iv]=sort(y1);
    for j=1:length(permlist)
        if (abs(permlist(j,:)-iv))==0
            c(j) = c(j) + 1 ;
            break;
        end
    end
end
    
p = c/(ly-(m-1)*L);

%figure
%plot(p)
%title(num2str(sum(p)))

e=0;
for i=1:length(p)
    if p(i)~=0
        e=e+(p(i)*log(p(i)));
    end
end  
e=-e;
e=e/log(factorial(m));