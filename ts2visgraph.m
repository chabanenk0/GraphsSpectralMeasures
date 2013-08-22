function Adj=ts2visgraph(y);
% Transforming the time series to the visibility graphs 
n=length(y);
Adj=zeros(n,n);
for i=1:n-1
    for j=i+1:n
        Adj(i,j)=visib_left_right(y,i,j);
        Adj(j,i)=Adj(i,j);
        %Adj(j,i)=visib_left_right(y,j,i);
    end
end

function vis=visib_left_right(y,i,j);
% return whether the node j is visible from i (from left to right)
% calculating the line equation, that crosses two points:
% y(i)=a0+a1*i
% y(j)=a0+a1*j
D=j-i;
Da0=y(i)*j-y(j)*i;
a0=Da0/D;
Da1=y(j)-y(i);
a1=Da1/D;
vis=1;
for k=i+1:j-1
    if (y(k)>a0+a1*k)
        vis=0;
        break;
    end
end

function vis=visib_right_left(y,i,j);
% return whether the node j is visible from i (from left to right)
% calculating the line equation, that crosses two points:
% y(i)=a0+a1*i
% y(j)=a0+a1*j
D=j-i;
Da0=y(i)*j-y(j)*i;
a0=Da0/D;
Da1=y(j)-y(i);
a1=Da1/D;
vis=1;
for k=i+1:j-1
    if (y(k)>a0+a1*k)
        vis=0;
        break;
    end
end

            