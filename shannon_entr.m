function Sh=shannon_entr(TS);

minTS = min(TS);
maxTS = max(TS);
n = length(TS);
n_ver=ceil(sqrt(n));
if (0)
X = [minTS:(maxTS-minTS)/n_ver:maxTS];
H = hist(TS,X);
else
    H = hist(TS,n_ver);
end    

P = H/n;

% shannon Sh=sum(pi*ln(pi))
Sh=0;
for i=1:n_ver
    if (P(i)>0)
        Sh=Sh+P(i)*log(P(i));
    end
end
Sh=-Sh;
