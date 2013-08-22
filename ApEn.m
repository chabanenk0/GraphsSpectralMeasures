function en=ApEn(x,l,st)
% Calculate Approximate Entropy (апроксимаційна ентропія):
%
%    en=ApEn(x,l,st)
%    Вхідні дані:
%       x      - вектор-рядок чи вектор-стовпчик
%       l      - довжина паттерна (вектора) (default 5)
%       st     - значення критерія подібності (default 1)
%    Вихідні дані:
%       en     - значення ентропії
%
%    Джерела:
%       Moody G. Approximate Entropy (ApEn) // http://www.physionet.org/
%
% (c) Сердюк О., 19.11.2006

% обробка параметрів
if (nargin==0)                          % обидва параметри не задано
    disp('  Не вистачає параметрів!');
    help ApEn
    return
elseif (nargin==1)                      % задано лише вектор вхідних даних
    l=5; st=1;
elseif (nargin==2)                      % не задано критерій подібності
    st=1;
end

% перевірка вхідного вектора
[m,n]=size(x);                          % якщо вектор-рядок - транспонуємо
if (m<n) x=x'; [m,n]=size(x); end
if (n>1) x=x(:,1); end                  % беремо лише перший стовпчик

% обраховуємо коефіцієнти подібності векторів довжиною l та l+1
c_curr=ApEnSimCoeff(x,l,st);
c_next=ApEnSimCoeff(x,l+1,st);

% обраховуємо значення ентропії
en=log(c_curr/c_next);

clear m n c_curr c_next;