function en=ApEn(x,l,st)
% Calculate Approximate Entropy (�������������� �������):
%
%    en=ApEn(x,l,st)
%    ����� ���:
%       x      - ������-����� �� ������-��������
%       l      - ������� �������� (�������) (default 5)
%       st     - �������� ������� �������� (default 1)
%    ������ ���:
%       en     - �������� �����ﳿ
%
%    �������:
%       Moody G. Approximate Entropy (ApEn) // http://www.physionet.org/
%
% (c) ������ �., 19.11.2006

% ������� ���������
if (nargin==0)                          % ������ ��������� �� ������
    disp('  �� ������� ���������!');
    help ApEn
    return
elseif (nargin==1)                      % ������ ���� ������ ������� �����
    l=5; st=1;
elseif (nargin==2)                      % �� ������ ������� ��������
    st=1;
end

% �������� �������� �������
[m,n]=size(x);                          % ���� ������-����� - �����������
if (m<n) x=x'; [m,n]=size(x); end
if (n>1) x=x(:,1); end                  % ������ ���� ������ ��������

% ���������� ����������� �������� ������� �������� l �� l+1
c_curr=ApEnSimCoeff(x,l,st);
c_next=ApEnSimCoeff(x,l+1,st);

% ���������� �������� �����ﳿ
en=log(c_curr/c_next);

clear m n c_curr c_next;