For matlab version 9.4 (R2018a) or higher

English part:
load ('MD_10.mat')
load ('MI_10.mat')

[T, p] = TCheb1 (5,0,0.5);
function [T, p] = TCheb1 (n, a, b)
%% The function of computing the Chebyshev polynomials of the first kind
%% to enter:
% n - number of decomposition members
% a and b - segment coordinates
%% at the exit:
% T - column vector of Chebyshev orthonormal polynomials
% p - weight of polynomials

Cy = TChebC (cos (t), 5,0,0.5);
Cy = TChebC (cos (t), 5,0,0.5, T, p);
function Cy = TChebC (f, n, a, b, T, p)
%% The function of computing the Chebyshev polynomials of the first kind
%% to enter:
% f - function of the variable t
% n - number of decomposition members
% a and b - segment coordinates
%% optional parameter:
% T - column vector of Chebyshev orthonormal polynomials
% p - weight of polynomials
%% at the exit:
% �y - column vector of expansion coefficients

fi = TChebF (T, Cy)
fi = TChebF (T, Cy, Mi)
fi = TChebF (T, Cy, Mi, 1)
t0 = 0: 0.01: 0.5;
plot (t0, subs (f, t0), t0, subs (fi, t0), 'r')
function f = TChebF (T, Cy, A, F0)
%% polynomial function function
%% to enter:
% T - column vector of Chebyshev orthonormal polynomials
% �y - column vector of expansion coefficients
%% optional parameter:
% A is the integration / differentiation matrix.
% F0 - displacement of the integrated function (antiderivative value at point 0)
%% at the exit:
% f - total function

Mi = TChebMI (5,0,0.5)
Mi = TChebMI (5,0,0.5, T, p)
function Axi = TChebMI (n, a, b, T, p)
%% The function of computing the integration matrix for the expansion of the function in Chebyshev polynomials
%% to enter:
% n - number of decomposition members
% a - coordinate of the beginning of the segment
% b - coordinate of the end of the segment
%% optional parameter:
% T - column vector of Chebyshev orthonormal polynomials
% p - weight of polynomials
%% at the exit:
% Axi - integration matrix

Md = TChebMI (5,0,0.5)
Md = TChebMI (5,0,0.5, T, p)
function Axd = TChebMD (n, a, b, T, p)
%% The function of computing the differentiation matrix for the expansion of the function in Chebyshev polynomials
%% to enter:
% n - number of decomposition members
% a - coordinate of the beginning of the segment
% b - coordinate of the end of the segment
%% optional parameter:
% T - column vector of Chebyshev orthonormal polynomials
% p - weight of polynomials
%% at the exit:
% Axd - differentiation matrix

Mi = TChebMI (5,0,0.5, t)
Mi = TChebMI (5,0,0.5, t, T, p)
function Axk = TChebMK (n, a, b, K, T, p)
%% The function of computing the integration matrix for the expansion of the function in Chebyshev polynomials
%% to enter:
% n - number of decomposition members
% a - coordinate of the beginning of the segment
% b - coordinate of the end of the segment
% K - function for multiplication
%% optional parameter:
% T - column vector of Chebyshev orthonormal polynomials
% p - weight of polynomials
%% at the exit:
% Axk - multiplication matrix


Check
>> syms t;
>> I = eye (10);
>> A = (2 * I + 5 * MI + MI ^ 2) ^ (- 1) * MI ^ 2;
>> C = TChebC (1 + 0 * t, 10.0.1);
>> Cx = A * C;
>> x = TChebF (T, Cx);

>> syms xy (t); eqn = diff (xy, t, 2) * 2 + 5 * diff (xy, t, 1) + xy == 1; cond = xy (0) == 0; c = diff (xy, t, 1); cond1 = c (0) == 0; x1 = dsolve (eqn, [cond cond1])
>> plot (t0, subs (x1, t0), t0, subs (x, t0), 'r')



������� �����:
load('MD_10.mat')
load('MI_10.mat')

	[ T,p ] =TCheb1(5,0,0.5);
function [ T,p] = TCheb1( n,a,b )
%% ������� ���������� ����������� �������� ������� ����
%% �� ���� ��������:
% n - ���������� ������ ����������
% a � b - ���������� �������
%% �� ������:
% T - ������-������� ����������������� ����������� ��������
% p - ��� �����������

	Cy=TChebC(cos(t),5,0,0.5);
	Cy=TChebC(cos(t),5,0,0.5,T,p);
function Cy = TChebC(f,n,a,b  ,T,p )
%% ������� ���������� ����������� �������� ������� ����
%% �� ���� ��������:
% f - ������� �� ���������� t
% n - ���������� ������ ����������
% a � b - ���������� �������
%% �������������� ��������:
% T - ������-������� ����������������� ����������� ��������
% p - ��� �����������
%% �� ������:
% �y - ������-������� ������������� ���������� � ���

	fi=TChebF(T,Cy)
	fi=TChebF(T,Cy,Mi)
	fi=TChebF(T,Cy,Mi,1)
	   t0=0:0.01:0.5;
	   plot(t0,subs(f,t0),t0,subs(fi,t0),'r')
function f= TChebF( T,Cy ,A ,F0 )
%% ������� ���������� ������� �� �����������
%% �� ���� ��������:
% T - ������-������� ����������������� ����������� ��������
% �y - ������-������� ������������� ���������� � ���
%% �������������� ��������:
% A - ������� ��������������/�����������������. 
% F0 - �������� ������������������ ������� (�������� ������������� � ����� 0)
%% �� ������:
% f - �������� �������

	Mi=TChebMI(5,0,0.5)
	Mi=TChebMI(5,0,0.5,T,p)
function Axi = TChebMI( n,a,b  ,T,p )
%% ������� ���������� ������� �������������� ��� ���������� ������� �� ����������� ��������
%% �� ���� ��������:
% n - ���������� ������ ����������
% � - ���������� ������ �������
% b - ���������� ����� �������
%% �������������� ��������:
% T - ������-������� ����������������� ����������� ��������
% p - ��� �����������
%% �� ������:
% ��i - ������� ��������������

	Md=TChebMI(5,0,0.5)
	Md=TChebMI(5,0,0.5,T,p)
function Axd = TChebMD( n,a,b  ,T,p )
%% ������� ���������� ������� ����������������� ��� ���������� ������� �� ����������� ��������
%% �� ���� ��������:
% n - ���������� ������ ����������
% � - ���������� ������ �������
% b - ���������� ����� �������
%% �������������� ��������:
% T - ������-������� ����������������� ����������� ��������
% p - ��� �����������
%% �� ������:
% ��d - ������� �����������������

	Mi=TChebMI(5,0,0.5,t)
	Mi=TChebMI(5,0,0.5,t,T,p)
function Axk = TChebMK( n,a,b,K  ,T,p )
%% ������� ���������� ������� �������������� ��� ���������� ������� �� ����������� ��������
%% �� ���� ��������:
% n - ���������� ������ ����������
% � - ���������� ������ �������
% b - ���������� ����� �������
% K - ������� ��� ���������
%% �������������� ��������:
% T - ������-������� ����������������� ����������� ��������
% p - ��� �����������
%% �� ������:
% ��k - ������� ���������


��������
>> syms t;
>> I=eye(10);
>> A=(2*I+5*MI+MI^2)^(-1)*MI^2;
>> C=TChebC(1+0*t,10,0,1);
>> Cx=A*C;
>> x=TChebF(T,Cx);

>> syms xy(t); eqn=diff(xy,t,2)*2+5*diff(xy,t,1)+xy==1; cond=xy(0)==0; c=diff(xy,t,1); cond1=c(0)==0; x1=dsolve(eqn,[cond cond1])
>> plot(t0,subs(x1,t0),t0,subs(x,t0),'r')