syms t;
n=7;
a=0;
b=10;
 
load('MI_10.mat')
MI=MI(1:n,1:n);
MI=MI*(b-a);
I=eye(n);
 
[T ~]=TCheb1(n,a,b);
A=inv(I+11*MI+35*MI^2+37*MI^3+30*MI^4+6*MI^5)*(80*MI^2+200*MI^3+40*MI^4);
Cy=TChebC(1,n,a,b);
f=TChebF(T,Cy,A);
