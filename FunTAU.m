function F = FunTAU(x)
%
%   x=[Kp Kd Ki2 Ki1] (Kp)-1000 +(Kd)100s +(Ki2)14900/(s+(Ki1)14)
n=3;
a=0.01;
b=10;

load('MD_10.mat')
MD=MD(1:n,1:n);
MD=MD/(b-a);
load('MI_10.mat')
MI=MI(1:n,1:n);
MI=MI*(b-a);
I=eye(n);
A0=inv(I+4*MI+9*MI^2)*MI^3;
Kp=x(1);
Kd=x(2);
Ki2=x(3);
Ki1=x(4);
Aky=I*Kp+MD*Kd+inv(I+MI*Ki1)*MI*Ki2;
A=inv(I+A0*Aky)*A0*Aky;

Ae=inv(I+14*MI)*100*MI^2; % Разомкнутая
Ae=inv(I+Ae)*Ae; % Замкнутая
%Cy=TChebC(1,n,a,b);
F= reshape (abs(A-Ae), [numel(A) 1]);%*Cy;% для развернутых в столбец матриц, без умножения на Су
end

