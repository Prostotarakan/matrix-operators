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
%% 
syms t;
%% ����� ������� ���������� ����������� �������� ������� ����:
if nargin==4
[T p]=TCheb1( n,a,b );
end
%% �������� �y - ������� ������������� ���������� � ���:
%[x1 x2]=size(f); % ���������� ����� � ��������
%if length(f)==1 %���������� �������
    Cy=double(vpa(int(p*T*f,t,a+0.00001,b-0.00001),3));
%else
%     Cy=zeros(n,1);
%     if x1==1 || x2==1 % �������� ����� ������ ����������
%         h=(b-a)/(length(f)-1); %���������� �� ���������� ��������
%         for i=1:length(f)
%             b=a+h;
%             Cy=Cy+double(vpa(int(p*T,t,a+0.00001,b-0.00001)*f(i),3)); %����������������
%             a=b;
%         end
%     else
%         if x1==2
%             f=f';
%         end
%         t0=f(:,1);
%         f=f(:,2);
%         for i=1:length(f)
%            Cy=Cy+double(vpa(int(p*T,t,a+0.00001,t0(i)-0.00001)*f(i),3)); %����������������
%            a=t0(i);
%         end
%     end
%     Cy=real(Cy)
% end

end

