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
%%
syms t tau;
%if n<=10 %!!!!!  
    if nargin==3
        %% ����� ������� ���������� ����������� �������� ������� ����:
        [T p]=TCheb1( n,a,b );
        %% �������� ������� ��������������:
        Axi=zeros(n);
        Axi(1,:)=double(vpa( int(  int( subs(T',tau),a,t )*T(1)*p   ,t,(a+0.0000001),(b-0.0000001) ),3))'; %������ �������
        for f=2:n
            d=f-1;
            Axi(f,d)=double(vpa( int(  int( subs(T(d),tau),a,t )*T(f)*p   ,t,(a+0.0000001),(b-0.0000001) ),3));
            if f<n
                c=f+1;
                Axi(f,c)=double(vpa( int(  int( subs(T(c),tau),a,t )*T(f)*p   ,t,(a+0.0000001),(b-0.0000001) ),3));
            end
        end
        %Axi=double(vpa( int(  (int( subs(T,tau),(a),t )*subs(T',t)*p)   ,t,(a+0.0000001),(b-0.0000001) ),3))';  %%�� ������� ������� ��� n>6;
    elseif nargin==5
        %% �������� ������� ��������������:
        Axi=zeros(n);
        Axi(1,:)=double(vpa( int(  int( subs(T',tau),a,t )*T(1)*p   ,t,(a+0.0000001),(b-0.0000001) ),3))'; %������ �������
        for f=2:n
            if f==7
                f
            end
            d=f-1;
            Axi(f,d)=double(vpa( int(  int( subs(T(d),tau),a,t )*T(f)*p   ,t,(a+0.0000001),(b-0.0000001) ),3));
            if f<n
                c=f+1;
                Axi(f,c)=double(vpa( int(  int( subs(T(c),tau),a,t )*T(f)*p   ,t,(a+0.0000001),(b-0.0000001) ),3));
            end
        end
    end
%else
%    Axi=0;
%end
end

