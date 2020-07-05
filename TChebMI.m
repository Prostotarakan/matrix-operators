function Axi = TChebMI( n,a,b  ,T,p )
%% Функция вычисления матрицы интегрирования для разложения функции по многочленам Чебышева
%% на вход подаются:
% n - количество членов разложения
% а - координата начала отрезка
% b - координата конца отрезка
%% необязательный параметр:
% T - вектор-столбец ортонормированных многочленов Чебышева
% p - вес многочленов
%% на выходе:
% Ахi - матрица интегрирования
%%
syms t tau;
%if n<=10 %!!!!!  
    if nargin==3
        %% Вызов функции вычисления многочленов Чебышева первого рода:
        [T p]=TCheb1( n,a,b );
        %% Создание матрицы интегрирования:
        Axi=zeros(n);
        Axi(1,:)=double(vpa( int(  int( subs(T',tau),a,t )*T(1)*p   ,t,(a+0.0000001),(b-0.0000001) ),3))'; %первая строчка
        for f=2:n
            d=f-1;
            Axi(f,d)=double(vpa( int(  int( subs(T(d),tau),a,t )*T(f)*p   ,t,(a+0.0000001),(b-0.0000001) ),3));
            if f<n
                c=f+1;
                Axi(f,c)=double(vpa( int(  int( subs(T(c),tau),a,t )*T(f)*p   ,t,(a+0.0000001),(b-0.0000001) ),3));
            end
        end
        %Axi=double(vpa( int(  (int( subs(T,tau),(a),t )*subs(T',t)*p)   ,t,(a+0.0000001),(b-0.0000001) ),3))';  %%не считает матрицу для n>6;
    elseif nargin==5
        %% Создание матрицы интегрирования:
        Axi=zeros(n);
        Axi(1,:)=double(vpa( int(  int( subs(T',tau),a,t )*T(1)*p   ,t,(a+0.0000001),(b-0.0000001) ),3))'; %первая строчка
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

