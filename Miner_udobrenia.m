%Исходные данные
n=5;
m=11;
delta=0.5;
Q=0.35;
ro=8;

Mat=1.57;
Sigm=0.5;

tmax=(m-1)*delta;
for j=0:m-1
    t(j+1)=j*delta;
end

%реализация случайного процесса
P=normrnd(Mat, Sigm, [m 5]).*Q/ro;

d=[];
P=P';
for i=1:m
    d(i)=mean(P(:, i));
end
M=d;

%График реализации случайного процесса
figure 
plot(t, P, '-r');
title('График реализации случайного процесса');
grid on
hold on;
plot(t, M,  '--b', 'LineWidth', 2);
hold off;

%Вычисление дисперсии
d=[];
for i=1:m
    d(i)=var(P(:, i));
end
D=d';

%Вычисление корреляционных моментов случайных процессов

n=0;
d=[];
for j=1:m
    for i=1:n+1
        d(i, j)=M(i)*M(j);
    end
    n=n+1;
end

Mx=d;

n=5;
d=zeros(m, m);
for s=1:n
   d=d+KS(s, P, m);
end
ks=d;
Kx=(ks./n-Mx).*(n/(n-1));

%Среднеквадратическое отклонение случайного процесса
Sigm=sqrt(D);

%Вычисление нормированной корреляционной функции

n=0;
d=[];
for j=1:m
    for i=1:n+1
        d(i, j)=Sigm(i)*Sigm(j);
    end
    n=n+1;
end
S=d;

n=0;
d=[];
for j=1:m
    for i=1:n+1
        d(i, j)=Kx(i,j)/S(i,j);
    end
    n=n+1;
end
Px=d;

%Осреднение математического ожидания, дисперсии и среднеквадратического отклонения

MN=mean(M);
N=mean(D);
SigmN=sqrt(sum(D))/m;

%Осреднение нормированной корреляционной функции по ансамблю реализаций

n=0;
d=[];
for j=1:m
    d(j)=0;
    s=1;
    for i=n+1:m
        d(j)=d(j)+Px(s,i);
        s=s+1;
    end
    Sr(j)=d(j)/s;
    n=n+1;
end
PxN=Sr;

for i=0:m-1
    tau(i+1)=i*delta;
end

figure
plot(tau, PxN)

%Вычисление нормированной спектральной плотности
%Добавление к исходной матрице нормированной корреляционной функции времени
%нулей

Q1=[0 0 0 0 0];

Q=[PxN Q1];
T=[t Q1];

%Число элементов новой матрицы
N=m+length(Q1);

%Выполнение быстрого преобразования Фурье
W=fft(Q);
W1=abs(W);

%Вычисление частот
i=1:length(W);
Omega=(i+1).*delta/m;

Omega0=delta/m;
OmegaN=(N/2+1)*delta/m;


figure
plot(Omega, W1);

function [d]=KS(s, P, m)
d=[];
%n=5;
%s=n-1;
P=P';
for k=1:s
    n=0;
    for j=1:m
        for i=1:n+1
            d(i,j)=P(i,k)*P(j,k);
        end
        n=n+1;
    end
end
end
