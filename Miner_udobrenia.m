%�������� ������
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

%���������� ���������� ��������
P=normrnd(Mat, Sigm, [m 5]).*Q/ro;

d=[];
P=P';
for i=1:m
    d(i)=mean(P(:, i));
end
M=d;

%������ ���������� ���������� ��������
figure 
plot(t, P, '-r');
title('������ ���������� ���������� ��������');
grid on
hold on;
plot(t, M,  '--b', 'LineWidth', 2);
hold off;

%���������� ���������
d=[];
for i=1:m
    d(i)=var(P(:, i));
end
D=d';

%���������� �������������� �������� ��������� ���������

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

%�������������������� ���������� ���������� ��������
Sigm=sqrt(D);

%���������� ������������� �������������� �������

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

%���������� ��������������� ��������, ��������� � ��������������������� ����������

MN=mean(M);
N=mean(D);
SigmN=sqrt(sum(D))/m;

%���������� ������������� �������������� ������� �� �������� ����������

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

%���������� ������������� ������������ ���������
%���������� � �������� ������� ������������� �������������� ������� �������
%�����

Q1=[0 0 0 0 0];

Q=[PxN Q1];
T=[t Q1];

%����� ��������� ����� �������
N=m+length(Q1);

%���������� �������� �������������� �����
W=fft(Q);
W1=abs(W);

%���������� ������
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
