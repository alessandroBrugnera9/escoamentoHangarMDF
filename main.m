close all;clc; clear;

%Constantes
V=100/3.6;
h=3;
d=5*h;
L=2*h;
H=8*h;
ro=1.25;
gama=1.4;
k=0.026;
cp=1002;

Tdentro=40;
Tfora=20;


%variaveis computacionais
dx=.01
dy=.01
lambda=1.85
eps=0.01



xi=2;
xf=(2*d+L)/dx;
yi=2;
yf=H/dy;


corr1=zeros(yf,xf);
corr2=zeros(yf,xf);

%popula com 1 o que tem dentro
for i = yi:yf
	for j = xi:xf
		corr1(i,j) = 1;
	end
end

gxi=d/dx+1; %x esquerda
gxf=gxi+(L)/dx; %x direita
gyi=yi; %
gyf=gyi+h/dy-1; %

%popula com 2 o que tem dentro do galpao
for i = gyi:gyf
	for j = gxi:gxf
		corr1(i,j) = 2;
	end
end

%popula com 3 o que tem no topo do galpao
for i = gyf:(gyf+(L/2)/dy-1)
	for j = gxi:gxf
		if (i*dy) <= (sqrt((L/2)^2 - (j*dx-d-L/2)^2)+h)
			corr1(i,j) = 3;
		end
	end
end






