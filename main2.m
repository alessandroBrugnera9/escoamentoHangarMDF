close all;clc; clear;


load('corr.mat','corr')
load('deltap.mat','deltap')
load('pressao.mat','pressao')
load('telhadoP.mat','telhadoP')
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
patm=10000;
comprimento=60; %fundura

Tdentro=40;
Tfora=20;



%2a aprte
dx=h/8;
dy=h/8;
lambda=1.15;
eps=0.01;

xi=2;
xf=(2*d+L)/dx-1;
yi=2;
yf=H/dy-1;


T=zeros(yf+1,xf+1);

%condicao de dischlet a esquerda
for i=1:yf+1
	T(i,1) = Tfora;
end

gxi=d/dx+1; %x esquerda
gxf=gxi+(L)/dx; %x direita
gyi=yi; %
gyf=gyi+h/dy-1; %

%popula com Tdentro o que tem dentro do galpao
for i = gyi:gyf
	for j = gxi:gxf
		T(i,j) = Tdentro;
	end
end

%popula com 3 o que tem no topo do galpao
for i = gyf:(gyf+(L/2)/dy-1)
	for j = gxi:gxf
		if (i*dy) <= (sqrt((L/2)^2 - ((j-1)*dx-d-L/2)^2)+h)
			T(i,j) = Tdentro;
		end
	end
end




contador=0;
convergiu=false;

while ~convergiu
	contador = contador+1;
	convergindo=true;

	%aplicando condicoes nas bordas utilizando equacionamento de taylor
	%para considerar condicoes de Neumann


	%primeiro nas pontas
	%cima
	noAtual=(T(yf+1,xf)+T(yf,xf+1))/2;
	noAntigo =T(i,j);
	T(yf+1,xf+1) = lambda*noAtual + (1-lambda)*noAntigo; %sobrerrelaxacao
	%baixo
	noAtual=(T(1,xf+1)+T(2,xf))/2;
	noAntigo =T(i,j);
	T(1,xf+1) = lambda*noAtual + (1-lambda)*noAntigo; %sobrerrelaxacao
	%segundo no topo
	u=V;
	for j = xi:xf
		noAtual=(k/dx^2*(2*T(yf,j) +  T(yf+1,j+1)+T(yf+1,j-1)) + ro*cp*u/dx*(T(yf+1,j-1)))/(4*k/dx^2 + ro*cp*u/dx);
		noAntigo =T(i,j);
		T(yf+1,j) = lambda*noAtual + (1-lambda)*noAntigo; %sobrerrelaxacao
	end

	%segundo emabvixo
		%antes do galpao
	for j = xi:gxi-1
		u=(corr(i+1,j)-corr(i-1,j))/(2*dy);
		noAtual=(k/dx^2*(2*T(2,j) +  T(1,j+1)+T(1,j-1)) + ro*cp*u/dx*(T(1,j-1)))/(4*k/dx^2 + ro*cp*u/dx);
		noAntigo =T(i,j);	
		T(1,j) = lambda*noAtual + (1-lambda)*noAntigo; %sobrerrelaxacao
	end
		%depois do galpao
	for j = gxf+1:xf
		u=(corr(i+1,j)-corr(i-1,j))/(2*dy);
		noAtual=(k/dx^2*(2*T(2,j) +  T(1,j+1)+T(1,j-1)) + ro*cp*u/dx*(T(1,j-1)))/(4*k/dx^2 + ro*cp*u/dx);
		noAntigo =T(i,j);	
		T(1,j) = lambda*noAtual + (1-lambda)*noAntigo; %sobrerrelaxacao
	end

	%terceito a direita
	for i = yi:yf
		u=(corr(i+1,xf+1)-corr(i-1,xf+1))/(2*dy);
		v=0;
		noAtual=(k/dx^2*(2*T(i,xf) +  T(i+1,xf+1)+T(i-1,xf+1)) + ro*cp*u/dx*(T(i-1,xf+1)))/(4*k/dx^2 + ro*cp*u/dx);
		noAntigo =T(i,j);
		T(i,xf+1) = lambda*noAtual + (1-lambda)*noAntigo; %sobrerrelaxacao
	end


	%Aplicando o metrodo nos nos internos
	%primeiro na altura inferior ao topo do telhado
	for i = gyi:gyf
		for j = xi:gxi-1 %antes do galpao
			u=(corr(i+1,j)-corr(i-1,j))/(2*dy);
			v=(corr(i,j+1)-corr(i,j-1))/(2*dx);
			if u<0 & v<0
				noAtual = (k/dx^2*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)) - ro*cp*u/dx*(2*T(i,j+1)))/(4*k/dx^2 + 2*ro*cp*u/dx);
			elseif u>0 & v>0
				noAtual = (k/dx^2*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)) + ro*cp*u/dx*(2*T(i,j-1)))/(4*k/dx^2 - 2*ro*cp*u/dx);
			else
				noAtual = (k/dx^2*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)) + ro*cp*u/dx*(T(i,j+1) - T(i,j-1)))/(4*k/dx^2);
			end

			noAntigo =T(i,j);
			T(i,j) = lambda*noAtual + (1-lambda)*noAntigo; %sobrerrelaxacao

			if convergindo %para melhorar desempenho e nao fazer contas desnecessarias
				if abs((T(i,j)-noAntigo))/T(i,j) > eps
					convergindo = false; %fazer mais interacoes
				end
			end
		end

		for j = gxf+1:xf %depois do galpao
			u=(corr(i+1,j)-corr(i-1,j))/(2*dy);
			v=(corr(i,j+1)-corr(i,j-1))/(2*dx);
			if u<0 & v<0
				noAtual = (k/dx^2*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)) - ro*cp*u/dx*(2*T(i,j+1)))/(4*k/dx^2 + 2*ro*cp*u/dx);
			elseif u>0 & v>0
				noAtual = (k/dx^2*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)) + ro*cp*u/dx*(2*T(i,j-1)))/(4*k/dx^2 - 2*ro*cp*u/dx);
			else
				noAtual = (k/dx^2*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)) + ro*cp*u/dx*(T(i,j+1) - T(i,j-1)))/(4*k/dx^2);
			end

			noAntigo =T(i,j);
			T(i,j) = lambda*noAtual + (1-lambda)*noAntigo; %sobrerrelaxacao

			if convergindo %para melhorar desempenho e nao fazer contas desnecessarias
				if abs((T(i,j)-noAntigo))/T(i,j) > eps
					convergindo = false; %fazer mais interacoes
				end
			end
		end
	end


	%segundo na altura superior ao topo do telhado
	for i = gyf+1:yf
		for j = xi:xf
			if  (i*dy) > (sqrt((L/2)^2 - ((j-1)*dx-d-L/2)^2)+h) %checa se esta nos limites externos do telhado
				u=(corr(i+1,j)-corr(i-1,j))/(2*dy);
				v=(corr(i,j+1)-corr(i,j-1))/(2*dx);
				if u<0 & v<0
					noAtual = (k/dx^2*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)) - ro*cp*u/dx*(2*T(i,j+1)))/(4*k/dx^2 + 2*ro*cp*u/dx);
				elseif u>0 & v>0
					noAtual = (k/dx^2*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)) + ro*cp*u/dx*(2*T(i,j-1)))/(4*k/dx^2 - 2*ro*cp*u/dx);
				else
					noAtual = (k/dx^2*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)) + ro*cp*u/dx*(T(i,j+1) - T(i,j-1)))/(4*k/dx^2);
				end

				noAntigo =T(i,j);
				T(i,j) = lambda*noAtual + (1-lambda)*noAntigo; %sobrerrelaxacao

				if convergindo %para melhorar desempenho e nao fazer contas desnecessarias
					if abs((T(i,j)-noAntigo))/T(i,j) > eps
						convergindo = false; %fazer mais interacoes
					end
				end
			end
		end
	end
	if convergindo %checa se pode parar de iterar
		convergiu=true;
	end
end

contador








%calculando Q
Q=0;

%percorrendo da cima para baixo
for j = 1:xf+1
	for i = fliplr(1:yf+1)
		if T(i,j)==Tdentro
			%prmeira diferente progressiva
			Q = Q -k*(T(i+1,j)-T(i,j))/dy *(dx*comprimento); %pressao*area
			break % pula para proximo x
		end
	end
end


%percorrendo da esquerda para direita
for i = 1:yf+1
	for j = 1:xf+1
		if T(i,j)==Tdentro
			%prmeira diferente progressiva
			Q = Q -k*(T(i,j-1)-T(i,j))/dx *(dy*comprimento); %pressao*area
			break % pula para proximo y
		end
	end
end


%percorrendo da direita para esquerda
for i = 1:yf+1
	for j = fliplr(1:xf+1)
		if T(i,j)==Tdentro
			%prmeira diferente progressiva
			Q = Q -k*(T(i,j+1)-T(i,j))/dx *(dy*comprimento); %pressao*area
			break % pula para proximo y
		end
	end
end

