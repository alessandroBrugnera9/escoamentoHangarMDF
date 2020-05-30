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
dx=0.1
dy=0.1
lambda=1.85
eps=1



xi=2;
xf=(2*d+L)/dx-1;
yi=2;
yf=H/dy-1;


corr=zeros(yf+1,xf+1);

%popula com 1 o que tem dentro
for i = yi:yf
	for j = xi:xf
		corr(i,j) = 50*rand();
	end
end

gxi=d/dx+1; %x esquerda
gxf=gxi+(L)/dx; %x direita
gyi=yi; %
gyf=gyi+h/dy-1; %

%popula com 2 o que tem dentro do galpao
for i = gyi:gyf
	for j = gxi:gxf
		corr(i,j) = 0;
	end
end

%popula com 3 o que tem no topo do galpao
for i = gyf:(gyf+(L/2)/dy-1)
	for j = gxi:gxf
		if (i*dy) <= (sqrt((L/2)^2 - (j*dx-d-L/2)^2)+h)
			corr(i,j) = 0;
		end
	end
end


contador=0;
convergiu=false;

while ~convergiu
	contador = contador+1
	convergindo=true;

	%aplicando condicoes nas bordas utilizando equacionamento de taylor
	%para considerar condicoes de Neumann


	%primeiro nas pontas superiores
	corr(yf+1,1)=(corr(yf+1,2)+corr(yf,1)+V*dy)/2
	corr(yf+1,xf+1)=(corr(yf+1,xf-1)+corr(yf,xf+1)+V*dy)/2

	%segundo no topo
	for j = xi:xf
		corr(yf+1,j)=(dx*V + corr(yf,j) + (corr(yf+1,j-1)+corr(yf+1,j+1))/2)/2;
	end

	%segundo nos lados
	for i = yi:yf
		corr(i,1)=(corr(i,2) + (corr(i-1,1)+corr(i+1,1))/2)/2; %esquerdo
		corr(i,xf+1)=(corr(i,xf) + (corr(i-1,xf+1)+corr(i+1,xf+1))/2)/2; %direito
	end


	%Aplicando o metrodo nos nos internos
	%primeiro na altura inferior ao topo do telhado
	for i = gyi:gyf
		for j = xi:gxi-1 %antes do galpao
			noAntigo =corr(i,j);
			noAtual = (corr(i+1,j)+corr(i-1,j)+corr(i,j+1)+corr(i,j-1))/4;
			corr(i,j) = lambda*noAtual + (1-lambda)*noAntigo; %sobrerrelaxacao

			if convergindo %para melhorar desempenho e nao fazer contas desnecessarias
				if abs((corr(i,j)-noAntigo))/corr(i,j) > eps
					convergindo = false; %fazer mais interacoes
				end
			end
		end

		for j = gxf+1:xf %depois do galpao
			noAntigo =corr(i,j);
			noAtual = (corr(i+1,j)+corr(i-1,j)+corr(i,j+1)+corr(i,j-1))/4;
			corr(i,j) = lambda*noAtual + (1-lambda)*noAntigo; %sobrerrelaxacao
			
			if convergindo %para melhorar desempenho e nao fazer contas desnecessarias
				if abs((corr(i,j)-noAntigo))/corr(i,j) > eps
					convergindo = false; %fazer mais interacoes
				end
			end
		end
	end


	%segundo na altura superior ao topo do telhado
	for i = gyf+1:yf
		for j = xi:xf
			if  (i*dy) > (sqrt((L/2)^2 - (j*dx-d-L/2)^2)+h) %checa se esta nos limites externos do telhado
				noAntigo =corr(i,j);
				noAtual = (corr(i+1,j)+corr(i-1,j)+corr(i,j+1)+corr(i,j-1))/4;
				corr(i,j) = lambda*noAtual + (1-lambda)*noAntigo; %sobrerrelaxacao

				if convergindo %para melhorar desempenho e nao fazer contas desnecessarias
					if abs((corr(i,j)-noAntigo))/corr(i,j) > eps
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
