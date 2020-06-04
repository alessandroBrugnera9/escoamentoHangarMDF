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
patm=10000;
comprimento=60; %fundura

Tdentro=40;
Tfora=20;


%variaveis computacionais
dx=0.1
dy=0.1
lambda=1.85
eps=0.01



xi=2;
xf=(2*d+L)/dx-1;
yi=2;
yf=H/dy-1;


corr=zeros(yf+1,xf+1);
pressao=zeros(yf+1,xf+1);
deltap=zeros(yf+1,xf+1);


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
		if ((i)*dy) <= (sqrt((L/2)^2 - ((j-1)*dx-d-L/2)^2)+h)
			corr(i,j) = 0;
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


	%primeiro nas pontas superiores
	corr(yf+1,1)=(corr(yf+1,2)+corr(yf,1)+V*dy)/2;
	corr(yf+1,xf+1)=(corr(yf+1,xf-1)+corr(yf,xf+1)+V*dy)/2;

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
			if  (i*dy) > (sqrt((L/2)^2 - ((j-1)*dx-d-L/2)^2)+h) %checa se esta nos limites externos do telhado
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

contador








%calculando pressao utilizando primeira diferenca central
for i = yi:yf
	for j = xi:xf %antes do galpao
		u=(corr(i+1,j)-corr(i-1,j))/(2*dy);
		v=(corr(i,j+1)-corr(i,j-1))/(2*dx);
		pressao(i,j) = -ro*(gama-1)/gama*(u^2+v^2)/2 + patm;
	end
end

%calculando nas bordas com diferenças progressivas e regressivas
for j = xi:xf
	%em cima
	u=(corr(yf+1,j)-corr(yf,j))/(dy);
	v=(corr(yf+1,j+1)-corr(yf+1,j-1))/(2*dx);
	pressao(yf+1,j)=-ro*(gama-1)/gama*(u^2+v^2)/2 + patm;

	%embaixo
	u=(corr(2,j)-corr(1,j))/(dy);
	v=(corr(1,j+1)-corr(1,j-1))/(2*dx);
	pressao(1,j)=-ro*(gama-1)/gama*(u^2+v^2)/2 + patm;
end

for i = yi:yf
	%lado esquerdo
	u=(corr(i+1,1)-corr(i-1,1))/(2*dy);
	v=(corr(i,2)-corr(i,1))/(dx);
	pressao(i,1)=-ro*(gama-1)/gama*(u^2+v^2)/2 + patm;

	%lado direito
	u=(corr(i+1,xf+1)-corr(i-1,xf+1))/(2*dy);
	v=(corr(i,xf+1)-corr(i,xf))/(dx);
	pressao(i,xf+1)=-ro*(gama-1)/gama*(u^2+v^2)/2 + patm;
end

%calculando nas pontas
%embaixo esquerda
u=(corr(2,1)-corr(1,1))/(dy);
v=(corr(1,2)-corr(1,1))/(dx);
pressao(1,1)=-ro*(gama-1)/gama*(u^2+v^2)/2 + patm;

%embaixo direita
u=(corr(2,xf+1)-corr(1,xf+1))/(dy);
v=(corr(1,xf+1)-corr(1,xf))/(dx);
pressao(1,xf+1)=-ro*(gama-1)/gama*(u^2+v^2)/2 + patm;

%em cima esquerda
u=(corr(yf+1,1)-corr(yf,1))/(dy);
v=(corr(yf+1,2)-corr(yf+1,1))/(dx);
pressao(yf+1,1)=-ro*(gama-1)/gama*(u^2+v^2)/2 + patm;

%em cima direita
u=(corr(yf+1,xf+1)-corr(yf,xf+1))/(dy);
v=(corr(yf+1,xf+1)-corr(yf+1,xf))/(dx);
pressao(yf+1,xf+1)=-ro*(gama-1)/gama*(u^2+v^2)/2 + patm;









%variacao de pressao
for i = 1:yf+1
	for j = 1:xf+1
		deltap(i,j) = pressao(i,j) - patm;
	end
end



maiory=0;
maiorx=0;
menorPressaoTelhado=Inf;
telhadoP=zeros(yf+1,xf+1);
for i = gyf:yf+1
	for j = 1:xf+1
		if deltap(i,j)==0 %checa se é um ponto interno do telhado
			%testa se alguma do vizinhanca é o telhado vendo se é diferente de 0 e atribui a matriz do telhado
			%algum valores serao subscritos neles mesmos
			if deltap(i-1,j)==0
				telhadoP(i-1,j) = deltap(i-1,j);
			end
			if deltap(i,j-1)
				telhadoP(i,j-1) = deltap(i,j-1);
			end
			if deltap(i+1,j)
				telhadoP(i+1,j) = deltap(i+1,j);
			end
			if deltap(i,j+1)
				telhadoP(i,j+1) = deltap(i,j+1);
			end
		end
	end
end







%como o telhado eh simetrico só há uma parcela em y de forca
forca=0;

%percorre em x de cima para baixo ate encontrar o telhado
for j = 1:xf+1
	for i = fliplr(1:yf+1)
		if telhadoP(i,j)~=0
			seno=(i-gyi)/sqrt(((i-gyi))^2+(j-xf/2)^2);
			forca = forca + (telhadoP(i,j)*seno)*(dx*comprimento*seno); %pressao*area
			break % pula para proximo x
		end
	end
end

forca










































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
	noAtual=(T(1,xf+1)+T(2,xf))/2
	noAntigo =T(i,j);
	T(1,xf+1) = lambda*noAtual + (1-lambda)*noAntigo; %sobrerrelaxacao
	%segundo no topo
	u=V;
	for j = xi:xf
		noAtual=(k/dx*(2*T(yf,j) +  T(yf+1,j+1)+T(yf+1,j-1)) + ro*cp*u/dx*(T(yf+1,j-1)))/(4*k/dx + ro*cp*u/dx);
		noAntigo =T(i,j);
		T(yf+1,j) = lambda*noAtual + (1-lambda)*noAntigo; %sobrerrelaxacao
	end

	%segundo emabvixo
	for j = xi:xf
		u=(corr(i+1,j)-corr(i-1,j))/(2*dy);
		noAtual=(k/dx*(2*T(2,j) +  T(1,j+1)+T(1,j-1)) + ro*cp*u/dx*(T(1,j-1)))/(4*k/dx + ro*cp*u/dx);
		noAntigo =T(i,j);	
		T(1,j) = lambda*noAtual + (1-lambda)*noAntigo; %sobrerrelaxacao
	end

	%terceito a direita
	for i = yi:yf
		u=(corr(i+1,xf+1)-corr(i-1,xf+1))/(2*dy);
		v=0;
		noAtual=(k/dx*(2*T(i,xf) +  T(i+1,xf+1)+T(i-1,xf+1)) + ro*cp*u/dx*(T(i-1,xf+1)))/(4*k/dx + ro*cp*u/dx);
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
				noAtual = (k/dx*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)) - ro*cp*u/dx*(2*T(i,j+1)))/(4*k/dx + 2*ro*cp*u/dx);
			elseif u>0 & v>0
				noAtual = (k/dx*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)) + ro*cp*u/dx*(2*T(i,j-1)))/(4*k/dx - 2*ro*cp*u/dx);
			else
				noAtual = (k/dx*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)) + ro*cp*u/dx*(T(i,j+1) - T(i,j-1)))/(4*k/dx);
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
				noAtual = (k/dx*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)) - ro*cp*u/dx*(2*T(i,j+1)))/(4*k/dx + 2*ro*cp*u/dx);
			elseif u>0 & v>0
				noAtual = (k/dx*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)) + ro*cp*u/dx*(2*T(i,j-1)))/(4*k/dx - 2*ro*cp*u/dx);
			else
				noAtual = (k/dx*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)) + ro*cp*u/dx*(T(i,j+1) - T(i,j-1)))/(4*k/dx);
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
					noAtual = (k/dx*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)) - ro*cp*u/dx*(2*T(i,j+1)))/(4*k/dx + 2*ro*cp*u/dx);
				elseif u>0 & v>0
					noAtual = (k/dx*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)) + ro*cp*u/dx*(2*T(i,j-1)))/(4*k/dx - 2*ro*cp*u/dx);
				else
					noAtual = (k/dx*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)) + ro*cp*u/dx*(T(i,j+1) - T(i,j-1)))/(4*k/dx);
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