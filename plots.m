 xt = get(gca, 'XTick');
 set(gca, 'XTick', xt, 'XTickLabel', xt*dx)
 yt = get(gca, 'YTick');
set(gca, 'YTick', yt, 'YTickLabel', yt*dy)
xlabel('Eixo x [m]')
ylabel('Eixo y [m]')

contourf(corr,30)
title('Função Corrente','FontSize',20);
c = colorbar('FontSize',20);
 c.Label.String='\Psi';
xt = get(gca, 'XTick');
 set(gca, 'XTick', xt, 'XTickLabel', xt*dx)
 yt = get(gca, 'YTick');
set(gca, 'YTick', yt, 'YTickLabel', yt*dy)
xlabel('Eixo x [m]')
ylabel('Eixo y [m]')
 saveas(gcf,'corr1.jpg')


surfc(deltap)
title('Variação de pressão','FontSize',20)
c = colorbar('FontSize',20);
c.Label.String='\Delta p [Pa]';
 xt = get(gca, 'XTick');
 set(gca, 'XTick', xt, 'XTickLabel', xt*dx)
 yt = get(gca, 'YTick');
set(gca, 'YTick', yt, 'YTickLabel', yt*dy)
xlabel('Eixo x [m]')
ylabel('Eixo y [m]')
saveas(gcf,'1c1.png')

 
 
 
 surfc(-telhadoP)
 title('Variação de pressão no telhado','FontSize',20)
 c = colorbar('FontSize',20);
 c.Label.String='\Delta p [Pa]';
  xt = get(gca, 'XTick');
 set(gca, 'XTick', xt, 'XTickLabel', xt*dx)
 yt = get(gca, 'YTick');
set(gca, 'YTick', yt, 'YTickLabel', yt*dy)
xlabel('Eixo x [m]')
ylabel('Eixo y [m]')

 saveas(gcf,'1d1.jpg')
















  xt = get(gca, 'XTick');
 set(gca, 'XTick', xt, 'XTickLabel', xt*dx)
 yt = get(gca, 'YTick');
set(gca, 'YTick', yt, 'YTickLabel', yt*dy)
xlabel('Eixo x [m]')
ylabel('Eixo y [m]')

contourf(corr,30)
title('Função Corrente','FontSize',20);
c = colorbar('FontSize',20);
 c.Label.String='\Psi';
xt = get(gca, 'XTick');
 set(gca, 'XTick', xt, 'XTickLabel', xt*dx)
 yt = get(gca, 'YTick');
set(gca, 'YTick', yt, 'YTickLabel', yt*dy)
xlabel('Eixo x [m]')
ylabel('Eixo y [m]')
 saveas(gcf,'corr0.5.jpg')


surfc(deltap)
title('Variação de pressão','FontSize',20)
c = colorbar('FontSize',20);
c.Label.String='\Delta p [Pa]';
 xt = get(gca, 'XTick');
 set(gca, 'XTick', xt, 'XTickLabel', xt*dx)
 yt = get(gca, 'YTick');
set(gca, 'YTick', yt, 'YTickLabel', yt*dy)
xlabel('Eixo x [m]')
ylabel('Eixo y [m]')
saveas(gcf,'1c05.png')

 
 
 
 surfc(-telhadoP)
 title('Variação de pressão no telhado','FontSize',20)
 c = colorbar('FontSize',20);
 c.Label.String='\Delta p [Pa]';
  xt = get(gca, 'XTick');
 set(gca, 'XTick', xt, 'XTickLabel', xt*dx)
 yt = get(gca, 'YTick');
set(gca, 'YTick', yt, 'YTickLabel', yt*dy)
xlabel('Eixo x [m]')
ylabel('Eixo y [m]')

 saveas(gcf,'1d05.jpg')