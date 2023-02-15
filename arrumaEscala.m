function arrumaEscala()
global dx dy
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel', xt*dx)
yt = get(gca, 'YTick');
set(gca, 'YTick', yt, 'YTickLabel', yt*dy)
xlabel('Eixo x [m]')
ylabel('Eixo y [m]')
end