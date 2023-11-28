function [] = plot_eigenvalues(eigenvalues)

fig=figure;
plot(real(eigenvalues),imag(eigenvalues),".",'LineWidth',1)
hold on;
xL = xlim;
yL = ylim;
line([0 0], yL,'LineWidth',2,'Color','k');  %x-axis
line(xL, [0 0],'LineWidth',2,'Color','k');  %y-axis
viscircles([0,0],1,'Color','b');
plot(real(eigenvalues),imag(eigenvalues),"x",'LineWidth',4,'Color','r')
axis equal
grid on
xlabel("Real")
ylabel("Imag")
saveas(fig,'eigenvalues_w20.pdf')
end