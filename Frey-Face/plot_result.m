K = 8:2:30;
marker = {'d-','o--','x:','-+','-*','-s','-d','->'};
figure(10)
for i = 1:7
    semilogy(K, result(i+1,:,1),marker{i},'Linewidth',2,'MarkerSize',4)
    hold on
end
legend({'RQMF-E','RQMF-K','Local-PCA','KDE','LOG-KDE','Mfit','Moving-LS','SPH-PCA'})
xlabel('K')
ylabel('MSE')
set(gca,'FontSize',10)
axis([8 30 0.20 0.44])
