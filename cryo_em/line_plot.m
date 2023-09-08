figure
for i = 5:8
    plot(result(i,:,1),'o-','Linewidth',1.5)
    hold on
end
for i = 9:13
    plot(result(i,:,1),'d-.','Linewidth',1.5)
    hold on
end

str = {'K=28','K=32','K=36','K=40','K=44','K=48','K=52','K=56','K=60'};
h = legend(str);
set(h,'box','off')
set(gca,'FontSize',18)
xlabel('\delta')
ylabel('MSE')
set(gca,'XTickLabel',{'0.1','1','5','10','50','100','300'});