figure
%t = tiledlayout(1,2,'TileSpacing','Compact');
%nexttile
%s = [1,3,5,7,9,11];
s = 2:9;
h = bar3(squeeze(mse(s,:,2)));
for n=1:numel(h)
    cdata=get(h(n),'zdata');
    set(h(n),'cdata',cdata,'facecolor','interp')
end
set(gca,'YTickLabel',{'7','10','13','16','19','22','25','28'});
set(gca,'XTickLabel',{'0.1','1','10','100','500'});
axis([0 6 0 10 0 0.04])
ylabel('K','Interpreter','latex')
xlabel('\delta')
zlabel('MSE')
set(gca,'FontSize',18)
%nexttile
% figure
% h = bar3(squeeze(mse(s,:,1)));
% for n=1:numel(h)
%     cdata=get(h(n),'zdata');
%     set(h(n),'cdata',cdata,'facecolor','interp')
% end
% set(gca,'YTickLabel',{'3','5','7','9','11','13'});
% set(gca,'XTickLabel',{'0.1','1','10','100','500'});
% axis([0 6 0 7 0 0.04])
% ylabel('K','Interpreter','latex')
% xlabel('\delta')
% zlabel('MSE')
% set(gca,'FontSize',18)