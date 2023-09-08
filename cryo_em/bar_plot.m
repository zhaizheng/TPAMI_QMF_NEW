figure
%t = tiledlayout(1,2,'TileSpacing','Compact');
%nexttile
mse = result;
%s = [1,3,5,7,9,11];
s = 4:13;
h = bar3(squeeze(mse(s,:,1))-1);
for n=1:numel(h)
    cdata=get(h(n),'zdata');
    set(h(n),'cdata',cdata,'facecolor','interp')
end
%Neig = [12,16,20,24,28,32,36,40,44,48,52,56,60];
%set(gca,'YTickLabel',{'7','10','13','16','19','22','25','28'});
set(gca,'YTickLabel',{'24','28','32','36','40','44','48','52','56','60'});
%set(gca,'YTickLabel',{'28','32','36','40','44','48','52','56','60'});

% %set(gca,'XTickLabel',{'0.1','1','10','100','500'});
% % [0.1, 1, 5,10, 50, 100, 300]
  set(gca,'XTickLabel',{'0.1','1','5','10','50','100','300'});
%  axis([0 8 3 14 0 1.8])
%  ylabel('K','Interpreter','latex')
%  xlabel('\delta')
%  zlabel('MSE')
%  set(gca,'FontSize',18)
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