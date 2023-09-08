%% this file need to run after the empliment of adaptive_fitting.m

for j = 1:3
    subplot(1,3,j)
    for i = 1:size(data,2)
        if mod(i,30)==0
            hold on
            plot(data_new{2*j-1,i}(1,:),data_new{2*j-1,i}(2,:),'.','markersize',8);
            hold on
            text(P{2*j-1}(1,i)+0.1,P{2*j-1}(2,i),num2str(Rho(2*j-1,i),'%.3f'),'FontSize', 14)
        end
    end
    axis([-1.5 1.5 -1.5 1.5])
end

figure
for j = 1:3
    subplot(1,3,j)
    %hold on
    plot(data_t(1,:), data_t(2,:),'-','linewidth',2)
    hold on
    plot(P{2*j-1}(1,:), P{2*j-1}(2,:),'.','markersize',8)
    %axis([-1 1 -1 1])
end

%%
figure
plot(data(1,:),data(2,:),'*')
hold on
plot(data_t(1,:),data_t(2,:),'-','linewidth',2)

box off
axis off

%%
%%figure                 
%plot(Delta,recover_result(3,:),'-o')

