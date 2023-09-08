0%clear
cd('/home/zzh/Documents/TPAMI_QMF_new/code/AutoEncoder')
addpath('/home/zzh/Downloads/Pytorch-VAE-tutorial-master')
addpath('/home/zzh/Downloads/Pytorch-VAE-tutorial-master/npy-matlab/')
X1 = readNPY('Noisy.npy');
X2 = readNPY('Signal.npy');
y = [];
for i = 1:10
    y = [y,i*ones(1,200)];
end


d = 1;
%k = 15;
delta = 5;
Delta = [0.1,1,2,5,10];
K = [5,10,15,20,25];
re = zeros(5,5);
%X3 = cell(5,5);

% for i = 1:5
%     for j = 1:5
%         X3{i,j} = quadratic_kernel(X1, X1, K(i), d, Delta(j));
%         %%
%         re(i,j) = norm(X2-X3{i,j},'fro');
%     end
% end
% [r,c] = find(re==min(re(:)));
X3 = quadratic_kernel(X1, X1, K(2), d, Delta(2));

%%
t = tiledlayout(1,4,'TileSpacing','compact');
nexttile
scatter3(X2(1,:),X2(2,:),X2(3,:),[],y)
nexttile
scatter3(X1(1,:),X1(2,:),X1(3,:),[],y,'filled')
axis([-3 3 -3 3 -3 4])
nexttile
%scatter3(X3{r,c}(1,:),X3{r,c}(2,:),X3{r,c}(3,:),[],y,'filled')
scatter3(X3(1,:),X3(2,:),X3(3,:),[],y,'filled')
axis([-3 3 -3 3 -3 4])
%nexttile
%scatter3(X3{4,3}(1,:),X3{4,3}(2,:),X3{4,3}(3,:),[],y,'filled')
%axis([-3 3 -3 3 -3 4])

%%
%W = X3{r,c};
%save 'After_process.mat' W
