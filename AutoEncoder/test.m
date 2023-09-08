clear
cd('/home/zzh/Documents/TPAMI_QMF_new/code/AutoEncoder')
addpath('/home/zzh/Documents/TPAMI_QMF_new/code/methods')
addpath('/home/zzh/Downloads/Pytorch-VAE-tutorial-master')
addpath('/home/zzh/Downloads/Pytorch-VAE-tutorial-master/npy-matlab/')
X1 = readNPY('Noisy.npy');
X2 = readNPY('Signal.npy');

d = 1;
delta = 2;
k = 5;
ind = 1:100:2000;
X1S = X1(:,ind);
Quadratic_K = quadratic_kernel(X1,X1S,k,d,delta);

%%
norm(Quadratic_K-X2(:,ind),'fro')^2/length(ind)

% label = [];
% for i = 1:10
%     label = [label,i*ones(1,200)];
% end
% %c = label(ind);
% 
% scatter3(Quadratic_K(1,:),Quadratic_K(2,:),Quadratic_K(3,:),[],label)

