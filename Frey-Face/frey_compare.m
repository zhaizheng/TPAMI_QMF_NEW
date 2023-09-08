%%
clear
cd('/home/zzh/Documents/TPAMI_QMF_new/code/Frey-Face')
addpath('/home/zzh/Documents/TPAMI_QMF_new/code/methods')
addpath('/home/zzh/Downloads/Pytorch-VAE-tutorial-master/Fray_Face/')
addpath('/home/zzh/Downloads/Pytorch-VAE-tutorial-master/npy-matlab/')
T = readNPY('Signal_Frey.npy');
X1 = readNPY('Noise_Frey.npy');

kk = 8:2:30;
d = 1;
D = 3;
delta = 2;
n_sample = size(X1,2);
n_interest = 100;
ind  =  randperm(n_sample,n_interest);
target = X1(:,ind);
result = zeros(9,length(kk),2);
XT = T(:,ind);

%%
for i = 1:length(kk)
    k = kk(i);
    %PCA DENOSING
    fprintf('PCA in process:\n')
    PCA = PCA_refine_centralized(double(target), double(X1), k, d);

    %%QUADRATIC APPROACH
    fprintf('Quadratic in Process:\n')
    Quadratic_R = quadratic(X1, target, k, k, d, delta);
    Quadratic_K = quadratic_kernel(X1, target,k,d,delta);

    %%KDE APPROACH
    fprintf('KDE in Process:\n')
    KDE_Fitting = linear_KDE(k, X1, target, D, d);

    %%LOG-KDE APPROACH
    fprintf('log_KDE in Process:\n')
    log_KDE_Fitting = linear_log_KDE(k, X1, target, D, d);

    %%MFIT APPROACH
    fprintf('mfint in Process:\n')
    mfit_Fitting = linear_mfit(X1, target, d, k);

     %MovingLS
     fprintf('MLS in Process:\n')
     MLS_Fitting = MovingLS(X1, target, k/2-2, d);

     %Spherelet
     fprintf('Sphere in Process:\n')
     Sphere_fitting = Spherelet(target,X1,k,d);

    [result(1,i,1),result(1,i,2)] = measure_distance(XT, target);
    [result(2,i,1),result(2,i,2)] = measure_distance(XT, Quadratic_R);
    [result(3,i,1),result(3,i,2)] = measure_distance(XT, Quadratic_K);
    [result(4,i,1),result(4,i,2)] = measure_distance(XT, PCA);
    [result(5,i,1),result(5,i,2)] = measure_distance(XT, KDE_Fitting);
    [result(6,i,1),result(6,i,2)] = measure_distance(XT, log_KDE_Fitting);
    [result(7,i,1),result(7,i,2)] = measure_distance(XT, mfit_Fitting);
    [result(8,i,1),result(8,i,2)] = measure_distance(XT, MLS_Fitting.result);
    [result(9,i,1),result(9,i,2)] = measure_distance(XT, Sphere_fitting);
end
%%
[val, ind2] = sort(ind,'ascend');
XTh = XT(:,ind2);
y = [2*ones(1,50),3*ones(1,50)];
t = tiledlayout(1,3,'TileSpacing','Compact');
nexttile
scatter3(XTh(1,:),XTh(2,:),XTh(3,:),[],XTh(1,:),'filled');
axis([-3 3 -4 4 -2 2])
nexttile
scatter3(X1(1,val),X1(2,val),X1(3,val),[],XTh(1,:),'filled');
axis([-3 3 -4 4 -2 2])
%hold on
%scatter3(Quadratic_K(1,:),Quadratic_K(2,:),Quadratic_K(3,:),[],'*');
nexttile
scatter3(Quadratic_K(1,ind2),Quadratic_K(2,ind2),Quadratic_K(3,ind2),[],XTh(1,:),'filled');
axis([-3 3 -4 4 -2 2])
% subplot(2,2,3)
% scatter3(Quadratic_K(1,:),Quadratic_K(2,:),Quadratic_K(3,:),[],'*');
% subplot(2,2,4)
% scatter3(Quadratic_R(1,:),Quadratic_R(2,:),Quadratic_R(3,:),[],'*');
%%
[val, ind2] = sort(ind,'ascend');
W_After = Quadratic_K(:,ind2);
W_True = XT(:,ind2);
W_Orig = target(:,ind2);
save Data_to_show.mat W_Orig W_True W_After 

%%
function [dis,sd] = measure_distance(A, T)
    dis = norm((A(:)-T(:)),'fro')^2/size(A,2);
    s = [];
    for i = 1:size(A,2)
        s(i) = norm(squeeze(A(:,i))-squeeze(T(:,i)),'fro')^2;
    end
    sd = std(s,1); 
end


