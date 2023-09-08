%%
clear
cd('/home/zzh/Documents/TPAMI_QMF_new/code/AutoEncoder')
addpath('/home/zzh/Documents/TPAMI_QMF_new/code/methods')
addpath('/home/zzh/Downloads/Pytorch-VAE-tutorial-master')
addpath('/home/zzh/Downloads/Pytorch-VAE-tutorial-master/npy-matlab/')
X1 = readNPY('Noisy.npy');
X2 = readNPY('Signal.npy');


%%
kk = [8 10 12 14 16 18 20];
n_sample = size(X1,2);
D = size(X1,1);
n_interest = 100;
ind  =  randperm(n_sample,n_interest);
target = X1(:,ind);
XT = X2(:,ind);
result = zeros(9,5,2);
d = 1;
delta = 2;
%%
for i = 1:7
    k = kk(i);
    %PCA DENOSING
    fprintf('PCA in process:\n')
    PCA = PCA_refine_centralized(double(target), double(X1), k, d);

    %%QUADRATIC APPROACH
    fprintf('Quadratic in Process:\n')
    Quadratic_R = quadratic(X1,target,k, k, d, delta);
    Quadratic_K = quadratic_kernel(X1,target,k,d,delta);

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

     fprintf('Spherelet in Process:\n')
     Sphere = Spherelet(target, X1, k, d);
    % 
    % 
    [result(1,i,1),result(1,i,2)] = measure_distance(XT, target);
    [result(2,i,1),result(2,i,2)] = measure_distance(XT, Quadratic_R);
    [result(3,i,1),result(3,i,2)] = measure_distance(XT, Quadratic_K);
    [result(4,i,1),result(4,i,2)] = measure_distance(XT, PCA);
    [result(5,i,1),result(5,i,2)] = measure_distance(XT, KDE_Fitting);
    [result(6,i,1),result(6,i,2)] = measure_distance(XT, log_KDE_Fitting);
    [result(7,i,1),result(7,i,2)] = measure_distance(XT, mfit_Fitting);
    [result(8,i,1),result(8,i,2)] = measure_distance(XT, MLS_Fitting.result);
    [result(9,i,1),result(9,i,2)] = measure_distance(XT, Sphere);
end
%%
Output = [];
for i = 1:9
    Output = [Output,[result(i,:,1)';result(i,:,2)']];
end
format_print(Output')
%%
label = [];
for i = 1:10
    label = [label,i*ones(1,200)];
end
c = label(ind);
subplot(1,3,1)
%scatter3(Quadratic_R(1,:),Quadratic_R(2,:),Quadratic_R(3,:),[],c)
scatter3(PCA(1,:),PCA(2,:),PCA(3,:),[],c)

subplot(1,3,2)
scatter3(XT(1,:),XT(2,:),XT(3,:),[],c)

subplot(1,3,3)
scatter3(X1(1,:),X1(2,:),X1(3,:),[],label)

function [dis,sd] = measure_distance(A, T)
    dis = norm((A(:)-T(:)),'fro')^2/size(A,2);
    s = [];
    for i = 1:size(A,2)
        s(i) = norm(squeeze(A(:,i))-squeeze(T(:,i)),'fro')^2;
    end
    sd = std(s,1); 
end