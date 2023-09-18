 %% DATA CONSTRUCTION
addpath('./../methods');

%K = linspace(18,30,5);
K = 4:3:30;

data = build_sphere(0.2, 240);
d = 2; D = 3;

result = zeros(9,6);
sta_div = zeros(9,6);

data_ = data;

for i =1:length(K)
    k = K(i);

    fprintf('PCA in process:\n')
    PCA_result = PCA_refine(data_, data, k, d);
    
    fprintf('Quadratic in Process:\n'); 
    delta = max(0.1,8*k-125);
    QuadraticE = quadratic(data, data_, k, k,  d, delta);
    delta = 500;
    kernel_k = floor(k/3+3);
    QuadraticK = quadratic_kernel(data, data_, kernel_k,  d, delta);

    %%KDE APPROACH
    fprintf('KDE in Process:\n')
    KDE = linear_KDE(k, data, data_, D, d);

    %%LOG-KDE APPROACH
    fprintf('log_KDE in Process:\n')
    log_KDE = linear_log_KDE(k, data, data_, D, d);

    %%MFIT APPROACH
    fprintf('mfint in Process:\n')
    Mfit = linear_mfit(data, data_, d, k);

    %%MovingLS
    fprintf('MLS in Process:\n')
    MLS = MovingLS(data, data_, k, d);

    %%Spherelet
    fprintf('Spherelet:\n')
    Sph = Spherelet(data_, data, k, d);

    [result(1,i),sta_div(1,i)] = measure_distance(data_, projection(data_));
    [result(2,i),sta_div(2,i)] = measure_distance(QuadraticE, projection(QuadraticE));
    [result(3,i),sta_div(3,i)] = measure_distance(QuadraticK, projection(QuadraticK));
    [result(4,i),sta_div(4,i)] = measure_distance(PCA_result,projection(PCA_result));
    [result(5,i),sta_div(5,i)] = measure_distance(KDE, projection(KDE));
    [result(6,i),sta_div(6,i)] = measure_distance(log_KDE,projection(log_KDE));
    [result(7,i),sta_div(7,i)] = measure_distance(Mfit, projection(Mfit));
    [result(8,i),sta_div(8,i)] = measure_distance(MLS.result, projection(MLS.result));
    [result(9,i),sta_div(9,i)] = measure_distance(Sph, projection(Sph));
end
%% output
RE = zeros(9,16);
for i = 1:8
    RE(:,2*i-1) = result(:,i+1);
    RE(:,2*i) = sta_div(:,i+1);
end
format_print(RE)



function re = projection(A)
    re = bsxfun(@rdivide,A,sqrt(sum(A.^2,1)));
end


function [dis, s] = measure_distance(A, T)
    dis = norm(A-T,'fro')^2/size(A,2);
    S = A-T;
    s = std(sum(S.^2,1),1);
end


function data = build_sphere(sigma, num)

    data = randn(3,num);
    data = bsxfun(@rdivide,data, sqrt(sum(data.^2,1)));
    data = data + sigma*randn(size(data));
end
