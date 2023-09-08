addpath('./methods')
clear
data = build_sphere(0.2, 240);
data_ = data; 


Delta = [0.1,1,10,100,500];
%Delta = 1;
%K = linspace(3,8,6);
%%
%K = 3:2:15;%linspace(3,15,13); 
K = 4:3:30;
d = 2;
%mse = zeros(length(K), length(Delta),2);
%sd = mse;
%%
for k = 1:length(K)
    for j = 1:length(Delta)
        %Quadratic_K = quadratic_kernel(data, data_, K(k),  d, Delta(j));
        Quadratic_E = quadratic(data, data_, K(k), K(k),  d, Delta(j));
        %[mse(k,j,1), sd(k,j,1)] = measure_distance(Quadratic_K, projection(Quadratic_K));
        [mse(k,j,2), sd(k,j,2)] = measure_distance(Quadratic_E, projection(Quadratic_E));
        %fprintf('MSE:%.4f: ',mse(:,:,1));
        fprintf('MSE:%.4f: ',mse(:,:,2));
        fprintf('\n')
    end
end
% fprintf('MSE:%.4f: ',mse(:));
% fprintf('\n')
% fprintf('STD:%.4f ',sd);
% fprintf('\n')


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