 %% DATA CONSTRUCTION
addpath('./../methods');

%K = linspace(18,30,5);
K = 4:3:30;

data = build_sphere(0.2, 240);
d = 2; D = 3;

result = zeros(1,length(K));
sta_div = zeros(1,length(K));

data_ = data;

for i =1:length(K)
    k = K(i);
    %%MovingLS
    fprintf('MLS in Process:\n')
    MLS = MovingLS(data, data_, k, d);
    [result(i),sta_div(i)] = measure_distance(MLS.result, projection(MLS.result));
end
%scatter3(MLS.result(1,:),MLS.result(2,:),MLS.result(3,:))




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
