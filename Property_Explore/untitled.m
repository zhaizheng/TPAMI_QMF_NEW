K = 6:3:30;
data = cell(1,5);
data_ = cell(1,5);
for i = 1:5
    data{i} = build_sphere(0.1*i, 240);
    data_{i} = data{i};
end
%%
QuadraticE = cell(length(K),5);
Rho = [0, 0.01, 0.05, 0.25, 0.5, 1];
Error = cell(1,5);%zeros(length(K),5);
Time = cell(1,5);%zeros(length(K),5);

d = 2; D = 3; delta = 0; rho = 0;
for s = 1:length(Rho)
    rho = Rho(s);
    Error{s} = zeros(length(K),5);
    Time{s} = zeros(length(K),5);
    for i = 1:length(K)
        k = K(i);
        for j = 1:5
            t0 = clock;
            QuadraticE{i,j} = quadratic(data{j}, data_{j}(:,1:20), k, k,  d, delta, rho);
            eltime = etime(clock, t0);
            Time{s}(i,j) = eltime;
            Error{s}(i,j) = measure_distance(QuadraticE{i,j}, projection(QuadraticE{i,j}));
        end
    end
end
%%
Error(i,j) = measure_distance(QuadraticE{i,j}, projection(QuadraticE{i,j}));

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