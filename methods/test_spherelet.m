X = build_sphere(0.1, 240);
d = 2; D = 3;
K = 50;
for i = 1:length(K)
    k = K(i);
    [Output, R, C] = Spherelet(X, X, k, d);
    [result(i),~] = measure_distance(Output, projection(Output));
end
result


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