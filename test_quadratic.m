data = build_sphere(0.1, 240);
m = 20;
data_ = data(:,1:m);


delta = 100; d = 2; k = 20; n = 14;
P = quadratic(data, data_, k, n, d, delta);


A = P-P*diag(1./sqrt(sum(P.^2,1)));
fprintf('error=%f\n',norm(A,'fro')/size(A,2));
%fprintf('measure2=%f\n',mean(A(:).^2));


function data = build_sphere(sigma, num)

    data = randn(3, num);
    %data = data*diag(1./sqrt(sum(data.^2.1)));
    data = bsxfun(@rdivide,data, sqrt(sum(data.^2,1)));
    data = data + sigma*randn(size(data));
end