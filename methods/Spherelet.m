function Output = Spherelet(DATA_TOBE, Data, k, d)
    Output = zeros(size(DATA_TOBE));
    for i = 1:size(Output,2)
        [~, ind] = sort(sum((Data-DATA_TOBE(:,i)).^2,1),'ascend');
        S = Data(:,ind(2:k+1));
        c1 = mean(S,2);
        [U,~,~] = svds(double(S-c1), d+1);
        Y = c1+U*U'*(S-c1);
        [c,r] = compute_spherelet_center(Y);
        PU = U*U'*(DATA_TOBE(:,i)-c);
        res = c + r*PU./(norm(PU)+eps);
        if isnan(res(1))
            pause;
        end
        Output(:,i) = res;
    end
end

function [c,r] = compute_spherelet_center(Y)
    center = mean(Y,2);
    H = (Y-center)*(Y-center)';
    center2 = mean(diag(Y'*Y));
    xi = zeros(size(Y,1),1);
    for j = 1:size(Y,2)
        xi = xi + (Y(:,j)'*Y(:,j)-center2)*(Y(:,j)-center);
    end
    c = pinv(H)*xi/2;
    r = 0;
    n = size(Y,2);
    for j = 1:size(Y,2)
        r = r+norm(Y(:,j)-c)/n;
    end
    if isnan(r)
        pause;
    end
end