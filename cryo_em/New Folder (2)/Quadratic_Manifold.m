

function points_projection = Quadratic_Manifold(k, Data, points, D, d, tol)
       
        points_projection = zeros(size(points));
        for i = 1:size(points,2)
            fprintf('processing number %d\n',i);
            data_old = zeros(size(points(:,i)));
            data_move = points(:,i);
            sk = 1;
            while norm(data_old-data_move)>tol && sk<3
                data_old = data_move;
                sigma = find_sigma(data_old, Data, k);
                %x = shift_mean(data_move, Data, sigma, k);% try another approach
                [U, center] = principal(Data, sigma, data_move, d);
                result = Least_Fitting(Data, U, center, sigma); 
%                 
                U = search_principal_space(Data, result, sigma, d);
                result = Least_Fitting(Data, U, result, sigma); 
                
                U = search_principal_space(Data, result, sigma, d);
                result = Least_Fitting(Data, U, result, sigma); 
                
%                [result, ~, ~, U] = find_origin(data_move, Data, sigma, d);
%                [Theta, ~, ~, U_h, U_v] = fitting2(result, Data, sigma, D, d, U);
                [Theta, ~, ~, U_h, U_v] = fitting(result, Data, sigma, D, d);
                [~, data_move, ~] = nonlinear_solution(U_h, U_v, points(:,i), result, Theta);
                sk = sk+1;
            end
            points_projection(:,i) = data_move;
        end
end


function [q, k, Q, C] = find_origin(x, Data, h, d)
    q = x;
    U = principal(Data, h, q, d);
    k = 0;
    Q = [];
    while k<100
        Data_c = Data-repmat(q,1,size(Data,2));
        Theta = build_theta(Data, h, q);
        R = Data_c*Theta;
        Tau = U(:,1:d)'*R;
        Tau1 = [ones(1,size(Tau,2));Tau];
        A = R*Tau1'/(Tau1*Tau1');
        q_temp = q + A(:,1);
        [U,~] = qr(A(:,2:end));
        q_new = q_temp+U(:,1:d)*U(:,1:d)'*(x-q_temp);
        if norm(q_new-q)<0.001
            break;
        end
        q = q_new;
        Q = [Q,q];
        k = k+1;
    end
    C = U(:,1:d);
end


function [U,center] = principal(Data, h, q, d)
    Theta = (build_theta(Data, h, q)).^2;
    center = sum(Data*Theta, 2)/sum(diag(Theta));
    [V,~,~] = svd((Data-center)*Theta*(Data-center)');
    U = V(:,1:d);
end


function U = search_principal_space(Data, center, h, d)
    Theta = (build_theta(Data, h, center)).^2;
    [V,~,~] = svd((Data-center)*Theta*(Data-center)');
    U = V(:,1:d);
end


function result = Least_Fitting(Data, U, q, h)
    Tau = U'*(Data-q);
    T = Construct_Higher_Order(Tau);
    Theta = (build_theta(Data, h, q)).^2;
    A = Data*Theta*T'/(T*Theta*T');
    result = A(:,1);
end


function sigma = find_sigma(x, Data, k)
    s_distance = sum((Data-x).^2, 1);
    [~,ind] = sort(s_distance,'ascend');
    Neig = Data(:,ind(2:k+1)); 
    sigma = max(sqrt(sum((Neig-x).^2,1)));
end



function [theta, Theta] = least_square(Tau, Co_v, W)

    d = size(Tau, 1);
    ind = triu(true(size(Tau, 1)));
    d2 = d*(d+1)/2;
    Theta = cell(1,size(Co_v,1));
    for j = 1:size(Co_v, 1)
        G = zeros(d2, size(Tau,2));
        for i = 1:size(Tau,2)
            A = Tau(:,i)*Tau(:,i)';
            G(:,i) = A(ind);
        end
        theta = (G*W*G')\G*W*Co_v(j,:)';

        Theta_temp = zeros(d,d);
        Theta_temp(ind) = theta/2;
        Theta{j} = Theta_temp+Theta_temp';
    end
end



function [Theta, Co_h, Co_v, U_h, U_v] = fitting(x, Data, sigma, D, d)

    [Co_h, Co_v, U_h, U_v] = coordinate(x, Data, sigma, D, d);
    W = build_W(x, Data, sigma);
    [~, Theta] = least_square(Co_h, Co_v, W);
end


function [Theta, Co_h, Co_v, U_h, U_v] = fitting2(x, Data, sigma, D, d, U)

    [Co_h, Co_v, U_h, U_v] = coordinate2(x, Data, sigma, D, d, U);
    W = build_W(x, Data, sigma);
    [~, Theta] = least_square(Co_h, Co_v, W);
end


function mean = shift_mean(x, Data, h, k)

    [~,ind] = sort(sum((Data-x).^2,1),'ascend');
    mean = zeros(size(x));
    s_weight = 0;
    for i = 1:k
        w = exp(-norm(Data(:,ind(i))-x)^2/(h^2));
        %w = 1;
        s_weight = s_weight+w;
        mean = mean+w*Data(:,ind(i));
    end
    mean = mean/s_weight;
    
end


function [Co_h, Co_v, U_h, U_v] = coordinate(x, A, h, D, d)

    C = zeros(size(x,1),size(x,1));
    for i = 1: size(A,2)
        a = (x-A(:,i));
        C = C + exp(-norm(a)^2/(h^2))*a*a';
    end
    [V,~,~] = svd(C);
    U_h = V(:,1:d);
    U_v = V(:,d+1:D);
    Co_h = U_h'*(A-repmat(x,1,size(A,2)));
    Co_v = U_v'*(A-repmat(x,1,size(A,2)));
    
end


function [Co_h, Co_v, U_h, U_v] = coordinate2(x, A, h, D, d, U)

%     C = zeros(size(x,1),size(x,1));
%     for i = 1: size(A,2)
%         a = (x-A(:,i));
%         C = C + exp(-norm(a)^2/(h^2))*a*a';
%     end
%     [V,~,~] = svd(C);
    U_h = U;
    [V,~,~] = svd(eye(D)-U*U');
    U_v = V(:,1:D-d);
    Co_h = U_h'*(A-repmat(x,1,size(A,2)));
    Co_v = U_v'*(A-repmat(x,1,size(A,2)));
    
end


function W = build_W(x, A, h)

    centered = A-repmat(x,[1,size(A,2)]);
    W = zeros(size(centered,2));
    for k = 1:size(centered,2)
        W(k,k) = exp(-norm(centered(:,k))^2/(h^2));
    end
    
end



function [tau, x_new, error] = nonlinear_solution(U_h, U_v, x, x0, Theta)

    s = U_h'*(x-x0);
    c = U_v'*(x-x0);
    d = size(U_h, 2);
    D = size(x,1);
    A = zeros(d,D-d);
    tau = s;
    tau_old = zeros(size(s));
    error = [];
    while norm(tau-tau_old)> 1.0e-14
        tau_old = tau;
        for i = 1:D-d
            A(:,i) = Theta{i}*tau;
        end
        tau = (A*A'+eye(d)/2)\(s/2+A*c);
        for i = 1:D-d
            A(:,i) = Theta{i}*tau;
        end
        tau = (A*A'+eye(d)/2)\(s/2+A*c);
        error = [error,norm(tau-tau_old)];
    end
    iota = zeros(D-d,1);
    for j = 1:D-d
        iota(j) = tau'*Theta{j}*tau;
    end
    x_new = U_h*tau+U_v*iota+x0;
    
end


function theta = build_theta(Data, h, q)  
    theta = diag(sqrt(exp(-sum((Data-q).^2,1)/h^2)));
end


function T = Construct_Higher_Order(Tau) 
    d = size(Tau, 1);
    T = zeros(1+d+d*(d+1)/2, size(Tau,2));
    ind = triu(true(size(Tau, 1)));
    for i = 1:size(Tau,2)
        T(1:1+d,i) = [1; Tau(:,i)];
        temp = Tau(:,i)*Tau(:,i)';
        T(d+2:end,i) = temp(ind);
    end
end