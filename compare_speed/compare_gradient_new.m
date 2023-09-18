sigma = 0.15; num = 80;  
approach = {'gradient_optimal','gradient','alter'}; %'gradient';%'alter';
choice = {'least_square','gradient'};
%[data_true, data] = generate_data(sigma, num);
sigma = 0.02; num = 600;
[data_n, data_t] = build_sphere(sigma, num);
x = data_n(:,1);%[1;0;0];
cdata = data_n-x;
dis = sum(cdata.^2,1) ;
[~,ind] = sort(dis,'ascend');
data = data_n(:,ind(2:20));
d = 2;
[CS, err, P, Time] = compare_gradient(data,d);


function [CS, Err, R, Time] = compare_gradient(Data, d)
    Tau = initial(Data, d);
    Err = [];
    iter = 1;
    R = ones(size(Data,1),1+d+d*(d+1)/2);
    Time = [];
    while true
        t1 = clock;
        gR = gradient_R(Data, Tau, R);
        s = 1;
        while evaluate(Data, R-s*gR, Tau) > evaluate(Data, R, Tau)
           s = s/2;
        end
        R_new = R-s*gR;
        time1 = etime(clock, t1);
        d = size(Tau, 1);
        Tau_old = Tau;
        t2 = clock;
        [c, A, B] = R2abc(R, d);
        for i = 1:size(Data, 2)
            s = 1;
            g = gradient(Data(:,i), A, B, c, Tau(:,i));
            while f_value(Data(:,i), A, B, c, Tau(:,i)-s*g) > f_value(Data(:,i), A, B, c, Tau(:,i)) 
                s = s/2;
            end
            Tau(:,i) = Tau(:,i) - s*g;
        end
        time2 = etime(clock, t2);
        Tau_ee = qrs(Tau);
        Tau = Tau_ee(1:d,:);
        R = R_new;
        CS = R *Construct_Higher_Order(Tau);
        error = norm(Data-CS,'fro').^2;
        Err = [Err,error];
        Tau_error = norm(Tau'*Tau- Tau_old'*Tau_old,'fro');
        if Tau_error < 1.e-8 || iter>400000
            break;
        end
        Time = [Time;[time1,time2]];
        iter = iter+1;
    end
end

function re = ini(a,b,num,d)
    l1 = linspace(a,b,num);
    l2 = linspace(a,b,num);
    [A, B] = meshgrid(l1,l2);
    re = zeros(d,num*num);
    re(1,:) = A(:);
    re(2,:) = B(:);
end



function Tau = qrs(Tau)
    d = size(Tau, 1);
    [Q,~] = qr([ones(size(Tau, 2),1),Tau']);
    %Tau = size(Tau,2)*Q(:,2:d+1)';
    Tau = Q(:,2:d+1)';
end

function Tau = initial(Data, d)
    center = sum(Data,2)/size(Data,2);
    [V,~,~] = svd(Data-center);
    Tau = (V(:,1:d))'*(Data-center);
end


function R = Least(Data, Tau)
    T = Construct_Higher_Order(Tau);
    R =  Data*T'/(T*T');
    % c = P(:,1);
    % A = P(:,2:d+1);
end


function re = reshape_s(a)
    n = sqrt(length(a));
    re = reshape(a,[n,n]);
end


function gR = gradient_R(Data, Tau, R)
    T = Construct_Higher_Order(Tau);
    gR = -2*(Data-R*T)*T';
    % R = R - step*gR;
    % c = R(:,1);
    % A = R(:,2:d+1);
end

function re = evaluate(X, R, Tau)
    T = Construct_Higher_Order(Tau);
    re = norm(X-R*T,'fro').^2;
end

function [c,A,B] = R2abc(R,d)
    c = R(:,1);
    A = R(:,2:d+1);
    B = build_tensor(R, d);
end


%function result = tensor_fold(B, tau)
%     result = zeros(size(B,1),size(B,2));
%     for i = 1:size(B,1)
%         result(i,:) = squeeze(B(i,:,:))*tau;
%     end
% end


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


function [tau,Tau] = gradient_descent(x, A, B, c, tau)
    Tau = [];
    iter = 0;
    s = 1;
    while true
        Bm = tensor_fold(B, tau);
        v = x- (A*tau+ Bm*tau+c);
        g = (-A-2*Bm)'*v;
        while f_value(x, A, B, c,tau-s*g) > f_value(x, A, B, c, tau) 
            s = s/2;
        end
        tau_new = tau - s * g;
        if norm(tau_new-tau)<1.e-6 || iter>3000
            break;
        end
        tau = tau_new;
        Tau = [Tau, tau];
        iter = iter+1;
    end
end

function a = f_value(x, A, B, c, tau)
    Bm = tensor_fold(B, tau);
    a = norm(x-(A*tau+ Bm*tau+c))^2;
end

% function result = tensor_fold(B, tau)
%     result = zeros(size(B,1),size(B,2));
%     for i = 1:size(B,1)
%         result(i,:) = squeeze(B(i,:,:))*tau;
%     end
% end

function g = gradient(x, A, B, c, tau)
    Bm = tensor_fold(B, tau);
    v = x- (A*tau+ Bm*tau+c);
    g = (-A-2*Bm)'*v;
end


function tau = projection(x, A, B, c, tau) %project x onto f(tau) = A tau+ B(tau,tau)+c
    iter = 0; 
    while true
        Bm = tensor_fold(B, tau);
        tau_new = (2*Bm'*Bm+Bm'*A+A'*A+A'*Bm)\((2*Bm'+A')*(x-c)-Bm'*A*tau);
        if norm(tau_new-tau)<1.e-6 || iter>3000
%             if iter>300
%                 fprintf('diverge projecting tau\n');
%             end
            break;
        end
        tau = tau_new;
        iter = iter+1;
    end 
end



function B = build_tensor(para, d)
    B = zeros(size(para,1),d,d);
    ind = triu(true(d));
    for i = 1:size(para,1)
        temp = zeros(d, d);
        temp(ind) = para(i,d+2:end);
        B(i,:,:) = (temp+temp')/2;
    end
end


function result = tensor_fold(B, tau)
    result = zeros(size(B,1),size(B,2));
    for i = 1:size(B,1)
        result(i,:) = squeeze(B(i,:,:))*tau;
    end
end

function [data, data_t] = build_sphere(sigma, num)
    nor = randn(3,num);
    data_t = nor*diag(1./sqrt(sum(nor.^2,1)));
    data = data_t + sigma*randn(size(data_t));
end