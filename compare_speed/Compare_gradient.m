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

for ss = 1:2
    figure
    t = tiledlayout(2,3,'TileSpacing','Compact');
    Err = cell(1,3);
    Time = cell(1,3);
    d = 2;
    Result = zeros(3,4);
    for i = 1:3
        t0 = clock;
        [CS, err, P, Time{i}] = compare(data, d, approach{i}, choice{ss});
        Err{i} = err;
        fprintf('Method:%s,time:%f steps:%d, err:%f, ratio:%f \n',approach{i}, etime(clock,t0), length(Err{i}),Err{i}(end),sum(Time{i}(:,2)/sum(Time{i}(:,1))));
        Result(i,:) = [etime(clock,t0), length(Err{i}), Err{i}(end), sum(Time{i}(:,2)/sum(Time{i}(:,1)))];
        nexttile
        plot3(data(1,:),data(2,:),data(3,:),'b*');
        hold on
        plot3(CS(1,:), CS(2,:),CS(3,:),'r*');
        
        hold on
        S = ini(-0.30,0.30,20,2);
        data_new = P*Construct_Higher_Order(S);
        %plot(data_new(1,:),data_new(2,:),'-','linewidth',2);
        A = reshape_s(data_new(1,:));
        B = reshape_s(data_new(2,:));
        C = reshape_s(data_new(3,:));
        mesh(A,B,C)
    end
    for i = 1:3
        nexttile
        semilogy(Err{i}(1:10))
    end
endalter








function [data_true, data] = generate_data(sigma, num)
    theta = linspace(-pi/2, pi/2, num);%pi/4:0.1:3*pi/4;
    data_true = [cos(theta);sin(theta)];
    data = data_true+sigma*randn(2,length(theta));
end



function [CS, Err, R, Time] = compare(Data, d, alg, alg2)
    Tau = initial(Data, d);
    Err = [];
    iter = 1;
    R = ones(size(Data,1),1+d+d*(d+1)/2);
    Time = [];
    while true
        t1 = clock;
        switch alg2 
            case 'least_square'
               R = Least(Data, Tau);   
               [c, A, B] = R2abc(R, d);
            case 'gradient'
               s = 1; 
               gR = gradient_R(Data, Tau, R);
               while evaluate(Data, R-s*gR, Tau) > evaluate(Data, R, Tau)
                   s = s/2;
               end
               R = R-s*gR;
               [c, A, B] = R2abc(R, d);
        end
        time1 = etime(clock, t1);
        d = size(Tau, 1);
        Tau_old = Tau;
        t2 = clock;
        for i = 1:size(Data, 2)
            switch alg
                case 'alter'
                    Tau(:,i) = projection(Data(:,i), A, B, c, Tau(:,i));
                case 'gradient_optimal'
                    [temp,~] = gradient_descent(Data(:,i), A, B, c, Tau(:,i));
                    Tau(:,i) = temp;
                case 'gradient'
                    s = 1;
                    g = gradient(Data(:,i), A, B, c, Tau(:,i));
                    while f_value(Data(:,i), A, B, c, Tau(:,i)-s*g) > f_value(Data(:,i), A, B, c, Tau(:,i)) 
                        s = s/2;
                    end
                    Tau(:,i) = Tau(:,i) - s*g;
            end
          
        end
        time2 = etime(clock, t2);
        Tau_ee = qrs(Tau);
        Tau = Tau_ee(1:d,:);

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