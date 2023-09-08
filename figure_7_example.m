clear
t = linspace(pi/3,2*pi/3,20);
truth = [t; sin(t)];
train = truth+0.03*randn(2,20);
%test  = test_truth+0.03*randn(2,20);
figure
t = tiledlayout(1,3,'TileSpacing','Compact');
c = mean(train,2);
[~,~,V] = svd(train-c);
Tau = V(:,1)';
Rho = 0:0.002:0.1; 
Delta = []; S = [];
P = cell(1,50);
for j = 1:length(Rho)
    rho = Rho(j);
    [f, resl] = fit_nonlinear(train, Tau, rho, 0, 0);
    [g, t1, ~] = Eval_Rho(train, Tau, rho);
    Delta = [Delta, -g];
    S = [S, t1];
    P{j} = [];
    for i = 1:size(train,2)
        Tau_q = projection(train(:,i), f.A, f.B, f.c, 0);
        proj_q = f.Parm*Construct_Higher_Order(Tau_q);
        P{j} = [P{j}, proj_q];
    end
    error(j) = norm(P{j}-truth,'fro')^2/size(truth,2);
    A = 1:50;
    if j==1 || j==length(Rho) 
        nexttile
        plot(train(1,:),train(2,:),'>','MarkerSize',5,'MarkerFaceColor','r');
        hold on
        plot(P{j}(1,:),P{j}(2,:),'o-','MarkerSize',5,'Linewidth',2);
        hold on
        plot(truth(1,:),truth(2,:),'d-','MarkerSize',5,'Linewidth',2);
        title(['\lambda=',num2str(Rho(j))])
        set(gca,'FontSize',18)
    end
end

k = min(A(error==min(error)));
nexttile
plot(train(1,:),train(2,:),'>','MarkerSize',5,'MarkerFaceColor','r');
hold on
plot(P{k}(1,:),P{k}(2,:),'o-','MarkerSize',5,'Linewidth',2);
hold on
plot(truth(1,:),truth(2,:),'d-','MarkerSize',5,'Linewidth',2);
title(['\lambda=',num2str(Rho(k))])
legend({'Noisy Data','Fitted Curve','Truth Curve'})
set(gca,'FontSize',18)
figure
t2 = tiledlayout(1,3,'TileSpacing','Compact');
nexttile
plot(Rho, error,'*-')
ylabel('error')
xlabel('\lambda')
axis([min(Rho),max(Rho),min(error),max(error)])
set(gca,'FontSize',18)


nexttile
plot(Delta, error,'*-')
ylabel('error')
xlabel('-s''(\lambda)')
axis([min(Delta),max(Delta),min(error),max(error)])
set(gca,'FontSize',18)



nexttile
plot(S, error,'*-')
ylabel('error')
xlabel('s(\lambda)')
axis([min(S),max(S),min(error),max(error)])
set(gca,'FontSize',18)


function [f, resl] = fit_nonlinear(Data, Tau, rho, delta, adaptive, W)
    if ~exist('W','var')
        W = diag(ones(1,size(Data, 2))); %equal weight
    end
    iter = 1;
    f.q_sequence = [];
    while true
        % Parameter for regression Using Tau
        
        if adaptive == 1
            inter_l = 0.0001; inter_r = 1000; 
            resl = search_lambda(Data, Tau, inter_l, inter_r, inter_l, delta, W);
            rho = resl.rho;
        else
            resl = [];
        end
        
        [c, A, P] = Least_Fitting(Data, Tau, rho, W);
        d = size(Tau, 1);
        B = build_tensor(P, d);
        
        f.Data_Constructed = P*Construct_Higher_Order(Tau);
        f.A = A; f.B = B; f.c = c;

        Tau_old = Tau;
        for i = 1:size(Data, 2)
            Tau(:,i) = projection(Data(:,i), A, B, c, Tau(:,i));
        end
        Tau_ee = qrs(Tau);
        Tau = Tau_ee(1:d,:);
        f.Taus{iter} = Tau;
        f.Parm = P;
        f.Tau = Tau;
        f.Data_new_Constructed = P*Construct_Higher_Order(Tau);
        f.data_error = norm(f.Data_new_Constructed- f.Data_Constructed,'fro');
        f.Tau_error(iter) = norm(Tau'*Tau- Tau_old'*Tau_old,'fro');
        if f.Tau_error(iter) < 1.e-4 || iter>800
            break;
        end
        iter = iter+1;
    end   
end


function [g, t1, t2] = Eval_Rho(Data, Tau, rho, W)
        if ~exist('W','var')
            W = diag(ones(1,size(Data, 2))); 
        end

        T= Construct_Higher_Order(Tau);
        J = Construct_Regularization(Tau);
        R = Data*W*T'/(T*W*T'+rho*J);
        g = - trace(R*J/(T*W*T'+rho*J)*J*R'); 
        t1 = norm(R*J,'fro');
        t2 = norm(R*T-Data,'fro');
end



function theta = build_theta(Data, h, q)  
    theta = diag(sqrt(exp(-sum((Data-q).^2,1)/h^2)));
    %theta = diag(ones(1,size(Data, 2)));
end


function [c, A, P] = Least_Fitting(Data, Tau, rho, W)

    T= Construct_Higher_Order(Tau);
    d = size(Tau, 1);
    Theta = W.^2;
    R = Construct_Regularization(Tau);
    %R = Construct_Regularization2(d, T*Theta*T');
    P = Data*Theta*T'/(T*Theta*T'+rho*R);
    c = P(:,1);
    A = P(:,2:d+1);
end



function sigma = find_sigma(x, Data, k)
    s_distance = sum((Data-x).^2, 1);
    [~,ind] = sort(s_distance,'ascend');
    Neig = Data(:,ind(2:k+1)); 
    sigma = max(sqrt(sum((Neig-x).^2,1)));
end


function Tau = qrs(Tau)
    d = size(Tau, 1);
    [Q,~] = qr([ones(size(Tau, 2),1),Tau']);
    %Tau = size(Tau,2)*Q(:,2:d+1)';
    Tau = Q(:,2:d+1)';
end


function [U,center] = principal(Data, h, q, d)
    Theta = (build_theta(Data, h, q)).^2;
    center = sum(Data*Theta, 2)/sum(diag(Theta));
    [V,~,~] = svd((Data-center)*Theta*(Data-center)');
    U = V(:,1:d);
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


function tau = projection(x, A, B, c, tau) %project x onto f(tau) = A tau+ B(tau,tau)+c
    iter = 0; 
    while true
        Bm = tensor_fold(B, tau);
        tau_new = (2*Bm'*Bm+Bm'*A+A'*A+A'*Bm)\((2*Bm'+A')*(x-c)-Bm'*A*tau);
        if norm(tau_new-tau)<1.e-6 || iter>300
            if iter>300
                fprintf('diverge projecting tau\n');
            end
            break;
        end
        tau = tau_new;
        iter = iter+1;
    end 
end


function result = tensor_fold(B, tau)
    result = zeros(size(B,1),size(B,2));
    for i = 1:size(B,1)
        result(i,:) = squeeze(B(i,:,:))*tau;
    end
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


function R = Construct_Regularization(Tau)

    d = size(Tau, 1);
    R = zeros(1+d+d*(d+1)/2);
    R(d+2:end,d+2:end) = eye(d*(d+1)/2);
    
    %R = eye(1+d+d*(d+1)/2);
end


function R = Construct_Regularization2(d, A)
    
    [U,~, ~] = svd(A);
    R = U(:,d+1:end)*U(:,d+1:end)';
    %R = eye(1+d+d*(d+1)/2);
end


% function theta = build_theta(Data, h, q)  
%     %theta = diag(sqrt(exp(-sum((Data-q).^2,1)/h^2)));
%     theta = diag(ones(1,size(Data, 2)));
% end


function [data_true, data] = generate_data(sigma, num)
    theta = linspace(-pi/2, pi/2, num);%pi/4:0.1:3*pi/4;
    data_true = [cos(theta);sin(theta)];
    data = data_true+sigma*randn(2,length(theta));
end


function data = build_circle(sigma, num)
    theta = linspace(-pi, pi, num);
    data = [cos(theta);sin(theta)];
    data = data + sigma*randn(size(data));
end

function data = build_sphere(sigma, num)

    data = randn(3,num);
    %data = data*diag(1./sqrt(sum(data.^2.1)));
    data = bsxfun(@rdivide,data, sqrt(sum(data.^2,1)));
    data = data + sigma*randn(size(data));
end



function [data,Orgi] = build_sin(sigma, num)
    theta = linspace(-2*pi, 2*pi, num);
    Orgi = [theta;cos(theta)];
    data = Orgi + sigma*randn(size(Orgi));
end

function data = build_data()

    [data_true, data1] = generate_data(0.03, 50);

    [data_true, data2] = generate_data(0.03, 50);

    data = [[data1(1,:)+0.6*ones(size(data1(1,:)));data1(2,:)], data2];
end

