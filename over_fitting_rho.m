%data = build_data();
num = 300;
[data, data_t] = build_circle(0.05, num);

%plot(data_true(1,:),data_true(2,:),'-','Linewidth',4);


%plot(q(1),q(2),'d','MarkerSize',10,'MarkerFaceColor','b');
Rho = [0,0.003:0.003:0.02,1];
k = 20;
d = 1;
n = size(data, 2)-1;
n = 40;
Data = data;

P_ = [];
recover_result = [];
for j = 1:length(Rho)
    P = [];
   % subplot('position', [0.02*(j+1)+0.22*j 0.12 0.22 0.8]);
   % axis([-1.5 1.5 -1.5 1.5])
    rho = Rho(j);
    for i = 1:size(data,2)
        q = data(:,i);
        % k determine the sigma and n determines the sample size
        [Tau_p, Data_p, h, center] = initial_Tau(q, d, Data, k, n);  
        [data_new, Taus, iter, Error, q_p] = fit_nonlinear(Data_p, q, Data_p, Tau_p, h, center, rho);
        P = [P, q_p(:,end)]; 
%         if mod(i,35)==0
%             hold on
%             plot(data_new(1,:),data_new(2,:),'.');
%         end
        i
    end
    recover_result = [recover_result, [sqrt(norm(data-data_t,'fro')^2/300);sqrt(norm(P-data_t,'fro')^2/300)]];
    %title(['\lambda=',num2str(rho)])
    %box on
    subplot(1,length(Rho),j)
    plot(data(1,:),data(2,:),'*')
    hold on
    plot(data_t(1,:), data_t(2,:),'-')
    hold on
    plot(P(1,:), P(2,:),'.-')
end
%%
%subplot('position', [0.02 0.12 0.22 0.8]);
%plot(data(1,:),data(2,:),'*');
%plot(P(1,:),P(2,:),'o');
%hold on
%plot(P(1,:),P(2,:),'.','MarkerSize',10)


% hold on
% plot(q_p(1,end), q_p(2,end),'o','MarkerSize',10,'MarkerFaceColor','r');
% 
% subplot(1,3,2)
% semilogy(Error(:,1),'-');
% hold on
% semilogy(Error(:,2),'-');
% legend({'Data Error','Tau Error'});
% 
% subplot(1,3,3)
% % error = zeros(1,length(Taus));
% % for i = 1:length(Taus)-1
% %     %error = [error, norm(Taus{i}'*Taus{i}-Taus{i+1}'*Taus{i+1})];
% % end
% semilogy(Error(:,3));


function [Tau, Data_selection, h, center] = initial_Tau(q, d, Data, k, n)
    [~,ind] = sort(sum((Data-q).^2,1),'ascend');
    Data_selection = Data(:,ind(2:n+1));
    
    h = find_sigma(q, Data, k);
    [U, center] = principal(Data_selection, h, q, d);
    Tau = U'*(Data-center);
   

    Tau_p = Tau(:,ind(2:n+1));
    Tau = qrs(Tau_p);
end


function [Data_new_Constructed, Taus, iter, Error, Q_projection] = fit_nonlinear(Data, q, Data_p, Tau, h1, center, rho)
    Error = [];
    iter = 1;
%    center = q;
    Q_projection = [];
    while true
        % Parameter for regression Using Tau
        [c, A, R, Theta] = Least_Fitting(Data, Tau, h1, center, rho);
        Rho = 0:0.01:0.1;
        g = zeros(1,length(Rho));
        t1 = zeros(1,length(Rho));
        t2 = zeros(1,length(Rho));
        for i = 1:length(Rho)
            [g(i),t1(i),t2(i)] = Eval_Rho(Data, Tau, h1, center, Rho(i));
        end
        set(gca,'LooseInset',max(get(gca,'TightInset'),0.01))
        subplot(1,3,1)
        plot(Rho,g(:),'o-');
        subplot(1,3,2)
        plot(Rho,t1(1,:),'*-');
        subplot(1,3,3)
        plot(Rho,t2(1,:),'>-');
        d = size(Tau, 1);
        B = build_tensor(R, d);
        
        Data_Constructed = R*Construct_Higher_Order(Tau);
        
        Tau_q = projection(q, A, B, c, 0);
        q_projection = R*Construct_Higher_Order(Tau_q);
        Q_projection = [Q_projection, q_projection];
%        center = q_projection; 
        %h1 = find_sigma(q_projection, Data, k);
        % project data_p and find new tau using A B, c
        Tau_old = Tau;
        for i = 1:size(Data_p, 2)
            %fprintf('iter:%f\n',i/size(Data,2));
            Tau(:,i) = projection(Data_p(:,i), A, B, c, Tau(:,i));%projection(Data_p(:,i), A, B, c, Tau(:,i));
        end
        Tau_ee = qrs(Tau);
        Tau = Tau_ee(1:d,:);
        Taus{iter} = Tau;
        
        Data_new_Constructed = R*Construct_Higher_Order(Tau);
        error = norm(Data_new_Constructed- Data_Constructed,'fro');
        error2 = norm(Tau'*Tau- Tau_old'*Tau_old,'fro');
        
       
        error3 = norm((Data_Constructed-Data)*Theta,'fro');
        Error = [Error; [error, error2, error3]];
        if error < 1.e-4 || iter>800
            break;
        end
        iter = iter+1;
    end   
end


function result = search_lambda(Data, Tau, h1, center,interval_l, interval_r, epsilon, delta)
    while (interval_r-interval_l) > epsilon
           middle = (interval_r-interal_l)/2;
           [g,~,~] = Eval_Rho(Data, Tau, h1, center, middle);
           if abs(g)>deltacompare
               interval_l = middle;
           elseif abs(g)<delta
               interval_r = middle;    
           end
    end
    result = (interval_r+interval_l)/2;
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


function B = build_tensor(para, d)compare
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
        if norm(tau_new-tau)<1.e-6 || iter>300compare
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



function [c, A, P, Theta_h] = Least_Fitting(Data, Tau, h, q, rho)
%   Tau = U'*(Data-q);
    T= Construct_Higher_Order(Tau);
    d = size(Tau, 1);
    Theta_h = (build_theta(Data, h, q));
    Theta = Theta_h.^2;
    R = Construct_Regularization(Tau);
    %R = Construct_Regularization2(d, T*Theta*T');
    P = Data*Theta*T'/(T*Theta*T'+rho*R);
    c = P(:,1);
    A = P(:,2:d+1);
end


function [g, t1, t2] = Eval_Rho(Data, Tau, h,q, rho)
    T= Construct_Higher_Order(Tau);
    d = size(Tau, 1);
    Theta_h = (build_theta(Data, h, q));
    Theta = Theta_h.^2;
    J = Construct_Regularization(Tau);
    R = Data*Theta*T'/(T*Theta*T'+rho*J);
    g = - trace(R*J/(T*Theta*T'+rho*J)*J*R'); 
    t1 = norm(R*J,'fro');
    t2 = norm(R*T-Data,'fro');
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


function theta = build_theta(Data, h, q)  
    %theta = diag(sqrt(exp(-sum((Data-q).^2,1)/h^2)));
    theta = diag(ones(1,size(Data, 2)));
end


function [data_true, data] = generate_data(sigma, num)
    theta = linspace(-pi/2, pi/2, num);%pi/4:0.1:3*pi/4;
    data_true = [cos(theta);sin(theta)];
    data = data_true+sigma*randn(2,length(theta));
end


function [data, data_t] = build_circle(sigma, num)
    theta = linspace(-pi, pi, num);
    data_t = [cos(theta);sin(theta)];
    data = data_t + sigma*randn(size(data_t));
end


function [data,data_t] = build_sin(sigma, num)
    theta = linspace(-pi, pi, num);
    data_t = [cos(theta);(theta)];
    data = data_t + sigma*randn(size(data_t));
end

function data = build_data()

    [data_true, data1] = generate_data(0.03, 50);

    [data_true, data2] = generate_data(0.03, 50);

    data = [[data1(1,:)+0.6*ones(size(data1(1,:)));data1(2,:)], data2];
end