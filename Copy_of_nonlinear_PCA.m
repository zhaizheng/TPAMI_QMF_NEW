[data_true, data] = generate_data(0.15, 80);

%plot(data_true(1,:),data_true(2,:),'-','Linewidth',4);
%hold on
q = [0.1; 0];
%plot(q(1),q(2),'d','MarkerSize',10,'MarkerFaceColor','b');
k = 25;
d = 1;
n = 70;
Data = data;
[Tau_p, Data_p, h, center] = initial_Tau(q, d, Data, k, n);

%Tau_p=randn(1,size(Data_p,2));
[data_new, Taus, iter, Error, P, q_p] = fit_nonlinear(Data_p, q, Data_p, Tau_p, h, k, center);

% subplot(1,3,1)
% plot(data(1,:),data(2,:),'*','MarkerSize',7);
% %hold on
% %plot(Data_p(1,:),Data_p(2,:),'k*');
% hold on
% plot(data_new(1,:),data_new(2,:),'o','MarkerSize',5,'MarkerFaceColor','b');
%hold on
%plot(q_p(1,end), q_p(2,end),'o','MarkerSize',10,'MarkerFaceColor','r');


%subplot(1,2,1)
%plot(data(1,:),data(2,:),'*','MarkerSize',7);
data_g = linspace(-0.2,0.2,30);
data_g2 = linspace(-1.1,1.1,30);

tau_g = Construct_Higher_Order(data_g);
tau_g2 = Construct_Higher_Order(data_g2);

c = mean(data, 2);
[U,~,~] = svd(data-c);
%new = U(:,1)*U(:,1)'*(data-c)+c;
new = U(:,1)*data_g2+c;
%plot(new(1,:),new(2,:),'o','MarkerSize',5,'MarkerFaceColor','b');
plot(new(1,:),new(2,:),'b-','LineWidth',3);
hold on
plot(data(1,:),data(2,:),'k*','MarkerSize',9);
hold on
data = Data_p;
tau = U(:,1)'*(data-c);
Tau = Construct_Higher_Order(tau);
A = data*Tau'/(Tau*Tau');
Recover = A*tau_g2;

[~,ind] = sort(Recover(2,:),'descend');
plot(Recover(1,ind),Recover(2,ind),'b-','LineWidth',3);
%plot(Recover(1,ind),Recover(2,ind),'-o','MarkerSize',6,'MarkerFaceColor','b');

hold on
data_new = P*tau_g;
[~,ind2] = sort(data_new(2,:),'descend');
plot(data_new(1,ind2),data_new(2,ind2),'r-.','LineWidth',3);
%plot(data_new(1,ind2),data_new(2,ind2),'-.d','MarkerSize',6,'MarkerFaceColor','r');




function [Tau, Data_selection, h, center] = initial_Tau(q, d, Data, k, n)
    h = find_sigma(q, Data, k);
    [U, center] = principal(Data, h, q, d);
    Tau = U'*(Data-center);
   
    [~,ind] = sort(sum((Data-center).^2,1),'ascend');
    Data_selection = Data(:,ind(2:n+1));
    Tau_p = Tau(:,ind(2:n+1));
    Tau = qrs(Tau_p);
end


function [Data_new_Constructed, Taus, iter, Error, P, Q_projection] = fit_nonlinear(Data, q, Data_p, Tau, h1, k, center)
    Error = [];
    iter = 1;
%    center = q;
    Q_projection = [];
    while true
  %     
        % Parameter for regression Using Tau
        [c, A, P, Theta] = Least_Fitting(Data, Tau, h1, center);
        d = size(Tau, 1);
        B = build_tensor(P, d);
        
        Data_Constructed = P*Construct_Higher_Order(Tau);
        
        Tau_q = projection(q, A, B, c, 0);
        q_projection = P*Construct_Higher_Order(Tau_q);
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
        
        Data_new_Constructed = P*Construct_Higher_Order(Tau);
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



function sigma = find_sigma(x, Data, k)
    s_distance = sum((Data-x).^2, 1);
    [~,ind] = sort(s_distance,'ascend');
    Neig = Data(:,ind(2:k+1)); 
    sigma = max(sqrt(sum((Neig-x).^2,1)));
end


function Tau = qrs(Tau)
    d = size(Tau, 1);
    [Q,~] = qr([ones(size(Tau, 2),1),Tau']);
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



function [c, A, P, Theta_h] = Least_Fitting(Data, Tau, h, q)
%   Tau = U'*(Data-q);
    T = Construct_Higher_Order(Tau);
    d = size(Tau, 1);
    Theta_h = (build_theta(Data, h, q));
    Theta = Theta_h.^2;
    P = Data*Theta*T'/(T*Theta*T');
    c = P(:,1);
    A = P(:,2:d+1);
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


function theta = build_theta(Data, h, q)  
    %theta = diag(sqrt(exp(-sum((Data-q).^2,1)/h^2)));
    theta = diag(ones(1,size(Data, 2)));
end


function [data_true, data] = generate_data(sigma, num)
    theta = linspace(-pi/2, pi/2, num);%pi/4:0.1:3*pi/4;
    data_true = [cos(theta);sin(theta)];
    data = data_true+sigma*randn(2,length(theta));
end