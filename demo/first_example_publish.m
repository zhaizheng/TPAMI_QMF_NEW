%data = build_data();
%data = build_circle(sigma, num)
%%
%data = build_circle(0.1, 50);
%theta = linspace(1,3*pi,100);
theta = [];
current = 0;
r = linspace(1,4,100);
for i = 1:length(r)
    current = current+0.2/r(i);
    theta = [theta,current];
end
x = r.*cos(theta);
y = r.*sin(theta);
true = [x;y];
data = [x;y]+0.1*randn(2,100);
a = min(data(1,:));
b = max(data(1,:));
c = min(data(2,:));
e = max(data(2,:));

%%

figure
%subplot(1,3,1)
%plot(data_true(1,:),data_true(2,:),'-','Linewidth',4);
hold on
beta = linspace(0,2*pi,100);

%plot(q(1),q(2),'d','MarkerSize',10,'MarkerFaceColor','b');
Rho = [1000, 0.02, 0];
k = 10;
d = 1;
n = size(data, 2)-1;
n = 10;
Data = data;%
t = tiledlayout(1,2,'TileSpacing','Compact');
P_ = [];
TT = [];
str = cell(1,2);
str{1} = 'linear';
str{2} = 'quadratic';

%nexttile
%for i = 1:size(data,2)
%     set(gca,'FontSize',18)
%     axis([a b c e])
%     box on
%end

%title('data','FontSize',18)
for j = 1:2
    nexttile
    plot(data(1,:),data(2,:),'r.','markersize',12);
    hold on
    p2 = plot(true(1,:),true(2,:),'--','linewidth',1);
    %subplot('position', [0.02*(j+1)+0.22*j 0.12 0.22 0.8]);
    %plot(cos(beta),sin(beta),'--','linewidth',1)
    hold on
    axis([a b c e])
    rho = Rho(j);
    P{j} = [];
   
    for i = 1:size(data,2)
        q = data(:,i);
        % k determine the sigma and n determines the sample size
        [Tau_p, Data_p, h, center] = initial_Tau(q, d, Data, k, n);  
        f = fit_nonlinear(Data_p, Tau_p, rho);
        
        Tau_q = projection(data(:,i), f.A, f.B, f.c, 0);
        proj_q = f.Parm*Construct_Higher_Order(Tau_q);
        P{j} = [P{j}, proj_q]; 
        T = [];
    
        if mod(i,8)==0
%             plot(data(1,i),data(2,i),'d','markersize',8);
%             hold on
%             plot(proj_q(1,:),proj_q(2,:),'*','markersize',8);
%             hold on
            data_new = f.Parm*Construct_Higher_Order(linspace(-0.5,0.5,10));
            p3 = plot(data_new(1,:),data_new(2,:),'b-','linewidth',2);
            T = [T, Tau_q];
        end
        set(gca,'FontSize',18)
        TT = [TT;T];
        i
        axis([a b c e])
        %axis([-1.5 1.5 -1.5 1.5])
    end
    title(str{j},'FontSize',20);
    %title(['\lambda=',num2str(rho)])
    %legend([p2 p3],{'data','truth','fitted'})
    box on
end
%%
%subplot('position', [0.02 0.12 0.22 0.8]);
% subplot(1,3,1)
% for i = 1:size(data,2)
%     if mod(i,5)==0
%         [~, Data_p, ~, ~] = initial_Tau(data(:,i), d, Data, k, n);  
%         hold on
%         axis([-1.5 1.5 -1.5 1.5])
%         plot(Data_p(1,:),Data_p(2,:),'d');
%         box on
%     end
% end
%%
for i = 1:3
    result(i) = norm(P{i}-P{i}*diag(1./sqrt(sum(P{i}.^2,1))),'fro');
    fprintf('error%d=%.3f\n',i,result(i));
end

%%
% figure
% %plot(data_new(1,:),data_new(2,:),'.');
% subplot(1,2,1)
% plot(P{2}(1,:),P{2}(2,:),'o');
% subplot(1,2,2)
% plot(P{3}(1,:),P{3}(2,:),'o');
% axis([-1.5 1.5 -1.5 1.5])
%hold on
%plot(P(1,:),P(2,:),'.','MarkerSize',10)



function [Tau, Data_selection, h, center] = initial_Tau(q, d, Data, k, n)
    [~,ind] = sort(sum((Data-q).^2,1),'ascend');
    Data_selection = Data(:,ind(2:n+1));
    
    h = find_sigma(q, Data, k);
    [U, center] = principal(Data_selection, h, q, d);
    Tau = U'*(Data-center);
  
    Tau_p = Tau(:,ind(2:n+1));
    Tau = qrs(Tau_p);
end


function f = fit_nonlinear(Data, Tau, rho, W)
    if ~exist('W','var')
        W = diag(ones(1,size(Data, 2))); %equal weight
    end
    iter = 1;
    f.q_sequence = [];
    while true
        % Parameter for regression Using Tau
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


function data = build_sin(sigma, num)
    theta = linspace(-pi, pi, num);
    data = [cos(theta);(theta)];
    data = data + sigma*randn(size(data));
end

function data = build_data()

    [data_true, data1] = generate_data(0.03, 50);

    [data_true, data2] = generate_data(0.03, 50);

    data = [[data1(1,:)+0.6*ones(size(data1(1,:)));data1(2,:)], data2];
end