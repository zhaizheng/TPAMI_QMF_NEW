%%
datat = load('mfeat-pix');
data = datat';
cdata = data-mean(data,2);

[U,~,V] = svd(cdata);
label = [];
for i = 1:10
    label = [label,i*ones(1,200)]; 
end
subplot(1,2,1)
scatter(V(:,1),V(:,2),[],label')

%%
d = 2;
q = cdata(:,1);
k = 30;
n = 10;

data = V(:,1:3)';
NEW = zeros(size(data));
for i = 1:size(data,2)
    fprintf('iter:%d\n',i);
    [B, Tau] = select_neig(data, data(:,i), k);
    [~, ~, ~, ~, ~, NEW(:,i)] = fit_nonlinear(B, data(:,i), Tau);
end
%%
subplot(1,2,2)
scatter(NEW(1,:),NEW(2,:),[],label');

% [Tau, Data_selection, h1, center] = initial_Tau(q, d, cdata, k, n);
% 
% datas = Data_selection;
% 
% [Data_new_Constructed, Taus, iter, Error, P, Q_projection] = fit_nonlinear(datas, datas, Tau, h1,  center);


function [B, tau] = select_neig(A, a, k)
    [~,ind] = sort(sum((A-a).^2,1),'ascend');
    B = A(:,ind(2:k+1));
    [U,~,~] = svd(B-mean(B,2));
    tau = U(:,1:2)'*(B-mean(B,2));
end


function [Tau, Data_selection, h, center] = initial_Tau(q, d, Data, k, n)
    h = find_sigma(q, Data, k);
    [U, center] = principal(Data, h, q, d);
    Tau = U'*(Data-center);
   
    [~,ind] = sort(sum((Data-center).^2,1),'ascend');
    Data_selection = Data(:,ind(1:n));
    Tau_p = Tau(:,ind(1:n));
    Tau = qrs(Tau_p);
end


function [Data_new_Constructed, Taus, iter, Error, P, p] = fit_nonlinear(Data, point, Tau)
    Error = [];
    iter = 1;

    while true
  %     
        % Parameter for regression Using Tau
        [c, A, P, Theta] = Least_Fitting(Data, Tau);
        d = size(Tau, 1);
        B = build_tensor(P, d);
        
        Data_Constructed = P*Construct_Higher_Order(Tau);
        
        Tau_old = Tau;
        for i = 1:size(Data, 2)
           % fprintf('iter:%f\n',i/size(Data,2));
            Tau(:,i) = projection(Data(:,i), A, B, c, Tau(:,i));%projection(Data_p(:,i), A, B, c, Tau(:,i));
        end
        Tau_ee = qrs(Tau);
        Tau = Tau_ee(1:d,:);
        Taus{iter} = Tau;
        
        Data_new_Constructed = P*Construct_Higher_Order(Tau);
        error = norm(Data_new_Constructed- Data_Constructed,'fro');
        error2 = norm(Tau'*Tau- Tau_old'*Tau_old,'fro');
        
       
        error3 = norm((Data_Constructed-Data)*Theta,'fro');
        Error = [Error; [error, error2, error3]];
        fprintf('construction error=%f\n', error);
        if error < 1.e-4 || iter>800
            break;
        end
        iter = iter+1;
    end  
    tauu = projection(point, A, B, c, [0;0]);
    p = P*Construct_Higher_Order(tauu);
end



function sigma = find_sigma(x, Data, k)
    s_distance = sum((Data-x).^2, 1);
    [~,ind] = sort(s_distance,'ascend');
    Neig = Data(:,ind(1:k)); 
    sigma = max(sqrt(sum((Neig-x).^2,1)));
end


function Tau = qrs(Tau)
    d = size(Tau, 1);
    [Q,~] = qr([ones(size(Tau, 2),1),Tau']);
    Tau = sqrt(size(Tau,2))*Q(:,2:d+1)';
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



function [c, A, P, Theta_h] = Least_Fitting(Data, Tau)
%   Tau = U'*(Data-q);
    T = Construct_Higher_Order(Tau);
    d = size(Tau, 1);
    Theta_h = (build_theta(Data));
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


function theta = build_theta(Data)  
    %theta = diag(sqrt(exp(-sum((Data-q).^2,1)/h^2)));
    theta = diag(ones(1,size(Data, 2)));
end


function [data_true, data] = generate_data(sigma, num)
    theta = linspace(-pi/2, pi/2, num);%pi/4:0.1:3*pi/4;
    data_true = [cos(theta);sin(theta)];
    data = data_true+sigma*randn(2,length(theta));
end