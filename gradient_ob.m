rho = linspace(0,0.1,100);
for i = 1:100
    [g(i), t1(i), ~] = Eval_Rho(Data, Tau, rho(i), W);
end
%t = tiledlayout(1,2,'TileSpacing','Compact');
plot(rho,-g,'linewidth',2);
hold on
plot(rho,t1.^2,'.-','linewidth',1);
hold on
plot(0.01*ones(1,10000),linspace(-40,300,10000),'--','linewidth',2);
hold on
plot(linspace(0,0.1,1000),ones(1,1000)*100,'--','linewidth',2);
hold on
xlabel('\lambda','FontSize',18)
ylabel('function value','FontSize',18)
legend({'-s`(\lambda)','s(\lambda)'})
axis([0 0.1 -40 300])
set(gca,'FontSize',18)
box on;

%title('Original','FontSize',18)



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
   function R = Construct_Regularization(Tau)
            d = size(Tau, 1);
            R = zeros(1+d+d*(d+1)/2);
            R(d+2:end,d+2:end) = eye(d*(d+1)/2);
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
end