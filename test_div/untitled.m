
clear

t = tiledlayout(1,2,'TileSpacing','Compact');

names = {'test_div.mat','div2.mat','test_conv.mat'};

for i = 1
    load(names{i})
    %x = f(A,B,c,[0.4,-0.3]);
    nexttile

    T1 = conv_div(x,A,B,c);
%     nexttile
% 
%     T2 = conv_div(x,A,2*B,c);
    nexttile
    T3 = conv_div(x,A,10*B,c);
    legend({'\tau_n','\eta_n'})
end


   function T=conv_div(x,A,B,c)
        [tau, T] = projection(x,A,B,c,[0;0]);
        q = zeros(1,size(T,2));

        a = min(T(1,:));
        b = max(T(1,:));
        e = min(T(2,:));
        d = max(T(2,:));
        a = -(b-a)/10+a;
        b = (b-a)/10+b;
        e = -(d-e)/10+e;
        d = (d-e)/10+d;
        z1 = linspace(a,b,10);
        z2 = linspace(e,d,10);
        %z1 = min(T(1,:))-0.02:0.01:max(T(1,:))+0.02;
        %z2 = min(T(2,:))-0.02:0.01:max(T(2,:))+0.02;
        [X, Y] = meshgrid(z1,z2);
        for i = 1:size(T,2)
            q(i) = evalz(x,A,B,c,T(:,i)');
        end
        T = [T;q];
        Ind_odd = 1:2:size(T,2);
        Ind_even = 2:2:size(T,2);
        plot3(T(1,Ind_odd),T(2,Ind_odd),q(Ind_odd),'bo','MarkerSize',5,'Linewidth',2)
        hold on
        plot3(T(1,Ind_even),T(2,Ind_even),q(Ind_even),'ro','MarkerSize',5,'Linewidth',2)
        
        for i = 1:length(X(:))
            Z(i) = evalz(x,A,B,c,[X(i),Y(i)]);
            %text(X(i),Y(i),num2str(z(i),'%.2f'));
        end
        hold on
        surf(X,Y,reshape(Z,length(z2),length(z1)),'FaceAlpha',0.3)
        set(gca,'FontSize',14)
        axis([a b e d])
        box on
    end
 


function re = evalz(x, A, B, c, tau)
    for i = 1:size(B,1)
        q(i) = tau*squeeze(B(i,:,:))*tau';
    end
    re = norm(x-A*tau'-c-q')^2;
end

function re = f(A,B,c,tau)
    for i = 1:size(B,1)
        q(i) = tau*squeeze(B(i,:,:))*tau';
    end
    re = A*tau'+c+q';

end


function [tau, Tau] = projection(x, A, B, c, tau) %project x onto f(tau) = A tau+ B(tau,tau)+c
    Tau = [];
    iter = 0; 
    while true
        Bm = tensor_fold(B, tau);
        tau_new = (2*Bm'*Bm+Bm'*A+A'*A+A'*Bm)\((2*Bm'+A')*(x-c)-Bm'*A*tau);
        if norm(tau_new-tau)<1.e-6 || iter>1000
            if iter>1000
                fprintf('diverge projecting tau\n');
            end
            break;
        end
        tau = tau_new;
        Tau = [Tau, tau];
        iter = iter+1;
    end 
end

function result = tensor_fold(B, tau)
    result = zeros(size(B,1),size(B,2));
    for i = 1:size(B,1)
        result(i,:) = squeeze(B(i,:,:))*tau;
    end
end
