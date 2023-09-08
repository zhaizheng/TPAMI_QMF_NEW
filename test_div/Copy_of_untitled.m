
clear
figure
t = tiledlayout(2,2,'TileSpacing','Compact');

names = {'test_div.mat','div2.mat','test_conv.mat'};

for i = 1
    load(names{i})
    %x = f(A,B,c,[0.4,-0.3]);
    A1 = A(:,1);
    B1 = B(:,1,1);
%     nexttile
% 
%     T1 = conv_div(x,A1,B1,c);
%     nexttile
% 
%     T2 = conv_div(x,A,2*B,c);
%     nexttile
%     T3 = conv_div(x,A1,10*B1,c);
    
    
    [T, Ind_odd, W,Q, q, X ,Y, Z, z1, z2, Qq] = conv_div(x,A1,20*B1,c);
    
    

    [T2, Ind_odd2, W2,Q2, q2, X2, Y2, Z2,z12,z22,Qq2] = conv_div(x,A1,30*B1,c);
    
    draw_line1(T, Ind_odd, W,Q, q, X, Y,Z,z1, z2)
    draw_line1(T2, Ind_odd2, W2,Q2, q2, X2, Y2, Z2,z12, z22)
    
    draw_line2(T, Ind_odd, W,Q, Qq)
    draw_line2(T2, Ind_odd2, W2,Q2, Qq2)
    %legend({'\tau_n','\eta_n'})
end


   function [T, Ind_odd, W,Q,q, X, Y, Z, z1, z2, Qq] =conv_div(x,A,B,c)
        [tau, T] = projection(x,A,B,c,0);
   
        Ind_odd = 1:2:size(T,2)-1;
        q = zeros(1,length(Ind_odd));
        Qq = zeros(1,length(Ind_odd));
        Ind_even = 2:2:size(T,2);
 
         a = min(T(Ind_odd));
         b = max(T(Ind_odd));
         e = min(T(Ind_even));
         d = max(T(Ind_even));
         a = -(b-a)/30+a;
         b = (b-a)/30+b;
         e = -(d-e)/30+e;
         d = (d-e)/30+d;
         z1 = linspace(min(T),max(T),10);
         z2 = linspace(min(T),max(T),10);
        %z1 = min(T(1,:))-0.02:0.01:max(T(1,:))+0.02;
        %z2 = min(T(2,:))-0.02:0.01:max(T(2,:))+0.02;
        [X, Y] = meshgrid(z1,z2);
        k = 0;
        for i = 1:length(Ind_odd)
            q(i) = evalz2(x,A,B,c,T(Ind_odd(i)),T(Ind_odd(i)+1));
            Qq(i) = evalz(x,A,B,c,T(Ind_odd(i)));
        end
        
        W = linspace(min(T),max(T),20);
  
        for i = 1:length(W)
            Q(i) = evalz2(x,A,B,c,W(i),W(i));
        end
        
        for i = 1:length(X(:))
            Z(i) = evalz2(x,A,B,c,X(i),Y(i));
            %text(X(i),Y(i),num2str(z(i),'%.2f'));
        end
        
       
        
   end
 
    function draw_line1(T, Ind_odd, W,Q, q, X, Y, Z, z1, z2)
        
        nexttile

        plot3(T(1,Ind_odd),T(1,Ind_odd+1),q,'bo-','MarkerSize',5,'Linewidth',2)
        hold on
        

        plot3(W,W,Q,'--','Linewidth',2)
       % plot(T(1,Ind_even),'ro-','MarkerSize',5,'Linewidth',2)
        

        hold on
        surf(X,Y,reshape(Z,length(z2),length(z1)),'FaceAlpha',0.3)
  
        set(gca,'FontSize',22)
        axis([min(T) max(T) min(T) max(T)])
        grid off
        box on
        %axis off
        xlabel('\tau');
        ylabel('\eta');
        zlabel('g(\tau,\eta)')
       end
        
       function draw_line2(T, Ind_odd, W, Q, Qq)
            nexttile
            semilogy(T(1,Ind_odd),Qq,'*');
            hold on
            semilogy(W, Q,'-','linewidth',0.5)
            xlabel('\tau')
            ylabel('h(\tau)')
            set(gca,'FontSize',22)
       end


function re = evalz(x, A, B, c, tau)
    for i = 1:size(B,1)
        q(i) = tau*squeeze(B(i,:,:))*tau';
    end
    re = norm(x-A*tau'-c-q')^2;
end

function re = evalz2(x, A, B, c, tau, eta)
    for i = 1:size(B,1)
        q(i) = tau*squeeze(B(i,:,:))*eta';
    end
    re =0.5*norm(x-A*tau'-c-q')^2+0.5*norm(x-A*eta'-c-q')^2;
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
