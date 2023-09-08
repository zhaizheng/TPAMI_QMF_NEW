%% DATA CONSTRUCTION
cd('/home/zzh/Documents/TPAMI_QMF_new/code/cryo_em')
addpath('./../methods');
show_or_not = 0;
[clean, noisy] = cryo_em_data(show_or_not);
vec_clean = matrix_to_vec(clean);
[U,L,~] = svd(vec_clean);
Data = U(:,1:20)*U(:,1:20)'* vec_clean;
image_Data = vec_to_matrix(Data);
%PARAMETERS SETTING
%D = 20; d = 3;  tol=1.e-4; kk = [5,10,15,20,25]; NUM = 30; sigma = 4;


D = 20; d = 5;  tol=1.e-4; 
NUM = 20; sigma = 0.5;

low_dimension_noisy_data = U(:,1:D)'* vec_clean+sigma*randn(D,size(vec_clean,2));
Data_noisy = U(:,1:D)*low_dimension_noisy_data;
image_Data2 = vec_to_matrix(Data_noisy);


result = zeros(6,5,2);
%%

kk = [40,44,48,52,56,60,64,68];
for i = 1:length(kk)
    k = kk(i);
    %%PCA DENOSING
    fprintf('PCA in process:\n')
    PCA_Data_Low_dimension = PCA_refine_centralized(low_dimension_noisy_data(:,1:NUM), low_dimension_noisy_data, k, d);
    PCA_Data_High = U(:,1:D)*PCA_Data_Low_dimension;
    PCA_img = vec_to_matrix(PCA_Data_High);
    
    %%QUADRATIC APPROACH
    %Quadratic_Low_dimension = Quadratic_Manifold(k, low_dimension_noisy_data, low_dimension_noisy_data(:,1:NUM), D, d, tol);
    delta = 50;
    fprintf('Quadratic in Process:\n')
    W = diag(ones(1,size(low_dimension_noisy_data,2))); %n = 3*k;
    Quadratic_Low_dimension = quadratic(low_dimension_noisy_data,low_dimension_noisy_data(:,1:NUM),k, k, d, delta);
    
    Quadratic_Kernel = quadratic_kernel(low_dimension_noisy_data,low_dimension_noisy_data(:,1:NUM),k,d,delta);
    
    %Quadratic_Low_dimension = adaptive_fit(low_dimension_noisy_data, low_dimension_noisy_data(:,1:NUM), W, k, n, d, delta);

    %Projection_MovingLS(low_dimension_noisy_data, low_dimension_noisy_data(:,1:NUM), k, d);
    Quadratic_High = U(:,1:D)*Quadratic_Low_dimension;
    Quadratic_kernel = U(:,1:D)*Quadratic_Kernel;
    
    Quadratic_img = vec_to_matrix(Quadratic_High);
    
    Quadratic_img_kernel = vec_to_matrix(Quadratic_kernel);

    %%KDE APPROACH
    fprintf('KDE in Process:\n')
    KDE_Fitting = linear_KDE(k, low_dimension_noisy_data, low_dimension_noisy_data(:,1:NUM), D, d);
    KDE_Fitting_High = U(:,1:D)*KDE_Fitting;
    KDE_img = vec_to_matrix(KDE_Fitting_High);

    %%LOG-KDE APPROACH
    fprintf('log_KDE in Process:\n')
    log_KDE_Fitting = linear_log_KDE(k, low_dimension_noisy_data, low_dimension_noisy_data(:,1:NUM), D, d);
    log_KDE_Fitting_High = U(:,1:D)*log_KDE_Fitting;
    log_KDE_img = vec_to_matrix(log_KDE_Fitting_High);


    
    %%MFIT APPROACH
    fprintf('mfint in Process:\n')
    mfit_Fitting = linear_mfit(low_dimension_noisy_data, low_dimension_noisy_data(:,1:NUM), d, k);
    mfit_Fitting_High = U(:,1:D)*mfit_Fitting;
    mfit_img = vec_to_matrix(mfit_Fitting_High);
    
    
    %%MovingLS
    fprintf('MLS in Process:\n')
    MLS_Fitting = MovingLS(low_dimension_noisy_data, low_dimension_noisy_data(:,1:NUM), k, d);
    MLS_Fitting_High = U(:,1:D)*MLS_Fitting.result;
    MLS_img = vec_to_matrix(MLS_Fitting_High);


    fprintf('Spherelet in Process:\n')
    Sphere_Fitting = Spherelet(low_dimension_noisy_data(:,1:NUM), low_dimension_noisy_data, k, d);
    Sphere_Fitting_High = U(:,1:D)*Sphere_Fitting;
    Sphere_img = vec_to_matrix(Sphere_Fitting_High);

    [result(1,i,1),result(1,i,2)] = measure_distance(image_Data(1:NUM,:,:),image_Data2(1:NUM,:,:));
    [result(2,i,1),result(2,i,2)] = measure_distance(image_Data(1:NUM,:,:),Quadratic_img);
    [result(3,i,1),result(3,i,2)] = measure_distance(image_Data(1:NUM,:,:),Quadratic_img_kernel);
    [result(4,i,1),result(4,i,2)] = measure_distance(image_Data(1:NUM,:,:),PCA_img);
    [result(5,i,1),result(5,i,2)] = measure_distance(image_Data(1:NUM,:,:),KDE_img);
    [result(6,i,1),result(6,i,2)] = measure_distance(image_Data(1:NUM,:,:),log_KDE_img);
    [result(7,i,1),result(7,i,2)] = measure_distance(image_Data(1:NUM,:,:), mfit_img);
    [result(8,i,1),result(8,i,2)] = measure_distance(image_Data(1:NUM,:,:), MLS_img);
    [result(9,i,1),result(9,i,2)] = measure_distance(image_Data(1:NUM,:,:), Sphere_img);
end


%%
% figure(1)
% subplot('position',[0.02 0.05 0.22 0.9])
% show_im(4,4,image_Data(1:20,:,:));
% title('Noiseless Images')
% set(gca,'XTick',[], 'YTick', [])
% subplot('position',[0.27 0.05 0.22 0.9])
% show_im(4,4,image_Data2(1:20,:,:));
% title('Noisy Images')
% set(gca,'XTick',[], 'YTick', [])
% subplot('position',[0.52 0.05 0.22 0.9])
% show_im(4, 4, PCA_img)
% title('Linear Recovery')
% set(gca,'XTick',[], 'YTick', [])
% subplot('position',[0.77 0.05 0.22 0.9])
% show_im(4, 4, Quadratic_img)
% title('Quadratic Recovery')
% set(gca,'XTick',[], 'YTick', [])


% figure(2)
% subplot('position',[0.02 0.05 0.22 0.9])
% show_im(4,4,KDE_img(1:20,:,:));
% title('KDE Recovery')
% set(gca,'XTick',[], 'YTick', [])
% subplot('position',[0.27 0.05 0.22 0.9])
% show_im(4,4,log_KDE_img(1:20,:,:));
% title('LOG-KDE Recovery')
% set(gca,'XTick',[], 'YTick', [])
% subplot('position',[0.52 0.05 0.22 0.9])
% show_im(4,4, mfit_img)
% title('Mfit Recovery')
% set(gca,'XTick',[], 'YTick', [])
% subplot('position',[0.77 0.05 0.22 0.9])
% show_im(4,4, MLS_img)
% title('Moving-LS Recovery')
% set(gca,'XTick',[], 'YTick', [])

%fprintf('Original average distance:%f, PCA=%f ,Quadritic distance=nexttile
show_im(4,4, MLS_img)
title('Moving LS','FontSize',18)
set(gca,'XTick',[], 'YTick', [])%f, KDE=%f, log_KDE=%f,mfit=%f\n',r1,r2,r3,r4,r5,r6)

%%
figure
t = tiledlayout(2,5,'TileSpacing','compact');
nexttile
show_im(4,5,image_Data(1:20,:,:));
title('Noiseless Images','FontSize',18)
set(gca,'XTick',[], 'YTick', [])

nexttile
show_im(4,5,image_Data2(1:20,:,:));
title('Noisy Images','FontSize',18)
set(gca,'XTick',[], 'YTick', [])

nexttile
show_im(4, 5, Quadratic_img)
title('RQMF-E','FontSize',18)
set(gca,'XTick',[], 'YTick', [])

nexttile
show_im(4, 5, Quadratic_img_kernel)
title('RQMF-K','FontSize',18)
set(gca,'XTick',[], 'YTick', [])

nexttile
show_im(4, 5, PCA_img)
title('Local PCA','FontSize',18)
set(gca,'XTick',[], 'YTick', [])


nexttile
show_im(4,5, MLS_img)
title('Moving LS','FontSize',18)
set(gca,'XTick',[], 'YTick', [])
show_im(4,5,KDE_img);
title('KDE','FontSize',18)
set(gca,'XTick',[], 'YTick', [])

nexttile
show_im(4,5,log_KDE_img);
title('LOG-KDE','FontSize',18)
set(gca,'XTick',[], 'YTick', [])


nexttile
show_im(4,5, mfit_img)
title('Mfit','FontSize',18)
set(gca,'XTick',[], 'YTick', [])

nexttile
show_im(4,5, MLS_img)
title('Moving LS','FontSize',18)
set(gca,'XTick',[], 'YTick', [])

nexttile
show_im(4,5, Sphere_img)
title('SPH-PCA','FontSize',18)
set(gca,'XTick',[], 'YTick', [])


%%
final_result = [];
for i = 1:8
    final_result = [final_result,[result(9,i,1),result(9,i,2)]];
end
format_print(final_result)

% figure(2)
% subplot('position',[0.02 0.05 0.22 0.9])
% show_im(4,4,KDE_img(1:20,:,:));
% title('KDE Recovery')
% set(gca,'XTick',[], 'YTick', [])
% 
% 
% subplot('position',[0.77 0.05 0.22 0.9])
% show_im(4,4, MLS_img)
% title('Moving-LS Recovery')
% set(gca,'XTick',[], 'YTick', [])


%%

function [dis,sd] = measure_distance(A, T)
    dis = norm((A(:)-T(:)),'fro')^2/size(A,1);
    s = [];
    for i = 1:size(A)
        s(i) = norm(squeeze(A(i,:,:))-squeeze(T(i,:,:)),'fro')^2;
    end
    sd = std(s,1); 
end
% 
% 
% function new_points_linear = linear_KDE(k, Data, points, D, d)
%         new_points_linear = zeros(size(points));
%         for i = 1:size(points,2)
%             data_old = zeros(size(points(:,i)));
%             data_move = points(:,i);
%             sk = 1;
%             while norm(data_old-data_move)>10.e-3 && sk<50
%                 data_old = data_move;
%                 sigma = find_sigma(data_move, Data, k);
%                 x = shift_mean(data_move, Data, sigma, k);
%                 [~, ~, ~, ~, U_v] = fitting(data_move, Data, sigma, D, d);
%                 data_move = data_old + U_v*U_v'*(x-data_old);
%                 sk = sk+1;
%             end
%             new_points_linear(:,i) = data_move;
%         end
% end
%    
% 
% function new_points_linear = linear_log_KDE(k, Data, points, D, d)
%         new_points_linear = zeros(size(points));
%         for i = 1:size(points,2)
%             data_old = zeros(size(points(:,i)));
%             data_move = points(:,i);
%             sk = 1;
%             while norm(data_old-data_move)>10.e-3 && sk<50
%                 data_old = data_move;
%                 sigma = find_sigma(data_move, Data, k);
%                 x = shift_mean(data_move, Data, sigma, k);
%                 [~, ~, ~, ~, U_v] = fitting(x, Data, sigma, D, d);
%                 data_move = data_old + U_v*U_v'*(x-data_old);
%                 sk = sk+1;
%             end
%             new_points_linear(:,i) = data_move;
%         end
% end
% 
% 
% function new_points_linear = linear_mfit(Data, points, d, k)
%         new_points_linear = zeros(size(points));
%         for i = 1:size(points,2)
%             data_old = zeros(size(points(:,i)));
%             data_move = points(:,i);
%             sk = 1;
%             while norm(data_old-data_move)>10.e-3 && sk<50
%                 data_old = data_move;               
% %                 x = shift_mean(data_move, Data, sigma, k);
% %                 [~, ~, ~, ~, U_v] = fitting(data_move, Data, sigma, D, d);
% %                 data_move = data_old + U_v*U_v'*(x-data_old);
%                 step = 0.2; x = data_old; beta = 0.1;  data = Data;
%                 r = find_sigma(data_move, Data, k);
%                 direction = mfit(x, data, r, beta, d, step);
%     %            direction = xia(x, data, r, beta, d, step);              
%                 data_move = data_old + step*direction;              
%                 sk = sk+1;
%             end
%             new_points_linear(:,i) = data_move;
%         end
% end
% 
% 
% function [Theta, Co_h, Co_v, U_h, U_v] = fitting(x, Data, sigma, D, d)
% 
%     %mean = shift_mean(x, Data, sigma, k);
%     [Co_h, Co_v, U_h, U_v] = coordinate(x, Data, sigma, D, d);
%     W = build_W(x, Data, sigma);
%     [~, Theta] = least_square(Co_h, Co_v, W);
% end
% 
% 
% function mean = shift_mean(x, Data, h, k)
%     [~,ind] = sort(sum((Data-x).^2,1),'ascend');
%     mean = zeros(size(x));
%     s_weight = 0;
%     for i = 1:k
%         w = exp(-norm(Data(:,ind(i))-x)^2/(h^2));
%         %w = 1;
%         s_weight = s_weight+w;
%         mean = mean+w*Data(:,ind(i));
%     end
%     mean = mean/s_weight;
% end
% 
% 
% function Y_out = PCA_refine(Data_in, Data, k, d)
%     Y_out = zeros(size(Data_in));
%     for i = 1:size(Y_out,2)
%         [~,ind] = sort(sum((Data-Data_in(:,i)).^2,1),'ascend');
%         [U,~] = svds(Data(:,ind(2:k+1)), d);
%         Y_out(:,i) = U*U'*Data_in(:,i);
%     end
% end
% 
% 
% %%
function vec = matrix_to_vec(Data)
    vec = [];
    for i = 1:size(Data,1)
        temp = Data(i,:,:);
        vec = [vec, temp(:)];
    end
end
% 
% 
function mat = vec_to_matrix(Data)
    mat = zeros(size(Data,2),64,64);
    for i = 1:size(Data,2)
        mat(i,:,:) = reshape(Data(:,i),64,[]);
    end
end
% 
% 
% 
function show_im(m, n, data2)
    im = zeros(64*m,64*n);
    for i = 1:m
        for j = 1:n
            im((i-1)*64+1:i*64,(j-1)*64+1:j*64) = data2(((i-1)*n+j),:,:);
        end
    end
    image(im*200);
end
% 
% 
% 
% 
% 
% function sigma = find_sigma(x, Data, k)
%     s_distance = sum((Data-x).^2, 1);
%     [~,ind] = sort(s_distance,'ascend');
%     Neig = Data(:,ind(2:k+1)); 
%     sigma = mean(sqrt(sum((Neig-x).^2,1)));
% end
% 
% 
% 
% function [Co_h, Co_v, U_h, U_v] = coordinate(x, A, h, D, d)
%         C = zeros(size(x,1),size(x,1));
%         for i = 1: size(A,2)
%             a = (x-A(:,i));
%             C = C + exp(-norm(a)^2/(h^2))*a*a';
%         end
%         [V,~,~] = svd(C);
%         U_h = V(:,1:d);
%         U_v = V(:,d+1:D);
%         Co_h = U_h'*(A-repmat(x,1,size(A,2)));
%         Co_v = U_v'*(A-repmat(x,1,size(A,2)));
% end
% 
% 
% function W = build_W(x, A, h)
%         centered = A-repmat(x,[1,size(A,2)]);
%         W = zeros(size(centered,2));
%         for k = 1:size(centered,2)
%             W(k,k) = exp(-norm(centered(:,k))^2/(h^2));
%         end
% end
% 
% 
% function [theta, Theta] = least_square(Tau, Co_v, W)
%         d = size(Tau, 1);
%         ind = triu(true(size(Tau, 1)));
%         d2 = d*(d+1)/2;
%         Theta = cell(1,size(Co_v,1));
%         for j = 1:size(Co_v, 1)
%             G = zeros(d2, size(Tau,2));
%             for i = 1:size(Tau,2)
%                 A = Tau(:,i)*Tau(:,i)';
%                 G(:,i) = A(ind);
%             end
%             theta = (G*W*G')\G*W*Co_v(j,:)';
% 
%             Theta_temp = zeros(d,d);
%             Theta_temp(ind) = theta/2;
%             Theta{j} = Theta_temp+Theta_temp';
%         end
% end
% 
% 
% function direction = mfit(x, data, r, beta, dim, step)
%         d = size(data, 1);
%         n = size(data, 2);
%         d2 = 1 - sum((data-x).^2, 1)./(r^2);
%         d2(d2<0) = 0;
%         alpha_tilde = d2.^beta;
%         alpha_sum = sum(alpha_tilde(:));
%         alpha = alpha_tilde./alpha_sum;
%         %cx = sum(data.* alpha, 2);
%         Ns = zeros(d);
%         c_vec = zeros(d,1);
%         for i = 1:n
%             if d2(i)>0
%                 ns = normal_space(data, i, r, dim);
%                 Ns = Ns+ns*alpha(i);
%                 c_vec = c_vec+alpha(i)*ns*(data(:,i)-x);
%             end    
%         end
%         [U,~,~] = svd(Ns);
%         P = U(:,1:d-dim)*U(:,1:d-dim)';
%         direction = step*P*c_vec;
% end
% 
% 
% function P = normal_space(data, i, r, dim)
%     ds = sum((data-data(:,i)).^2, 1);
%     ds(ds > r^2) = 0;
%     n = size(data, 2);
%     d = size(data, 1);
%     indicator = zeros([1,n]);
%     indicator(ds>0) = 1;
%     select = (data-data(:,i)).*indicator;
%     cor = select*select';
%     [U,~,~] = svd(cor);
%     P = eye(d)-U(:,1:dim)*U(:,1:dim)';     
% end
% 
% function direction = xia(x, data, r, beta, dim, step)
%     d = size(data, 1);
%     n = size(data, 2);
%     d2 = 1 - sum((data-x).^2, 1)./(r^2);
%     d2(d2<0) = 0;
%     alpha_tilde = d2.^beta;
%     alpha_sum = sum(alpha_tilde(:));
%     alpha = alpha_tilde./alpha_sum;
%     cx = sum(data.* alpha, 2);
%     Ns = zeros(d);
%     for i = 1:n
%         if d2(i)>0
%             ns = normal_space(data, i, r, dim);
%             Ns = Ns+ns*alpha(i);
%         end    
%     end
%     [U,~,~] = svd(Ns);
%     P = U(:,1:d-dim)*U(:,1:d-dim)';
%     direction = step*P*(cx - x );
% end

