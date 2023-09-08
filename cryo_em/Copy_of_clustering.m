% %% build affinity matrix
% show_or_not = 0;
% [clean, noisy] = cryo_em_data(show_or_not);
% data2 = clean+0.6*randn(2000,64,64);
% %data2 = noisy;
% n = size(data2,1);
% S = zeros(n);
% for i = 1:n
%     S(i,i) = 1;
%     for j = i+1:n
%         a = data2(i,:,:);
%         b = data2(j,:,:);
%         S(i,j) = a(:)'*b(:)/norm(a(:))/norm(b(:));
%         S(j,i) = S(i,j);
%     end
% end
%%
show_or_not = 0;
[clean, noisy] = cryo_em_data(show_or_not);
vec_clean = matrix_to_vec(clean);
[U,L,~] = svd(vec_clean);
Data = U(:,1:20)*U(:,1:20)'* vec_clean;
image_Data = vec_to_matrix(Data);



Data_noisy = U(:,1:20)*(U(:,1:20)'* vec_clean+randn(20,size(vec_clean,2)));
image_Data2 = vec_to_matrix(Data_noisy);


low_dimension_noisy_data = U(:,1:20)'* vec_clean+randn(20,size(vec_clean,2));

k = 15;   d = 3; 
PCA_Data_Low_dimension = PCA_refine(low_dimension_noisy_data, low_dimension_noisy_data, k, d);
PCA_Data_High = U(:,1:20)*PCA_Data_Low_dimension;
PCA_img = vec_to_matrix(PCA_Data_High);



D = 20; d = 3;  tol=1.e-4; k = 15;
Quadratic_Low_dimension = Quadratic_Manifold(k, low_dimension_noisy_data, low_dimension_noisy_data(:,1:16), D, d, tol);
Quadratic_High = U(:,1:20)*Quadratic_Low_dimension;
Quadratic_img = vec_to_matrix(Quadratic_High);


KDE_Fitting = linear_KDE(sigma, k, Data, inputs, D, d);
log_KDE_Fitting = linear_log_KDE(sigma, k, Data, inputs, D, d);
mfit_Fitting = linear_mfit(Data, inputs, d, sigma, 0.4);


%%
subplot('position',[0.02 0.05 0.22 0.9])
show_im(4,4,image_Data(1:20,:,:));
xlabel('Noiseless Images')
set(gca,'XTick',[], 'YTick', [])
subplot('position',[0.27 0.05 0.22 0.9])
show_im(4,4,image_Data2(1:20,:,:));
xlabel('Noisy Images')
set(gca,'XTick',[], 'YTick', [])
subplot('position',[0.52 0.05 0.22 0.9])
show_im(4, 4, PCA_img)
xlabel('Linear Recovery')
set(gca,'XTick',[], 'YTick', [])
subplot('position',[0.77 0.05 0.22 0.9])
show_im(4, 4, Quadratic_img)
xlabel('Quadratic Recovery')
set(gca,'XTick',[], 'YTick', [])

r1 = mean(mean(mean((image_Data(1:16,:,:)-image_Data2(1:16,:,:)).^2)));
r2 = mean(mean(mean((image_Data(1:16,:,:)-PCA_img(1:16,:,:)).^2)));
r3 = mean(mean(mean((image_Data(1:16,:,:)-Quadratic_img(1:16,:,:)).^2)));


fprintf('Original average distance:%f, PCA=%f ,Quadritic distance=%f\n',r1,r2,r3)


%%
[V, ~, ~] = svds(S,20);
%[~,indx] = sort(diag(L),'descend');
%V = U(:,indx);
idx = kmeans(V,20);

%%
%hist(idx);

subplot('position',[0.05 0.55 0.43 0.43])
class = 14;
ind = find(idx == class);
m = 4;n = 4;
show_im(m, n, data2(ind,:,:))

 

subplot('position',[0.55 0.05 0.43 0.43])
Data = data2(idx == class,:,:);

Data_vec = matrix_to_vec(Data);
[U,~,~] = svds(Data_vec, 10);

ambient_dim = 10; k = 15;   d = 5; tol=1.e-4;

%Data_new = U*U'*Data_vec;
Data_new = PCA_refine(Data_vec, Data_vec, k, d);
PCA_img = vec_to_matrix(Data_new);
show_im(m, n, PCA_img)



subplot('position',[0.05 0.05 0.43 0.43])
Data_input = U(:,1:10)'*Data_vec;



D = size(Data_vec, 1); 
points_projection = Quadratic_Manifold(k, Data_vec, Data_vec, D, d, tol);

Data_output = points_projection;
Quadratic_img = vec_to_matrix(Data_output);
show_im(m, n, Quadratic_img);

subplot('position',[0.55 0.55 0.43 0.43])
show_im(m, n, clean(ind,:,:))

r1 = mean(mean(mean((data2(ind,:,:)-clean(ind,:,:)).^2)));
r2 = mean(mean(mean((PCA_img-clean(ind,:,:)).^2)));
r3 = mean(mean(mean((Quadratic_img-clean(ind,:,:)).^2)));


fprintf('Original average distance:%f, PCA=%f ,Quadritic distance=%f\n',r1,r2,r3)


%%
result_each_class = zeros(3,20);

%

% for class = 10
%     
%     ind = find(idx == class);
%     
%     Data = data2(idx == class,:,:);
%     Data_vec = matrix_to_vec(Data);
%     [U,~,~] = svds(Data_vec,10);
%     
%     %PCA
%     %Data_new = U(:,1:10)*U(:,1:10)'*Data_vec;
%     %PCA_image = vec_to_matrix(Data_new);
%     PCA_image_vec = PCA_refine(Data_vec, Data_vec, k, D);
%     PCA_image = vec_to_matrix(PCA_image_vec);
%     
%     %Quadratic
%     Data_input = U'*Data_vec;
%     points_projection = Quadratic_Manifold(k, Data_input, Data_input, D, d, tol);
%     Data_output = U*points_projection;
%     mat = vec_to_matrix(Data_output);
%     
%     result_each_class(1,class) = mean(mean(mean((data2(ind,:,:)-clean(ind,:,:)).^2)));
%     result_each_class(2,class) = mean(mean(mean((PCA_image-clean(ind,:,:)).^2)));
%     result_each_class(3,class) = mean(mean(mean((mat-clean(ind,:,:)).^2))); 
% end


%%


function new_points_linear = linear_KDE(sigma, k, Data, points, D, d)
        new_points_linear = zeros(size(points));
        for i = 1:size(points,2)
            data_old = zeros(size(points(:,i)));
            data_move = points(:,i);
            sk = 1;
            while norm(data_old-data_move)>10.e-3 && sk<50
                data_old = data_move;
                x = shift_mean(data_move, Data, sigma, k);
                [~, ~, ~, ~, U_v] = fitting(data_move, Data, sigma, D, d);
                data_move = data_old + U_v*U_v'*(x-data_old);
                sk = sk+1;
            end
            new_points_linear(:,i) = data_move;
        end
end
   

function new_points_linear = linear_log_KDE(sigma, k, Data, points, D, d)
        new_points_linear = zeros(size(points));
        for i = 1:size(points,2)
            data_old = zeros(size(points(:,i)));
            data_move = points(:,i);
            sk = 1;
            while norm(data_old-data_move)>10.e-3 && sk<50
                data_old = data_move;
                x = shift_mean(data_move, Data, sigma, k);
                [~, ~, ~, ~, U_v] = fitting(x, Data, sigma, D, d);
                data_move = data_old + U_v*U_v'*(x-data_old);
                sk = sk+1;
            end
            new_points_linear(:,i) = data_move;
        end
end


function new_points_linear = linear_mfit(Data, points, d, beta)
        new_points_linear = zeros(size(points));
        for i = 1:size(points,2)
            data_old = zeros(size(points(:,i)));
            data_move = points(:,i);
            sk = 1;
            while norm(data_old-data_move)>10.e-3 && sk<50
                data_old = data_move;               
%                 x = shift_mean(data_move, Data, sigma, k);
%                 [~, ~, ~, ~, U_v] = fitting(data_move, Data, sigma, D, d);
%                 data_move = data_old + U_v*U_v'*(x-data_old);
                step = 1; x = data_old; r = 0.3;  data = Data;
                direction = mfit(x, data, r, beta, d, step);
                data_move = data_old + step*direction;              
                sk = sk+1;
            end
            new_points_linear(:,i) = data_move;
        end
end


function [Theta, Co_h, Co_v, U_h, U_v] = fitting(x, Data, sigma, D, d)

    %mean = shift_mean(x, Data, sigma, k);
    [Co_h, Co_v, U_h, U_v] = coordinate(x, Data, sigma, D, d);
    W = build_W(x, Data, sigma);
    [~, Theta] = least_square(Co_h, Co_v, W);
end


function mean = shift_mean(x, Data, h, k)
    [~,ind] = sort(sum((Data-x).^2,1),'ascend');
    mean = zeros(size(x));
    s_weight = 0;
    for i = 1:k
        w = exp(-norm(Data(:,ind(i))-x)^2/(h^2));
        %w = 1;
        s_weight = s_weight+w;
        mean = mean+w*Data(:,ind(i));
    end
    mean = mean/s_weight;
end


function Y_out = PCA_refine(Data_in, Data, k, d)
    Y_out = zeros(size(Data_in));
    for i = 1:size(Y_out,2)
        [~,ind] = sort(sum((Data-Data_in(:,i)).^2,1),'ascend');
        [U,~] = svds(Data(:,ind(2:k+1)), d);
        Y_out(:,i) = U*U'*Data_in(:,i);
    end
end


%%
function vec = matrix_to_vec(Data)
    vec = [];
    for i = 1:size(Data,1)
        temp = Data(i,:,:);
        vec = [vec, temp(:)];
    end
end


function mat = vec_to_matrix(Data)
    mat = zeros(size(Data,2),64,64);
    for i = 1:size(Data,2)
        mat(i,:,:) = reshape(Data(:,i),64,[]);
    end
end



function show_im(m, n, data2)
    im = zeros(64*m,64*n);
    for i = 1:m
        for j = 1:n
            im((i-1)*64+1:i*64,(j-1)*64+1:j*64) = data2(((i-1)*n+j),:,:);
        end
    end
    image(im*200);
end







