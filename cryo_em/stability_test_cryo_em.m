%% DATA CONSTRUCTION
addpath('./../methods');
cd('/Users/zhengzhai/Dropbox (Personal)/Work/MF/code/cryo_em/');
show_or_not = 0;
[clean, noisy] = cryo_em_data(show_or_not);
vec_clean = matrix_to_vec(clean);
[U,L,~] = svd(vec_clean);
Data = U(:,1:20)*U(:,1:20)'* vec_clean;
image_Data = vec_to_matrix(Data);


D = 20; d = 5;  tol=1.e-4; NUM = 20; sigma = 0.5;
Neig = [12,16,20,24,28,32,36,40,44,48,52,56,60];
Delta = [0.1, 1, 5,10, 50, 100, 300];


low_dimension_noisy_data = U(:,1:D)'*vec_clean+sigma*randn(D,size(vec_clean,2));
Data_noisy = U(:,1:D)*low_dimension_noisy_data;
image_Data2 = vec_to_matrix(Data_noisy);

result = zeros(length(Neig),length(Delta),2);
%%
Neig = [12,16,20,24,28,32,36,40,44,48,52,56,60];
for i = 12:13
    for j = 1:length(Delta)
        delta = Delta(j);
        k = Neig(i);
        Quadratic_Low_dimension = quadratic(low_dimension_noisy_data,low_dimension_noisy_data(:,1:NUM),k, k, d, delta);
        Quadratic_Kernel = quadratic_kernel(low_dimension_noisy_data,low_dimension_noisy_data(:,1:NUM),k, d, delta);
        
        Quadratic_High = U(:,1:D)*Quadratic_Low_dimension;
        Quadratic_kernel = U(:,1:D)*Quadratic_Kernel;
    
        Quadratic_img = vec_to_matrix(Quadratic_High);
        Quadratic_img_kernel = vec_to_matrix(Quadratic_kernel);
        [result(i,j,1),~] = measure_distance(image_Data(1:NUM,:,:),Quadratic_img);
        [result(i,j,2),~] = measure_distance(image_Data(1:NUM,:,:),Quadratic_img_kernel);
        fprintf('\n')
        format_print(squeeze(result(:,:,1)))
        fprintf('\n')
        format_print(squeeze(result(:,:,2)))
    end
end



function [dis,sd] = measure_distance(A, T)
    dis = norm((A(:)-T(:)),'fro')^2/size(A,1);
    s = [];
    for i = 1:size(A)
        s(i) = norm(squeeze(A(i,:,:))-squeeze(T(i,:,:)),'fro')^2;
    end
    sd = std(s,1); 
end


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