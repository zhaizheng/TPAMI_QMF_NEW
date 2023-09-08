show_or_not = 0;
[clean, noisy] = cryo_em_data(show_or_not);
vec_clean = matrix_to_vec(clean);
[U,L,~] = svd(vec_clean);
Data = U(:,1:20)*U(:,1:20)'* vec_clean;
image_Data = vec_to_matrix(Data);
%PARAMETERS SETTING
D = 20; d = 3;  tol=1.e-4; kk = [10,15,20,25,30]; NUM = 30; 

for i = 1:5
    subplot(1,5,i)
    sigma = 1;
    low_dimension_noisy_data = U(:,1:D)'* vec_clean+sigma*randn(D,size(vec_clean,2));
    Data_noisy = U(:,1:D)*low_dimension_noisy_data;
    image_Data2 = vec_to_matrix(Data_noisy);
    if i~=1
        image(squeeze(image_Data2(1,:,:)*200));
    else
        image(squeeze(clean(1,:,:)*200));
    end
end

function mat = vec_to_matrix(Data)
    mat = zeros(size(Data,2),64,64);
    for i = 1:size(Data,2)
        mat(i,:,:) = reshape(Data(:,i),64,[]);
    end
end

function vec = matrix_to_vec(Data)
    vec = [];
    for i = 1:size(Data,1)
        temp = Data(i,:,:);
        vec = [vec, temp(:)];
    end
end