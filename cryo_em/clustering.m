%% build affinity matrix
show_or_not = 0;
[clean, noisy] = cryo_em_data(show_or_not);
data2 = clean+0.6*randn(2000,64,64);
%data2 = noisy;
n = size(data2,1);
S = zeros(n);
for i = 1:n
    S(i,i) = 1;
    for j = i+1:n
        a = data2(i,:,:);
        b = data2(j,:,:);
        S(i,j) = a(:)'*b(:)/norm(a(:))/norm(b(:));
        S(j,i) = S(i,j);
    end
end

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







