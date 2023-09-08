data_struc = load('mnist.mat');
addpath('./../methods/')
DataX = data_struc.testX;
DataY = data_struc.testY;
Part_image = DataX(DataY==2,:);
row_size = 2;
column_size = 5;
dim = 20;
large_image_1 = zeros(28*row_size,28*column_size);
large_image_2 = zeros(28*row_size,28*column_size);
large_image_3 = zeros(28*row_size,28*column_size);
[~,~,U] = svd(double(Part_image));
V = U(:,1:dim);
low_rank_data = V'*double(Part_image)'+0.1*randn(dim, size(Part_image,1));
interest_data = low_rank_data(:,1:row_size*column_size);

P = U(:,1:40)*U(:,1:40)';

subplot(1,3,1)
for i = 0:row_size-1
    for j = 0:column_size-1
        large_image_1(i*28+1:(i+1)*28, j*28+1:(j+1)*28) = reshape(double(V*interest_data(:,i*column_size+j+1)), [28, 28])';
    end
end
image(uint8(large_image_1))

subplot(1,3,2)
for i = 0:row_size-1
    for j = 0:column_size-1
        large_image_2(i*28+1:(i+1)*28, j*28+1:(j+1)*28) = reshape(Part_image(i*column_size+j+1,:), [28, 28])';
    end
end
image(uint8(large_image_2))


subplot(1,3,3)
delta = 2;
k = 20;
d = 10;

Quadratic_Low_dimension = adaptive_fit(low_rank_data, interest_data, k, d, delta);
for i = 0:row_size-1
    for j = 0:column_size-1
        large_image_3(i*28+1:(i+1)*28, j*28+1:(j+1)*28) = reshape(V*Quadratic_Low_dimension(:,i*column_size+j+1), [28, 28])';
    end
end
image(uint8(large_image_3))

