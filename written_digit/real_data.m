clear
n_image = 20;
n_column = 5;
n_row = 4;
data_struc = load('mnist.mat');
addpath('./../methods/')
DataX = data_struc.testX;
DataY = data_struc.testY;
Part_image = DataX(DataY==2,:);
figure(1)
t = tiledlayout(n_row, n_column,'TileSpacing','Compact');
for i = 1:n_image
    nexttile
    image_row = DataX(DataY==mod(i,10), :);
    image_v{i} = reshape(image_row(floor(i/10)*10+i,:),[28, 28])';
    [rows,columns] = find(image_v{i}>0);
    Data_val{i} = image_v{i}(image_v{i}>0);
    Data{i} = [rows'; columns'];
    image(image_v{i});
end

figure(2)
t = tiledlayout(n_row, n_column,'TileSpacing','Compact');
for i = 1:n_image
    nexttile
    image(image_v{i});
    axis off
end

figure(3)
t = tiledlayout(n_row, n_column,'TileSpacing','Compact');
d = 1; D = 2;
neigs = 20; 
Delta = 100; 
Projection = cell(8,n_image);
Result = cell(1,n_image);

for j = 1:n_image
    nexttile
    rho = 1;
    n = floor(rho*neigs);
    W = PixValue_Weight(Data_val{j});
    Projection{1,j} = quadratic(Data{j}, Data{j}, neigs, n, d, Delta);
    Projection{2,j} = quadratic_kernel(Data{j}, Data{j}, neigs,  d, Delta, W);
    Projection{3,j} = PCA_refine(Data{j}, Data{j}, neigs, d);
    Projection{4,j} = linear_KDE(neigs, Data{j}, Data{j}, D, d);
    Projection{5,j} = linear_log_KDE(neigs, Data{j}, Data{j}, D, d);
    Projection{6,j} = linear_mfit(Data{j}, Data{j}, d, neigs);
    Projection{7,j} = MovingLS(Data{j}, Data{j}, neigs, d);
    Projection{8,j} = Spherelet(Data{j}, Data{j}, neigs, d);
end


function A =  New_Image(C, m, n)
    A =  sparse(round(C(1,:)),round(C(2,:)),255*ones(1, size(C,2)),m, n);
    B = full(A);
end
% 
%     image(full(Result{j}));
%     axis off

function W = PixValue_Weight(Data_val)
    W = diag(double(Data_val)/255);
end