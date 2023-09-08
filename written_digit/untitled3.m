clear
data_struc = load('mnist.mat');
addpath('./../methods/')
DataX = data_struc.testX;
DataY = data_struc.testY;
Part_image = DataX(DataY==2,:);
for i = 1:9%1:size(Part_image,2)
    image_row = DataX(DataY==i, :);
    image_v = reshape(image_row(i,:),[28, 28])';
    [rows,columns] = find(image_v>0);
    Data_val = image_v(image_v>0);
    Data = [columns';-rows'];
end

plot(columns,-rows,'o');


d = 1;
neigs = 20; 
neigs = [15, 20, 25];
Delta = 0.5; 
Delta = [0.125, 0.25, 0.5, 1];
Projection = cell(5,5);
for k = 1:length(neigs)
    for s = 1:length(Delta)
        rho = 1.4;
        n = floor(rho*neigs(k));
        W = PixValue_Weight(Data_val);
        Projection{k,s} = adaptive_fit(Data, Data, W, neigs(k), n, d, Delta(s));
        %plot(Quadratic_Low_dimension(1,:),Quadratic_Low_dimension(2,:),'*');
        result(k,s) = norm(Projection{k,s}-Data,'fro');
    end
end

%%

for k = 1:length(neigs)
    for s = 1:length(Delta)
        subplot(5,5,(k-1)*5+s);
        plot(Data(1,:),Data(2,:),'o');
        hold on
        plot(Projection{k,s}(1,:),Projection{k,s}(2,:),'.');
    end
end
figure
        
for k = 1:length(neigs)
    plot(result(k,:),'-');
    hold on
end


function W = PixValue_Weight(Data_val)
    W = diag(double(Data_val)/255);
end