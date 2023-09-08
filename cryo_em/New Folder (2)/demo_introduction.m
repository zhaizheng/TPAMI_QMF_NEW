show_or_not = 0;
[clean, noisy] = cryo_em_data(show_or_not);
vec_clean = matrix_to_vec(clean);
[U,L,~] = svd(vec_clean);
Data = U(:,1:20)*U(:,1:20)'* vec_clean;
image_Data = vec_to_matrix(Data);
%PARAMETERS SETTING
D = 20; d = 3;  tol=1.e-4; kk = [10,15,20,25,30]; NUM = 30; sigma = 0;

low_dimension_noisy_data = U(:,1:D)'* vec_clean+sigma*randn(D,size(vec_clean,2));
Data_noisy = U(:,1:D)*low_dimension_noisy_data;
image_Data2 = vec_to_matrix(Data_noisy);
[mappedX, mapping] = lle(low_dimension_noisy_data', 2, 8);

%%
plot(mappedX(:,1),mappedX(:,2),'.');


Location = mappedX';
Xmin = min(Location(1,:));
Xmax = max(Location(1,:));
Ymin = min(Location(2,:));
Ymax = max(Location(2,:));
X_intervel = 8; Y_intervel = 8;

Interest_Point = [];
Xstep = (Xmax-Xmin)/X_intervel;
for i = 1:X_intervel
    Temp = [];
    for k = 1:size(Location,2)
        if Location(1,k)>= Xmin+(i-1)*Xstep && Location(1,k)<Xmin+i*Xstep
            Temp = [Temp, [Location(:,k);k]];
        end
    end
    lmax = find(Temp(2,:)==max(Temp(2,:)));
    lmin = find(Temp(2,:)==min(Temp(2,:)));
    Interest_Point = [Interest_Point, Temp(:,lmax), Temp(:,lmin)];
end

Ystep = (Ymax-Ymin)/Y_intervel;
for i = 1:Y_intervel
    Temp = [];
    for k = 1:size(Location,2)
        if Location(2,k)>= Ymin+(i-1)*Ystep && Location(2,k)<Ymin+i*Ystep
            Temp = [Temp, [Location(:,k);k]];
        end
    end
    lmax = find(Temp(1,:)==max(Temp(1,:)));
    lmin = find(Temp(1,:)==min(Temp(1,:)));
    Interest_Point = [Interest_Point, Temp(:,lmax), Temp(:,lmin)];
end

for k = 1:size(Interest_Point,2)
    hold on
    image(squeeze(image_Data2(Interest_Point(3,k),:,:))*200,'XData', ...
    [Interest_Point(1,k), Interest_Point(1,k)+0.01], 'YData',[Interest_Point(2,k), Interest_Point(2,k)+0.01]);
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