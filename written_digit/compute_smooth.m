n_alg = 8;
re_z = zeros(n_alg+1,n_image);
for j = 1:n_image
    re_z(1,j) = Calcu(double(image_v{j}));
end


for i = 1:n_alg
    for j = 1:n_image
        if i == 7
            temp = Projection{i,j}.result;
        else
            temp = Projection{i,j};
        end  
        temp = temp(:,any(~isnan(temp),1));
        re_z(i+1,j) = Calcu(New_Image(temp,28,28));
    end
end

function abs_re = Calcu(img_input)
    scrImg = img_input;
    rows = size(img_input, 1);
    cols = size(img_input, 2);
    mask = [0 1 0;...
        1 -4 1;...
        0 1 0];
      % mask = [1/4 1 1/4;...
      %   1 -5 1;...
      %   1/4 1 1/4];
    lapImg = zeros(rows,cols);
    for i = 1:rows-2
        for j = 1:cols-2
            temp = mask.*scrImg(i:i+2, j:j+2);
            lapImg(i+1,j+1) = sum(temp(:));
        end
    end
    Abs_Img = abs(lapImg);
    abs_re = sum(Abs_Img(Abs_Img>0));
end


function B =  New_Image(C, m, n)
    if max(round(C(1,:))) > m
        m = max(round(C(1,:)));
    end
    if max(round(C(2,:))) > n
        n = max(round(C(2,:)));
    end
    A =  sparse(round(C(1,:)),round(C(2,:)),255*ones(1, size(C,2)), m, n);
    B = full(A);
end
%lapImg_unscaling = uint8(abs(lapImg));
%figure
%imshow(lapImg_unscaling)