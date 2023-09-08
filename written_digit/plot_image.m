figure
n_alg = 8;
t = tiledlayout(n_alg+1,n_image,'TileSpacing','none');

for i = 1:n_image
    nexttile
    image(image_v{i});
    axis off
end

for i = 1:n_alg
    for j = 1:n_image
        nexttile
        if i == 7
            temp = Projection{i,j}.result;
        else
            temp = Projection{i,j};
        end
        temp = temp(:,any(~isnan(temp),1));
        re = New_Image(temp,28,28);
        image(re)
        axis off
    end
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