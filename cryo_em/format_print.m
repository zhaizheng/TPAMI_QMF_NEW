function format_print(result)
    [m,n] = size(result);
    for i = 1:m
        for j = 1:n
            if mod(j,2)== 1
                fprintf('%s %.4f','&',result(i,j));
            else
                fprintf('(%.4f) ',result(i,j));
            end
        end
        fprintf('%s','\\');
        fprintf('\n')
    end
end