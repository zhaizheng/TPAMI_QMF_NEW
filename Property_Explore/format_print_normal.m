function format_print_normal(result,digit)
    [m,n] = size(result);
    for i = 1:m
        for j = 1:n
            fprintf(['%s %.',num2str(digit),'f'],'&',result(i,j));
        end
        fprintf('%s','\\');
        fprintf('\n')
    end
end