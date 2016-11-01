function digit = double_digit(num)

    if(num < 10)
        digit = ['0', num2str(num)];
    else
        digit = num2str(num);
    end

end
