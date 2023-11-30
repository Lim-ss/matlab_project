%{
    for 循环变量 = 初值 : 步长 : 终值（省略步长则默认1）
        循环语句1
        ……
        循环语句n
    end
%}  
    sum = 0;
    for n = 1:5
        sum = sum + n^2;
    end

%{
    while 条件表达式
        循环语句1
        ……
        循环语句n
    end
%}

%{
    if 条件表达式
        ……
        语句体
        ……
    end
%}

%{
    if 条件表达式
        ……
    else
        ……
    end
%}

%{
    switch 表达式(数值或字符串)
        case 数值或字符串1
            ……；
        case 数值或字符串2
            ……；
        case 数值或字符串3
            ……；
        ……
        otherwise
            ……；
    end
%}