%林文浩PB21050974
clear;
clc;
freedown = input("请输入系统的自由度\n");
m = zeros(freedown);%预创建质量矩阵
k = zeros(freedown);%预创建刚度矩阵
P = sym(zeros(freedown,1));%激励向量
x0 = sym(zeros(freedown,1));%c初始位移
v0 = sym(zeros(freedown,1));%c初始速度
disp("请按先行后列的顺序输入质量矩阵m");
disp("注意每输一个值都要按一下回车");
for i = 1:freedown
    for j = 1:freedown
        fprintf("第%d行第%d列:",i,j);
        m(i,j) = input("");
    end
end
disp("请按先行后列的顺序输入刚度矩阵k");
disp("注意每输一个值都要按一下回车");
for i = 1:freedown
    for j = 1:freedown
        fprintf("第%d行第%d列:",i,j);
        k(i,j) = input("");
    end
end
disp("请输入激励列向量P");
disp("注意表达式以t为自变量,且符合matlab语法,比如2t是错的,2*t是对的");
disp("注意每输入一个激励表达式都要按一下回车");
for i = 1:freedown
    fprintf("第%d列:",i);
    P_str = input("","s");
    P(i,1) = str2sym(P_str);
end
disp("请输入初始位移列向量x(0)");
disp("注意每输入一个初始位移都要按一下回车");
for i = 1:freedown
    fprintf("第%d列:",i);
    x0_str = input("","s");
    x0(i,1) = str2sym(x0_str);
end
disp("请输入初始速度列向量x(0)");
disp("注意每输入一个初始速度都要按一下回车");
for i = 1:freedown
    fprintf("第%d列:",i);
    v0_str = input("","s");
    v0(i,1) = str2sym(v0_str);
end
%m = [1,0,0;0,1,0;0,0,1];%test
%k = [3,-1,0;-1,2,-1;0,-1,3];%test
%解特征值问题得到freedown个固有频率
A = sym(zeros(freedown));% 预创建符号矩阵，用于计算特征值问题
syms Omiga_square;
for i = 1:freedown
    for j = 1:freedown
        A(i,j) = k(i,j) - Omiga_square*m(i,j);
    end
end
determinant = det(A);%求解行列式的结果
Omiga_n_square = solve(determinant==0,Omiga_square);%令行列式为0，得到freedown个特征值，也就是固有频率
%由已经得到的freedown个特征值，求出freedown个特征向量
phi = zeros(freedown);% 预创建符号模态矩阵
for i = 1:freedown
    %B = double(A);%这里用不了，除非符号矩阵中的元素全部是常量
    B = double(subs(A,'Omiga_square',double(Omiga_n_square(i,1))));
    phi(1:freedown,i) = null(B,'r');
end
M = phi'*m*phi;%模态质量矩阵
K = phi'*k*phi;%模态刚度矩阵
% 利用方程组Mq"+Kq = Q = phi'*P已经解耦的性质，分别求广义坐标q1~qn
phiT = phi';
phi_inv = inv(phi);
for i = 1:freedown
    Qi = sym(0);
    qi0 = sym(0);
    dqi0 = sym(0);
    Mi = sym(M(i,i));
    Ki = sym(K(i,i));
    for j = 1:freedown
        Qi = Qi + phiT(i,j)*P(j,1);
    end
    for j = 1:freedown
        qi0 = qi0 + phi_inv(i,j)*x0(j,1);
    end
    for j = 1:freedown
        dqi0 = dqi0 + phi_inv(i,j)*v0(j,1);
    end
    syms qi(t);
    dqi = diff(qi,1);
    q(i,1) = dsolve(Mi*diff(qi,2)+Ki*qi == Qi,dqi(0)==dqi0,qi(0)==qi0);
    %simplify_result(i,1) = simplify(result(i,1));
end
%将广义坐标转为原本的坐标
for i = 1:freedown
    xi = sym(0);
    for j = 1:freedown
        xi = xi + phi(i,j)*q(j,1);
    end
    fprintf("x%d的解为\n",i);
    disp(xi);
end