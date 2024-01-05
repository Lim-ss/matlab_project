%林文浩PB21050974
clear;
clc;
m_str = input("请输入(等效)质量 m\n","s");
c_str = input("请输入(等效)阻尼系数 c\n","s");
k_str = input("请输入(等效)刚度 k\n","s");
P_str = input("请输入激励 P(t)\n","s");
x0_str = input("请输入初始位移 x(0)的值\n","s");
v0_str = input("请输入初始速度 x'(0)的值\n","s");

m = str2sym(m_str);
c = str2sym(c_str);
k = str2sym(k_str);
P = str2sym(P_str);
x0 = str2sym(x0_str);
v0 = str2sym(v0_str);

syms x(t);
dx = diff(x,1);
result = dsolve(m*diff(x,2)+c*diff(x,1)+k*x == P,dx(0)==v0,x(0)==x0);
simplify_result = simplify(result);
disp("振动方程的解为:");
disp(result);
disp("解的简化结果为:");
disp(result);

