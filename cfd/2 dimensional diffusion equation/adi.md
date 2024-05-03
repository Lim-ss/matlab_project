二维扩散方程表示为

$\frac{\partial u}{\partial t}=\sigma_x\frac{\partial^2u}{\partial^2x}+\sigma_y\frac{\partial^2u}{\partial^2y}$

本题中

$\sigma_x=\sigma_y=\sigma$

因此方程简化成

$\frac{\partial u}{\partial t}=\sigma(\frac{\partial^2u}{\partial^2x}+\frac{\partial^2u}{\partial^2y})$

交替方向隐式方法（Alternating Diretion Implicit, ADI）可以表示为以下两个式子

1式：$\frac{u^*_{ij}-u^n_{ij}}{\Delta t/2} = \sigma_x\frac{\delta_x^2u_{ij}^*}{\Delta^2x}+\sigma_y\frac{\delta_y^2u_{ij}^n}{\Delta^2y}$，$u^*_{ij}$为未知量

2式：$\frac{u^{n+1}_{ij}-u^*_{ij}}{\Delta t/2} = \sigma_x\frac{\delta_x^2u_{ij}^*}{\Delta^2x}+\sigma_y\frac{\delta_y^2u_{ij}^{n+1}}{\Delta^2y}$，$u^{n+1}_{ij}$为未知量

其中$\delta_x^2=E^1_x-2I+E^{-1}_x$，$\delta_y^2=E^1_y-2I+E^{-1}_y$

在本题中$\Delta x=\Delta y,\sigma_x=\sigma_y=\sigma$，令$\frac{\sigma\Delta t}{2\Delta^2x}=\frac{\sigma\Delta t}{2\Delta^2y}=\lambda$，两式子分别简化为

3式：$(1+2\lambda)u^*_{ij}-\lambda u^*_{i+1,j}-\lambda u^*_{i-1,j}=(1-2\lambda)u^n_{ij}+\lambda u^n_{i,j+1}+\lambda u^n_{i,j-1}$

4式：$(1+2\lambda)u^{n+1}_{ij}-\lambda u^{n+1}_{i,j+1}-\lambda u^{n+1}_{i,j-1}=(1-2\lambda)u^*_{ij}+\lambda u^*_{i+1,j}+\lambda u^*_{i-1,j}$

对于边界条件，以$u^*_{ij}$为例，如果简单地令$u^*_{ij}=u^{n+\frac12}_{ij}$，则只能得到一阶精度，用上方1式减去2式，得到表达式5式

5式：

$u^*_{ij}=\frac{u^n_{ij}+u^{n+1}_{ij}}{2}-\frac12\lambda\delta^2_y(u^n_{ij}-u^{n+1}_{ij})$

$=\frac12[(1-2\lambda)u^n_{ij}+(1+2\lambda)u^{n+1}_{ij}+\lambda(u^n_{i,j-1}+u^n_{i,j+1}-u^{n+1}_{i,j-1}+u^{n+1}_{i,j+1})]$

利用该式子作为边界条件则能够获得二阶精度

因此，