clear;
clc;

% In this code, we use TVD-minmod to reconstruct U_R and U_L, and use ROS
% to discrete space, and use RK3 as time-progress solver.

global gamma p1 rho1 u1 p2 rho2 u2;
gamma = 1.4;
u1 = 0;p1 = 1;rho1 = 1;
u2 = 0;p2 = 0.1;rho2 = 0.125;

%% 初始化
x = -1:0.01:1; % 空间区域
dx = 0.01;Nx = length(x);
t = 0:0.01:0.60; % 传播时间
dt = 0.01;Nt = length(t);
y = zeros(length(t),length(x),3); % 速度1 压强2 密度3
y(:,x<0,1) = u1;y(:,x<0,2) = p1;y(:,x<0,3) = rho1;
y(:,x>=0,1) = u2;y(:,x>=0,2) = p2;y(:,x>=0,3) = rho2;
u = y(1,:,1);p = y(1,:,2);rho = y(1,:,3);E = p/(gamma-1) + 0.5*rho.*u.^2;
U = [rho;rho.*u;E];

%% 迭代计算
for i = 2:Nt  % 计算的当前时刻i
    U_old = U;
    for j = 3:Nx-2 % 计算的当前节点j 
        f_jplus = flux(U_old,j);
        f_jmiu = flux(U_old,j-1);
        temp = -(f_jplus-f_jmiu)/dx;
        U(:,j) = temp*dt + U_old(:,j);
        [rho_new,u_new,p_new] = cons2prim(U(:,j));
        y(i,j,1)=u_new;y(i,j,2)=p_new;y(i,j,3)=rho_new;
    end
end

%% 限制器函数
function s = minmod(a,b)
    if(a*b>0)
        s = sign(a)*min(abs(a),abs(b));
    else
        s = 0;
    end
end

%% 从U反推密度 速度 压强
function [rho,u,p] = cons2prim(U)
    global gamma
    rho = U(1);
    u = U(2) / rho;
    E = U(3);
    p = (gamma - 1) * (E - 0.5*rho*u^2);
end

%% 计算通量
function [f_jplus] = flux(U,j)
    global gamma
    % 计算j+1/2 U_R,U_L
    U_j = U(:,j);U_jplus1 = U(:,j+1);U_jmiu1 = U(:,j-1);U_jplus2 = U(:,j+2);
    U_R = zeros(3,1);U_L = zeros(3,1);
    for temp = 1:3
        U_L(temp) = U_j(temp) + 0.5*minmod(U_j(temp)-U_jmiu1(temp),U_jplus1(temp)-U_j(temp));
        U_R(temp) = U_jplus1(temp) - 0.5*minmod(U_jplus1(temp) - U_j(temp),U_jplus2(temp)-U_jplus1(temp));
    end
    % 计算j+1/2 aveU
    [rho_l,u_l,p_l] = cons2prim(U_L);[rho_r,u_r,p_r] = cons2prim(U_R);
    h_l = U_L(3)/rho_l + 0.5*u_l^2 + p_l/rho_l;
    h_r = U_R(3)/rho_r + 0.5*u_r^2 + p_r/rho_r;
    averho = ((sqrt(rho_l) + sqrt(rho_r))/2)^2;
    aveu = (sqrt(rho_l)*u_l + sqrt(rho_r)*u_r)/(sqrt(rho_l) + sqrt(rho_r));
    aveh = (sqrt(rho_l)*h_l + sqrt(rho_r)*h_r)/(sqrt(rho_l) + sqrt(rho_r));
    aveU = zeros(3,1);
    avep = (gamma-1)/gamma*(averho*aveh - 0.5*averho*aveu^2);
    aveU(1) = averho;aveU(2) = averho*aveu;aveU(3) = averho*(aveh -0.5*aveu^2 - avep/averho);
    % 计算Jacobi矩阵
    aveA = [0,1,0;
        -(3-gamma)/2*aveu^2,(3-gamma)*aveu,gamma-1;
        (gamma-1)*aveu^3-gamma*aveU(3)*aveu/averho,gamma*aveU(3)/averho-3*(gamma-1)*aveu^2/2,gamma*aveu];
    [S,Lambda] = eig(aveA);
    % 计算f_j+1/2
    f_u_r = zeros(3,1);f_u_l = zeros(3,1);
    f_u_l(1) = rho_l*u_l;f_u_l(2) = rho_l*u_l^2+p_l;f_u_l(3) = u_l*(U_L(3)+p_l);
    f_u_r(1) = rho_r*u_r;f_u_r(2) = rho_r*u_r^2+p_r;f_u_r(3) = u_r*(U_R(3)+p_r);
    f_jplus = 0.5*(f_u_l+f_u_r) - 0.5*S*abs(Lambda)*(S\(U_R - U_L));
end