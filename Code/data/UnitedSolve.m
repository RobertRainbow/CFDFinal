clear;clc;

%% 读入精确解
u_exact = load('u.mat');
u_exact = u_exact.temp_u;
p_exact = load('p.mat');
p_exact = p_exact.temp_p;
rho_exact = load('rho.mat');
rho_exact = rho_exact.temp_rho;

%% 读入NND格式近似解
u_nnd = load('u_nnd.mat');
u_nnd = u_nnd.u;
p_nnd = load('p_nnd.mat');
p_nnd = p_nnd.p;
rho_nnd = load('rho_nnd.mat');
rho_nnd = rho_nnd.rho;

%% 读入WENO格式近似解
u_weno = load('u_weno.mat');
u_weno = u_weno.u;
p_weno = load('p_weno.mat');
p_weno = p_weno.p;
rho_weno = load('rho_weno.mat');
rho_weno = rho_weno.rho;

%% 读入TVD格式近似解
u_tvd = load('u_tvd.mat');
u_tvd = u_tvd.u;
p_tvd = load('p_tvd.mat');
p_tvd = p_tvd.p;
rho_tvd = load('rho_tvd.mat');
rho_tvd = rho_tvd.rho;

%% 合并到一张图
x = linspace(-1,1,length(u_exact));

figure;
hold on
plot(x,u_exact,x,u_nnd,x,u_weno,x,u_tvd,'LineWidth',2);
legend('exact','nnd','weno','tvd');
xlabel('x'); ylabel('u');
title('u-solve with methods');

figure;
hold on
plot(x,p_exact,x,p_nnd,x,p_weno,x,p_tvd,'LineWidth',2);
legend('exact','nnd','weno','tvd');
xlabel('x'); ylabel('p');
title('p-solve with methods');

figure;
hold on
plot(x,rho_exact,x,rho_nnd,x,rho_weno,x,rho_tvd,'LineWidth',2);
legend('exact','nnd','weno','tvd');
xlabel('x'); ylabel('rho');
title('rho-solve with methods');