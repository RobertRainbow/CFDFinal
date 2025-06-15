clear;
clc;

% In this code file, we use WENO to calculate U_R and U_L to discrete space
% and use Roe to calculate flux with RK3 to progress time-solver.
%% 网格划分
global gamma n_r n_l dx dt
dx = 0.001;
dt = 0.0001;
time_step = 0.40/dt;
length_l = 1;
length_r = 1;
n_l = length_l/dx;
n_r = length_r/dx;

%% 流场初始化
gamma=1.4;
p_new = [ones(1,n_l),0.1*ones(1,n_r)];
rho_new = [ones(1,n_l),0.125*ones(1,n_r)];
u_new = zeros(1,n_l+n_r);
mass_new = rho_new .* u_new;
e_new = p_new ./(gamma-1) + 0.5 * rho_new.*u_new.^2;

%% 开始计算
v = VideoWriter('approximate2.mp4','MPEG-4');
v.FrameRate = 30;
open(v);
for i=1:time_step
    p = p_new;rho=rho_new;u=u_new;mass=mass_new;e=e_new;
    for j=4:n_l+n_r-3
        % 初始化RK3所用变量
        rho1=rho;mass1=mass;e1=e;
        rho2=rho;mass2=mass;e2=e;

        % RK算法
        output = calculate_f(rho,mass,e,j);
        rho1(j) = rho(j) + 1/3*dt*output(1);
        mass1(j) = mass(j) + 1/3*dt*output(2);
        e1(j) = e(j) + 1/3*dt*output(3);

        output = calculate_f(rho1,mass1,e1,j);
        rho2(j) = rho(j) + 2/3*dt*output(1);
        mass2(j) = mass(j) + 2/3*dt*output(2);
        e2(j) = e(j) + 2/3*dt*output(3);

        output = calculate_f(rho2,mass2,e2,j);
        rho_new(j) = 0.25*(rho(j) + 3*rho1(j)) + 0.75*dt*output(1);
        mass_new(j) = 0.25*(mass(j) + 3*mass1(j)) + 0.75*dt*output(2);
        e_new(j) = 0.25*(e(j) + 3*e1(j)) + 0.75*dt*output(3);
    end
    % 更新其他数据
    u_new = mass_new ./ rho_new;
    p_new = (e_new - 0.5.*rho_new.*u_new.^2)*(gamma-1);

    % 画图
    x=linspace(-1,1,n_l+n_r);
    plot(x,rho_new,x,u_new,x,p_new,'LineWidth', 2);
    legend('rho', 'u', 'p');
    title({'Sod Problem Solution';['t=',num2str(dt*i),'s']});
    pause(0.0001);

    frame = getframe(gcf);
    writeVideo(v,frame);

end
close(v);

%% WENO重构U_R,U_L
function [output] = weno5(U)
    a = U(:,1);
    b = U(:,2);
    c = U(:,3);
    d = U(:,4);
    e = U(:,5);
    epsilon = 1e-6;

    h0 = (2*a - 7*b + 11*c) / 6;
    h1 = (-b + 5*c + 2*d) / 6;
    h2 = (2*c + 5*d - e) / 6;

    beta0 = 13/12 * (a-2*b+c) .^2 + 1/4 * (a-4*b+3*c).^2;
    beta1 = 13/12 * (b-2*c+d) .^2 + 1/4 * (b-d).^2;
    beta2 = 13/12 * (c-2*d+e).^2 + 1/4 * (3*c - 4*d + e).^2;

    w0 = 1./ (epsilon + beta0).^2;
    w1 = 6./ (epsilon + beta1).^2;
    w2 = 1./ (epsilon + beta2).^2;
    S = w0 + w1 + w2;
    w0 = w0./S;
    w1 = w1./S;
    w2 = w2./S;

    output = w0 .* h0 + w1 .* h1 + w2 .* h2;

end

%% 从U获取原始变量rho,u,p
function [rho,u,p] = getOringin(U)
    global gamma
    rho = U(1);
    u = U(2) / rho;
    p = (U(3) - 0.5*rho*u^2)*(gamma-1);
end

%% 计算通量
function f = Flux(U)
    global gamma
    rho=U(1);
    u = U(2)/rho;
    e = U(3);
    p = (gamma-1)*(e-0.5*rho*u^2);
    f = [rho*u;rho*u^2+p;u*(e+p)];
end

%% 计算重构量
function [output] = calculate_f(rho,mass,e,j)
    global gamma dx
    %% 计算f_plus_mao
    % 重构
    U = [rho;mass;e];
    U_L = weno5(U(:,j-2:j+2));
    U_R = weno5(flip(U(:,j-1:j+3),2));
    
    % 计算平均值
    [rhol,ul,pl] = getOringin(U_L);
    al = sqrt(gamma*pl/rhol);
    [rhor,ur,pr] = getOringin(U_R);
    ar = sqrt(gamma*pr/rhor);
    sl = min(ul,ur) - max(al,ar);
    sr = min(ul,ur) + max(al,ar);
    sm = (pr-pl+rhol*ul*(sl-ul)-rhor*ur*(sr-ur))/(rhol*(sl-ul)-rhor*(sr-ur));
    plr = 0.5*(pl+pr+rhol*(sl-ul)*(sm-ul)+rhor*(sr-ur)*(sm-ur));
    D = zeros(3,1);
    D(1)=0;D(2)=1;D(3)=sm;
    % 该段代码为尝试Roe处理通量，但是效果不佳故加以注释
    % hl = (U_L(3) + pl) / rhol;
    % hr = (U_R(3) + pr) / rhor;
    % 
    % averho = sqrt(rhol * rhor);
    % aveu = (sqrt(rhol)*ul + sqrt(rhor)*ur)/(sqrt(rhol) + sqrt(rhor));
    % aveh = (sqrt(rhol)*hl + sqrt(rhor)*hr)/(sqrt(rhol) + sqrt(rhor));
    % avea = sqrt((gamma-1)*aveh-0.5*aveu^2);
    % avep = avea^2*averho/gamma;
    % avee = avep / (gamma-1) + 0.5*averho*aveu^2;
    % 
    % A = [0,1,0;
    %     -(3-gamma)/2*aveu^2,(3-gamma)*aveu,gamma-1;
    %     (gamma-1)*aveu^3-gamma*avee*aveu/averho,gamma*avee/averho-3*(gamma-1)*aveu^2/2,gamma*aveu];
    % [S,Lambda] = eig(A);
    
    f_plus = Flux(U_R);
    f_miu = Flux(U_L);
    
    % alpha = max(abs(ul)+al,abs(ur)+ar);
    % f_plus_mao = 0.5*(f_plus+f_miu) - 0.5*alpha*(U_R-U_L);
    if sl>=0
        f_plus_mao = f_miu;
    elseif sr<=0
        f_plus_mao = f_plus;
    elseif sl<=0 && sm>=0
        f_plus_mao = (sm*(sl*U_L-f_miu)+sl*plr*D)/(sl-sm);
    elseif sr>=0 && sm<=0
        f_plus_mao = (sm*(sr*U_R-f_plus)+sr*plr*D)/(sr-sm);
    end
    
    %% 计算f_miu_mao
    % 重构
    U = [rho;mass;e];
    U_L = weno5(U(:,j-3:j+1));
    U_R = weno5(flip(U(:,j-2:j+2),2));
    
    % 计算平均值
    [rhol,ul,pl] = getOringin(U_L);
    al = sqrt(gamma*pl/rhol);
    [rhor,ur,pr] = getOringin(U_R);
    ar = sqrt(gamma*pr/rhor);
    sl = min(ul,ur) - max(al,ar);
    sr = min(ul,ur) + max(al,ar);
    sm = (pr-pl+rhol*ul*(sl-ul)-rhor*ur*(sr-ur))/(rhol*(sl-ul)-rhor*(sr-ur));
    plr = 0.5*(pl+pr+rhol*(sl-ul)*(sm-ul)+rhor*(sr-ur)*(sm-ur));
    D = zeros(3,1);
    D(1)=0;D(2)=1;D(3)=sm;
    % 该段代码为尝试Roe处理通量，但是效果不佳故加以注释
    % hl = (U_L(3) + pl) / rhol;
    % hr = (U_R(3) + pr) / rhor;
    % 
    % averho = sqrt(rhol * rhor);
    % aveu = (sqrt(rhol)*ul + sqrt(rhor)*ur)/(sqrt(rhol) + sqrt(rhor));
    % aveh = (sqrt(rhol)*hl + sqrt(rhor)*hr)/(sqrt(rhol) + sqrt(rhor));
    % avea = sqrt((gamma-1)*aveh-0.5*aveu^2);
    % avep = avea^2*averho/gamma;
    % avee = avep / (gamma-1) + 0.5*averho*aveu^2;
    % 
    % A = [0,1,0;
    %     -(3-gamma)/2*aveu^2,(3-gamma)*aveu,gamma-1;
    %     (gamma-1)*aveu^3-gamma*avee*aveu/averho,gamma*avee/averho-3*(gamma-1)*aveu^2/2,gamma*aveu];
    % [S,Lambda] = eig(A);
    
    f_plus = Flux(U_R);
    f_miu = Flux(U_L);
    
    % alpha = max(abs(ul)+al,abs(ur)+ar);
    % f_miu_mao = 0.5*(f_plus+f_miu) - 0.5*alpha*(U_R-U_L);
    if sl>=0
        f_miu_mao = f_miu;
    elseif sr<=0
        f_miu_mao = f_plus;
    elseif sl<=0 && sm>=0
        f_miu_mao = (sm*(sl*U_L-f_miu)+sl*plr*D)/(sl-sm);
    elseif sr>=0 && sm<=0
        f_miu_mao = (sm*(sr*U_R-f_plus)+sr*plr*D)/(sr-sm);
    end

    output = -1/dx*(f_plus_mao - f_miu_mao);
end