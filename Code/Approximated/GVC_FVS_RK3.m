clear;
clc;

% In this code file, we use NND combining with FVS(Steger_Warming method)
% to discrete space and solve with RK3 to progress.
%% 网格划分
global gamma n_r n_l dx dt
dx = 0.001;
dt = 0.0002;
time_step = 0.40/dt;
length_l = 1;
length_r = 1;
n_l = length_l/dx;
n_r = length_r/dx;

%% 流场初始化
gamma=1.4;
p = [ones(1,n_l),0.1*ones(1,n_r)];
rho = [ones(1,n_l),0.125*ones(1,n_r)];
u = zeros(1,n_l+n_r);
mass = rho .* u;
e = p ./(gamma-1) + 0.5 * rho.*u.^2;

%% 计算守恒量
v = VideoWriter('approximate1.mp4','MPEG-4');
v.FrameRate = 30;
open(v);
for i=1:time_step
    [f] = solver(p,rho,u,mass,e);
    for m=3:n_l+n_r-2
        rho1(m) = rho(m) + 1/3*dt*f(1,m);
        mass1(m) = mass(m) + 1/3*dt*f(2,m);
        u1(m) = mass1(m)/rho1(m);
        e1(m) = e(m) + 1/3*dt*f(3,m);
        p1(m) = (e1(m) - 0.5*rho1(m)*u1(m)^2)*(gamma-1);
    end
    % 数据对齐
    for m=n_l+n_r-1:n_l+n_r
        rho1(m) = rho(m);
        mass1(m) = mass(m);
        u1(m) = u(m);
        e1(m) = e(m);
        p1(m) = p(m);
    end
    for m=1:2
        rho1(m) = rho(m);
        mass1(m) = mass(m);
        u1(m) = u(m);
        e1(m) = e(m);
        p1(m) = p(m);
    end
    
    [f] = solver(p1,rho1,u1,mass1,e1);
    for m=3:n_l+n_r-2
        rho2(m) = rho(m) + 2/3*dt*f(1,m);
        mass2(m) = mass(m) + 2/3*dt*f(2,m);
        u2(m) = mass2(m)/rho2(m);
        e2(m) = e(m) + 2/3*dt*f(3,m);
        p2(m) = (e2(m) - 0.5*rho2(m)*u2(m)^2)*(gamma-1);
    end
    % 数据对齐
    for m=n_l+n_r-1:n_l+n_r
        rho2(m) = rho(m);
        mass2(m) = mass(m);
        u2(m) = u(m);
        e2(m) = e(m);
        p2(m) = p(m);
    end
    for m=1:2
        rho2(m) = rho(m);
        mass2(m) = mass(m);
        u2(m) = u(m);
        e2(m) = e(m);
        p2(m) = p(m);
    end

    [f] = solver(p2,rho2,u2,mass2,e2);
    for m=3:n_l+n_r-2
        rho(m) = 0.25*(rho(m)+3*rho1(m)) + 0.75*dt*f(1,m);
        mass(m) = 0.25*(mass(m)+3*mass1(m)) + 0.75*dt*f(2,m);
        u(m) = mass(m)/rho(m);
        e(m) = 0.25*(e(m)+3*e1(m)) + 0.75*dt*f(3,m);
        p(m) = (e(m) - 0.5*rho(m)*u(m)^2)*(gamma-1);
    end

    % 每个时间步作图
    x=linspace(-1,1,n_l+n_r);
    % axis([0 1 0 1]);
    plot(x,rho,x,u,x,p,'LineWidth', 2);
    legend('rho', 'u', 'p');
    title({'Sod Problem Solution';['t=',num2str(dt*i),'s']});
    pause(0.0001);

    frame = getframe(gcf);
    writeVideo(v,frame);

end
close(v);

function [f] = solver(p,rho,u,mass,e)
    global gamma n_l n_r dx dt
    % 计算声速和特征值
    c = sqrt(gamma*p./rho);
    lambda1 = u;
    lambda2 = u+c;
    lambda3 = u-c;
    
    % 根据FVS计算守恒量通量项
    F_plus(1,:)=rho./(2*gamma).*(2*(gamma-1)*max(lambda1,0)+max(lambda2,0)+max(lambda3,0));
    F_plus(2,:)=rho./(2*gamma).*(2*(gamma-1)*max(lambda1,0).*u+max(lambda2,0).*(u+c)+max(lambda3,0).*(u-c));
    F_plus(3,:)=rho./(2*gamma).*((gamma-1)*max(lambda1,0).*u.^2+(3-gamma)/(2*(gamma-1))*(max(lambda2,0)+max(lambda3,0)).*c.^2+0.5*max(lambda2,0).*(u+c).^2+0.5*max(lambda3,0).*(u-c).^2);
    
    F_minus(1,:)=rho./(2*gamma).*(2*(gamma-1)*min(lambda1,0)+min(lambda2,0)+min(lambda3,0));
    F_minus(2,:)=rho./(2*gamma).*(2*(gamma-1)*min(lambda1,0).*u+min(lambda2,0).*(u+c)+min(lambda3,0).*(u-c));
    F_minus(3,:)=rho./(2*gamma).*((gamma-1)*min(lambda1,0).*u.^2+(3-gamma)/(2*(gamma-1))*(min(lambda2,0)+min(lambda3,0)).*c.^2+0.5*min(lambda2,0).*(u+c).^2+0.5*min(lambda3,0).*(u-c).^2);

    % 采用NND格式（群速度控制）计算通量
    for m=3:n_l+n_r-2 % 空间位置x
        for j=1:3 % 分量
            F_plus_minus_1(j,m)=F_plus(j,m)-F_plus(j,m-1);
            F_minus_minus_1(j,m)=F_minus(j,m)-F_minus(j,m-1);
            
            F_plus_plus_1(j,m)=F_plus(j,m+1)-F_plus(j,m);
            F_minus_plus_1(j,m)=F_minus(j,m+1)-F_minus(j,m);
            
            F_plus_minus_3(j,m)=F_plus(j,m-1)-F_plus(j,m-2);
            F_minus_minus_3(j,m)=F_minus(j,m-1)-F_minus(j,m-2);
            
            F_plus_plus_3(j,m)=F_plus(j,m+2)-F_plus(j,m+1);
            F_minus_plus_3(j,m)=F_minus(j,m+2)-F_minus(j,m+1);
        end
    end
    % 计算minmod系数以简化计算
    x1=(sign(F_plus_minus_1)==sign(F_plus_plus_1));
    x2=(sign(F_minus_plus_1)==sign(F_minus_plus_3));
    x3=(sign(F_plus_minus_3)==sign(F_plus_minus_1));
    x4=(sign(F_minus_minus_1)==sign(F_minus_plus_1));
    
    % 迭代计算守恒量
    for m=3:n_l+n_r-2

        f(1,m) = -1/dx*(F_plus(1,m)+0.5*x1(1,m)*sign(F_plus_minus_1(1,m))*min(abs(F_plus_minus_1(1,m)),abs(F_plus_plus_1(1,m)))...
            +F_minus(1,m+1)-0.5*x2(1,m)*sign(F_minus_plus_1(1,m))*min(abs(F_minus_plus_1(1,m)),abs(F_minus_plus_3(1,m)))...
            -F_plus(1,m-1)-0.5*x3(1,m)*sign(F_plus_minus_3(1,m))*min(abs(F_plus_minus_3(1,m)),abs(F_plus_minus_1(1,m)))...
            -F_minus(1,m)+0.5*x4(1,m)*sign(F_minus_minus_1(1,m))*min(abs(F_minus_minus_1(1,m)),abs(F_minus_plus_1(1,m))));
        
        f(2,m) = -1/dx*(F_plus(2,m)+0.5*x1(2,m)*sign(F_plus_minus_1(2,m))*min(abs(F_plus_minus_1(2,m)),abs(F_plus_plus_1(2,m)))...
            +F_minus(2,m+1)-0.5*x2(2,m)*sign(F_minus_plus_1(2,m))*min(abs(F_minus_plus_1(2,m)),abs(F_minus_plus_3(2,m)))...
            -F_plus(2,m-1)-0.5*x3(2,m)*sign(F_plus_minus_3(2,m))*min(abs(F_plus_minus_3(2,m)),abs(F_plus_minus_1(2,m)))...
            -F_minus(2,m)+0.5*x4(2,m)*sign(F_minus_minus_1(2,m))*min(abs(F_minus_minus_1(2,m)),abs(F_minus_plus_1(2,m))));
        
        f(3,m) = -1/dx*(F_plus(3,m)+0.5*x1(3,m)*sign(F_plus_minus_1(3,m))*min(abs(F_plus_minus_1(3,m)),abs(F_plus_plus_1(3,m)))...
            +F_minus(3,m+1)-0.5*x2(3,m)*sign(F_minus_plus_1(3,m))*min(abs(F_minus_plus_1(3,m)),abs(F_minus_plus_3(3,m)))...
            -F_plus(3,m-1)-0.5*x3(3,m)*sign(F_plus_minus_3(3,m))*min(abs(F_plus_minus_3(3,m)),abs(F_plus_minus_1(3,m)))...
            -F_minus(3,m)+0.5*x4(3,m)*sign(F_minus_minus_1(3,m))*min(abs(F_minus_minus_1(3,m)),abs(F_minus_plus_1(3,m))));
        
    end
end