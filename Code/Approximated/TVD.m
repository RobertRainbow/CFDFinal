clear;
clc;

% In this code file, we use TVD to solve Sod problem.
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
p_new = [ones(1,n_l),0.1*ones(1,n_r)];
rho_new = [ones(1,n_l),0.125*ones(1,n_r)];
u_new = zeros(1,n_l+n_r);
mass_new = rho_new .* u_new;
e_new = p_new ./(gamma-1) + 0.5 * rho_new.*u_new.^2;

%% 开始计算
v = VideoWriter('approximate3.mp4','MPEG-4');
v.FrameRate=30;
open(v);
for i=1:time_step
    p=p_new;rho=rho_new;u=u_new;mass=mass_new;e=e_new;
    U = [rho;mass;e];
    for j = 3:n_l+n_r-1
        f_low_plus1 = 0.5 * (flux(U(:,j)) + flux(U(:,j+1))) - dx/2/dt*(U(:,j+1) - U(:,j));
        f_high_plus1 = 0.5 * (flux(U(:,j)) + flux(U(:,j+1)));
        pre_down_r = U(:,j+1) - U(:,j);
        pre_up_r = U(:,j) - U(:,j-1);
        r = cal_r(pre_up_r,pre_down_r);
        phi = constrain(r);
        f_plus_1 = f_low_plus1 - phi.*(f_low_plus1 - f_high_plus1);

        f_low_plus1 = 0.5 * (flux(U(:,j-1)) + flux(U(:,j))) - dx/2/dt*(U(:,j) - U(:,j-1));
        f_high_plus1 = 0.5 * (flux(U(:,j-1)) + flux(U(:,j)));
        pre_down_r = U(:,j) - U(:,j-1);
        pre_up_r = U(:,j-1) - U(:,j-2);
        r = cal_r(pre_up_r,pre_down_r);
        phi = constrain(r);
        f_miu_1 = f_low_plus1 - phi.*(f_low_plus1 - f_high_plus1);

        output = -(f_plus_1-f_miu_1)./dx;

        rho_new(j) = output(1) * dt + rho(j);
        mass_new(j) = output(2) * dt + mass(j);
        e_new(j) = output(3) * dt + e(j);

    end
    u_new = mass_new ./ rho_new;
    p_new = (e_new - 0.5.*rho_new.*u_new.^2)*(gamma-1);

    %% 画图
    x=linspace(-1,1,n_l+n_r);
    plot(x,rho_new,x,u_new,x,p_new,'LineWidth', 2);
    legend('rho', 'u', 'p');
    title({'Sod Problem Solution';['t=',num2str(dt*i),'s']});
    pause(0.0001);

    frame = getframe(gcf);
    writeVideo(v,frame);

end
close(v);

%% 计算通量
function f = flux(U)
    global gamma
    rho=U(1);
    u = U(2)/rho;
    e = U(3);
    p = (gamma-1)*(e-0.5*rho*u^2);
    f = [rho*u;rho*u^2+p;u*(e+p)];
end

%% 限制器函数
function phi = constrain(r)
    if(r(1)<=0)
        phi(1) = 0;
    else
        phi(1) = 2*r(1)/(1+r(1));
    end

    if(r(2)<=0)
        phi(2) = 0;
    else
        phi(2) = 2*r(2)/(1+r(2));
    end

    if(r(3)<=0)
        phi(3) = 0;
    else
        phi(3) = 2*r(3)/(1+r(3));
    end
end

%% 计算r
function r = cal_r(pre_up,pre_down)
    if(pre_down(1)==0)
        r(1)=0;
    else
        r(1) = pre_up(1)/pre_down(1);
    end

    if(pre_down(2)==0)
        r(2)=0;
    else
        r(2) = pre_up(2)/pre_down(2);
    end
    if(pre_down(3)==0)
        r(3)=0;
    else
        r(3) = pre_up(3)/pre_down(3);
    end
end