clear;
clc;

%% 参数设置
global gamma p1 rho1 c1 p2 rho2 c2 u1 u2;
global p_star
gamma = 1.4;
u1 = 0;u2 = 0;
p1 = 1;rho1 = 1;c1 = sqrt(gamma*p1/rho1);
p2 = 0.1;rho2 = 0.125;c2 = sqrt(gamma*p2/rho2);

%% 函数F(p_star)展示
p_star = 0:0.001:2;
y = zeros(1,length(p_star));
for i = 1:length(p_star)
    y(i) = F(p_star(i));
end
figure(1);
plot(p_star,y,'LineWidth',2);
legend("F")

%% 求解中心区域p_star
%%% 设定迭代控制条件
residual_newton_max = 0.0001;
residual_newton = 1E5;
newton_count_max = 20;
%%% 迭代公式 p_star_new=p_satr-F(p_star)/F_diff(p_star)
p_star0 = (p1+p2)/2;p_star = p_star0;% 首先给定初值
delta_p_star = 0;newton_count = 0;
while(residual_newton > residual_newton_max && newton_count < newton_count_max)
    delta_p_star = -F(p_star)/F_diff(p_star);
    p_star = p_star + delta_p_star;
    if (newton_count == 0)
        p_star0 = p_star;
    end
    newton_count=newton_count+1 ;
    residual_newton = abs(delta_p_star)/p_star0;
end

%% 求解中心区速度u_star
global u_star % 储存计算得到的全局变量p_star 
global rho_star_L rho_star_R cL cR Z1_head Z1_tail Z2_head Z2_tail; 
if(p_star >= p1)
    f1 = (p_star - p1)/rho1/c1/((gamma+1)/2/gamma *(p_star/p1) + (gamma-1)/2/gamma)^0.5;
end
if(p_star <p1)
    f1 = 2*c1/(gamma-1)*((p_star/p1)^((gamma-1)/2/gamma) - 1);
end
if(p_star >= p2)
    f2 = (p_star - p2)/rho2/c2/((gamma+1)/2/gamma *(p_star/p2) + (gamma-1)/2/gamma)^0.5;
end
if(p_star <p2)
    f2 = 2*c2/(gamma-1)*((p_star/p2)^((gamma-1)/2/gamma) - 1);
end
u_star=(u1+u2-f1+f2)/2;

%% 区域初始化（未扰动区）
dx = 0.001;dt = 0.0002;
x = linspace(-1,1,2/dx); % 空间区域
t = linspace(0,0.40,0.4/dt); % 传播时间
y = zeros(length(t),length(x),3);
y(:,x<0,1) = u1;y(:,x<0,2) = p1;y(:,x<0,3) = rho1;
y(:,x>=0,1) = u2;y(:,x>=0,2) = p2;y(:,x>=0,3) = rho2;

%% 开始传播
for i =1:length(t)
    for j =1:length(x)
      %% 左行波
      if(p1<p_star) 
          rho_star_L = rho1*(p1-p_star)/(rho1*(u_star-u1)^2+p1-p_star);
          Z_L= (rho1*u1-rho_star_L*u_star)/(rho1-rho_star_L);
          cL = sqrt(gamma*p_star/rho_star_L); 
          %中间区域-接触间断左侧
          if(x(j) <u_star*t(i) && x(j) >= Z_L*t(i))
              y(i,j,1) = u_star;y(i,j,2) = p_star;y(i,j,3) = rho_star_L;
          end
      end
      if(p1>p_star) 
          rho_star_L = (p_star*rho1^gamma/p1)^(1/gamma);
          cL = sqrt(gamma*p_star/rho_star_L);
          Z1_head=u1-c1;Z1_tail=u_star-cL;
          %中间区域-接触间断左侧
          if(x(j) <u_star*t(i) && x(j) >= Z1_tail*t(i))
              y(i,j,1) = u_star;y(i,j,2) = p_star;y(i,j,3) = rho_star_L;
          end
          % 在稀疏波内 稀疏波存在厚度
          if(x(j) < Z1_tail*t(i) && x(j) >= Z1_head*t(i))
              c_expansion = (gamma-1)/(gamma+1)*(u1-x(j)/t(i))+2/(gamma+1)*c1;
              u_expansion = c_expansion + x(j)/t(i);
              p_expansion = p1*(c_expansion/c1)^(2*gamma/(gamma-1));
              rho_expansion = gamma*p_expansion/c_expansion^2;
              y(i,j,1) = u_expansion;y(i,j,2) = p_expansion;y(i,j,3) = rho_expansion;
          end
      end
      %% 右行波
      if(p2<p_star) 
          rho_star_R = rho2*(p2-p_star)/(rho2*(u_star-u2)^2+p2-p_star);
          Z_R= (rho2*u2-rho_star_R*u_star)/(rho2-rho_star_R);
          cR = sqrt(gamma*p_star/rho_star_R); 
          %中间区域-接触间断右侧
          if(x(j) >=u_star*t(i) && x(j) <= Z_R*t(i))
              y(i,j,1) = u_star;y(i,j,2) = p_star;y(i,j,3) = rho_star_R;
          end
      end
      if(p2>p_star) 
          rho_star_R = (p_star*rho2^gamma/p2)^(1/gamma);
          cR = sqrt(gamma*p_star/rho_star_R); 
          Z2_head=u2+c2;Z2_tail=u_star+cR; 
          %中间区域-接触间断右侧
          if(x(j) >u_star*t(i) && x(j) <= Z2_tail*t(i))
              y(i,j,1) = u_star;y(i,j,2) = p_star;y(i,j,3) = rho_star_L;
          end
          % 在稀疏波内 稀疏波存在厚度
          if(x(j) > Z2_tail*t(i) && x(j) <= Z2_head*t(i))
              c_expansion =(u2+x(j)/t(i))*(gamma-1)/(gamma+1)+2*c2/(gamma+1);
              u_expansion = x(j)/t(i)-c_expansion  ;
              p_expansion = p2*((c_expansion/c2)^(2*gamma/(gamma-1)));
              rho_expansion = gamma*p_expansion/(c_expansion^2);
              y(i,j,1) = u_expansion;y(i,j,2) = p_expansion;y(i,j,3) = rho_expansion;
          end
      end
    end
end

%% 结果可视化
x = linspace(-1,1,2/0.001);
v = VideoWriter('RiemannExactSolve.mp4','MPEG-4');
v.FrameRate = 30;
open(v);

% 速度为1 压强为2 密度为3
u = y(:,:,1);
p = y(:,:,2);
rho = y(:,:,3);

for i=1:length(t)
    temp_u = u(i,:);
    temp_p = p(i,:);
    temp_rho = rho(i,:);
    plot(x,temp_rho,x,temp_u,x,temp_p,'LineWidth', 2);
    legend('rho', 'u', 'p');
    title({'Sod Problem Solution';['t=',num2str(dt*i),'s']});
    pause(0.0001);

    frame = getframe(gcf);
    writeVideo(v,frame);

end
close(v);


% %% 画曲线图
% x = -1:0.001:1;
% 
% u_final   = squeeze(y(end, :, 1));  % 速度
% p_final   = squeeze(y(end, :, 2));  % 压强
% rho_final = squeeze(y(end, :, 3));  % 密度
% 
% % 绘图
% figure;
% subplot(3,1,1);
% plot(x, u_final, 'b', 'LineWidth', 1.5);
% ylabel('u');
% title('t = 最后时间步');
% 
% subplot(3,1,2);
% plot(x, p_final, 'r', 'LineWidth', 1.5);
% ylabel('p');
% 
% subplot(3,1,3);
% plot(x, rho_final, 'k', 'LineWidth', 1.5);
% ylabel('\rho');
% xlabel('x');
% 
% sgtitle('最终时刻 Sod 激波管 Riemann 精确解');

%% F_diff(p_star) 导数
function[y] = F_diff(p_star)
global gamma p1 rho1 c1 p2 rho2 c2;
if(p_star >= p1)
    f1_diff = (( (gamma+1)/2 *(p_star/p1) + (gamma-1)/2/gamma )^0.5 - ...
        (p_star-p1)/2*( (gamma+1)/2/gamma *(p_star/p1) + (gamma-1)/2/gamma)^(-0.5) * (gamma+1)/2/gamma/p1)...
        /(rho1*c1*((gamma+1)/2/gamma *p_star/p1  + (gamma-1)/2/gamma));
end
if(p_star <p1)
    f1_diff = c1/gamma*(p_star/p1)^((gamma-1)/2/gamma)*(p_star)^(-1);
end
if(p_star>=p2)
    f2_diff = (( (gamma+1)/2 *(p_star/p2) + (gamma-1)/2/gamma )^0.5 - ...
        (p_star-p2)/2*( (gamma+1)/2/gamma *(p_star/p2) + (gamma-1)/2/gamma)^(-0.5) * (gamma+1)/2/gamma/p2)...
        /(rho2*c2*((gamma+1)/2/gamma *p_star/p2  + (gamma-1)/2/gamma));
end
if(p_star <p2)
    f2_diff = c2/gamma*(p_star/p2)^((gamma-1)/2/gamma)*(p_star)^(-1);
end
y = f1_diff+f2_diff;
end

%% F(p_star)
function[y] = F(p_star)
global gamma p1 rho1 c1 p2 rho2 c2 u1 u2;
if(p_star >= p1)
    f1 = (p_star - p1)/rho1/c1/( (gamma+1)/2/gamma *(p_star/p1) + (gamma-1)/2/gamma)^0.5;
end
if(p_star < p1)
    f1 = 2*c1/(gamma-1)*((p_star/p1)^((gamma-1)/2/gamma) - 1);
end
if(p_star >= p2)
    f2 = (p_star - p2)/rho2/c2/( (gamma+1)/2/gamma *(p_star/p2) + (gamma-1)/2/gamma)^0.5;
end
if(p_star < p2)
    f2 = 2*c2/(gamma-1)*((p_star/p2)^((gamma-1)/2/gamma) - 1);
end
y = f1+f2-u1+u2 ;
end

