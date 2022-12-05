function cancf(m_units_num,m_t)
%%
%计算常量定义
L = 1.8;            %单摆长度，m
A = 2.5*10^-4;      %单摆横截面积，m^2
Ro = 2766.67;       %密度，kg/m^3
g = 9.81;           %重力加速度，m/s^2
units_num = m_units_num;      %划分单元数量
t_end = m_t;        %计算截止时间
initial_angle = 0*pi/180;  %柔性梁初始角度，水平起始

%%
%离散
tic;
unit_l = L/units_num; %单元长度
%计算节点坐标列阵
points_e = zeros(4*(units_num+1),1);
for index = 1:units_num+1
    step = 4*(index-1);
    l = unit_l*(index-1);
    points_e(1+step:2+step) = [l*cos(initial_angle) l*sin(initial_angle)]';
    points_e(3+step:4+step) = [cos(initial_angle) sin(initial_angle)]';
end
%计算广义单元质量阵
fun = @(x) Ro*A*sfunction(x,unit_l)'*sfunction(x,unit_l);
Me = integral(fun,0,unit_l,'ArrayValued',true);
%计算广义单元重力阵
G = [0 -Ro*g]';
fun = @(x) A*sfunction(x,unit_l)';
Qge = integral(fun,0,unit_l,'ArrayValued',true)*G;
%组装
Ma = zeros(4*(units_num+1),4*(units_num+1));
Qga = zeros(4*(units_num+1),1);
for index = 1:units_num
    step = 4*(index-1);
    Ma(1+step:8+step,1+step:8+step) = Ma(1+step:8+step,1+step:8+step) + Me;
    Qga(1+step:8+step,1) = Qga(1+step:8+step,1) + Qge;
end
Ma(1:2,:) = [];
Ma(:,1:2) = [];
Qga(1:2,:) = [];
%赋初值
y0 = zeros(1,4*(units_num+1)*2);
y0(1:2:end) = points_e(:);
yy0 = y0(:,5:end);
Qfa = zeros(4*(units_num+1),1);
toc;
tic;
options = odeset('OutputFcn',@myOutputFcn);
[t,y] = ode15s(@(t,y) odefun(t,y,units_num,Qga,@myfile,Qfa,Ma),0:0.001:t_end,yy0,options);
toc;

end

%%
%求解迭代函数
function dudt = odefun(t,u,units_num,Qga,Qfe_fun,Qfa,Ma)
dudt = zeros(4*(units_num+1)*2 - 4,1);
%提取总体节点坐标列阵
Qe = u(1:2:end);
%单元弹性力更新
arg = [[0 0]';Qe(1:6,1)];
Qfa(1:8,1) = Qfe_fun(arg(1),arg(2),arg(3),arg(4),arg(5),arg(6),arg(7),arg(8));
for index = 2:units_num
    step = 4*(index-1);
    arg = Qe(4*(index-1)-1:4*(index+1)-2,1);
    Qfa(1+step:8+step,1) = Qfa(1+step:8+step,1) + Qfe_fun(arg(1),arg(2),arg(3),arg(4),arg(5),arg(6),arg(7),arg(8));
end
Qfa(1:2,:) = [];
expression = Ma\(Qga + Qfa);
dudt(1:2:end) = u(2:2:end);
dudt(2:2:end) = expression(:);
end
%%
%函数定义
function s = sfunction(x,m_l)
e = x/m_l;
s = zeros(2,8);
s1 = 1 - 3*e^2 + 2*e^3;
s2 = e - 2*e^2 + e^3;
s3 = 3*e^2 - 2*e^3;
s4 = -e^2 + e^3;
s(1,1) = s1;
s(2,2) = s1;
s(1,3) = s2*m_l;
s(2,4) = s2*m_l;
s(1,5) = s3;
s(2,6) = s3;
s(1,7) = s4*m_l;
s(2,8) = s4*m_l;
end %返回给定点的形函数矩阵
function status = myOutputFcn(t,y,flag)
if flag == "init"
elseif flag == "done"
else
    if mod(t,0.1) == 0
        disp(t);
    end
    status = 0;
end
end