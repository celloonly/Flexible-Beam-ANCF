clc;
clear;
%%
%计算常量定义
L = 1.8;            %单摆长度，m
b = 3.125*10^-2;    %单摆宽度，m
h = 0.8*10^-2;      %单摆高，m
A = 2.5*10^-4;      %单摆横截面积，m^2
Iz = 0.133*10^-8;   %截面惯性矩，m^4
Ro = 2766.67;       %密度，kg/m^3
E = 6.8952*10^7;    %弹性模量，Pa
g = 9.81;           %重力加速度，m/s^2
units_num = 8;      %划分单元数量
t_end = 3;        %计算截止时间
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
%计算常质量阵
bin = zeros(4*(units_num+1)*2-4,1);
bin(1:2:end) = 1;
M = diag(bin,0);
M(2:2:end,2:2:end) = Ma;
%广义弹性力表达式计算
syms x;
syms qe [8 1];
sp_square = qe.' * sfunction_p(x,unit_l).' * sfunction_p(x,unit_l) * qe;
f_rootMeanSquare = sqrt(int(sp_square,x,[0 unit_l]) / unit_l);
GreenLagrangian = (qe.'*sfunction_p(x,unit_l).'*sfunction_p(x,unit_l)*qe - 1)/2;
fun1 = E*A*GreenLagrangian * sfunction_p(x,unit_l).' * sfunction_p(x,unit_l)*qe;
fun2 = -E*Iz*sfunction_pp(x,unit_l).' * sfunction_pp(x,unit_l) * qe / f_rootMeanSquare^4;
% fun3 = 2*E*Iz*qe.'*sfunction_pp(x,unit_l).'*sfunction_pp(x,unit_l)*qe*sfunction_p(x,unit_l).'*sfunction_p(x,unit_l)*qe/f_rootMeanSquare^6;
% fun4 = E*Iz*qe.'*sfunction_p(x,unit_l).'*sfunction_pp(x,unit_l)*qe*(sfunction_p(x,unit_l).'*sfunction_pp(x,unit_l)+sfunction_pp(x,unit_l).'*sfunction_p(x,unit_l))*qe/f_rootMeanSquare^6;
% fun5 = -3*E*Iz*(qe.'*sfunction_p(x,unit_l).'*sfunction_pp(x,unit_l)*qe)^2*sfunction_p(x,unit_l).'*sfunction_p(x,unit_l)*qe/f_rootMeanSquare^8;
fun_accurate = fun2 - fun1;
Qfe = int(fun_accurate,x,[0 unit_l]);
Qfe_fun = matlabFunction(Qfe);%,'File','myfile','Optimize',false);
Qfa = zeros(4*(units_num+1),1);
Ma_inv = inv(Ma);
toc;
tic;
% options = odeset('Mass',M,'Stats','on','OutputFcn',@myOutputFcn,'RelTol',1e-10,'AbsTol',1e-12);
options = odeset('Stats','on','OutputFcn',@myOutputFcn,'RelTol',1e-6,'AbsTol',1e-8);
[t,y] = ode113(@(t,y) odefun(t,y,units_num,Qga,Ma_inv,Qfe_fun,Qfa),0:0.001:t_end,yy0,options);
toc;
tic;
data_e = [zeros(size(t,1),2) y(:,1:2:end)];
data_ep = [zeros(size(t,1),2) y(:,2:2:end)];
%组装
Ma = zeros(4*(units_num+1),4*(units_num+1));
Qga = zeros(4*(units_num+1),1);
for index = 1:units_num
    step = 4*(index-1);
    Ma(1+step:8+step,1+step:8+step) = Ma(1+step:8+step,1+step:8+step) + Me;
    Qga(1+step:8+step,1) = Qga(1+step:8+step,1) + Qge;
end
ke = diag(data_ep*Ma*data_ep')/2;
kg = data_e*Qga;
% f_rootMeanSquare = (qe.') * (sfunction_p(x,unit_l).') * sfunction_p(x,unit_l) * qe;
kf1 = int(E*A*GreenLagrangian^2,x,[0 unit_l])/2;
kf2 = int(E*Iz*cross([sfunction_p(x,unit_l)*qe;0],[sfunction_pp(x,unit_l)*qe;0]).'*cross([sfunction_p(x,unit_l)*qe;0],[sfunction_pp(x,unit_l)*qe;0])/f_rootMeanSquare^6,x,[0 unit_l])/2;
kf_fun = matlabFunction(kf1+kf2);
kf = 0;
for i = 1:units_num
    step = 4*(i-1);
    arg = num2cell(data_e(:,1+step:8+step)',2);
    kf = kf + kf_fun(arg{:});
end
kf = kf';
toc;

%%
%求解迭代函数
function dudt = odefun(~,u,units_num,Qga,Ma_inv,Qfe_fun,Qfa)
dudt = u;
%提取总体节点坐标列阵
Qe = u(1:2:end);
%单元弹性力更新
arg = [[0 0]';Qe(1:6)];
Qfa(1:8) = Qfe_fun(arg(1),arg(2),arg(3),arg(4),arg(5),arg(6),arg(7),arg(8));
for index = 2:units_num
    step = 4*(index-1);
    arg = Qe(4*(index-1)-1:4*(index+1)-2);
    Qfa(1+step:8+step) = Qfa(1+step:8+step) + Qfe_fun(arg(1),arg(2),arg(3),arg(4),arg(5),arg(6),arg(7),arg(8));
end
dudt(1:2:end) = u(2:2:end);
dudt(2:2:end) = Ma_inv*(Qga + Qfa(3:end,:));
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
function s = sfunction_p(x,m_l)
e = x/m_l;
s1 = -6*e + 6*e^2;
s2 = 1 - 4*e + 3*e^2;
s3 = 6*e - 6*e^2;
s4 = -2*e + 3*e^2;
s(1,1) = s1;
s(2,2) = s1;
s(1,3) = s2*m_l;
s(2,4) = s2*m_l;
s(1,5) = s3;
s(2,6) = s3;
s(1,7) = s4*m_l;
s(2,8) = s4*m_l;

s = s/m_l;
end %形函数一阶导数矩阵
function s = sfunction_pp(x,m_l)
e = x/m_l;
s1 = -6 + 12*e;
s2 = -4 + 6*e;
s3 = 6 - 12*e;
s4 = -2 + 6*e;
s(1,1) = s1;
s(2,2) = s1;
s(1,3) = s2*m_l;
s(2,4) = s2*m_l;
s(1,5) = s3;
s(2,6) = s3;
s(1,7) = s4*m_l;
s(2,8) = s4*m_l;

s = s/m_l/m_l;
end %形函数二阶导数矩阵
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






