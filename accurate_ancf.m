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
units_num = 20;      %划分单元数量
t_end = 10;         %计算截止时间
initial_angle = 0*pi/180;  %柔性梁初始角度，水平起始
AbsTol = 1e-9;
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
syms x;
fun = Ro*A*sfunction(x,unit_l)'*sfunction(x,unit_l);
Me = int(fun,x,0,unit_l);
%计算广义单元重力阵
G = [0 -Ro*g]';
fun = A*sfunction(x,unit_l)';
Qge = int(fun,x,0,unit_l)*G;
%组装
Ma = zeros(4*(units_num+1),4*(units_num+1));
Qga = zeros(4*(units_num+1),1);
for index = 1:units_num
    step = 4*(index-1);
    Ma(1+step:8+step,1+step:8+step) = Ma(1+step:8+step,1+step:8+step) + Me;
    Qga(1+step:8+step,1) = Qga(1+step:8+step,1) + Qge;
end
Ma(abs(Ma) < 1e-10) = 0;
Qga(abs(Qga) < 1e-10) = 0;
Ma(1:2,:) = [];
Ma(:,1:2) = [];
Qga(1:2,:) = [];
Ma_inv = inv(Ma);
%赋初值
y0 = zeros(1,4*(units_num+1)*2);
y0(1:2:end) = points_e(:);
yy0 = y0(:,5:end);
% %计算常质量阵
% bin = zeros(4*(units_num+1)*2-4,1);
% bin(1:2:end) = 1;
% M = diag(bin,0);
% M(2:2:end,2:2:end) = Ma;
Qfa = zeros(4*(units_num+1),1);
%计算广义弹性力表达式
syms qe [8 1];
f = sqrt(qe.' * sfunction_p(x,unit_l).' * sfunction_p(x,unit_l) * qe);
GreenLagrangian = (qe.' * sfunction_p(x,unit_l).' * sfunction_p(x,unit_l) * qe - 1)/2;
fun1 = E*A*GreenLagrangian * sfunction_p(x,unit_l).' * sfunction_p(x,unit_l) * qe;
fun2 = -E*Iz*sfunction_pp(x,unit_l).' * sfunction_pp(x,unit_l) * qe / f^4;
fun3 = 2*E*Iz*qe.'*sfunction_pp(x,unit_l).'*sfunction_pp(x,unit_l)*qe*sfunction_p(x,unit_l).'*sfunction_p(x,unit_l)*qe/f^6;
fun4 = E*Iz*qe.'*sfunction_p(x,unit_l).'*sfunction_pp(x,unit_l)*qe*(sfunction_p(x,unit_l).'*sfunction_pp(x,unit_l)+sfunction_pp(x,unit_l).'*sfunction_p(x,unit_l))*qe/f^6;
fun5 = -3*E*Iz*(qe.'*sfunction_p(x,unit_l).'*sfunction_pp(x,unit_l)*qe)^2*sfunction_p(x,unit_l).'*sfunction_p(x,unit_l)*qe/f^8;
qfe_fun = matlabFunction(-fun1+fun2+fun3+fun4+fun5);
%生成弹性力列阵更新函数
Cal_Qfe = @(qe,unit_l) integral(@(x) qfe_fun(qe(1),qe(2),qe(3),qe(4),qe(5),qe(6),qe(7),qe(8),x),0,unit_l,'ArrayValued',true,'RelTol',AbsTol,'AbsTol',AbsTol);
toc;
%%
%求解及后处理
tic;
% options = odeset('Mass',M,'Stats','on','OutputFcn',@myOutputFcn,'RelTol',1e-11,'AbsTol',1e-12);
options = odeset('Stats','on','OutputFcn',@myOutputFcn,'RelTol',AbsTol,'AbsTol',AbsTol);
[t,y] = ode113(@(t,y) odefun(t,y,units_num,unit_l,Qga,Ma_inv,Qfa,zeros(4*(units_num+1)*2-4,1),Cal_Qfe),0:0.001:t_end,yy0,options);
toc;
tic;
data_e = [zeros(size(t,1),2) y(:,1:2:end)];
data_ep = [zeros(size(t,1),2) y(:,2:2:end)];
%重新组装
Ma = zeros(4*(units_num+1),4*(units_num+1));
Qga = zeros(4*(units_num+1),1);
for index = 1:units_num
    step = 4*(index-1);
    Ma(1+step:8+step,1+step:8+step) = Ma(1+step:8+step,1+step:8+step) + Me;
    Qga(1+step:8+step,1) = Qga(1+step:8+step,1) + Qge;
end
%计算能量（小规模）
% ke = diag(data_ep*Ma*data_ep')/2;
% kg = data_e*Qga;
% kf = Cal_kf(data_e,E,A,Iz,units_num,unit_l);
toc;

%%
%求解迭代函数
function dudt = odefun(~,u,units_num,unit_l,Qga,Ma_inv,Qfa,zero,Cal_Qfe)
dudt = zero;
%提取总体节点坐标列阵
Qe = u(1:2:end);
%单元弹性力更新
Qfa(1:8) = Cal_Qfe([[0 0]';Qe(1:6)],unit_l);
for index = 2:units_num
    step = 4*(index-1);
    Qfa(1+step:8+step) = Qfa(1+step:8+step) + Cal_Qfe(Qe(4*(index-1)-1:4*(index+1)-2),unit_l);
end
dudt(1:2:end) = u(2:2:end);
dudt(2:2:end) = Ma_inv * (Qga + Qfa(3:end,:));
end
%%
%函数定义
%计算能量（大规模）
function [ke,kg,kf] = Cal_kf(data_e,data_ep,Ma,Qga,E,A,Iz,units_num,unit_l)
syms x;
syms qe [8 1];
f = sqrt(qe.' * sfunction_p(x,unit_l).' * sfunction_p(x,unit_l) * qe);
GreenLagrangian = qe.' * sfunction_p(x,unit_l).' * sfunction_p(x,unit_l) * qe/2 - 1/2;
kf1 = E*A*GreenLagrangian^2/2;
kf2 = E*Iz*cross([sfunction_p(x,unit_l)*qe;0],[sfunction_pp(x,unit_l)*qe;0]).'*cross([sfunction_p(x,unit_l)*qe;0],[sfunction_pp(x,unit_l)*qe;0])/f^6/2;
kf_fun = matlabFunction(kf1+kf2);
ke = zeros(size(data_e,1),1);
kg = zeros(size(data_e,1),1);
kf = zeros(size(data_e,1),1);
parfor t = 1:size(data_e,1)
    ke(t) = ke(t) + data_ep(t,:)*Ma*data_ep(t,:)'/2;
    kg(t) = kg(t) + data_e(t,:)*Qga;
    for i = 1:units_num
        step = 4*(i-1);
        arg = data_e(t,1+step:8+step);
        kf(t) = kf(t) + integral(@(x) kf_fun(arg(1),arg(2),arg(3),arg(4),arg(5),arg(6),arg(7),arg(8),x),0,unit_l,'RelTol',1e-10,'AbsTol',1e-12);
    end
    if mod(t,1000) == 0
        disp(t);
    end
end
end
%返回给定点的形函数矩阵
function s = sfunction(x,m_l)
e = x/m_l;
% s = zeros(2,8);
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
end 
%形函数一阶导数矩阵
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
end 
%形函数二阶导数矩阵
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
end

function status = myOutputFcn(t,~,flag)
if flag == "init"
elseif flag == "done"
else
    disp(t);
%     if mod(t,0.1) == 0
%         disp(t);
%     end
    status = 0;
end
end






