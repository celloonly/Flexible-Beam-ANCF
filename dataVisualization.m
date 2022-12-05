%%
%获取当前工作区的结果数据
m_data = evalin('base','data_e');
m_t = evalin('base','t');
unit_resolving = 20;

%%
%整理数据
points_num = size(m_data,2)/4; %计算节点数
units_num = points_num - 1;
inital_position = m_data(1,1:4:end); %节点初始位置
units_length = zeros(1,units_num); %储存每个单元的长度
points_r = zeros(unit_resolving,units_num*2,size(m_data,1));
for index = 1:units_num
    units_length(1,index) = inital_position(1,index+1) - inital_position(1,index);
end

row_num = 1;
for unit_x = 0:units_length(1,1)/(unit_resolving-1):units_length(1,1)
    column_num = 1;
    for index = 1:units_num
        points_r(row_num,column_num:column_num+1,:) = sfunction(unit_x,units_length(1,1))*m_data(:,1+4*(index-1):4*(index+1))';
        column_num = column_num + 2;
    end
    row_num = row_num + 1;
end
global_r1 = reshape(points_r(:,1:2:end,:),unit_resolving*units_num,1,size(t,1));
global_r2 = reshape(points_r(:,2:2:end,:),unit_resolving*units_num,1,size(t,1));
for index = flip(unit_resolving:unit_resolving:unit_resolving*units_num-1)
    global_r1(index,:,:) = [];
    global_r2(index,:,:) = [];
end

%%
%绘制结果
%F(size(t,1)) = struct('cdata',[],'colormap',[]);
temp = plot(NaN);
for m_time = 1:size(t,1)
  set(temp, 'XData', global_r1(:,1,m_time), 'YData', global_r2(:,1,m_time),'LineWidth',2);
  ylim([-2 0.5]);
  xlim([-2 2]);
  drawnow limitrate
%   pause(0.0001)
  %F(m_time) = getframe;
end
%movie(F,1,size(t,1)/2.5);


%%
%形函数定义
function s = sfunction(unit_x,unit_l)
e = unit_x/unit_l;
s = zeros(2,8);
s1 = 1 - 3*e^2 + 2*e^3;
s2 = e - 2*e^2 + e^3;
s3 = 3*e^2 - 2*e^3;
s4 = -e^2 + e^3;
s(1,1) = s1;
s(2,2) = s1;
s(1,3) = s2*unit_l;
s(2,4) = s2*unit_l;
s(1,5) = s3;
s(2,6) = s3;
s(1,7) = s4*unit_l;
s(2,8) = s4*unit_l;
end