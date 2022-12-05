subplot(3,1,1:2);
plot(t,ke,t,-kg,t,kf,t,ke-kg+kf);
h1 = legend("动能","重力势能","弹性势能","总能量",'Location','north');
set(h1,'Orientation','horizon');
subplot(3,1,3);
plot(t,ke-kg+kf);
h2 = legend("总能量",'Location','north');
set(h2,'Orientation','horizon');%,'Box','off');