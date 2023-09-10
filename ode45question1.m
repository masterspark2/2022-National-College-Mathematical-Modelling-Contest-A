close all
clc;
clear;
global m1 m2 k omega f m3
m3 = 1335.535;
m1 = 4866 + m3;
m2 = 2433;
k = 80000;
omega = 1.4005;
f = 6250;
Time = 2*pi/omega;
Timetotal = 40*Time;
dt = 0.2;
timespan = 0:dt:Timetotal;
y0 = [0;0;0;0];
[timespan1,solution1] = ode45(@constantode1,timespan,y0);
[timespan2,solution2] = ode45(@noconstantode2,timespan,y0);
figure(1) 
set (gcf,'Position',[100 100 800 400])
subplot(211)
plot(timespan1,solution1(:,1));
xlabel('time/s')
ylabel('displacement/m')
title('定常阻尼系数,浮子位移')
subplot(212)
plot(timespan1,solution1(:,3));
xlabel('time/s')
ylabel('displacement/m')
title('非定常阻尼系数,振子位移')
figure(2) 
set (gcf,'Position',[100 100 800 400])
subplot(211)
plot(timespan1,solution1(:,2));
xlabel('time/s')
ylabel('displacement/m')
title('定常阻尼系数,浮子速度')
subplot(212)
plot(timespan1,solution1(:,4));
xlabel('time/s')
ylabel('displacement/m')
title('非定常阻尼系数,振子速度')
for time =[10 20 40 60 100]
    n = time/dt +1;
    d1_fuzi = solution1(n,1);
    v1_fuzi = solution1(n,2);
    d1_zhenzi = solution1(n,3);
    v1_zhenzi = solution1(n,4);
    fprintf('\t时间为%ds时，浮子位移为%f/m,速度为%f/(m/s);振子位移为%f/m，速度为%f/(m/s).\n',time,d1_fuzi,v1_fuzi,d1_zhenzi,v1_zhenzi)
end
disp('第二种情况：非定常阻尼系数')
for time =[10 20 40 60 100]
    n = time/dt +1;
    d1_fuzi = solution2(n,1);
    v1_fuzi = solution2(n,2);
    d1_zhenzi = solution2(n,3);
    v1_zhenzi = solution2(n,4);
    fprintf('\t时间为%ds时，浮子位移为%f/m,速度为%f/(m/s);振子位移为%f/m，速度为%f/(m/s).\n',time,d1_fuzi,v1_fuzi,d1_zhenzi,v1_zhenzi)
end
a = max(solution2(:,1));
b = max(solution2(:,2));
c = max(solution2(:,3));
d = max(solution2(:,4));
fprintf('非常阻尼浮子最大位移为%f/m,最大速度为%f/(m/s);振子最大位移为%f/m，最大速度为%f/(m/s).\n',a,b,c,d)
a = max(solution1(:,1));
b = max(solution1(:,2));
c = max(solution1(:,3));
d = max(solution1(:,4));
fprintf('常阻尼浮子最大位移为%f/m,最大速度为%f/(m/s);振子最大位移为%f/m，最大速度为%f/(m/s).\n',a,b,c,d)
function differentitaly =constantode1(t,y)
global m1 m2 k omega f m3
c2=10000;
a=0.8;
c1=656.3616;
h=1.9447;
differentitaly = zeros(4,1);
differentitaly(1) = y(2);
if y(1)> h
    F=1025*9.8/(3*(a^2))*pi*(a+h-y(1))^3;
else 
    F=1025*9.8*(pi*a/3+pi*(h-y(1))); 
end
F=-F+(m1+m2-m3)*9.8;
differentitaly(2) = -(c2+c1)/m1*y(2)+c2/m1*y(4)+(k*(y(3)-y(1))-F)/m1+f/m1*cos(omega*t);
differentitaly(3) = y(4);
differentitaly(4) = -c2/m2*(y(4)-y(2))-k/m2*(y(3)-y(1));
end
function differentitaly =noconstantode2(t,y)
global m1 m2 k omega f m3
c2=10000*sqrt(abs(y(2)-y(4))); 
a=0.8;
c1=656.3616;
h=1.9447;
differentitaly = zeros(4,1);
differentitaly(1) = y(2);
if y(1)> h
    F=1025*9.8/(3*(a^2))*pi*(a+h-y(1))^3;
else 
    F=1025*9.8*(pi*a/3+pi*(h-y(1))); 
end
F=-F+(m1+m2-m3)*9.8;
differentitaly(2) = -(c2+c1)/m1*y(2)+c2/m1*y(4)+(k*(y(3)-y(1))-F)/m1+f/m1*cos(omega*t);
differentitaly(3) = y(4);
differentitaly(4) = -c2/m2*(y(4)-y(2))-k/m2*(y(3)-y(1));
end


