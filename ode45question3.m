close all
clc;
clear;
global m1 m2 mfujia c1 k omega f seadensity gravity radius J1 J_fu C_huifu C_xingbo k_niuzhuan k_zuni L
m1=4866;
mfujia=1028.876;
J1=zhuandongguanliangfuzi;
J_fu=7001.914;
c1=683.4558;
C_xingbo=654.3383;
C_huifu=8890.7;
m2=2433;
k=80000;
k_niuzhuan=250000;
k_zuni=1000;
omega=1.7152;
f=3640;
L=1690;
seadensity=1025;
gravity=9.8;
radius=1;
Time=2*pi/omega;
Timetotal=40*Time;
dt=0.2;
timespan=0:dt:Timetotal;
y0=[0;0;0;0;0;0;0;0];
[timespan1,solution1]=ode45(@constantode1,timespan,y0);
figure(1)
subplot(211)
set(gcf,'Position',[100 100 700 400])
plot(timespan1,solution1(:,1));
xlabel('时间（s）','FontName','宋体')
ylabel('振幅（m）','FontName','宋体')
title('浮子垂荡位移')
subplot(212)
plot(timespan1,solution1(:,3));
xlabel('时间（s）','FontName','宋体')
ylabel('振幅（m）','FontName','宋体')
title('振子垂荡位移')
figure(2)
set(gcf,'Position',[100 400 700 400])
subplot(211)
plot(timespan1,solution1(:,5));
xlabel('时间（s）','FontName','宋体')
ylabel('转动角度（rad）','FontName','宋体')
title('浮子纵摇角位移')
subplot(212)
plot(timespan1,solution1(:,7));
xlabel('时间（s）','FontName','宋体')
ylabel('转动角度（rad）','FontName','宋体')
title('振子纵摇角位移')
figure(3)
subplot(211)
set(gcf,'Position',[100 100 700 400])
plot(timespan1,solution1(:,2));
xlabel('时间（s）','FontName','宋体')
ylabel('振幅（m）','FontName','宋体')
title('浮子垂荡速度')
subplot(212)
plot(timespan1,solution1(:,4));
xlabel('时间（s）','FontName','宋体')
ylabel('振幅（m）','FontName','宋体')
title('振子垂荡速度')
figure(4)
set(gcf,'Position',[100 400 700 400])
subplot(211)
plot(timespan1,solution1(:,6));
xlabel('时间（s）','FontName','宋体')
ylabel('转动角度（rad）','FontName','宋体')
title('浮子纵摇角速度')
subplot(212)
plot(timespan1,solution1(:,8));
xlabel('时间（s）','FontName','宋体')
ylabel('转动角度（rad）','FontName','宋体')
title('振子纵摇角速度')
for time=[10 20 40 60 100]
 n=time/dt+1;
 d_fuzi=solution1(n,1);
 v_fuzi=solution1(n,2);
 d_zhenzi=solution1(n,3);
 v_zhenzi=solution1(n,4); 
 angle_fuzi=solution1(n,5);
 anglev_fuzi=solution1(n,6);
 angle_zhenzi=solution1(n,7);
 anglev_zhenzi=solution1(n,8); 
 fprintf('\t时间为%ds时，浮子位移为%f/m,速度为%f/(m/s);振子位移为%f/m,速度为%f/(m/s);\n\t\t浮子角位移为%f/rad,角速度为%f/(rad/s);振子角位移为%f/rad,角速度为%f/(rad/s).\n',.....
  time,d_fuzi,v_fuzi,d_zhenzi,v_zhenzi,angle_fuzi,anglev_fuzi,angle_zhenzi,anglev_zhenzi);
end
a = max(solution1(:,1));
b = max(solution1(:,2));
c = max(solution1(:,3));
d = max(solution1(:,4));
e = max(solution1(:,5));
f123 = max(solution1(:,6));
g = max(solution1(:,7));
h = max(solution1(:,8));
fprintf('常阻尼浮子最大位移为%f/m,最大速度为%f/(m/s);振子最大位移为%f/m,最大速度为%f/(m/s);\n\t\t浮子最大角位移为%f/rad,最大角速度为%f/(rad/s);振子最大角位移为%f/rad,最大角速度为%f/(rad/s).\n',a,b,c,d,e,f123,g,h)
function differentitaly=constantode1(t,y)
global m1 m2 mfujia c1 k omega f seadensity gravity radius J1 J_fu C_huifu C_xingbo k_niuzhuan k_zuni L
h0=2;
c2=10000; 
differentitaly=zeros(8,1);
differentitaly(1)=y(2);
differentitaly(2)=-c1/(m1+mfujia)*y(2)-(pi*radius^2*seadensity*gravity/(m1+mfujia)*y(1)*(y(1)>=h0-3)+pi*radius^2*seadensity*gravity/(m1+mfujia)*(h0-3)*(y(1)<h0-3))-c2/(m1+mfujia)*(y(2)-y(4))-k/(m1+mfujia)*(y(1)-y(3))+f/(m1+mfujia)*cos(omega*t);
differentitaly(3)=y(4);
differentitaly(4)=-c2/m2*(y(4)-y(2))-k/m2*(y(3)-y(1));
differentitaly(5)=y(6);
differentitaly(6)=(L*cos(omega*t)-(C_huifu*y(5)+C_xingbo*y(6)+k_niuzhuan*(y(5)-y(7))+k_zuni*(y(6)-y(8))))/(J1+J_fu); 
differentitaly(7)=y(8);
J2=(0.5-m2*gravity/k+y(3)-y(1))^2*m2; 
differentitaly(8)=-(k_niuzhuan*(y(7)-y(5))+k_zuni*(y(8)-y(6)))/J2;%+2*m2*(y(4)-y(2))*y(8)*(0.5-m2*gravity/k+y(3)-y(1)))/J2; 
end
function J1=zhuandongguanliangfuzi()
m1=4866;
hshangdinggai=2/3*0.8;
hyuanzhuti=2.3;
hxiachuiti=3.8;
Sshangdinggai=pi*sqrt(1+0.8^2);
Syuanzhuti=6*pi;
Sxiadinggai=pi;
Sfuti=Sshangdinggai+Syuanzhuti+Sxiadinggai;
zc=Sshangdinggai/Sfuti*hshangdinggai+Syuanzhuti/Sfuti*hyuanzhuti+Sxiadinggai/Sfuti*hxiachuiti;
seadensity=m1/Sfuti; 
Ja=seadensity*sqrt(1+0.8^2)*(1/2+0.8^2)/4*2*pi+seadensity*sqrt(1+0.8^2)*pi*zc^2-2*seadensity*sqrt(1+0.8^2)*zc*0.8/3*2*pi; 
Jb=seadensity*(pi*3+2*pi*zc^2*3+2*pi/3*54.36-2*pi*zc*13.8); 
Jc=seadensity*(3.8-zc)^2*pi+seadensity*2*pi/2/4; 
J1=Ja+Jb+Jc; 
end