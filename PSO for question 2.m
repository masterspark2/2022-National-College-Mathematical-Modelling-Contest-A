close all;
clc;
format long
tic;
c23 = 2;                            
c13 = 0.75;                         
w = 0.2;                           
xlimit = [0,100000;...
  0,1];
popsize = 50;                      
maxgen = 35;                       
vrange = 0.3;                                         
SolveTimes = 1;                    

popdimension = size(xlimit,1);
vlimit = vrange*[-xlimit(:,2)+xlimit(:,1),xlimit(:,2)-xlimit(:,1)];
xrecode = zeros(SolveTimes,popdimension); 
fitrecord = zeros(SolveTimes,maxgen);
 
for times = 1:1:SolveTimes
    
    pop= repmat(xlimit(:,1)',popsize,1)+repmat(diff(xlimit'),popsize,1).*rand(popsize,popdimension);
    v = repmat(vlimit(:,1)',popsize,1)+repmat(diff(vlimit'),popsize,1).*rand(popsize,popdimension);
    
    fitness = zeros(popsize,1);
    for i=1:popsize
        fitness(i)=fcn(pop(i,:));
    end
    
    
  
    [Firstbestfitness,Firstbestindex]=max(fitness,[],'all','linear');   
    gbestnum = pop(Firstbestindex,:);    
    gbesfitness = Firstbestfitness;       
    pbestnum = pop;        
    pbestfitness = fitness;        
    clear Firstbestfitness Firstbestindex
    
    %
    for m=1:maxgen  
        
        
        v = w*v + ...
            c13*rand(popsize,1).*(pbestnum -pop)+...
            c23*rand(popsize,1).*(repmat(gbestnum,popsize,1) - pop);
        
       
        vlimitbelow = repmat(vlimit(:,1).',popsize,1);
        vlimitup = repmat(vlimit(:,2).',popsize,1);
        bigger = v > vlimit(:,2).';
        v(bigger) = vlimitup(bigger);
        smaller = v < vlimit(:,1).';
        v(smaller) = vlimitbelow(smaller);
        
       
        pop = pop +v;
        
        
        poplimitbelow = repmat(xlimit(:,1).',popsize,1);
        poplimitup = repmat(xlimit(:,2).',popsize,1);
        bigger = pop > xlimit(:,2).';
        pop(bigger) = poplimitup(bigger);
        smaller = pop < xlimit(:,1).';
        pop(smaller) = poplimitbelow(smaller);
        
      
        for i=1:popsize
            fitness(i)=fcn(pop(i,:));
        end   
        
        
        toreplacepop = (fitness > pbestfitness);
        pbestfitness(toreplacepop) = fitness(toreplacepop);    
        pbestnum(toreplacepop,:) = pop(toreplacepop,:);
        
       
        [maxfitness,maxfitnessindex] = max(fitness,[],'all','linear');
        if maxfitness > gbesfitness
            gbesfitness = maxfitness;    
            gbestnum = pop(maxfitnessindex,:);
        end
        
       
        fitrecord(times,m)=gbesfitness;
        xrecode(times,:)=gbestnum;
    end
end
[PSObestfit,bestfitindex] = max(fitrecord(:,end),[],'all','linear');
PSObestnum = xrecode(bestfitindex,:);
disp('%%%%%%%(粒子群算法)%%%%%%%');
disp(['最优值为:',num2str(PSObestfit)]);
disp(['最优解为:',num2str(PSObestnum)]);
toc
for i=1:1
    Fontsize = 10.5;
    LineWidth = 1.2;
    f1 = figure;
    f1.Visible = 'on';
    f1.Units = 'centimeters';
    f1.Position = [20 8 16 9];
    f1.Color = 'w';
    f1.NumberTitle = 'off';
    f1.Name = '迭代收敛图';
    f1.Resize = 'off';
    
    ax1=axes(f1);
    ax1.Box = 'off';
    ax1.Visible = 'on';
    ax1.NextPlot = 'add';
    ax1.Color= 'w';
    ax1.TickDir = 'in';
    ax1.FontUnits =  'points';
    ax1.FontSize = Fontsize;
    ax1.LineWidth = LineWidth;
    ax1.FontName = 'Times New Roman';
    ax1.FontWeight = 'bold';
    ax1.TitleFontSizeMultiplier = 1;
    ax1.LabelFontSizeMultiplier = 1;
    ax1.XColor = 'k';
    ax1.YColor = 'k';
    ax1.YAxisLocation = 'left';
    ax1.XAxisLocation = 'bottom';
    ax1.TickLength = [0.01,0.035];
    ax1.XDir = 'normal';
    ax1.YDir = 'normal';
    grid(ax1,'on');
    ax1.GridLineStyle = '--';
    ax1.GridAlpha = 0.85;
    ax1.GridColor =  [0.85,0.85,0.85];
    
    x1 = 1:1:maxgen;
    y1 = fitrecord(bestfitindex,:);
    p1 = plot(ax1,x1,y1);
    p1.Visible = 'on';
    p1.LineWidth = LineWidth;
    p1.Color = 'k';
    p1.LineStyle = '-';
    p1.LineJoin = 'round';
    p1.Marker = 'none';
    p1.MarkerIndices = 1:round(length(x1)/10):length(x1);
    p1.MarkerSize = 6;
    p1.MarkerEdgeColor = 'auto';
    p1.MarkerFaceColor = 'none';
    ax1.XScale = 'linear';
    ax1.YScale = 'linear';
    axis(ax1,[min(x1,[],'all'),max(x1,[],'all'),min(y1,[],'all') - 1,max(y1,[],'all') + 1]);
    ax1.XTickMode = 'auto';
    ax1.YTickMode = 'auto';
    ax1.XTickLabelMode = 'auto';
    ax1.YTickLabelMode = 'auto';
    xlabel('迭代次数','FontSize',Fontsize,'FontName','宋体');
    ylabel('适应度','FontSize',Fontsize,'FontName','宋体');
end
%适应度函数
function y=fcn(x)
global m1 m2 k omega f1 m3 c2 bilixishu mizhishu
m3 = 1165.992;
m1 = 4866 + m3;
m2 = 2433;
k = 80000;
omega = 2.2413;
f1 = 4890;
T = 2*pi/omega;
Ttotal = 40*T;
dt = 0.1;
tspan = Ttotal/2:dt:Ttotal;
y0 = [0;0;0;0];
bilixishu=x(:,1);
mizhishu=x(:,2);
[tt2,yy2]=ode45(@noconstantode2,tspan,y0);
vr=yy2(:,2)-yy2(:,4);
c2=bilixishu*(abs(vr)).^mizhishu;
summ=2*trapz(tt2,c2.*vr.^2)/Ttotal;
y = summ;
end
function differentitaly =noconstantode2(t,y,bilixishu,mizhishu)
global m1 m2 k omega f1 m3 bilixishu mizhishu
a=0.8;
c1=167.8395;
h=1.9447;
differentitaly = zeros(4,1);
differentitaly(1) = y(2);
if y(1)> h
    F=1025*9.8/(3*(a^2))*pi*(a+h-y(1))^3;
else 
    F=1025*9.8*(pi*a/3+pi*(h-y(1))); 
end
F=-F+(m1+m2-m3)*9.8;
c2=bilixishu*abs(y(2)-y(4))^mizhishu; 
differentitaly(2) = -(c2+c1)/m1*y(2)+c2/m1*y(4)+(k*(y(3)-y(1))-F)/m1+f1/m1*cos(omega*t);
differentitaly(3) = y(4);
differentitaly(4) = -c2/m2*(y(4)-y(2))-k/m2*(y(3)-y(1));
end
