clc;
close all;
clear all;

%绘制函数图象在画布1
figure(1);
lbx=-1;ubx=1;
lby=-1;uby=1;
ezmesh('20*exp(-0.2*sqrt(0.5*(X^2+Y^2)))+exp(0.5*(cos(2*pi*X)+cos(2*pi*Y)))',[lbx,ubx,lby,uby],100);
hold on;

% 参数初始化
NIND=500;        %种群大小
MAXGEN=300;      %最大遗传代数
PRECI=32;       %个体编码的二进制位数
GGAP=1;      %代沟值
SUBPOP=1;     %子种群的数量
px=0.7;         %交叉概率
pm=0.01;        %变异概率
trace=zeros(3,MAXGEN);                        %寻优结果的初始值
FieldD=[PRECI,PRECI;lbx,lby;ubx,uby;0,0;0,0;1,1;1,1];                      %译码矩阵FieldD
Chrom=crtbp(NIND,PRECI*2);                      %初始种群

% 求解
gen=0;                                  %遗传代数计数器
XY=bs2rv(Chrom,FieldD);                 %计算初始种群的二进制转换成十进制
X=XY(:,1);Y=XY(:,2);
ObjVCH=20*exp(-0.2*sqrt(0.5*(power(X,2)+power(Y,2))))+exp(0.5*(cos(2*pi.*X)+cos(2*pi.*Y)));        %计算初始目标函数值
while gen<MAXGEN
    FitnV=ranking(-ObjVCH);                              %分配适应度值(ranking由小到大排序)
    SelCh=select('rws',Chrom,FitnV,GGAP,SUBPOP);              %从种群中选择出个体
    SelCh=recombin('xovsp',SelCh,px);                  %对个体进行交叉操作
    SelCh=mut(SelCh,pm);                               %对个体进行变异操作
    XY=bs2rv(SelCh,FieldD);               %子代个体的二进制转换成十进制
    X=XY(:,1);Y=XY(:,2);
    ObjVSel=20*exp(-0.2*sqrt(0.5*(power(X,2)+power(Y,2))))+exp(0.5*(cos(2*pi.*X)+cos(2*pi.*Y)));             %计算子代的目标函数值
    [Chrom,ObjVCH]=reins(Chrom,SelCh,1,1,ObjVCH,ObjVSel); %将子代插入到父代，得到新种群
    XY=bs2rv(Chrom,FieldD);
    gen=gen+1;                                             %遗传代数计数器计数值+1

    %获取每代的最优解及其序号，Y为最优解,I为个体的序号
    [Y,I]=max(ObjVCH);
    trace(1:2,gen)=XY(I,:);                       %记下每代的最优值
    trace(3,gen)=Y;                               %记下每代的最优值
end
plot3(trace(1,:),trace(2,:),trace(3,:),'bo');                            %画出每代的最优点
grid on;
plot3(XY(:,1),XY(:,2),ObjVCH,'bo');  %画出最后一代的种群
hold off

% 画进化过程图在画布2
figure(2);
plot(1:MAXGEN,trace(3,:));
grid on
xlabel('代数')
ylabel('目标函数值')
title('过程')
bestZ=trace(3,end);
bestX=trace(1,end);
bestY=trace(2,end);
fprintf(['最优解:\nX=',num2str(bestX),'\nY=',num2str(bestY),'\nZ=',num2str(bestZ),'\n'])
