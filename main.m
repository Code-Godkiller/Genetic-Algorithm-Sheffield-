clc;
close all;
clear all;

%���ƺ���ͼ���ڻ���1
figure(1);
lbx=-1;ubx=1;
lby=-1;uby=1;
ezmesh('20*exp(-0.2*sqrt(0.5*(X^2+Y^2)))+exp(0.5*(cos(2*pi*X)+cos(2*pi*Y)))',[lbx,ubx,lby,uby],100);
hold on;

% ������ʼ��
NIND=500;        %��Ⱥ��С
MAXGEN=300;      %����Ŵ�����
PRECI=32;       %�������Ķ�����λ��
GGAP=1;      %����ֵ
SUBPOP=1;     %����Ⱥ������
px=0.7;         %�������
pm=0.01;        %�������
trace=zeros(3,MAXGEN);                        %Ѱ�Ž���ĳ�ʼֵ
FieldD=[PRECI,PRECI;lbx,lby;ubx,uby;0,0;0,0;1,1;1,1];                      %�������FieldD
Chrom=crtbp(NIND,PRECI*2);                      %��ʼ��Ⱥ

% ���
gen=0;                                  %�Ŵ�����������
XY=bs2rv(Chrom,FieldD);                 %�����ʼ��Ⱥ�Ķ�����ת����ʮ����
X=XY(:,1);Y=XY(:,2);
ObjVCH=20*exp(-0.2*sqrt(0.5*(power(X,2)+power(Y,2))))+exp(0.5*(cos(2*pi.*X)+cos(2*pi.*Y)));        %�����ʼĿ�꺯��ֵ
while gen<MAXGEN
    FitnV=ranking(-ObjVCH);                              %������Ӧ��ֵ(ranking��С��������)
    SelCh=select('rws',Chrom,FitnV,GGAP,SUBPOP);              %����Ⱥ��ѡ�������
    SelCh=recombin('xovsp',SelCh,px);                  %�Ը�����н������
    SelCh=mut(SelCh,pm);                               %�Ը�����б������
    XY=bs2rv(SelCh,FieldD);               %�Ӵ�����Ķ�����ת����ʮ����
    X=XY(:,1);Y=XY(:,2);
    ObjVSel=20*exp(-0.2*sqrt(0.5*(power(X,2)+power(Y,2))))+exp(0.5*(cos(2*pi.*X)+cos(2*pi.*Y)));             %�����Ӵ���Ŀ�꺯��ֵ
    [Chrom,ObjVCH]=reins(Chrom,SelCh,1,1,ObjVCH,ObjVSel); %���Ӵ����뵽�������õ�����Ⱥ
    XY=bs2rv(Chrom,FieldD);
    gen=gen+1;                                             %�Ŵ���������������ֵ+1

    %��ȡÿ�������Ž⼰����ţ�YΪ���Ž�,IΪ��������
    [Y,I]=max(ObjVCH);
    trace(1:2,gen)=XY(I,:);                       %����ÿ��������ֵ
    trace(3,gen)=Y;                               %����ÿ��������ֵ
end
plot3(trace(1,:),trace(2,:),trace(3,:),'bo');                            %����ÿ�������ŵ�
grid on;
plot3(XY(:,1),XY(:,2),ObjVCH,'bo');  %�������һ������Ⱥ
hold off

% ����������ͼ�ڻ���2
figure(2);
plot(1:MAXGEN,trace(3,:));
grid on
xlabel('����')
ylabel('Ŀ�꺯��ֵ')
title('����')
bestZ=trace(3,end);
bestX=trace(1,end);
bestY=trace(2,end);
fprintf(['���Ž�:\nX=',num2str(bestX),'\nY=',num2str(bestY),'\nZ=',num2str(bestZ),'\n'])
