%修正Hough变换
%Author:Shen Baoyin
%Time:2018/8/1
close all
clear all

target=2;%目标数
n=3;%起始拍数
k=90;%sig分的个数
m=500;%p分的个数
Monte_Carlo=100;%Monte_Carlo仿真次数
L=150;%雷达量测距离 单位（km）
Vmax=0.7;%目标最大速度 单位（km/s）
Vmin=0.05;%目标最小速度

%目标起始坐标及速度
x1=40;y1=20;vx1=0.3;vy1=0.18;%单位km,km/s
x2=20;y2=80;vx2=0.3;vy2=-0.18;
Ts=4;%采样周期，单位s

success=zeros(Monte_Carlo,target);%目标航迹成功起始矩阵
fake(1:Monte_Carlo)=0;%目标航迹虚假起始矩阵
track_number(1:Monte_Carlo)=0;%总航迹起始数

N=0:n-1;
X1_init=x1+Ts*N*vx1;%真实航迹1
Y1_init=y1+Ts*N*vy1;
Y1_0=y1-vy1*x1/vx1;
offset(1)=Y1_0*cos(atan(abs(vy1)/abs(vx1)));%航迹1真实垂距
X2_init=x2+Ts*N*vx2;%真实航迹2
Y2_init=y2+Ts*N*vy2;
Y2_0=y2-vy2*x2/vx2;
offset(2)=Y2_0*cos(atan(abs(vy2)/abs(vx2)));%航迹2真实垂距

Np=1:k;
dNp=pi/k;%参数空间角度间隔
angle=(Np-1/2)*dNp;

dMp=6*0.1;%参数空间垂距间隔
 
R = poissrnd(20,1,n);%每拍杂波个数，服从泊松分布
Rn=R(1);
X_za=unifrnd (0, 100, 1, R(1));%随机产生（x,y）坐标，服从0-100的均匀分布
Y_za=unifrnd (0, 100, 1, R(1));
for i=2:n
    X_za(Rn+1:Rn+R(i))=unifrnd (0, 100, 1, R(i));
    Y_za(Rn+1:Rn+R(i))=unifrnd (0, 100, 1, R(i));
    Rn=Rn+R(i);
end

noisex=normrnd(0,0.05,1,n);%x量测噪声
noisey=normrnd(0,0.05,1,n);%y量测噪声

X1=X1_init+noisex;X2=X2_init+noisex;%实际量测
Y1=Y1_init+noisey;Y2=Y2_init+noisey;
    
Z1(1,1:2)=[X1(1),Y1(1)];
Z1(2,1:2)=[X2(1),Y2(1)];
Z1(3:2+R(1),1:2)=[(X_za(1:R(1)))',(Y_za(1:R(1)))'];
Z2(1,1:2)=[X1(2),Y1(2)];
Z2(2,1:2)=[X2(2),Y2(2)];
Z2(3:2+R(2),1:2)=[(X_za(1:R(2)))',(Y_za(1:R(2)))'];
Z3(1,1:2)=[X1(3),Y1(3)];
Z3(2,1:2)=[X2(3),Y2(3)];
Z3(3:2+R(3),1:2)=[(X_za(1:R(3)))',(Y_za(1:R(3)))'];

%遍历式关联
number=1;
for i=1:2+R(1)
    for j=1:2+R(2)
        for h=1:2+R(3)
            x11=Z1(i,1);y11=Z1(i,2);
            x21=Z2(j,1);y21=Z2(j,2);
            x31=Z3(h,1);y31=Z3(h,2);
            V12=sqrt((x11-x21)*(x11-x21)+(y11-y21)*(y11-y21))/Ts;
            V23=sqrt((x31-x21)*(x31-x21)+(y31-y21)*(y31-y21))/Ts;
            if Vmin<=V12 && V12<=Vmax && Vmin<=V23 && V23<=Vmax
                Zc1(number,1:2)=[x11,y11];            
                Zc2(number,1:2)=[x21,y21];
                Zc3(number,1:2)=[x31,y31];
                number=number+1;
            end
        end
    end
end

%计算垂距
for i=1:number-1
    x11=Zc1(i,1);y11=Zc1(i,2);
    x21=Zc2(i,1);y21=Zc2(i,2);
    x31=Zc3(i,1);y31=Zc3(i,2);
    for j=1:k 
        P1(i,j)=x11*cos(angle(j))+y11*sin(angle(j));
        P2(i,j)=x21*cos(angle(j))+y21*sin(angle(j));
        P3(i,j)=x31*cos(angle(j))+y31*sin(angle(j));
    end
end

%解过零处角度
for i=1:number-1
    x11=Zc1(i,1);y11=Zc1(i,2);
    x21=Zc2(i,1);y21=Zc2(i,2);
    x31=Zc3(i,1);y31=Zc3(i,2);
    sign12(i)=atan((x21-x11)/(y11-y21));
    sign23(i)=atan((x31-x21)/(y21-y31));
end

%寻找符合要求的量测点
count=1;
for i=1:number-1
    x11=Zc1(i,1);y11=Zc1(i,2);
    x21=Zc2(i,1);y21=Zc2(i,2);
    x31=Zc3(i,1);y31=Zc3(i,2);
    slope12=(y21-y11)/(x21-x11);
    slope23=(y31-y21)/(x31-x21);
    if abs(sign12(i)-sign23(i))<=3*dNp && (slope23/slope12)>0
        para(count,:)=[x11,y11,x21,y21,x31,y31];
        count=count+1;
    end
end

%输出符合要求的参数
for h=1:count-1
    for j=1:6
        fprintf('%f ',para(h,j));
    end
    fprintf('\n');
end
fprintf('*****************************************\n');

%绘图
figure
subplot(1,2,1)
scatter(X1,Y1,5,'r');
hold on
scatter(X2,Y2,6,'*','r');
hold on
scatter(X_za,Y_za,3,'filled','g')
hold on
xlabel('x(km)');
ylabel('y(km)');
legend('真实航迹1','真实航迹2','杂波');
title('量测图')
axis([0 100 0 100])
axis on

%figure
subplot(1,2,2)
for h=1:count-1
    X=[para(h,1),para(h,3),para(h,5)];
    YS=[para(h,2),para(h,4),para(h,6)];
    plot(X,YS,'b');
    xlabel('x(km)');
    ylabel('y(km)');
    legend('起始航迹');
    title('起始结果图')
    hold on
    axis([0 100 0 100])
    axis on
end
