%序列Hough变换
%Author:Shen Baoyin
%Time:2018/8/1
close all
clear all

target=2;%目标数
targets=2;%序列数
n=15;%起始拍数
k=90;%sig分的个数
m=200;%p分的个数
Monte_Carlo=10;%Monte_Carlo仿真次数
L=150;%雷达量测距离
Pd=1;%检测概率

%目标起始坐标及速度
x1=40;y1=20;vx1=0.3;vy1=0.18;%单位km,km/s
x2=20;y2=80;vx2=0.3;vy2=-0.18;
Ts=4;%采样周期，单位s

success=zeros(Monte_Carlo,target);%目标航迹成功起始矩阵
fake(1:Monte_Carlo)=0;%目标航迹虚假起始矩阵
track_number(1:Monte_Carlo)=0;%总航迹起始数

N=0:n-1;
X1_init=x1+Ts*N*vx1;%目标航迹
Y1_init=y1+Ts*N*vy1;
Y1_0=y1-vy1*x1/vx1;
offset(1)=Y1_0*cos(atan(abs(vy1)/abs(vx1)));%航迹1真实垂距
X2_init=x2+Ts*N*vx2;
Y2_init=y2+Ts*N*vy2;
Y2_0=y2-vy2*x2/vx2;
offset(2)=Y2_0*cos(atan(abs(vy2)/abs(vx2)));%航迹2真实垂距

Np=1:k;
dNp=pi/k;%参数空间角度间隔
angle=(Np-1/2)*dNp;

dMp=6*0.1;%参数空间垂距间隔

for monte=1:Monte_Carlo
    clear A0 P0 R X_za Y_za noisex noisey
    R = poissrnd(50,1,n);%每拍杂波个数
    Rn=R(1);
    X_za=unifrnd (0, 100, 1, R(1));
    Y_za=unifrnd (0, 100, 1, R(1));
    for i=2:n
        X_za(Rn+1:Rn+R(i))=unifrnd (0, 100, 1, R(i));
        Y_za(Rn+1:Rn+R(i))=unifrnd (0, 100, 1, R(i));
        Rn=Rn+R(i);
    end

    noisex=normrnd(0,0.1,1,n);%x量测噪声
    noisey=normrnd(0,0.1,1,n);%y量测噪声

    X1=X1_init+noisex;X2=X2_init+noisex;%真实量测
    Y1=Y1_init+noisey;Y2=Y2_init+noisey;

    A=zeros(k,2*m);%积累矩阵

    %航迹1 Hough变换
    for i=1:n
        for j=1:k
            P(i,j)=X1(i)*cos(angle(j))+Y1(i)*sin(angle(j));
        end
    end

    %航迹2 Hough变换
    for i=(n+1):(2*n)
        for j=1:k
            P(i,j)=X2(i-n)*cos(angle(j))+Y2(i-n)*sin(angle(j));
        end
    end

    %杂波的Hough变换
    for i=2*n+1:(2*n+Rn)
        for j=1:k
            P(i,j)=X_za(i-2*n)*cos(angle(j))+Y_za(i-2*n)*sin(angle(j));
        end
    end

    %对积累矩阵投票
    for i=1:k
        for j=1:2*m
            a=-L+(j-1)*dMp;
            b=-L+j*dMp;
           for h=1:2*n+Rn
               if (P(h,i)>=a && P(h,i)<b) 
                   A(i,j)=A(i,j)+1;
               end
           end
        end
    end
    
     %寻找投票数最大的两个参数
    for s=1:targets
        max=A(1,1);
        maxi=1;maxj=1;
        for i=1:k
            for j=1:2*m
                if A(i,j)>=max
                    max=A(i,j);
                    maxi=i;
                    maxj=j;
                end
            end
        end
        for i=maxi-5:maxi+5
            for j=maxj-5:maxj+5
                if i<=0
                    i=1;
                end
                if j<=0
                    j=1;
                end
                A(i,j)=0;
            end
        end
        P0(s)=-L+(maxj-1/2)*dMp;
        A0(s)=angle(maxi);
 
        %输出符合要求的参数
        flag=0;
        for din=1:target
            if abs(P0(s)-offset(din))<=3
                success(monte,din)=1;
                flag=1;
            end
        end
        if flag==0
            fake(monte)=fake(monte)+1;
        end
        fprintf('the value of P0 is %f;the value of A0 is %f\n',P0(s),A0(s)); 
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
        
    subplot(1,2,2)
    scatter(X1,Y1,5,'r');
    hold on
    scatter(X2,Y2,6,'*','r');
    hold on
    scatter(X_za,Y_za,3,'filled','g')
    hold on
    for s=1:targets     
        X=0:1:100;
        YS=(P0(s)-X*cos(A0(s)))/(sin(A0(s)));
        plot(X,YS,'b');
        hold on
        xlabel('x(km)');
        ylabel('y(km)');
        legend('真实航迹1','真实航迹2','杂波','起始航迹');
        title('起始结果图')
        hold on
        axis([0 100 0 100])
        axis on           
    end
end

%计算航迹起始成功率
success_number=0;
for i=1:Monte_Carlo
    for j=1:target
        success_number=success_number+success(i,j);
    end
end
success_rate=success_number/(Monte_Carlo*target);
fprintf('the rate of the successful Track initialization is %f%%\n',success_rate*100);

%计算航迹虚假率
fake_rate=fake./targets;
fake_mean=0;
for i=1:Monte_Carlo
    fprintf('the rate of the fake Track initialization is %f%%\n',fake_rate(i)*100);
    fake_mean=fake_mean+fake_rate(i);
end
fprintf('the mean of the fake Track initialization is %f%%\n',fake_mean/Monte_Carlo*100);
figure
plot(fake_rate.*100);
xlabel('Monte_Carlo');
ylabel('fake_rate(%)');
title('虚假航迹率');
