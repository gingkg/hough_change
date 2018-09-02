%����Hough�任
%Author:Shen Baoyin
%Time:2018/8/1
close all
clear all

target=2;%Ŀ����
n=3;%��ʼ����
k=90;%sig�ֵĸ���
m=500;%p�ֵĸ���
Monte_Carlo=100;%Monte_Carlo�������
L=150;%�״�������� ��λ��km��
Vmax=0.7;%Ŀ������ٶ� ��λ��km/s��
Vmin=0.05;%Ŀ����С�ٶ�

%Ŀ����ʼ���꼰�ٶ�
x1=40;y1=20;vx1=0.3;vy1=0.18;%��λkm,km/s
x2=20;y2=80;vx2=0.3;vy2=-0.18;
Ts=4;%�������ڣ���λs

success=zeros(Monte_Carlo,target);%Ŀ�꺽���ɹ���ʼ����
fake(1:Monte_Carlo)=0;%Ŀ�꺽�������ʼ����
track_number(1:Monte_Carlo)=0;%�ܺ�����ʼ��

N=0:n-1;
X1_init=x1+Ts*N*vx1;%��ʵ����1
Y1_init=y1+Ts*N*vy1;
Y1_0=y1-vy1*x1/vx1;
offset(1)=Y1_0*cos(atan(abs(vy1)/abs(vx1)));%����1��ʵ����
X2_init=x2+Ts*N*vx2;%��ʵ����2
Y2_init=y2+Ts*N*vy2;
Y2_0=y2-vy2*x2/vx2;
offset(2)=Y2_0*cos(atan(abs(vy2)/abs(vx2)));%����2��ʵ����

Np=1:k;
dNp=pi/k;%�����ռ�Ƕȼ��
angle=(Np-1/2)*dNp;

dMp=6*0.1;%�����ռ䴹����
 
R = poissrnd(20,1,n);%ÿ���Ӳ����������Ӳ��ɷֲ�
Rn=R(1);
X_za=unifrnd (0, 100, 1, R(1));%���������x,y�����꣬����0-100�ľ��ȷֲ�
Y_za=unifrnd (0, 100, 1, R(1));
for i=2:n
    X_za(Rn+1:Rn+R(i))=unifrnd (0, 100, 1, R(i));
    Y_za(Rn+1:Rn+R(i))=unifrnd (0, 100, 1, R(i));
    Rn=Rn+R(i);
end

noisex=normrnd(0,0.05,1,n);%x��������
noisey=normrnd(0,0.05,1,n);%y��������

X1=X1_init+noisex;X2=X2_init+noisex;%ʵ������
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

%����ʽ����
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

%���㴹��
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

%����㴦�Ƕ�
for i=1:number-1
    x11=Zc1(i,1);y11=Zc1(i,2);
    x21=Zc2(i,1);y21=Zc2(i,2);
    x31=Zc3(i,1);y31=Zc3(i,2);
    sign12(i)=atan((x21-x11)/(y11-y21));
    sign23(i)=atan((x31-x21)/(y21-y31));
end

%Ѱ�ҷ���Ҫ��������
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

%�������Ҫ��Ĳ���
for h=1:count-1
    for j=1:6
        fprintf('%f ',para(h,j));
    end
    fprintf('\n');
end
fprintf('*****************************************\n');

%��ͼ
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
legend('��ʵ����1','��ʵ����2','�Ӳ�');
title('����ͼ')
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
    legend('��ʼ����');
    title('��ʼ���ͼ')
    hold on
    axis([0 100 0 100])
    axis on
end
