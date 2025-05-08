clc
clear all;
close all;
%%
Noise_data=0.0;
sigma=1.5;
lamda=10;
lamda2=0.0;
GK=1.5;
I=1.5;
%% Problem setting 
fp=50;          % peak frequency of ricker in hertz
dt=0.001;       % sample time
Nt=2000;        % number of samples
t0=1.5/fp;      % timeshift of ricker
t=0:dt:(Nt-1)*dt;
fmax=4*fp;
fsample=1/dt;   % number of samples
df=fsample/Nt;
Srct=(1-0.5*(2*pi*fp*(t-t0)).^2).*exp(-0.25*(2*pi*fp*(t-t0)).^2); %雷克子波
Srcf=fft(Srct);

Q_define=[40,70,90,50,80,120];
c=[2500,2600,2800,2700,2900,3000];
p=[2.3,2.4,2.5,2.35,2.6,2.7];
c=[2500,2600,3500,3200,4500,5000];
p=[2.3,2.4,2.6,2.45,2.6,2.7];

dh=20;
h1=[20,40,60];
h2=[20,40,60,80];
h3=[20,40,60,80,100];
h4=[20,40,60,80,100,120];
h5=[20,40,60,80,100,120,140];
h6=[20,40,60,80,100,120,140,160];

Layer_number=length(Q_define);
for i=1:Layer_number
    eval(['Receiver_spacing(i)=length(h',num2str(i),');']);
end

f=0:fsample/Nt:fsample/2;
DATA=[Srct];
Recf=Srcf;
T=1;

for j=1:Layer_number
    Rect=zeros(Receiver_spacing(j),length(Srct));%介质1中时域接收信号矩阵，每一行代表一道2000采样点
    Recf_now=zeros(Receiver_spacing(j),length(Srcf));%介质1中频域接收信号矩阵，每一行代表一道2000采样点
    eval(['h=h',num2str(j),';']);
    TravelTime=h/c(j);
    for k=1:1:Receiver_spacing(j)
        Recf_now(k,1:Nt/2+1)=T*Recf(1:Nt/2+1).*exp(-pi*TravelTime(k).*f/Q_define(j)).*exp(-1i*2*pi.*f*TravelTime(k));%接收信号衰减
        Recf_now(k,Nt/2+2:end)=conj(Recf_now(k,Nt/2:-1:2));
        Rect(k,:)=real(ifft(Recf_now(k,:)));
    end
    if j<Layer_number
        T=(2*p(j)*c(j))/(p(j)*c(j)+p(j+1)*c(j+1)); %界面1的透射系数
    end
    Recf=Recf_now(Receiver_spacing(j),:);%介质1中频域接收信号矩阵，每一行代表一道2000采样点
    DATA=[DATA;Rect];
end

Noise=randn(length(DATA(:,1)),length(DATA(1,:)));
DATA=DATA+Noise_data*Noise;
NUMBER=length(DATA(:,1))-1;

iterative_number=500;
Noise_time1=0.3;
Noise_time2=0.5;
Noise_time3=0.8;
[loc1,loc_noise1,result1,Qk1,v_orignal1,v1]=IterativeMAP_v(dt,df,f,fp,Noise_time1,sigma,lamda,lamda2,GK,I,DATA,iterative_number);
[loc2,loc_noise2,result2,Qk2,v_orignal2,v2]=IterativeMAP_v(dt,df,f,fp,Noise_time2,sigma,lamda,lamda2,GK,I,DATA,iterative_number);
[loc3,loc_noise3,result3,Qk3,v_orignal3,v3]=IterativeMAP_v(dt,df,f,fp,Noise_time3,sigma,lamda,lamda2,GK,I,DATA,iterative_number);



Q_Given=zeros(NUMBER,1);
Sum_Receiver_spacing=1;
for i=1:Layer_number
   Q_Given(Sum_Receiver_spacing:Sum_Receiver_spacing+Receiver_spacing(i)-1)=Q_define(i);
   Sum_Receiver_spacing=Sum_Receiver_spacing+Receiver_spacing(i);
end

h=0:dh:(length(DATA(:,1))-1)*dh;

%%
v_ture1=[c(1),c(1),c(1)];
v_ture2=[c(2),c(2),c(2),c(2)];
v_ture3=[c(3),c(3),c(3),c(3),c(3)];
v_ture4=[c(4),c(4),c(4),c(4),c(4),c(4)];
v_ture5=[c(5),c(5),c(5),c(5),c(5),c(5),c(5)];
v_ture6=[c(6),c(6),c(6),c(6),c(6),c(6),c(6),c(6)];
v_ture=[v_ture1,v_ture2,v_ture3,v_ture4,v_ture5,v_ture6];
%%
B_T=DATA(:,1:300);
h=0:dh:(length(B_T(:,1))-1)*dh;
T=0:1:300-1;
F1=figure;
set(F1,'position',[20 50 600 280]);
subplot(1,2,1);
wigb_my(-B_T',1,h,T,1);
xlabel('Time(ms)','FontName','Times New Roman','FontSize',20);
% set(gca,'xaxislocation','top');
ylabel('Depth(m)','FontName','Times New Roman','FontSize',20);
set(gca,'FontName','Times New Roman','FontSize',12,'linewidth',1);
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.09 0.25 0.4 0.7]);
subplot(1,2,2);
hold on;
box on;
stem(h,loc1-loc_noise1,':bo','linewidth',1);
stem(h,loc2-loc_noise2,':kv','linewidth',1);
stem(h,loc3-loc_noise3,':gs','linewidth',1);
xlabel('Depth(m)','FontName','Times New Roman','FontSize',20);
ylabel('Time disturbance(ms)','FontName','Times New Roman','FontSize',20);
set(gca,'FontName','Times New Roman','FontSize',12,'linewidth',1);
set(gca,'TickLength',[0 0.001]);
set(gca,'position',[0.575 0.25 0.4 0.7]);
%%
h=h(2:end)';
F2=figure;
set(F2,'position',[20 50 600 800]);
subplot(3,2,1);
hold on;
grid on;
box on;
plot(h,Q_Given,'*k','linewidth',1);
plot(h,result1,':bo','linewidth',1);
plot(h,Qk1,'--r','linewidth',1.5);
xlabel('Depth(m)','FontName','Times New Roman','FontSize',20);
ylabel('Q value','FontName','Times New Roman','FontSize',20);
legend('Real Q','SRM','IMAP','Location','NorthWest');
set(gca,'FontName','Times New Roman','FontSize',12,'linewidth',1);
set(gca,'TickLength',[0 0.001]);
axis([0 660,0 200]);
set(gca,'position',[0.09 0.76 0.4 0.23]);
subplot(3,2,2);
hold on;
grid on;
box on;
plot(h,v_ture./1000,'*k','linewidth',1);
plot(h,v_orignal1./1000,':bo','linewidth',1);
plot(h,v1./1000,'--r','linewidth',1.5);
legend('Real velocity','Uncorrected velocity','Corrected velocity ','Location','NorthWest');
xlabel('Depth(m)','FontName','Times New Roman','FontSize',20);
ylabel('Velocity(km/s)','FontName','Times New Roman','FontSize',20);
set(gca,'FontName','Times New Roman','FontSize',12,'linewidth',1);
set(gca,'TickLength',[0 0.001]);
axis([0 660,2 8]);
set(gca,'position',[0.56 0.76 0.4 0.23]);
subplot(3,2,3);
hold on;
grid on;
box on;
plot(h,Q_Given,'*k','linewidth',1);
plot(h,result2,':kv','linewidth',1);
plot(h,Qk2,'--r','linewidth',1.5);
xlabel('Depth(m)','FontName','Times New Roman','FontSize',20);
ylabel('Q value','FontName','Times New Roman','FontSize',20);
legend('Real Q','SRM','IMAP','Location','NorthWest');
set(gca,'FontName','Times New Roman','FontSize',12,'linewidth',1);
set(gca,'TickLength',[0 0.001]);
axis([0 660,0 200]);
set(gca,'position',[0.09 0.43 0.4 0.23]);

subplot(3,2,4);
hold on;
grid on;
box on;
plot(h,v_ture./1000,'*k','linewidth',1);
plot(h,v_orignal2./1000,':kv','linewidth',1);
plot(h,v2./1000,'--r','linewidth',1.5);
legend('Real velocity','Uncorrected velocity','Corrected velocity ','Location','NorthWest');
xlabel('Depth(m)','FontName','Times New Roman','FontSize',20);
ylabel('Velocity(km/s)','FontName','Times New Roman','FontSize',20);
set(gca,'FontName','Times New Roman','FontSize',12,'linewidth',1);
set(gca,'TickLength',[0 0.001]);
axis([0 660,2 8]);
set(gca,'position',[0.56 0.43 0.4 0.23]);
subplot(3,2,5);
hold on;
grid on;
box on;
plot(h,Q_Given,'*k','linewidth',1);
plot(h,result3,':gs','linewidth',1);
plot(h,Qk3,'--r','linewidth',1.5);
xlabel('Depth(m)','FontName','Times New Roman','FontSize',20);
ylabel('Q value','FontName','Times New Roman','FontSize',20);
legend('Real Q','SRM','IMAP','Location','NorthWest');
set(gca,'FontName','Times New Roman','FontSize',12,'linewidth',1);
set(gca,'TickLength',[0 0.001]);
axis([0 660,0 200]);
set(gca,'position',[0.09 0.1 0.4 0.23]);
subplot(3,2,6);
hold on;
grid on;
box on;
plot(h,v_ture./1000,'*k','linewidth',1);
plot(h,v_orignal3./1000,':gs','linewidth',1);
plot(h,v3./1000,'--r','linewidth',1.5);
legend('Real velocity','Uncorrected velocity','Corrected velocity ','Location','NorthWest');
xlabel('Depth(m)','FontName','Times New Roman','FontSize',20);
ylabel('Velocity(km/s)','FontName','Times New Roman','FontSize',20);
set(gca,'FontName','Times New Roman','FontSize',12,'linewidth',1);
set(gca,'TickLength',[0 0.001]);
axis([0 660,2 8]);
set(gca,'position',[0.56 0.1 0.4 0.23]);
