function [loc,loc_noise,Q_SRM,Qk,v_orignal,v]=IterativeMAP_v(dt,df,f,fp,Noise_time,sigma,lamda,lamda2,GK,I,DATA,iterative_number)

%% 谱比法估计Q值
NUMBER=length(DATA(:,1))-1;
% frequency band parameter slection
Begf=3;
Endf=3*fp;
BegNum=floor(Begf/df)+1; %floor向下取整
EndNum=floor(Endf/df)+1;
Bandf=((BegNum-1):(EndNum-1))*df;
% 利用相邻两道计算
nn=length(Bandf);
big_G=zeros(NUMBER*nn,NUMBER*2);
for i=1:NUMBER
    big_G((i-1)*nn+1:i*nn,2*(i-1)+1)=ones(nn,1);
    big_G((i-1)*nn+1:i*nn,2*(i-1)+2)=f(BegNum:EndNum);
end
D=[];
Q_SRM=zeros(NUMBER,1);
TraT=[];
mid=abs(hilbert(DATA(1,:)));
mid2=max(mid);
mid3=find(mid==mid2);
loct=mid3(1);%初至波峰值到达时
loc=[loct];
loct=round(loct+Noise_time*randn(1));
loc_noise=[loct];
M=zeros(size(Q_SRM));
for first=1:1:NUMBER
    st=DATA(first,:);
    AbsSrcf=abs(fft(st));
    sr=DATA(first+1,:);
%     [peakr,locr]=max(sr);
%     loc= [loc;locr];
    mid=abs(hilbert(sr));
    mid2=max(mid);
    mid3=find(mid==mid2);
    locr=mid3(1);%初至波峰值到达时
    loc= [loc;locr];
    
    locr= round(locr+Noise_time*randn(1));
    TravelTime=(locr-loc_noise(end))*dt; %两个波的时间差
    AbsRecf=abs(fft(sr));
    BandAbsSrcf=AbsSrcf(BegNum:EndNum);
    BandAbsRecf=AbsRecf(BegNum:EndNum);
    %%  Log spectral ratio method
    Specratio=log(BandAbsRecf./BandAbsSrcf);
     % 法方程计算最小二乘解(与上面解一致)  上面采用polyfit自带函数 进行最小二乘拟合 下面利用最小二乘解 结果一致
    y=Specratio;
    d=y';%观测数据
    ff=f(BegNum:EndNum);
    L=length(y);
    G=ones(L,2);
    G(:,2)=ff';
    E=eye(length(G(1,:)));
    m=inv(G'*G+0*E)*G'*d;
    M(first)=m(2);
    Q=-pi*TravelTime/m(2);
    Q_SRM(first)=Q;
    TraT=[TraT;TravelTime];
    loc_noise= [loc_noise;locr];
    D=[D;d];
end
noise=zeros(size(loc_noise));
%%
TraT_orignal=TraT;
Qk=Q_SRM;
for j=1:iterative_number
    A=zeros(NUMBER*2,NUMBER*2);
    b=zeros(NUMBER*2,1);
     for i=1:NUMBER
        if i>=2&&i<=NUMBER-1
            hq1=(abs(Qk(i+1)-Qk(i)))^GK+sigma;
            hq2=(abs(Qk(i)-Qk(i-1)))^GK+sigma;
            hq1add=hq2/(hq1+hq2);
            hq2add=hq1/(hq2+hq1);
            hqi=I;
            A(2*i,2*i)=(hq2add*Qk(i-1)+hqi*Qk(i)+hq1add*Qk(i+1))/(hqi+hq1add+hq2add); 
        elseif  i==1
            hqi=I;
            A(2*i,2*i)=(hqi*Qk(i)+Qk(i+1))/(hqi+1);
        else 
            hqi=I;
            A(2*i,2*i)=(hqi*Qk(i)+Qk(i-1))/(hqi+1); 
        end
            b(2*i)=pi*TraT(i);
     end

    E=eye(length(big_G(1,:)));
    mk=inv((big_G'*big_G)+lamda*(A'*A)+lamda2*E)*(big_G'*D-lamda*A'*b);
    Qk_mid=zeros(size(Qk));
    for i=1:NUMBER  
    Qk_mid(i)=-pi*TraT(i)/mk(2*i); 
    TraT(i)=-M(i).*Qk_mid(i)/pi;
    end
    Qk=Qk_mid;
end

v=zeros(size(TraT));
v_orignal=zeros(size(TraT));
for i=1:NUMBER
    v(i)=20./TraT(i);
    v_orignal(i)=20/TraT_orignal(i);
end
    





