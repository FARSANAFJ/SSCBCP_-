clear all;
close all;
clc;
%%%%%%%%%%%%%%%%%%%%%speech signal compression
Input_speech=audioread('female600.wav','native');     %Input speech file for encryption 
d=Input_speech;
d1=d(4097:8192); 
subplot(2,3,1);
plot(d1);
d2_=fwht(double(d1));
d2=d2_(1:4096);%%%%%%% compression 

S=reshape(d2,64,64);
subplot(2,3,2);
plot(S);

A=S;
[M,N]=size(A);                  
data=zeros(M,N); 
 
a=1;     
b=1; 
c=1; 
d=1; 
  
for loop=1:5096
%     for z=1:K 
        for x=1:M 
            for y=1:N 
                x1=mod((1-a*x^2+b*y),M);%???? 
                if x1==0 
                    x1=M; 
                end 
         
                y1=mod((c*x+d*x1^2),N);%???? 
                if y1==0 
                    y1=N; 
                end 
         
                data(x1,y1)=A(x,y);%????? 
                %data(x1,y1,z)=A(x,y,z);  
            end 
        end 
    A=data; 
end 
 
subplot(2,3,3); 
d3=reshape(A,4096,1);
plot(d3);
title('permuted signal'); 
%%
%substitution 
%Second phase (substitution)
z1(1)=0.1;  %x(1)
z2(1)=0.1;  %y(1)
e=0.3;      %e=0.3
v=0.2;
r_=5;
m=(1-exp(-r_))/r_;
omega=100;
k=9;
a_=1.885;
% keystream generation 
% Computing the values for x(i),y(i)
% We will compute the first 4096 values. 
for i=2:4096
    z1(i)=mod(z1(i-1)+omega/(2*pi)+(a_*omega)/(2*pi*r_)*(1-exp(-r_))*z2(i-1)+...
    (k/r_)*(1-exp(-r_))*cos(2*pi*z1(i-1)),1);
    % Since x(i) is computed mod1 we always have  0<=x(i)<1
    z2(i)=exp(-r_)*(z2(i-1)+e*cos(2*pi*z1(i-1)));
end

%%
%subtitution phase
x_key=int64(fix(z1.*10^3));
Input_speech=audioread('male600.wav','native');     %Input speech file for encryption  
data_sub=int64(d3);
d4=bitxor(int64(d3),x_key');%encryption
d4_d=ifwht(double(d4));
subplot(2,3,4)
plot(d4);
d5=bitxor(d4,x_key');       %decryption
plot(d4)
%% 
d4_=reshape(d5,64,64);
A2=zeros(M,N); 
A3=zeros(M,N); 
for loop=1:5096 
%     for z1=1:K 
        for x1=1:M 
            for y1=1:N 
                x2=mod(y1-d*x1^2,M);%???? 
                if x2==0 
                    x2=M; 
                end 
                y2=mod(x1+a*x2^2-1,N);%???? 
                if y2==0 
                    y2=N; 
                end 
%                 z2=z1; 
%                 A2(x2,y2)=mod(A1(x1,y1,z1)-x2^2-y2^2,256);%????? 
                A2(x2,y2)=d4_(x1,y1); 
            end 
        end 
%     end 
    d4_=A2; 
end 
 subplot(2,3,5); 
 d6=reshape(A2,4096,1);
 d7=ifwht(double(int64(d6)));
 plot(d7);
sound(d7)
    
%%%%%%% analysis-------------------------
%% crosscorrelaton
Input=audioread('Male 1.wav');
ff=Input;
d1_d=ff(4096:8191);
XX=corrcoef(double(d1),d4_d);
%%
coX=cov(double(d1),d4_d);
std_d1=std(double(d1));
std_d4_d=std(d4_d);
std_deviation=std_d1.*std_d4_d;
xxx=coX/std_deviation;
%%
%SNR
S_p=sum(double(d1).^2);
n_p=sum((d4_d-double(d1)).^2);
ratio=S_p./n_p;
SNR=10*log(ratio);
%%
save('scatter','d7','d1');
%histogram
hist(double(d1),100)
d23=randn(4096:1)
 hist(double(d4),100);
 grid on