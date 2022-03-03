%%%%% generate the simulated data of Hydrological curve fitting
clc;
clear;
h0=25;a=0.012;b=1.55;
% % % num_t=300;
% % % num_p=randperm(num_t,4);
% % % num1=num_p(1,1);num2=num_p(1,2);num3=num_p(1,3);num4=num_p(1,4);

num=796;
temp_num=round((rand(1,3)*num/8-num/16)+num/4);

% % % for i=1:3
% % %     if temp_num(1,i)<2
% % %         temp_num(1,i)=round((rand(1,1)*num/3-num/6)+num/2);      
% % %     end
% % % end

 num1=temp_num(1,1);num2=temp_num(1,2);num3=temp_num(1,3);num4=num-num1-num2-num3;
% % num1=200;num2=189;num3=210;num4=197;

% % num11=120;num12=num1-num11;
num11=round(rand(1,1)*num1/6+num1/2);num12=num1-num11;

h11=rand(num11,1)*5+28;
h12=rand(num12,1)*5+33;
h1=[h11;h12];
for i=1:num1
   r1(i,1)=a*(h1(i,1)-h0)^b; 
end

% % % num21=107;num22=num2-num21;
num21=round(rand(1,1)*num2/6+num2/2);num22=num2-num21;
h21=rand(num21,1)*5+28;
h22=rand(num22,1)*5+33;
h2=[h21;h22];
for i=1:num2
   r2(i,1)=a*(h2(i,1)-h0)^b; 
end

% % num31=123;num32=num3-num31;
num31=round(rand(1,1)*num3/6+num3/2);num32=num3-num31;
h31=rand(num31,1)*5+28;
h32=rand(num32,1)*5+33;
h3=[h31;h32];
for i=1:num3
   r3(i,1)=a*(h3(i,1)-h0)^b; 
end


% % num41=116;num42=num4-num41;
num41=round(rand(1,1)*num4/6+num4/2);num42=num4-num41;
h41=rand(num41,1)*5+28;
h42=rand(num42,1)*5+33;
h4=[h41;h42];
% % % % h4=rand(num4,1)*20+26;
for i=1:num4
   r4(i,1)=a*(h4(i,1)-h0)^b; 
end
R=[r1;r2;r3;r4];
h=[h1;h2;h3;h4];

%%% generate the std of R and h with 1th to 4th periods data%%%%
% % std_R1=unifrnd(0.04,0.05,num1,1);
% % % std_h1=unifrnd(0.08,0.09,num1,1);
% % % 
% % std_R2=unifrnd(0.03,0.04,num2,1);
% % % std_h2=unifrnd(0.05,0.07,num2,1);
% % % 
% % % std_R3=unifrnd(0.02,0.03,num3,1);
% % % std_h3=unifrnd(0.02,0.04,num3,1);
% % % 
% % std_R4=unifrnd(0.01,0.02,num4,1);
% % % std_h4=unifrnd(0.009,0.02,num4,1);
std_R1=0.1*r1;
% % std_h1=unifrnd(0.02,0.05,num1,1);
% std_R1=rand(num1,1)*0.001+0.01;
std_h1=rand(num1,1)*0.75+0.05;

std_R2=0.1*r2;
% std_h2=unifrnd(0.02,0.05,num2,1);
% std_R2=rand(num2,1)*0.001+0.02;
std_h2=rand(num2,1)*0.75+0.05;

std_R3=0.1*r3;
% % std_h3=unifrnd(0.02,0.05,num3,1);
% std_R3=rand(num3,1)*0.001+0.03;
std_h3=rand(num3,1)*0.75+0.05;

std_R4=0.1*r4;
% % std_h4=unifrnd(0.02,0.05,num4,1);
% std_R4=rand(num4,1)*0.001+0.04;
std_h4=rand(num4,1)*0.75+0.05;
%%% generate the VCO of R and h %%%%
D_R1=diag(std_R1.^2);
D_h1=diag(std_h1.^2);

D_R2=diag(std_R2.^2);
D_h2=diag(std_h2.^2);

D_R3=diag(std_R3.^2);
D_h3=diag(std_h3.^2);

D_R4=diag(std_R4.^2);
D_h4=diag(std_h4.^2);

% % % chol_R1=chol(D_R1);
% % % chol_h1=chol(D_h1);
% % % 
% % % chol_R2=chol(D_R2);
% % % chol_h2=chol(D_h2);
% % % 
% % % chol_R3=chol(D_R3);
% % % chol_h3=chol(D_h3);
% % % 
% % % chol_R4=chol(D_R4);
% % % chol_h4=chol(D_h4);

num_s=1000;
for i=1:num_s
    
    for j=1:num1
      e_R1(1,j)=normrnd(0,std_R1(j,1),1,1);
      e_h1(1,j)=normrnd(0,std_h1(j,1),1,1);  
    end
    
    for j=1:num2
      e_R2(1,j)=normrnd(0,std_R2(j,1),1,1);
      e_h2(1,j)=normrnd(0,std_h2(j,1),1,1);  
    end
    for j=1:num3
      e_R3(1,j)=normrnd(0,std_R3(j,1),1,1);
      e_h3(1,j)=normrnd(0,std_h3(j,1),1,1);  
    end
    for j=1:num4
      e_R4(1,j)=normrnd(0,std_R4(j,1),1,1);
      e_h4(1,j)=normrnd(0,std_h4(j,1),1,1);  
    end

    
% % %     
% % %     e_R1=normrnd(0,1,1,num1)*chol_R1;
% % %     e_h1=normrnd(0,1,1,num1)*chol_h1;
% % %     
% % %     e_R2=normrnd(0,1,1,num2)*chol_R2;
% % %     e_h2=normrnd(0,1,1,num2)*chol_h2;
% % %     
% % %     e_R3=normrnd(0,1,1,num3)*chol_R3;
% % %     e_h3=normrnd(0,1,1,num3)*chol_h3;
% % %     
% % %     e_R4=normrnd(0,1,1,num4)*chol_R4;
% % %     e_h4=normrnd(0,1,1,num4)*chol_h4;
    
    R1(1:num1,i)=r1+e_R1';
    H1(1:num1,i)=h1+e_h1';
    
    R2(1:num2,i)=r2+e_R2';
    H2(1:num2,i)=h2+e_h2';
    
    R3(1:num3,i)=r3+e_R3';
    H3(1:num3,i)=h3+e_h3';
    
    R4(1:num4,i)=r4+e_R4';
    H4(1:num4,i)=h4+e_h4';    
    
    
end

plot(H1(:,1),R1(:,1),'r*');
hold on
plot(H2(:,1),R2(:,1),'g*');
hold on
plot(H3(:,1),R3(:,1),'b*');
hold on
plot(H4(:,1),R4(:,1),'k*');
legend('the 1th period of data','the 2th period of data','the 3th period of data','the 4th period of data');
xlabel('水位h(m)')
ylabel('水位流量变换Q(m3/s)')
% % % plot_point=[H1(:,1),R1(:,1),H2(:,1),R2(:,1),H3(:,1),R3(:,1),H4(:,1),R4(:,1)];
save R1.mat R1;
save R2.mat R2;
save R3.mat R3;
save R4.mat R4;

save H1.mat H1;
save H2.mat H2;
save H3.mat H3;
save H4.mat H4;

save D_R1.mat D_R1;
save D_R2.mat D_R2;
save D_R3.mat D_R3;
save D_R4.mat D_R4;

save D_h1.mat D_h1;
save D_h2.mat D_h2;
save D_h3.mat D_h3;
save D_h4.mat D_h4;

save h.mat h;
save R.mat R;
save a.mat a;
save b.mat b;
save h0.mat h0;

