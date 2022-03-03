clc;
clear;
load R1.mat;load R2.mat;load R3.mat;load R4.mat;
load H1.mat;load H2.mat;load H3.mat;load H4.mat;
load D_R1.mat;load D_R2.mat;load D_R3.mat;load D_R4.mat
load D_h1.mat;load D_h2.mat;load D_h3.mat;load D_h4.mat;
load R.mat;load h.mat;load a.mat;load b.mat;load h0.mat;

% % % h0=26;
tic;
[row1,com1]=size(R1);
[row2,com2]=size(R2);
[row3,com3]=size(R3);
[row4,com4]=size(R4);
kX0=[2,3];
DR=blkdiag(D_R1,D_R2,D_R3,D_R4);
DH=blkdiag(D_h1,D_h2,D_h3,D_h4);
QR=DR./0.1^2;
QH=DH./0.1^2;
num_r=[row1;row2;row3;row4];
com1=10;
for i=1:com1
    i
R=[R1(:,i);R2(:,i);R3(:,i);R4(:,i)];
H=[H1(:,i);H2(:,i);H3(:,i);H4(:,i)];

[pe1,t1,iter1]=LS_total(R,H,DR,DH,num_r,h0);
PE1(i,:)=pe1;runtime1(i,1)=t1;
Iter1(i,1)=iter1;
[pe2,t2,iter2]=LS_sequential(R,H,DR,DH,num_r,h0);
PE2(i,:)=pe2;runtime2(i,1)=t2;
Iter2(i,1)=iter2;
[pe3,t3,iter3]=TLS_total(R,H,DR,DH,num_r,h0);
PE3(i,:)=pe3;runtime3(i,1)=t3;
Iter3(i,1)=iter3;
[pe4,t4,iter4]=TLS_sequential(R,H,DR,DH,num_r,h0);
PE4(i,:)=pe4;runtime4(i,1)=t4;
Iter4(i,1)=iter4;

end

[rms1,da1,db1] = accuracy_calculation(PE1(:,1),PE1(:,2),a,b);
[rms2,da2,db2] = accuracy_calculation(PE2(:,1),PE2(:,2),a,b);
[rms3,da3,db3] = accuracy_calculation(PE3(:,1),PE3(:,2),a,b);
[rms4,da4,db4] = accuracy_calculation(PE4(:,1),PE4(:,2),a,b);

dab=[da1 db1 da2 db2 da3 db3 da4 db4];
x=1:1:com1;
plot(x,runtime1,'r',x,runtime2,'g',x,runtime3,'b',x,runtime4,'k');
legend('TSGM','SSGM','TSGH','SSGH');
xlabel('拟合次数');
ylabel('运行时间/s');
rms=[rms1;rms2;rms3;rms4;]
Max_ab=max(dab);
Ave_ab=mean(dab);
Min_ab=min(dab);
MAN_ab=[Max_ab;Ave_ab;Min_ab];
Max_t=[max(runtime1) max(runtime2) max(runtime3) max(runtime4)];
Ave_t=[mean(runtime1) mean(runtime2) mean(runtime3) mean(runtime4)];
Min_t=[min(runtime1) min(runtime2) min(runtime3) min(runtime4)];
MAN_t=[Max_t Ave_t Min_t];
