function [ls_s,run_time,ls_iter] = LS_sequential(R,H,DR,DH,num_r,h0)


num1=num_r(1,1);num2=num_r(2,1);num3=num_r(3,1);num4=num_r(4,1);

temp_R=mat2cell(R,[num1 num2 num3 num4],[1]);
R1=temp_R{1,1};R2=temp_R{2,1};R3=temp_R{3,1};R4=temp_R{4,1};

temp_H=mat2cell(H,[num1 num2 num3 num4],[1]);
H1=temp_H{1,1};H2=temp_H{2,1};H3=temp_H{3,1};H4=temp_H{4,1};
DR=DR./0.1^2;
temp_DR=mat2cell(DR,[num1 num2 num3 num4],[num1 num2 num3 num4]);
D_R1=temp_DR{1,1};D_R2=temp_DR{2,2};D_R3=temp_DR{3,3};D_R4=temp_DR{4,4};
DH=DH./0.1^2;
temp_DH=mat2cell(DH,[num1 num2 num3 num4],[num1 num2 num3 num4]);
D_h1=temp_DH{1,1};D_h2=temp_DH{2,2};D_h3=temp_DH{3,3};D_h4=temp_DH{4,4};




i=1; %%%

    tstart=tic;
    a0=1;b0=2;iter1=0;
 while 1
     iter1=iter1+1;
     
     %%% the first period data
  for j=1:num1
     L1(j,1)=R1(j,i)-a0*(H1(j,i)-h0)^b0;
     A1(j,:)=[(H1(j,i)-h0)^b0 a0*(H1(j,i)-h0)^b0*log((H1(j,i)-h0))];
  end
     ls1=inv(A1'*inv(D_R1)*A1)*(A1'*inv(D_R1)*L1);
     Q1=inv(A1'*inv(D_R1)*A1);
% %      a0=a0+ls1(1,1);b0=b0+ls1(2,1);
     %%% the second period data
  for j=1:num2
     L2(j,1)=R2(j,i)-a0*(H2(j,i)-h0)^b0;
     A2(j,:)=[(H2(j,i)-h0)^b0 a0*(H2(j,i)-h0)^b0*log((H2(j,i)-h0))];
  end
   J1=Q1*A2'*inv(D_R2+A2*Q1*A2');
   dL1=L2-A2*ls1;
   ls2=ls1+J1*dL1;
   Q2=inv(inv(Q1)+A2'*inv(D_R2)*A2);
% % %    a0=a0+ls2(1,1);b0=b0+ls2(2,1);
 %%% the third period data
   for j=1:num3
     L3(j,1)=R3(j,i)-a0*(H3(j,i)-h0)^b0;
     A3(j,:)=[(H3(j,i)-h0)^b0 a0*(H3(j,i)-h0)^b0*log((H3(j,i)-h0))];
   end
   J2=Q2*A3'*inv(D_R3+A3*Q2*A3');
   dL2=L3-A3*ls2;
   ls3=ls2+J2*dL2;
   Q3=inv(inv(Q2)+A3'*inv(D_R3)*A3);
% %    a0=a0+ls3(1,1);b0=b0+ls3(2,1);
   
   %%% the fouth period data
   for j=1:num4
     L4(j,1)=R4(j,i)-a0*(H4(j,i)-h0)^b0;
     A4(j,:)=[(H4(j,i)-h0)^b0 a0*(H4(j,i)-h0)^b0*log((H4(j,i)-h0))];
   end
   
   J4=Q3*A4'*inv(D_R4+A4*Q3*A4');
   dL3=L4-A4*ls3;
   ls4=ls3+J4*dL3;
   Q4=inv(inv(Q3)+A4'*inv(D_R4)*A4);
   a0=a0+ls4(1,1);b0=b0+ls4(2,1);
   if norm(ls4)<10^-10
     ls_a(i,1)=a0;ls_b(i,1)=b0;
       ls_iter(i,1)=iter1;
        break     
    end    
 end
 ls_s=[ls_a ls_b];
 toc;
 telapsed=toc(tstart);
 run_time(i,1)=telapsed;


end

