function [ls_t,run_time,ls_iter] = LS_total(R,H,DR,DH,num_r,h0)

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
D_R=blkdiag(D_R1,D_R2,D_R3,D_R4);
% % P=inv(D_R);
%%%% LS total estimation %%%%%
   
    i=1;
    tstart=tic;
    a0=1;b0=2;iter1=0;
 while 1
     iter1=iter1+1;
     
     %%% the first period data
  for j=1:num1
     L1(j,1)=R1(j,i)-a0*(H1(j,i)-h0)^b0;
     A1(j,:)=[(H1(j,i)-h0)^b0 a0*(H1(j,i)-h0)^b0*log((H1(j,i)-h0))];
  end
     %%% the second period data
  for j=1:num2
     L2(j,1)=R2(j,i)-a0*(H2(j,i)-h0)^b0;
     A2(j,:)=[(H2(j,i)-h0)^b0 a0*(H2(j,i)-h0)^b0*log((H2(j,i)-h0))];
  end
 %%% the third period data
   for j=1:num3
     L3(j,1)=R3(j,i)-a0*(H3(j,i)-h0)^b0;
     A3(j,:)=[(H3(j,i)-h0)^b0 a0*(H3(j,i)-h0)^b0*log((H3(j,i)-h0))];
   end
   %%% the fouth period data
   for j=1:num4
     L4(j,1)=R4(j,i)-a0*(H4(j,i)-h0)^b0;
     A4(j,:)=[(H4(j,i)-h0)^b0 a0*(H4(j,i)-h0)^b0*log((H4(j,i)-h0))];
  end
   L=[L1;L2;L3;L4];
   A=[A1;A2;A3;A4];
   
   
    ls=inv(A'*inv(D_R)*A)*(A'*inv(D_R)*L);
    a0=a0+ls(1,1);b0=b0+ls(2,1);
    if norm(ls)<10^-10 
     ls_a(i,1)=a0;ls_b(i,1)=b0; ls_iter(i,1)=iter1;
        break     
    end    
 end

 ls_t=[ls_a ls_b];
 toc;
 telapsed=toc(tstart);
 run_time(i,1)=telapsed;


% % % % % 
% % % % % d_ls_a=ls_a-a*ones(num_s,1);
% % % % % 
% % % % % 
% % % % % d_ls_b=ls_b-b*ones(num_s,1);
% % % % % 
% % % % % 
% % % % % rms_ls_a=sqrt(sum((ls_a-a*ones(num_s,1)).^2,1)/num_s);
% % % % % rms_ls_b=sqrt(sum((ls_b-b*ones(num_s,1)).^2,1)/num_s);
% % % % % rms_ab=[rms_ls_a;rms_ls_b];
% % % % % 
% % % % % rms_ls_total=rms_ab;
% % % % % time_ls_total=run_time;
end

