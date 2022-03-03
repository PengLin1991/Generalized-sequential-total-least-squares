function [tls_s,run_time,tls_iter] = TLS_sequential(R,H,DR,DH,num_r,h0)


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





%%%% single TLS estimation for the 1th observation %%%

    i=1;
    tstart=tic;
    a0=1;b0=2;iter2=0;
    e_h1=zeros(num1,1);e_h2=zeros(num2,1);e_h3=zeros(num3,1);e_h4=zeros(num4,1);
 while 1
     iter2=iter2+1;
  for j=1:num1
     
     A1(j,:)=[(H1(j,i)-e_h1(j,1)-h0)^b0 a0*(H1(j,i)-e_h1(j,1)-h0)^b0*log((H1(j,i)-e_h1(j,1)-h0))];
     B1(j,j)=-a0*b0*(H1(j,i)-e_h1(j,1)-h0)^(b0-1);
     L1(j,1)=R1(j,i)-a0*(H1(j,i)-e_h1(j,1)-h0)^b0+B1(j,j)*e_h1(j,1);
     
  end
     G1=[eye(num1) B1];
     D1=blkdiag(D_R1,D_h1);
     Qg1=G1*D1*G1';
     tls1=inv(A1'*inv(Qg1)*A1)*(A1'*inv(Qg1)*L1);
     Q1=inv(A1'*inv(Qg1)*A1);

     
  
  
  for j=1:num2
     
     A2(j,:)=[(H2(j,i)-e_h2(j,1)-h0)^b0 a0*(H2(j,i)-e_h2(j,1)-h0)^b0*log((H2(j,i)-e_h2(j,1)-h0))];
     B2(j,j)=-a0*b0*(H2(j,i)-e_h2(j,1)-h0)^(b0-1);
     L2(j,1)=R2(j,i)-a0*(H2(j,i)-e_h2(j,1)-h0)^b0+B2(j,j)*e_h2(j,1);
     
  end
     G2=[eye(num2) B2];
     D2=blkdiag(D_R2,D_h2);
     Qg2=G2*D2*G2';
     
     J1=Q1*A2'*inv(Qg2+A2*Q1*A2');
     dL1=L2-A2*tls1;
     tls2=tls1+J1*dL1;
     Q2=inv(inv(Q1)+A2'*inv(Qg2)*A2);

     
  for j=1:num3
     
     A3(j,:)=[(H3(j,i)-e_h3(j,1)-h0)^b0 a0*(H3(j,i)-e_h3(j,1)-h0)^b0*log((H3(j,i)-e_h3(j,1)-h0))];
     B3(j,j)=-a0*b0*(H3(j,i)-e_h3(j,1)-h0)^(b0-1);
     L3(j,1)=R3(j,i)-a0*(H3(j,i)-e_h3(j,1)-h0)^b0+B3(j,j)*e_h3(j,1);
     
  end
  
     G3=[eye(num3) B3];
     D3=blkdiag(D_R3,D_h3);
     Qg3=G3*D3*G3';
     
     J2=Q2*A3'*inv(Qg3+A3*Q2*A3');
     dL2=L3-A3*tls2;
     tls3=tls2+J2*dL2;
     Q3=inv(inv(Q2)+A3'*inv(Qg3)*A3);

  
  for j=1:num4
     
     A4(j,:)=[(H4(j,i)-e_h4(j,1)-h0)^b0 a0*(H4(j,i)-e_h4(j,1)-h0)^b0*log((H4(j,i)-e_h4(j,1)-h0))];
     B4(j,j)=-a0*b0*(H4(j,i)-e_h4(j,1)-h0)^(b0-1);
     L4(j,1)=R4(j,i)-a0*(H4(j,i)-e_h4(j,1)-h0)^b0+B4(j,j)*e_h4(j,1);
     
  end
  
     G4=[eye(num4) B4];
     D4=blkdiag(D_R4,D_h4);
     Qg4=G4*D4*G4';
     
     J3=Q3*A4'*inv(Qg4+A4*Q3*A4');
     dL3=L4-A4*tls3;
     tls4=tls3+J3*dL3;
     Q4=inv(inv(Q3)+A4'*inv(Qg4)*A4);
     
     e_rh1=D1*G1'*inv(Qg1)*(L1-A1*tls4);
     e_h1=e_rh1(num1+1:num1*2,1);
     
     e_rh2=D2*G2'*inv(Qg2)*(L2-A2*tls4);
     e_h2=e_rh2(num2+1:num2*2,1);     
 
     e_rh3=D3*G3'*inv(Qg3)*(L3-A3*tls4);
     e_h3=e_rh3(num3+1:num3*2,1);
     
     e_rh4=D4*G4'*inv(Qg4)*(L4-A4*tls4);
     e_h4=e_rh4(num4+1:num4*2,1); 
     
    
    a0=a0+tls4(1,1);b0=b0+tls4(2,1);
    if norm(tls4)<10^-10 
     tls_a(i,1)=a0;tls_b(i,1)=b0; tls_iter(i,1)=iter2;
        break     
    end    
 end
 tls_s=[tls_a tls_b];
 toc;
 telapsed=toc(tstart);
 run_time(i,1)=telapsed;


end

