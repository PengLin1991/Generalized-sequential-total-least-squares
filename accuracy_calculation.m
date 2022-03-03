function [rms,da,db] = accuracy_calculation(arg1,arg2,a,b)
num_s=size(arg1,1);

da=arg1-a*ones(num_s,1);


db=arg2-b*ones(num_s,1);


rms_a=sqrt(sum(da.^2,1)/num_s);
rms_b=sqrt(sum(db.^2,1)/num_s);
rms=[rms_a rms_b];

end

