% Tester

tic;
for n = 1:1e3
    z = afsH_varHsize_Addy_den(PI_single,Rg_H(i1,i2),Rh_H(i1,i2),0,1,1,1,'n-1');
end
toc;