params.NumPosts = 3;
params.distOP = 50;
params.latPost=60.065843;
params.lonPost=30.246163;
params.lonstr1KA1='1 43074U 17083E   20352.49900893  .00000110  00000-0  32273-4 0  9994';  
params.lonstr2KA1='2 43074  86.3979  81.0759 0001886 103.9895 256.1511 14.34217107156373';
params.distOA=350;  
params.tolx = 1e-4;  
params.tolf=1e-3;  
params.maxfun = 200*3;  
params.maxiter= 200*3;  
params.DT_=abs(30*1);  
params.DT_ALL=abs(1);  
params.dXYZ=[1 1 1] ; 
params.freqDownLink=7*10^9;  
params.RMS_DOPPLER=1;
params.OKSType=0;
params.min_elev_deg=10;
tabl_Znach = zeros (2,3);
k=1;
for distOP = 50:50:150;
    params.distOP = distOP;
    disp(distOP )
results=loshadka(params);
tabl_Znach(1,k) = distOP;
tabl_Znach(2,k) = results.SrZnachOKS;
disp(['d OKS = ' num2str(results.SrZnachOKS)  ' км ' ]) %
disp(['distOP = ' num2str(distOP)  ' км ' ])
k=k+1;
end
disp (tabl_Znach)

