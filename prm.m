params.RMS_DOPPLER=1; %% заменил
params.NumPosts = 3; %% заменил
params.distOP = 350; %% заменил
params.dX = 0.5;  %% заменил
params.dY = 0.5;  %% заменил
params.dZ = -0.5;  %% заменил
params.dXYZ=[1 1 -1] ; %% заменил
params.latPost=60.065843;
params.lonPost=30.246163;
params.lonstr1KA1='1 43074U 17083E   20352.49900893  .00000110  00000-0  32273-4 0  9994';  
params.lonstr2KA1='2 43074  86.3979  81.0759 0001886 103.9895 256.1511 14.34217107156373';
params.distOA=1000;  
params.tolx = 1e-4;  
params.tolf=1e-3;  
params.maxfun = 200*3;  
params.maxiter= 200*3;  
params.DT_=abs(30*1);  
params.DT_ALL=abs(1);  
params.freqDownLink=7*10^9;  
params.OKSType=0;
params.min_elev_deg=10;
tabl_Znach_RMS = zeros (4,12);
tabl_Znach_Numposts = zeros (2,3);
tabl_Znach_DistOP = zeros (2,3);

tic
k=1;
for RMS_DOPPLER = 0.8:0.2:1; %цикл р��?�?чета OKS от SKO, при кол-во по�?тов 3шт, ра�?�?тро�?нии по�?тов от центра 350км
    params.RMS_DOPPLER = RMS_DOPPLER;
    
results=loshadka(params);
tabl_Znach_RMS(1,k) = RMS_DOPPLER;
tabl_Znach_RMS(2,k) = results.SrZnachOKS;
disp(['d OKS = ' num2str(results.SrZnachOKS)  ' км ' ]) 
disp(['RMS_DOPPLER = ' num2str(RMS_DOPPLER)  ' Гц ' ])
k=k+1;
end

dX = 1;   
dY = 1;   
dZ = -1;   
params.dX = dX;  
params.dY = dY;  
params.dZ = dZ;   


for RMS_DOPPLER = 0.8:0.2:1; %цикл ра�?�?чета OKS от SKO, при кол-во по�?тов 3шт, ра�?�?тро�?нии по�?тов от центра 350км
    params.RMS_DOPPLER = RMS_DOPPLER;
    
results=loshadka(params);
% tabl_Znach(1,k) = RMS_DOPPLER;
tabl_Znach(3,k) = results.SrZnachOKS;
disp(['d OKS = ' num2str(results.SrZnachOKS)  ' км ' ]) 
disp(['RMS_DOPPLER = ' num2str(RMS_DOPPLER)  ' Гц ' ])
% k=k+1;
end

dX = 2;   
dY = 2;   
dZ = -2;   
params.dX = dX; 
params.dY = dY;  
params.dZ = dZ;  


for RMS_DOPPLER = 0.8:0.2:1; %цикл ра�?�?чета OKS от SKO, при кол-во по�?тов 3шт, ра�?�?тро�?нии по�?тов от центра 350км
    params.RMS_DOPPLER = RMS_DOPPLER;
    
results=loshadka(params);
tabl_Znach(1,k) = RMS_DOPPLER;
tabl_Znach(4,k) = results.SrZnachOKS;
disp(['d OKS = ' num2str(results.SrZnachOKS)  ' км ' ]) 
disp(['RMS_DOPPLER = ' num2str(RMS_DOPPLER)  ' Гц ' ])
k=k+1;
end

% plot (tabl_Znach_RMS(1,:),tabl_Znach_RMS(2,:), tabl_Znach_RMS(3,:),tabl_Znach_RMS(4,:));
X = tabl_Znach_RMS(1,:);
Y1 = tabl_Znach_RMS(2,:);
Y2 = tabl_Znach_RMS(3,:);
Y3 = tabl_Znach_RMS(4,:);

plot (X,Y1,X,Y2,X,Y3);

% dX = 0.5;  
% dY = 0.5;  
% dZ = -0.5;  
% params.dX = dX;  
% params.dY = dY;  
% params.dZ = dZ;  

% RMS_DOPPLER = 1;
%  params.RMS_DOPPLER = RMS_DOPPLER;
% for NumPosts = 3:1:20; %цикл ра�?�?чета OKS от кол-ва по�?тов, SKO = 1гц, ра�?�?тро�?нии по�?тов от центра 350км
%    
%     params.NumPosts = NumPosts;
%     
% results=loshadka(params);
% tabl_Znach(1,k) = NumPosts;
% tabl_Znach(2,k) = results.SrZnachOKS;
% disp(['d OKS = ' num2str(results.SrZnachOKS)  ' км ' ]) %
% disp(['NumPosts = ' num2str(NumPosts)  ' Штук' ])
% k=k+1;
% end
% 
% dX = 1;  
% dY = 1;  
% dZ = -1;  
% params.dX = dX;  
% params.dY = dY;  
% params.dZ = dZ;  
% 
% 
% for NumPosts = 3:1:10; %цикл ра�?�?чета OKS от кол-ва по�?тов, SKO = 1гц, ра�?�?тро�?нии по�?тов от центра 350км
%    
%     params.NumPosts = NumPosts;
%     
% results=loshadka(params);
% tabl_Znach(1,k) = NumPosts;
% tabl_Znach(2,k) = results.SrZnachOKS;
% disp(['d OKS = ' num2str(results.SrZnachOKS)  ' км ' ]) %
% disp(['NumPosts = ' num2str(NumPosts)  ' Штук' ])
% k=k+1;
% end

% dX = 2;  
% dY = 2;  
% dZ = -2;  
% params.dX = dX;  
% params.dY = dY;  
% params.dZ = dZ;  
% 
% 
% for NumPosts = 3:1:20; %цикл ра�?�?чета OKS от кол-ва по�?тов, SKO = 1гц, ра�?�?тро�?нии по�?тов от центра 350км
%    
%     params.NumPosts = NumPosts;
%     
% results=loshadka(params);
% tabl_Znach(1,k) = NumPosts;
% tabl_Znach(2,k) = results.SrZnachOKS;
% disp(['d OKS = ' num2str(results.SrZnachOKS)  ' км ' ]) %
% disp(['NumPosts = ' num2str(NumPosts)  ' Штук' ])
% k=k+1;
% end
% 
% dX = 0.5;  
% dY = 0.5;  
% dZ = -0.5;  
% params.dX = dX;  
% params.dY = dY;  
% params.dZ = dZ;
% 
%  NumPosts = 3;
%  params.NumPosts = NumPosts;
% for distOP = 200:50:600;  %цикл ра�?чета OKS от ра�?�?тро�?нии по�?тов от центра, SKO = 1гц, кол-во по�?тов 3
% 
%    params.distOP = distOP;
%    
% results=loshadka(params);
% tabl_Znach(1,k) = distOP;
% tabl_Znach(2,k) = results.SrZnachOKS;
% disp(['d OKS = ' num2str(results.SrZnachOKS)  ' км ' ]) 
% disp(['distOP = ' num2str(distOP)  ' км' ])
% k=k+1;
% end

% dX = 1;  
% dY = 1;  
% dZ = -1;  
% params.dX = dX;  
% params.dY = dY;  
% params.dZ = dZ;
% 
% 
% for distOP = 200:50:600;  %цикл ра�?чета OKS от ра�?�?тро�?нии по�?тов от центра, SKO = 1гц, кол-во по�?тов 3
% 
%    params.distOP = distOP;
%    
% results=loshadka(params);
% tabl_Znach(1,k) = distOP;
% tabl_Znach(2,k) = results.SrZnachOKS;
% disp(['d OKS = ' num2str(results.SrZnachOKS)  ' км ' ]) 
% disp(['distOP = ' num2str(distOP)  ' км' ])
% k=k+1;
% end
% 
% dX = 2;  
% dY = 2;  
% dZ = -2;  
% params.dX = dX;  
% params.dY = dY;  
% params.dZ = dZ;
% 
% 
% 
% for distOP = 200:50:600;  %цикл ра�?чета OKS от ра�?�?тро�?нии по�?тов от центра, SKO = 1гц, кол-во по�?тов 3
% 
%    params.distOP = distOP;
%    
% results=loshadka(params);
% tabl_Znach(1,k) = distOP;
% tabl_Znach(2,k) = results.SrZnachOKS;
% disp(['d OKS = ' num2str(results.SrZnachOKS)  ' км ' ]) 
% disp(['distOP = ' num2str(distOP)  ' км' ])
% k=k+1;
% end
 disp (tabl_Znach)

toc
