function [loshadka] =prm(params)
loshadka = GEOLOCATION_LEO_DOPPLER_N_POST_2023
function []=GEOLOCATION_LEO_DOPPLER_N_POST_2023(varargin)
   
%**************************************************************************
%
% ������� GEOLOCATION_LEO_DOPPLER_N_POST_2023 ���������� �������� ����������� 
% ��������� �������� (OKS) �� ��������� ��������� ������� ������� N �������
% (N_POST). ����� ������ � ������ ����� ���� �����.
%
% � �������� �������� ������������ ���������������� ����������� ������� (���, LEO). 
% ������ ����(�������� ����������), � ��������, ����� ������� K ��������� ������� ������� "����" �� ����� ��
% ����� ���. ��������� ���������� �� ��������� ��������� ��� �������. 
% �� ����� ��� ��������� ��������� ������ ������ �����, �, �������������
% ������� ��������� ���������� ���������.
%
% � ����������� � ����������� �������� � ������ ������������ ���� ���� 
% (������) ��������� �� ������ �����.
%
% ��� ������������� �������� (TRU) ���������� �������� �������� ������������ 
% ������ �� ��������� �� "������" ���������� TLE/SGP4. ��� ���� ���� 
% "������" ���������� ����������� � �������� ��������� ������� ���
% ������ �������� ����������� �������. ��� ���� ������ ������������� 
% ������� � ������� ��� �������� �������� ���������� �� ������ 
% (�� ���� ��������� �������� ������� �������� ��� �� ����������� �� ����).
%
% ������ ���������� 4 ��������������� �������: 
% - dopplerSAT1_downLink 
% - OKS_Doppler
% - SOMP_Doppler
% - UKFupdate
%
% ������ dopplerSAT1_downLink ���������� ��������� ������� ������� ���� (�
% ��������� ��������� ����������� ����� ������������� ������� dopplerSAT1)
%
% ������� OKS_Doppler - ��������� ������ �������� ���������� ��������. ���
% � ������ �������� ����������� ������� ������������ �������� �������-����.
%
% ������� SOMP_Doppler ������ ������ ���� ��� (����) �����, �.�. �� K ����������
% ������� ������� �� ����� ��� � �������� � ������ ���������� (TLE/SGP4) ������ 
% ������ ��������� ����� (�� �������� -�������). ��� ������ �������� ����������� 
% ������� ������������ �������� �������-����. �������, � ��������, ��������� 
% ������ ������� ������� ���������� (�������� � ������). ������� �� ������� 
% ���������� �������� ������ ���� ������������� ���������� meanErr_Doppler.
%
% � �������� ������������ ������� OKS_Doppler ������ ��������� ������������
% �������� �������. ��� ����� ������������ ������� UKFupdate, ������� ����������
% ������ ��������� ����������� �����-��������� ������� (���������� ������ 
% ������� �������). � �������� ��������� ��� ������� ������������ ��������� 
% ������� ������� ����, � ������ ��������� ����� 6 ���������: 3 ���������� 
% �������������� � 3 ���������� �������� ���. 
%
% �� ������ ������ ���������� ������ ������� ��� ��������� ���������� ��  ���������� 
% ������� �� ��� �������������� ����������. ��-��������, ������� ������� �
% ���, ��� ������ ��������� ���������� � ������ (� ���� ������ ������� ����
% ���� ��������� �������), � ������ ��������� ����� ������ 6. (����� � 
% �������� ������ ����� ��� ������ �������� �� �����������  �������������
% (������ ��� ����� � � �) � ��� ������ �������. ������ ��� ������ ���������
% (� ���� ������ �������) ���� ��� ����������: ������  � ���� �����.) 
%
% ���������
% � ���� ������ ��� �� �����. ������ ��� ������������� � �����������
% �������� � ������ ����� ���, �������������� ��������, ��������� � 
% ��������������� ������� �����.
%
% ��� ���� (��� �������� ������ � �������) 
% � �������� ����� ������������ ����� ��� ��� �� ������ ���������������-������������� 
% ������ (���, DDM, ��. ���������� DDM_OMP), ����������� ������ � ������� GEOLOCATION_LEO_DDM. 
%
%%*************************************************************************
% ��������� ������ ���������� �������� �� ���������� ����������� �� ������:
% https://yadi.sk/d/dNrJubwO3JFoVB/
%**************************************************************************
%
% var@mail.spbstu.ru, V.A. Vargauzin, �������-������ 2020, 
%
%**************************************************************************

if isempty(varargin)
    help(mfilename) % ������ help
    close all       % ������� ��� �������
    clc
end

%% GPU

n = gpuDeviceCount;
for ii = 1:n
    g = gpuDevice(ii);
    fprintf(1, 'Device %i has ComputeCapability %s \n', ...
            g.Index, g.ComputeCapability);
end

%% ���������� ����������

global min_per_day sec_per_day UnitDist...
       Xp Yp TAI_UTC LOD TT_TAI  GPS_UTC UT1_UTC...
       omega_Earth JD1858 JD1950 JD2000 MJD2000

%% ���������� ����������

GUIinput=boolean(1); % ������������ GUI (����� - ���� ����������� � ��������� ������)
Geoid=boolean(0); % ������������ ��� �� ������������ �����
DDM_OMP=boolean(0); % ��������� �������� DDM ��� ��� ���
rng('default')

%% ��������� ��� ��������� � ������ ������� ���������

UnitDist=10^3; % 1 ��  = 1000 �
Arcs      = 3600*180/pi;         % Arcseconds per radian
Rad       = pi/180;              % Radians per degree
UnitFreq=10^3; % 1 ��� =1000 ��
OneGHz=UnitFreq^3; % 1 ���

%% ���������� ���������

% Earth rotation (derivative of GMST at J2000; differs from inertial period by precession)
omega_Earth = 15.04106717866910/3600*Rad;  % [rad/s]; WGS-84
% omega_Earth=7.29211514670698e-05
c = 299792458; % �������� �����, �/c

%% ��������� �������

min_per_day=60*24;    % ����� ����� � ���� = 1440
sec_per_day=60*24*60; % ����� ������ � ���� = 8640
sec_per_mim=60;       % ����� ������ � ������

[year,mon,day,hr,minute,sec] = datevec('18-Nov-1858 00:00');
JD1858=jday(year,mon,day,hr,minute,sec); % % ��������� ���� (����� ���� �� 1 ������ 4713 �.�.�.�.)
% J1858=2400000.5; % �������� ��� ��������� �� ��������� ���� ���������������� ��������� ���� 

[year,mon,day,hr,minute,sec] = datevec('1-Jan-2000 12:00');
JD2000=jday(year,mon,day,hr,minute,sec);
% JD2000=2451545.0

[year,mon,day,hr,minute,sec] = datevec('1-Jan-1950 00:00');
JD1950=jday(year,mon,day,hr,minute,sec);
% JD1950=2433281.5

MJD2000=JD2000- JD1858;  % ���������������� ��������� ���� JD2000, � ������� ��������� ������������ ������� ������� ��������� IERC
% MJD_J2000 = 51544.5; 

TAI_UTC= +37; % ���, ���������� �������� ������� �� UTC �� 2017 ��� (Beginning 1 January 2017:)

TT_TAI  = +32.184; % ���, ������� ����� ������ � ������� �������� https://ru.wikipedia.org/wiki/%D0%AD%D1%84%D0%B5%D0%BC%D0%B5%D1%80%D0%B8%D0%B4%D0%BD%D0%BE%D0%B5_%D0%B2%D1%80%D0%B5%D0%BC%D1%8F

GPS_UTC= +8*0; % ���, ����� GPS �� ����������� �������� ����� UTC �� 16 ��� !!! 

% ���� 2017:     
UT1_UTC= 0.51363; % ���
% ������������ ��������� ������:
Xp= 0.0804/Arcs; % ��� 
Yp= 0.2640/Arcs; % ��� 
LOD= 0; %                                  
     
%% ������������� ���������� ������ ������
% ������ lat, ���� (-90...+90), ������� lon, ���� (-180...+180)

hPostSee=0; % ������������� ������ �������������� ����� P


               latPost=60.065843; 
               lonPost=30.246163; 
            
          
      
if Geoid       
% �������������� ������ ��� ������� ���� (������� 'egm96') � ������ ��� �����������
% ��� ����� ������������ ����� ����� �����,���������� � ���� geoidegm96grid.mat        

       hPostGeoid = geoidheight( latPost, lonPost); % ������ ������ ��� ����������� � ����� ������������ ����� 
       hPost=hPostSee+hPostGeoid; % ������������� ������ �������������� �����     
       
end
 
%% ����� ����� ������
 
NumPosts= 3; % ����� ������                          


 %% ������������ ������� ������ (OP) ������ ������

    % ����� ����������� ������� 3 x NumPosts (���������� x �����) - ���������
    % ���������� ������������� ������ P � (������) ��������������� 
    % ������� ��������� (����)
    % ����� ���� ��������� � ������� ���� �����
    % ��� OZ - ����� �� ����� ������� ���������� � ���������� �� �������� �����
    % ��� OX - ���������� � ����� ����������� �������� (����������) � ���������� ��������.     
    
    % ����� OP:
    latOP_0=latPost;
    lonOP_0=lonPost;
%     hOP_0=hPost;
    distOP=0; % ������ ������ �����

if NumPosts > 1
    
    distOP=350*UnitDist; % ��������� �� ������ ������, �
    
    
 end    
    
    xyzOP=nan(3,NumPosts);
    az=linspace(-180,+180,NumPosts+1);
    az=az(2:end);
    latOP=nan(size(az));
    lonOP=nan(size(az));                  
         
    hOP=zeros(size(az)); % ������ ������ ��� �������
    hSeeOP=zeros(size(az)); % ������ ��� ������� ���� ������
    hAntennaOP=zeros(size(az)); % ������ ������ ������          
         
    for j=1:length(az)         
        [latOP(j), lonOP(j)]=vreckon(latOP_0,lonOP_0, distOP, az(j)); 
        [x, y, z] = xyz( latOP(j), lonOP(j), hOP(j), Geoid, hSeeOP(j), hAntennaOP(j) ); 
        xyzOP(:,j)=[x, y, z]';              
    end

% ��� ������������� � ���������� �������:
latPost=latOP(1);
lonPost=lonPost(1);
  
%% ����� ���           

    
      %TLE:
      % https://www.space-track.org/#/tle
      %
      % ����� ����� �������� TLE ��� �����, �������� ������� �� ������� ���������.
      %
      % ��������� TLE, ������, ������ ����� �������� (��� �����������, � ������� �� ����,
      % ��� ��������� �� ����� space-track) �������� �� ������:
      %
      % https://www.n2yo.com
      %
      % ��� ����, �� ���� �� ����� ����������� � ������ �������� �� ���������� � �������� �������, 
      % ��� ������, �������, ������ �������� ����������� ����� ����������� (��� - ������ �����). 
      % � ���������, ���� �������� � �� UFO-10 (TLE ��������, ��� ����������, ����������� ����� �����),
      % ��.:
      %
      % https://www.n2yo.com/satellite/?s=25967    
    
           
               % IRIDIUM 151               
                lonstr1KA1='1 43074U 17083E   20352.49900893  .00000110  00000-0  32273-4 0  9994';
                lonstr2KA1='2 43074  86.3979  81.0759 0001886 103.9895 256.1511 14.34217107156373';               
                            

           
%% ������������ ��������� ��� 
% � ���� ������ ��� �� �����. ������ ��� ������������� � �����������
% �������� ������� �������������� ���, ����������� � ��������������� �������
% �����

hIRI=0;    % ������������� ������ ���
hIRIsee=0; % ������ ��� ������� ����

% Menu= { '��������� � ������ ������'...
%         '��������� � ������� ������'...
%         '������� �� ������ ������ �� �����'...
%         '������� �� ������ ������ �� ��'...
%         '������� �� ������ ������ �� ������'...
%         '������� �� ������ ������ �� �����'};
% 
%       
%       titl='�������������� ���';
%       [select,OK]=listdlg('ListString',Menu,...
%                             'Name',titl,...
%                             'SelectionMode','single',... 
%                             'ListSize',[350 100]); % [width height]
%        if ~OK
%            disp('���� ��������� �� �������')           
%            return                     
%        end       

     
               
           % ��������� � 1-� ������
               latIRI=latOP(1); 
               lonIRI=lonOP(1); 
               hIRIsee=0;
               az=0;           
              
                            
    

    % ����� ��������� �� �����
    distOA=350*UnitDist; % ��������� ������� ������� (OA) �� �����, �
   

    [latIRI, lonIRI]=vreckon(latOP_0,lonOP_0,distOA,az);

    % �������� ��������� ����-��� ����� ������������� ����������: 
    d = vdist(latIRI,lonIRI,latOP_0,lonOP_0);
    disp(['��������� ����� ������-��� recon = ' num2str(d/UnitDist) ' ��'])
    disp(' ')
       
     
% �������������� ������ ��� ������� ���� (������� 'egm96') � ������ ��� �����������
% ��� ����� ������������ ����� ����� ������ �� �����������,���������� � ���� geoidegm96grid.mat        

       hIRIGeoid = geoidheight( latIRI, lonIRI); % ������ ������ ��� ����������� � ����� ������������ ����� 
       hIRI=hIRIsee+hIRIGeoid; % ������������� ������ �������������� �����             


% % ��������������� ���������� xyzIRI: 
 [xIRI, yIRI, zIRI]=llh2xyz( latIRI, lonIRI, hIRI );
 xyzIRI=[xIRI, yIRI, zIRI];

    
%% ������������� ������������� (���������������) ���������� �������� ��������� �� ��������� SPG4

    [satrecKA1,~, utc0_KA1, NORAD_ID_KA1] = spg4_initVARG(lonstr1KA1, lonstr2KA1);   
    
     % ��������� satrecKA1 (���������� ������):
     
%          inclo: 1.7027 - inclination - ���������� ������ � ��������, (180/pi) * angleInRadians= �������
%           ecco: 0.0059 - eccentricity - ��������������� ������
%          nodeo: 0.5558 - right ascension of ascending node - �������                        
%                          ����������� ����,  (180/pi) * angleInRadians= ������� 
%          argpo: 0.5347 - argument of perigee (output if ds) - ��������
%                          �������
%             mo: 5.7566 - mean anomaly (output if ds)
%             no: 0.0647 - mean motion - ������� ��������
%              a: 1.0977 - ������� ������� ������� ? - ��������
%           alta: 0.1041 - ������ � ������� ?
%           altp: 0.0912 - ������ � �������-?

% ����� �������������� ( -> ) ������������ ������� ��������� (ECI), 
% ������������ � ��������� SGP4 (TEME), � ������������� 
% ������ ��������������� ������� (ECEF)
    TEME2ECEF='ECI SGP4 (TEME, True Equator Mean Equinox) -> PEF (Pseudo Earth Fixed) + TOD (True of Date)  = ECEF'; 
    % �������� ������� TEME2ECEF='ECI TEME -> PEF' ; 

%%  ����� ���������� ��������� ������ �������� ����������� ������� ������� �������-����
%
% ���������� ������ ������ �������� ������� ���������� ���������� 
% ����������� ������� �������-���� (Nelder-Mead) ������ ������� �����:
%
% http://www.jasoncantarella.com/downloads/NelderMeadProof.pdf
%
% � �����:
%
% http://www.mathworks.com/help/matlab/math/optimizing-nonlinear-functions.html#bsgpq6p-11
%
% ��� ���� ���������� (������ � �������) �������� - ��������������
% �����������.
% 
% ��� ��� ���������� (������, ������� � ������) �������� - ��������������
% ��������.
%
% �������� ���� ����������� ��������. ����������� ��������� ��������, � 
% ������� ���������� ����� ��������� �� ��������� �������� �������� tolx_TDOA,
% �, ��� ���� �������� ����� ����������� ��������� ����������� � ����� �� 
% ������ ��������� (� ����� ������ 1) ����� ���������� ����������� � ��������� 
% �������� ��������� �� ��������� �������� tolf_TDOA.
%
% ����������� ��������� ����� ��������: 
% - ������������ ����� ��������  maxiter ��� ������ ��������
% - ������������ ����� ���������� �������� ����������� maxfun ��� 
%   ��� ������ ��������
%
% ������������ ��������� �������� ���������� � ������������� ���������
% �����������, ��������� � �.�.
%
% �������� �������-���� �� �++:
% http://www.codeguru.com/cpp/article.php/c17505/Simplex-Optimization-Algorithm-and-Implemetation-in-C-Programming.htm   
       
tolx = 1e-4;
tolf=1e-3;    
maxfun = 200*3; % ������������ ����� ���������� �������� �����������
                % 3 - ��� ������ ����� ���������� ��� �����������
maxiter= 200*3; % ������������ ����� ��������                          



%% ��������� ������

timeEst=cell(1);
% ���� ������ ���������:
timeEst{1}=utc0_KA1; % ����� TLE

% DT_=30; % �������� ���������, ���
% DT_ALL= 1*min_per_day/2; % �������� ����������, ��� (1440 ��� = min_per_day = �����)  
% dX=1;  % ��
% dY=-1; % ��
% dZ=1;   % ��
% DT_=60*10; % �������� ���������, ���



DT_=abs(30*1); % ���
DT_ALL=abs(1*min_per_day/1); % ���
dXYZ=UnitDist*[1 1 -1]; % ����� 

freqDownLink=7*OneGHz; % ��
RMS_DOPPLER=1; % ��
OKSType=0;
min_elev_deg=10;

if OKSType > 3
    OKSType=2;
end   
    
%% ������������� ��������� �������� ������� ��������� ��������� 
% (�������� ��������� ����� ���� ����������� ��� ������� ������ ���������
%  ���-���-����)
        
        [year,mon,day,hr,minute,sec] = datevec(timeEst{1}); % ���� ��������
        jdEst = jday(year,mon,day,hr,minute,sec); % ��������� ���� ��������, ���    
        jd1 = jday(year,mon,day,hr,minute,sec); % ��������� ���� ������� ��������� ��������, ���                                            
        DT=0; % �������� ������� �������� ������������ ������� ���������, ���           
        i=1;
        DT_min=DT_/60; % �������� ���������, ���
        jdEst=(jdEst*min_per_day +0)/min_per_day;

        while DT < DT_ALL
          
          i=i+1;           
          jdEst=(jdEst*min_per_day +DT_min)/min_per_day;         
          [year,mon,day,hr,minute,sec] = invjday(jdEst);
          dtm = datenum(year,mon,day,hr,minute,sec);
          timeEst{i}=datestr(dtm);            
          DT=(jdEst-jd1)*min_per_day; % �������� ������� �������� ������������ ������� ���������, ���         

        end        
                
%% ������������� ���������� � ��������� ���

        Nexp=numel(timeEst); % ����� ���������      
        adrKA1=zeros(1,Nexp);
        
        % ���������� SGP4:        
        SGP4=zeros(Nexp,3);  % "������" ���������� (X,Y,Z), ���������� ������� ��������������� TLE/SGP4            
        SGP4_V=zeros(Nexp,3);%  ���������� ������� �������� (Vx,Vy,Vz)
        
        % ������ ��� ������� ������ (OP):
        SGP4_DOPPLER = zeros(Nexp,NumPosts); % ���������� ������������� ��������
        SGP4_V_PROECT=zeros(size(SGP4_DOPPLER)); % ���������� �������� ������� �������� � ����������� �� ����                       
        
        % �������� ���������� TRU(E):
        TRU=zeros(Nexp,3);  % �������� ���������� (X,Y,Z), ��������� �� TLE/SGP4            
        TRU_V=zeros(Nexp,3);%  ���������� ������� �������� (Vx,Vy,Vz)
        
         % ������ ��� ������� ������ (OP):
        TRU_DOPPLER = zeros(Nexp,NumPosts); % ���������� ������������� ��������
        TRU_V_PROECT=zeros(size(SGP4_DOPPLER)); % ���������� �������� ������� �������� � ����������� �� ���� 
        MEAS_DOPPLER=zeros(size(SGP4_DOPPLER)); % ��������� ������������� ��������
        MEAS_V_PROECT=zeros(size(SGP4_DOPPLER)); % ��������� �������� ������� �������� � ����������� �� ����         
                
       % ���������� dt - �������� ������� ��������� (� ���)
       % ������������ ������� ���������        
        dt=zeros(Nexp,1); 
        [year,mon,day,hr,minute,sec] = datevec(timeEst{1}); % ���� ���������
        jd1 = jday(year,mon,day,hr,minute,sec); % ��������� ���� ������� ��������� ��������, ��� 
         
        % ���� �� ���-�� ��������� ��� ���������� ���� ���������:       
        if Nexp >1
            for iexp = 2:Nexp
                % ���������� ������� ������� ���������������
                [year,mon,day,hr,minute,sec] = datevec(timeEst{iexp}); % ���� ��������
                jdEst = jday(year,mon,day,hr,minute,sec); % ��������� ���� ��������, ���                                  
                dt(iexp)=(jdEst-jd1)*min_per_day; % �������� ������� �������� ������������ ������� ����������, ���           
            end         
        end 
        
     % ���� �� ���-�� ��������� ��� ���������� ���������� ���:       
     for iexp = 1:Nexp             

                [ r_ecf, v_ecf, satrecKA1, adrKA1(iexp) ] = SPG4_ecef(satrecKA1, timeEst{iexp},TEME2ECEF );
                sat_TLE_xyz(1) = r_ecf(1); % �   
                sat_TLE_xyz(2) = r_ecf(2); % �
                sat_TLE_xyz(3) = r_ecf(3); % �                
                SGP4(iexp,:)=sat_TLE_xyz;                
                              
                sat_TLE_Vxyz(1) = v_ecf(1); %  �/�    
                sat_TLE_Vxyz(2) = v_ecf(2); %  �/�   
                sat_TLE_Vxyz(3) = v_ecf(3); %  �/�                
                
                SGP4_V(iexp,:)=sat_TLE_Vxyz;
                
               [SGP4_DOPPLER(iexp,:), SGP4_V_PROECT(iexp,:)] =...
                   dopplerSAT1_downLink(sat_TLE_xyz, sat_TLE_Vxyz, xyzOP, freqDownLink); 
               
                % �������� ����������:
                sat_TRUE_xyz=sat_TLE_xyz+dXYZ; % �������� ���������� TRUE ������������ SGP4
                TRU(iexp,:)=sat_TRUE_xyz;                
                               
                % �������, ��� �������� � ���������� TRU �� �� �������
                % ��������, ��� � � ���������� SGP4:                
                TRU_V(iexp,:)=SGP4_V(iexp,:);
               
                % ������ ��� �������� ����������:
                [TRU_DOPPLER(iexp,:), TRU_V_PROECT(iexp,:)] =...
                    dopplerSAT1_downLink(sat_TRUE_xyz, sat_TLE_Vxyz, xyzOP, freqDownLink);
                
                MEAS_DOPPLER(iexp,:)=TRU_DOPPLER(iexp,:)+ RMS_DOPPLER*randn(1,NumPosts);
                MEAS_V_PROECT(iexp,:)=c*MEAS_DOPPLER(iexp,:)/freqDownLink;                
              
    end       
    
%% ������ ���������� ������ ��������� (los) ��� ������ � ���

% min_elev_deg=10; % minimal elevation, deg

[lat,lon,h]=xyz2llh(TRU(:,1), TRU(:,2), TRU(:,3) ); % lat,lon,h - �������-������� 
losCUBIK=ones(size(lat));

losAnaliz=boolean(1); % ���� losAnaliz=0, �� ����� ������������� ��� 
% ���������� �������� ���������� �� ��� ��������� ������� � ��� (�����
% ���������� ��������� ��������� "����� ������")

if losAnaliz

    % ��� ����� ������ �� NumPosts ����:
    losPost=nan(length(lat), NumPosts);
    losCUBIK_Post=ones(size(lat)); % ���������� ��������� ����� �������

    for j=1:NumPosts

        % ENU �� �����:
        refLat=latOP(j);
        refLong=lonOP(j);
        refH=hOP(j);
        [~,~,losPost(:,j)]=peleng(refLat,refLong,refH,lat,lon,h,min_elev_deg);
        losCUBIK_Post=losCUBIK_Post & losPost(:,j);
    end

    % ������ ������ ��������� ��� - ��� 

    % ENU �� ���:
    refLat=latIRI;
    refLong=lonIRI;
    refH=hIRI;
    [~,~,losCUBIK_IRI]=peleng(refLat,refLong,refH,lat,lon,h,min_elev_deg);

    % ���������� ���������:
    losCUBIK=losCUBIK_Post & losCUBIK_IRI;
    % ��������� ���������� ��� ��� ���������� ������ ��������� ������ � ���:

end

% ���������� SGP4 - > sgp4:
sgp4=SGP4;
sgp4(losCUBIK==0,:)=nan; % ������� ���������� ���������� ���

sgp4_V=SGP4_V;
sgp4_V(losCUBIK==0,:)=nan; % ������� ���������� ���������� �������� ���

sgp4_DOPPLER=SGP4_DOPPLER;
sgp4_DOPPLER(losCUBIK==0,:)=nan; % ������� ���������� ���������� ������� ���

sgp4_V_PROECT=SGP4_V_PROECT;
sgp4_V_PROECT(losCUBIK==0,:)=nan; % ������� ���������� ���������� �������� �������� ���

% ���������� TRU - > tru:
tru=TRU;
tru(losCUBIK==0,:)=nan; % ������� ���������� ���������� ���

tru_V=TRU_V;
tru_V(losCUBIK==0,:)=nan; % ������� ���������� ���������� ������� �������� ���

tru_DOPPLER=TRU_DOPPLER;
tru_DOPPLER(losCUBIK==0,:)=nan; % ������� ���������� ������� �������

tru_V_PROECT=TRU_V_PROECT;
tru_V_PROECT(losCUBIK==0,:)=nan; % ������� ���������� ���������� �������� �������� ���

% ���������:
meas_DOPPLER=MEAS_DOPPLER;
meas_DOPPLER(losCUBIK==0,:)=nan; % ������� ���������� ������� �������

meas_V_PROECT=MEAS_V_PROECT;
meas_V_PROECT(losCUBIK==0,:)=nan; % ������� ���������� ���������� �������� �������� ���

latSAT1=lat;
latSAT1(losCUBIK==0,:)=nan;  % ������� ������ ���

lonSAT1=lon;
lonSAT1(losCUBIK==0,:)=nan; % ������� ������� ���

%% ���������� �������� ������� �� ������� �����

Vitok=cell(1);
ind=find(losCUBIK);

if ~isempty(ind) && length(ind) > 1
    
    v=1; % ������ ������� �����
    Vitok{v}=ind(1);

    for i=2: length(ind)
        if (ind(i)-ind(i-1))==1
           Vitok{v}=[Vitok{v} ind(i)] ; % ���������� �������� �������� �������� �����
        else
            v=v+1; % ��������� ������� �����
            Vitok{v}=ind(i);
        end
    end  
    
end
% 
%% �������� ������� ������ ���������� �� ����� �� ���������� �������� �������

Kalman=boolean(0);
if OKSType==1 || OKSType==2
  Kalman=boolean(1);  
  KalmanType=OKSType;
end

 if Kalman
     
     if RMS_DOPPLER ~=0
        RMS_KALMAN=RMS_DOPPLER*c/freqDownLink; % ����������� ��������� ������� �������     
     else         
        RMS_KALMAN=10; % ��        
     end
     
     % ������ ��������� ������� =
     % S=[X Y Z Vx Vy Vz]
     S=zeros(1,6);
     
     D=10^4; % ����������� ���������� ���������� �������� ������������� ��������� (� ������)      
     V=10^4; % �/� 
     
     accelX=10^5; % ��������� ��������� �� ���������� X
     accelY=10^5; % ��������� ��������� �� ���������� Y
     accelZ=10^5; % ��������� ��������� �� ���������� Z     
 
     dT=DT_; % ���
     % ������� �������� � ����� ��������� (���������� �������� ���� dx=x+v*dt):
     A = [1   0  0   dT  0   0
          0   1  0   0   dT  0 
          0   0  1   0   0   dT
          0   0  0   1   0   0
          0   0  0   0   1   0
          0   0  0   0   0   1];

    % ������������� ������� ���������� ���� �������� � ���������� �����������
    % ��������� �� ����������� X � Y:  
      U=[dT^4/4   dT^3/2
         dT^3/2   dT^2];

    % ������� ��������� ��� ����������� ���� ��������:   
      accel=[accelX     0        0
               0      accelY     0
               0        0      accelZ];

    % ������������� ������� ���� ����������� ��������, ������������ �����������, 
    % ���������� �������� accel:         
    Q=[accel*U(1,1)    accel*U(1,2);...
       accel*U(2,1)    accel*U(2,2)];       

    if KalmanType ==1
        % � �������� ��������� ������������ ������ �������� �������
         R=RMS_KALMAN .^ 2*diag(ones(1,NumPosts));                
    else 
        % � �������� ��������� ������������ �������� ������� � ������
        % TLE/SGP4
         Q=NumPosts + numel(S); % ����� ��������� � ������� ���������
         R=diag(ones(1,Q));
         for j=1:NumPosts
             R(j,:)=RMS_KALMAN .^ 2*R(j,:);
         end

         dxyz=10^2;
         for j=NumPosts+1:NumPosts+1+2
             R(j,:)=dxyz^2*R(j,:);
         end

         dVxyz=10^1;
         for j=NumPosts+1+2+1:Q 
             R(j,:)=dVxyz^2*R(j,:);
         end         
    end               
     
     NumVitkov=length(Vitok); % ����� ������� ������
     tru_est=sgp4;
     tru_V_est=sgp4_V;
     b=nan(1,NumVitkov);
     q=nan(1,NumVitkov);
     
   for v=1:NumVitkov
       
     ind=Vitok{v};
     M=length(ind); % ����� ��������� �� �����      
          
     % ��������� �������� ������� ��������� ������� =
     S=[sgp4(ind(1),:) sgp4_V(ind(1),:)];       
     
    % ��������� �������� �������������� ������� 
    % ����������� ����������  ��������� ������� ���������
     P=[D^2  0   0   0   0   0
        0  D^2  0   0   0   0
        0   0  D^2  0   0   0 
        0   0   0  V^2  0   0
        0   0   0   0  V^2  0
        0   0   0   0   0  V^2];         

    for t = 1:M
         % ���������� ������� ��������� ���������� �������� ������� (UTF)
         % � �������� ��������� ������������ �� ������������ ������� f, � 
         % ���������� �������� v (tru_V_PROECT), ����� � �������� �� ������� 
         % �������� ������� ���� freqDownlink. ��� ��������, ��������� 
         % f=(v/c)*freqDownlink, ��� c - �������� �����.

         if KalmanType==1         
            u=meas_V_PROECT(ind(t),:);     
         else
            u=[meas_V_PROECT(ind(t),:) sgp4(ind(t),:) sgp4_V(ind(t),:)];
         end

          [S, P] = UKFupdate(u,xyzOP,S,P, A, Q, R, KalmanType);

          tru_est(ind(t),:)=S(1:3);
          tru_V_est(ind(t),:)=S(4:6);      
    
    end     
    
      err=tru_est(ind,:)-tru(ind,:);
      err=err';
      err_=sum(err.^2); % ������������ �� ��������, � ������� ���������� x y z
      b(v)=sqrt(mean(err_));
      
      err=tru_V_est(ind,:)-tru_V(ind,:);
      err=err';
      err_=sum(err.^2); % ������������ �� ��������, � ������� ���������� x y z
      q(v)=sqrt(mean(err_)); 
      
   end

  pos_err_Kalman=mean(b);
  vel_err_Kalman=mean(q);
% 
 end

%% �������� ������ �������� �������� (OKS) �� ������� �� N ������ �� �������

OKS=boolean(0);
if OKSType==0
    OKS=boolean(1);
end

if OKS
 NumVitkov=length(Vitok); % ����� ������� ������
 % ��������� �KS:
 tru_est_OKS=nan(size(sgp4));
 tru_V_est_OKS=nan(size(sgp4));
 b=nan(1,NumVitkov);
 q=nan(1,NumVitkov); 

 for v=1:NumVitkov

     ind=Vitok{v};
     ind=ind(1); % ��������� ���������
     % ������� ���������� sgp4: 
     SAT=sgp4(ind,:); 
     SAT_V=sgp4_V(ind,:); 
     % ��������� ������� �� �������� ����������:
     Dop=meas_V_PROECT(ind,:); % � �������� ��������� ������������ ��      
    %  Dop=tru_V_PROECT(ind,:); % � �������� ��������� ������������ �� 
     % ������������ ������� f (tru_DOPPLER), �  ���������� �������� v, 
     % ����� �� ������� ������� ������� ���� freqDownLink. ��� ��������,
     % ��������� f=(v/c)*freqDownLink, ��� c - �������� �����
     [SAT_est,SAT_V_est]=OKS_Doppler(SAT,...
                                     SAT_V,...                          
                                     Dop,...
                                     xyzOP,...                                   
                                     tolx,...
                                     tolf,...
                                     maxiter,...
                                     maxfun);  
                                     
        tru_est_OKS(ind,:)=SAT_est; 
        tru_V_est_OKS(ind,:)=SAT_V_est;
        
      err=tru_est_OKS(ind,:)-tru(ind,:);
      err=err';
      err_=sum(err.^2); % ������������ �� ��������, � ������� ���������� x y z
      b(v)=sqrt(mean(err_));
      
      err=tru_V_est_OKS(ind,:)-tru_V(ind,:);
      err=err';
      err_=sum(err.^2); % ������������ �� ��������, � ������� ���������� x y z
      q(v)=sqrt(mean(err_));         
 end
 
   pos_err_OKS=mean(b);
   vel_err_OKS=mean(q);

   disp(['������� �������� ������ ��S = ' num2str(pos_err_OKS/UnitDist)  ' �� ' ])
   disp(' ')
 
end

%% ����� ����� ��� ���
PostNumber=1;

if DDM_OMP


  %% ������������ ��������� ��� ��� ��� ������� ���
  
  diffD=DDM_measure(tru,xyzOP(:,PostNumber),xyzIRI); % ������-������� ��������� ��������� (� ������� ���������)
  noise_Diff_D=0; % ���
  diffD_measure=diffD +noise_Diff_D; % ��������� � ����� 
  
%% �������� ��� ��� ������� ��� �� ������ ����������� ������������� ����������� �������
           
 opt=[latOP(PostNumber), lonOP(PostNumber)] ; % ��������� ������ ��� ������ �������� ������������� �����������                           
 h0=0; % ������������� ������ ��� ������ 
 
 NumVitkov=length(Vitok); % ����� ������� ������
 latIRI_DDM=nan(NumVitkov,1);
 lonIRI_DDM=nan(NumVitkov,1);
 NumdiffD=nan(NumVitkov,1); % ����� ��������� �� ������� ����� 

 for v=1:NumVitkov

     ind=Vitok{v};
     % ��� ��� ������������ �� �������� ���������� tru, � ������ sgp4
     % !!!!!!!!!!!!! ����� ������ ���������� ���� �� �������� sgp4 �� tru_est:
%          SAT=tru_est(ind,:); 
     SAT=sgp4(ind,:); % ������� ������� N x 3
     diffD=diffD_measure(ind(1:end-1)); % 
     NumdiffD(v)=length(diffD);

     [opt,~, exitflag]=DDM(SAT,...
                           diffD,...
                           xyzOP(:,PostNumber),...
                           opt,...
                           h0,...
                           tolx,...
                           tolf,...
                           maxiter,...
                           maxfun);                

     if ~ exitflag
          MinFlagStr='��� ��� ������� ���: min ������������� ����������� ������� �� ������ !!!';
          disp(MinFlagStr)
     end

    latIRI_DDM(v)=opt(1); % ������ ������ ���
    lonIRI_DDM(v)=opt(2); % ������ ������� ���    
 
end

%%  ������ ������ ��� ��� ������� ���

err_DDM=nan(NumVitkov,1);

if ~ isempty(Vitok{1})
    
    for v=1:NumVitkov

        % ����� ��������� ����������:
        h0=0;
        [xIRIest, yIRIest, zIRIest]=llh2xyz( latIRI_DDM(v), lonIRI_DDM(v), h0 );             
        err_DDM(v)=sqrt((xIRIest-xIRI).^2 + yIRIest-yIRI).^2+(zIRIest-zIRI).^2; 

%         if boolean(err_DDM(v))
%             % ����� ������������� ����������:                               
%              err_DDM(v)= vdist(latIRI_DDM(v),lonIRI_DDM(v),latIRI,lonIRI);
%         end
    end
end

meanErr=nanmean(err_DDM);
disp(['������� �������� ��� ��� ������� ���, ��,  = ' num2str(meanErr/UnitDist)])
disp(' ')

end % if DDM_OMP

  
%% �������� C��� ����� �� ������� �� ��S �� ������ ����������� ������������� ����������� �������

%  opt=[latPost, lonPost] ; % ��������� ������ ��� ������ �������� ������������� ����������� 
 opt=[latOP(PostNumber), lonOP(PostNumber)] ; % ��������� ������ ��� ������ �������� ������������� ����������� 
 h0=0; % ������������� ������ ��� ������ 
 
 NumVitkov=length(Vitok); % ����� ������� ������
 % ��������� ���:
 latPost_Doppler=nan(NumVitkov,1);
 lonPost_Doppler=nan(NumVitkov,1);

 for v=1:NumVitkov

     ind=Vitok{v};
     ind=ind(1); % ��������� ���������
     % ������� ���������� sgp4: 
     SAT=sgp4(ind,:); 
     SAT_V=sgp4_V(ind,:); 

%      SAT=tru_est_OKS(ind,:);
%      SAT_V=tru_V_est_OKS(ind,:);

     % ��������� ������� �� �������� ����������:
     %Dop=tru_V_PROECT(ind,PostNumber); % �������� ������ 
     Dop=meas_V_PROECT(ind,PostNumber); % ��������� (� �������)
     % � �������� ��������� ������������ �� 
     % ������������ ������� f (tru_DOPPLER), �  ���������� �������� v, 
     % ����� � �������� SOMP_Doppler �� ������� ������� ������� ����
     % freqDownLink. ��� ��������, ��������� f=(v/c)*freqDownLink, ��� c - 
     % �������� �����
          [opt,~, exitflag]=SOMP_Doppler(SAT,...
                                         SAT_V,...                          
                                         Dop,...
                                         opt,...
                                         h0,...
                                         tolx,...
                                         tolf,...
                                         maxiter,...
                                         maxfun);  

     if ~ exitflag
          MinFlagStr='���� �� �������: min ������������� ����������� ������� �� ������ !!!';
          disp(MinFlagStr)
     end

    latPost_Doppler(v)=opt(1); % ������ ������ �����
    lonPost_Doppler(v)=opt(2); % ������ ������� �����    
 
end

%%  ������ ������ ��� ����� �� ������� �� ��S

errPost_Doppler=nan(NumVitkov,1);

if ~ isempty(Vitok{1})
    
    for v=1:NumVitkov

        % ����� ��������� ����������:
        h0=0;
        [xIRIest, yIRIest, zIRIest]=llh2xyz( latPost_Doppler(v), lonPost_Doppler(v), h0 );             
        errPost_Doppler(v)=sqrt((xIRIest-xyzOP(1,PostNumber)).^2 + (yIRIest-xyzOP(2,PostNumber)).^2+(zIRIest-xyzOP(3,PostNumber)).^2); 

        if boolean(errPost_Doppler(v))
            % ����� ������������� ���������� ����������:                               
             errPost_Doppler(v)= vdist(latPost_Doppler(v),lonPost_Doppler(v),latOP(PostNumber),lonOP(PostNumber));
        end
    end
end

meanErr_Doppler=nanmean(errPost_Doppler);
disp(['������� �������� ������ ���� ����� �� ������� = ' num2str(meanErr_Doppler/UnitDist)  ' �� ' ])
disp(' ')


    
%% ������� ����� ���������



end % GEOLOCATION_LEO_DOPPLER_N_POST
end