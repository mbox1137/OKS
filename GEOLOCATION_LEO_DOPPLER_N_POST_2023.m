function []=GEOLOCATION_LEO_DOPPLER_N_POST_2023(varargin)
%**************************************************************************
%
% ������� GEOLOCATION_LEO_DOPPLER_N_POST_2023 ���������� �������� ����������� 
% ��������� �������� (OKS) �� ��������� ��������� ������� ������� N �������
% (N_POST). ����� ������ � ������ ����� ���� �����.

% � �������� �������� ������������ ���������������� ����������� ������� (���, LEO). 
% ������ ����, � ��������, ����� ������� K ��������� ������� ������� "����" �� ����� ��
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
% (� ���� ������ �������) ���� ��� ����������: ������ � ���� �����.) 
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

Menu= {'���-�����'...  
       '���������'...
       '�����'...
       '��������'...   
       '�������'...
       '�����������'...
       '����'...
       '������'...
       '�����-������'...
       '������� �� ������� ����������'};   
      
      titl='����� ������';
      [select,OK]=listdlg('ListString',Menu,...
                            'Name',titl,...
                            'SelectionMode','single',... 
                            'ListSize',[300 160]); % [width height]
       if ~OK
           disp('����� ������ �� ������')           
           return                     
       end
       
     titlePost=Menu{select};
  
       switch select
           case 1
               % ���
               latPost=60.065843; 
               lonPost=30.246163; 
            
           case 2
               % ���������
              latPost=45.026896; 
              lonPost=39.001564;              
           case 3
               % �����
              latPost=56.752483; 
              lonPost=37.202713; 
           case 4
               % ��������
               latPost=54.46; 
               lonPost=32.2;  
              
           case 5
               % �������
               latPost=44.59; 
               lonPost=41.7;
          case 6
               % �����������               
               latPost=44.58883; 
               lonPost=33.5224; 
               hPostSee=0;               
           case 7
               % ����
               latPost=44.36; 
               lonPost=39.44; 
               hPostSee=0;               
           case 8
               % ������
               latPost=33.30; 
               lonPost=36.17; 
               hPostSee=0;               
           case 9
               % �����-�����
               latPost=9.1; 
               lonPost=38.44; 
               hPostSee=0;
           case 10
              % ������� �� ������ ����������
              latPost=0.0; 
              lonPost=39.001564;             
       end       
      
if Geoid       
% �������������� ������ ��� ������� ���� (������� 'egm96') � ������ ��� �����������
% ��� ����� ������������ ����� ����� �����,���������� � ���� geoidegm96grid.mat        

       hPostGeoid = geoidheight( latPost, lonPost); % ������ ������ ��� ����������� � ����� ������������ ����� 
       hPost=hPostSee+hPostGeoid; % ������������� ������ �������������� �����     
       
end
 
%% ����� ����� ������
 
NumPosts= 1; % ����� ������                          

Default_Answer = {int2str(NumPosts)};  
Menu={'����� ������'}; 
select=dispInput(Menu,GUIinput, '����� ������', Default_Answer);
if isempty(select)
       disp('����� ������ �� �������')                
       return
end                                  

NumPosts=abs(select(1)); 
 
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
    Default_Answer = {num2str(distOP/UnitDist)};  
    Menu={'���������� �� ������ ������, ��'}; 
    select=dispInput(Menu, GUIinput, '', Default_Answer);
    if isempty(select)
     disp('���������� �� ������ ������ �� �������')                
     return
    end
    distOP=select(1)*UnitDist; % �����
    if distOP == 0
        NumPosts=1; % ��� ������� ��������� �� ������ ������ ������������� 
        % ������������� ���� ����, ����������� � ������� �����
    end
    
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

Menu= {'IRIDIUM 151' ...
       'IRIDIUM 171'...
       'IRIDIUM 173 '...
       'GLOBALSTAR M076'...
       'GLOBALSTAR M077'...
       'GLOBALSTAR M085'...
       'GLOBALSTAR M004'...
       '����� FUNcube-1 (AO-73)'...   
       '����� Cute-1.7+APD II'...            
       '����� SO-50 (SAUDISAT 1C)'...
       'MERIDIAN 2'...
       'MERIDIAN 4'...
       'MERIDIAN 6'...
       'MERIDIAN 7'...
       'MERIDIAN 8'...
       'MERIDIAN 9'}; 
             
      titl='���';
      [selectNKA,OK]=listdlg('ListString',Menu,...
                            'Name',titl,...
                            'SelectionMode','single',... 
                            'ListSize',[500 300]); % [width height]
       if ~OK
           disp('��� �� ������')           
           return                     
       end

switch selectNKA
    
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
    
           case 4
                % GLOBALSTAR M076 [+]     
                lonstr1KA1='1 37190U 10054C   20338.24947745 -.00000114  00000-0 -58496-4 0  9992';
                lonstr2KA1='2 37190  52.0006 143.3762 0000591 108.6725  64.0302 12.62272048467887';
           case 5
                % GLOBALSTAR M077 [+]     
                lonstr1KA1='1 37191U 10054D   20337.83272597 -.00000119  00000-0 -84904-4 0  9991';
                lonstr2KA1='2 37191  52.0066 101.0273 0000802 117.1670  19.7231 12.62267223468058';
           case 6
                % GLOBALSTAR M085 [+]     
                lonstr1KA1='1 37742U 11033D   20338.16851090 -.00000078  00000-0  14180-3 0  9990';
                lonstr2KA1='2 37742  51.9859 230.3343 0000981  75.2803  89.7754 12.62264295434924';
           case 7 % 4 GLOBALSTAR M004
                lonstr1KA1='1 25163U 98008B   21016.82236812 -.00000072  00000-0  17848-3 0  9995';
                lonstr2KA1='2 25163  52.0002 292.6049 0003193 328.4593  72.7602 12.63295617 57895';
           case 1
               % IRIDIUM 151               
                lonstr1KA1='1 43074U 17083E   20352.49900893  .00000110  00000-0  32273-4 0  9994';
                lonstr2KA1='2 43074  86.3979  81.0759 0001886 103.9895 256.1511 14.34217107156373';               
           case 2
               % IRIDIUM 171                
                lonstr1KA1='1 43929U 19002H   20351.87452659  .00000089  00000-0  24893-4 0  9993';
                lonstr2KA1='2 43929  86.3979 112.9788 0001936 102.7395 257.4017 14.34216896101106';               
           case 3
               % IRIDIUM 173
                lonstr1KA1='1 43925U 19002D   20351.83013171  .00000084  00000-0  23109-4 0  9991';
                lonstr2KA1='2 43925  86.3970 112.9124 0001898  98.7341 261.4069 14.34216573101091';                              
           case 8 % (�=97.5, e=0.006, h=630)
                 % ����������� ������
                 lonstr1KA1='1 39444U 13066AE  20023.53316916 +.00000248 +00000-0 +36725-4 0  9992';
                 lonstr2KA1='2 39444 097.5574 031.8478 0058731 030.6351 329.8271 14.82049256331758';             
           case 9 % (�=97.6, e=0.001, h=615)
                 % ����������� ������
                 lonstr1KA1='1 32785U 08021C   20040.99782450  .00000178  00000-0  24810-4 0  9991';
                 lonstr2KA1='2 32785  97.5638  34.7370 0012605 193.0078 167.0816 14.88414842638832';           
           case 10 % SO-50
                 % https://www.n2yo.com/satellite/?s=27607#results
                 % ������� � ���������, �����������������               
                 lonstr1KA1='1 27607U 02058C   20076.92790432 -.00000027 +00000-0 +16916-4 0  9995';
                 lonstr2KA1='2 27607 064.5547 297.2597 0059449 200.8117 159.0565 14.75615979927076';                                   
           case 11 % MERIDIAN 2               
                lonstr1KA1='1 35008U 09029A   20257.08987808  .00000306  00000-0  90674-4 0  9999';
                lonstr2KA1='2 35008  63.0139 195.0790 7283882 249.6108  23.5675  2.24593017 92652';
           case 12 % MERIDIAN 4                             
                lonstr1KA1='1 37398U 11018A   20257.85220613  .00000090  00000-0  10000-3 0  9997';
                lonstr2KA1='2 37398  63.9528  21.6889 7304490 253.3011 352.6865  2.00636009 68602';
           case 13 % % MERIDIAN 6
                lonstr1KA1='1 38995U 12063A   20256.95521867  .00000160  00000-0  00000-0 0  9991';
                lonstr2KA1='2 38995  65.3273 221.9150 6467353 260.9824  25.9153  2.00584980 57369';
           case 14 % MERIDIAN 7                      
                lonstr1KA1='1 40296U 14069A   20256.10795484  .00000086  00000-0  00000-0 0  9998';
                lonstr2KA1='2 40296  63.0034 112.1745 6935080 273.3063 345.0263  2.00620096 43011';
           case 15 % MERIDIAN 8                      
                lonstr1KA1='1 44453U 19046A   20257.47447177  .00000049  00000-0  00000+0 0  9999';
                lonstr2KA1='2 44453  62.7045 289.7770 7102336 272.7416 258.0061  2.00639516  8245';
           case 16
                lonstr1KA1='1 45254U 20015A   20257.32545561  .00000121  00000-0  00000-0 0  9999';
                lonstr2KA1='2 45254  62.8995 199.0723 7094530 295.4213 321.1634  2.00612051  4120';                 

end 
               
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

       select=1;

       switch select
           
           case 1
               
               % ��������� � 1-� ������
               latIRI=latOP(1); 
               lonIRI=lonOP(1); 
               hIRIsee=0;
               az=0;           
              
           case 2
               
               % ��������� � ������
               latIRI=latOP_0; 
               lonIRI=lonOP_0; 
               hIRIsee=0;
               az=0;

           case 3
               % ������� �� ������ ������ �� �����               
               hIRIsee=0;                             
               az=0;               
                                          
           case 4
               % ������� �� ������ ������ �� ��
               hIRIsee=0;                
               az=180;
               
           case 5
               % ������� �� ������ ������ �� ������               
               hIRIsee=0;                 
               az=90;
               
           case 6
               % ������� �� ������ ������ �� �����
               hIRIsee=0;                              
               az=-90;                              
       end

if ~OK
    disp('��� �� ������')           
    return                     
end

titleIRI=Menu{select};
 
if select > 2   

    % ����� ��������� �� �����
    distOA=1000*UnitDist; % ��������� ������� ������� (OA) �� �����, �
    Default_Answer = {num2str(distOA/UnitDist)};  
    Menu={'���������� �� �����, ��'}; 
    select=dispInput(Menu, GUIinput, '', Default_Answer);
    if isempty(select)
     disp('���������� �� ����� �� �������')                
     return
    end
    titleIRI=[titleIRI ' �� ' num2str(select(1)) ' ��'];
    distOA=select(1)*UnitDist; % �����

    [latIRI, lonIRI]=vreckon(latOP_0,lonOP_0,distOA,az);

    % �������� ��������� ����-��� ����� ������������� ����������: 
    d = vdist(latIRI,lonIRI,latOP_0,lonOP_0);
    disp(['��������� ����� ������-��� recon = ' num2str(d/UnitDist) ' ��'])
    disp(' ')
       
end
      
if Geoid       
% �������������� ������ ��� ������� ���� (������� 'egm96') � ������ ��� �����������
% ��� ����� ������������ ����� ����� ������ �� �����������,���������� � ���� geoidegm96grid.mat        

       hIRIGeoid = geoidheight( latIRI, lonIRI); % ������ ������ ��� ����������� � ����� ������������ ����� 
       hIRI=hIRIsee+hIRIGeoid; % ������������� ������ �������������� �����             
end 

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

Default_Answer = {num2str(tolx, '% 10.8f')...
                  num2str(tolf, '% 10.8f')...
                  int2str(maxiter)...
                  int2str(maxfun)};  
Menu={'�������� ��������� ���.������������ (���������)'...
      '�������� �������� ����������� � �������� ���.������������'...
      '����. ����� ��������'...
      '����. ����� ���������� �����������'}; 
select=dispInput(Menu,GUIinput, '�������� ������ ���.�����������', Default_Answer);
if isempty(select)
       disp('��������� �������� ������ ���.����������� �� �������')                
       return
end                                  

tolx=abs(select(1)); 
tolf=abs(select(2));
maxiter=abs(select(3));
maxfun=abs(select(4));

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
DT_=30*1; % �������� ���������, ���
DT_ALL= 1*min_per_day/1; % �������� ����������, ��� (1440 ��� = min_per_day = �����)  
dX=1;  % ��
dY=-1; % ��
dZ=1;   % ��
freqDownLink=7; % ���, 
                  % 7 ��� - C ��������, � GLOBALSTAR ������� 
                  % ������� �� ���� � ���� ��������� -
                  % downLink �����
                  % L - 1.6: ������� ������� �������� -
                  % downLink ��������
                  % S - 2.4: ����������� upLink ��������
% RMS_DOPPLER=100;   % ��
RMS_DOPPLER=1;   % ��
OKSType=0; % 0 - �������� ��� �� ������ �������
           % 1 - ������1 ( � �������� ��������� ������������ ������ ��������
           %               �������)
           % 2 - ������2 (� �������� ��������� ������������ ��������
           %               ������� � ������ TLE/SGP4)
min_elev_deg=10;

Default_Answer = {num2str(DT_)...
                  num2str(DT_ALL)...
                  num2str(dX)...
                  num2str(dY)...
                  num2str(dZ)...
                  num2str(freqDownLink)...                  
                  num2str(RMS_DOPPLER)...
                  int2str(OKSType)...
                  num2str(min_elev_deg)};
                    
Menu={'�������� ���������, ���'...
      '�������� ����������, ��� (1440 ��� = �����)'...
      '�������� �� X �� TLE/SGP4, ��'...
      '�������� �� Y �� TLE/SGP4, ��'...
      '�������� �� Z �� TLE/SGP4, ��'...
      '������� ����, ���'...
      '��� ���������, ��'...
      ['����� ������ ����������' sprintf('\n') '0-OKS ��������, 1 - ������1, 2 -������2']...
      '���. ���� ������� ���������'};  
      
select=dispInput(Menu,GUIinput, '��������� ������', Default_Answer);
if isempty(select)
       disp('��������� ������ �� �������')                
       return
end

DT_=abs(select(1)); % ���
DT_ALL=abs(select(2)); % ���
dXYZ=UnitDist*[select(3) select(4) select(5) ]; % ����� 

freqDownLink=select(6)*OneGHz; % ��
RMS_DOPPLER=select(7); % ��
OKSType=select(8);
min_elev_deg=select(9);

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
% %% �������� ������� ������ ���������� �� ����� �� ���������� �������� �������
% 
% Kalman=boolean(0);
% if OKSType==1 || OKSType==2
%   Kalman=boolean(1);  
%   KalmanType=OKSType;
% end
% 
%  if Kalman
%      
%      if RMS_DOPPLER ~=0
%         RMS_KALMAN=RMS_DOPPLER*c/freqDownLink; % ����������� ��������� ������� �������     
%      else         
%         RMS_KALMAN=10; % ��        
%      end
%      
%      % ������ ��������� ������� =
%      % S=[X Y Z Vx Vy Vz]
%      S=zeros(1,6);
%      
%      D=10^4; % ����������� ���������� ���������� �������� ������������� ��������� (� ������)      
%      V=10^4; % �/� 
%      
%      accelX=10^5; % ��������� ��������� �� ���������� X
%      accelY=10^5; % ��������� ��������� �� ���������� Y
%      accelZ=10^5; % ��������� ��������� �� ���������� Z     
%  
%      dT=DT_; % ���
%      % ������� �������� � ����� ��������� (���������� �������� ���� dx=x+v*dt):
%      A = [1   0  0   dT  0   0
%           0   1  0   0   dT  0 
%           0   0  1   0   0   dT
%           0   0  0   1   0   0
%           0   0  0   0   1   0
%           0   0  0   0   0   1];
% 
%     % ������������� ������� ���������� ���� �������� � ���������� �����������
%     % ��������� �� ����������� X � Y:  
%       U=[dT^4/4   dT^3/2
%          dT^3/2   dT^2];
% 
%     % ������� ��������� ��� ����������� ���� ��������:   
%       accel=[accelX     0        0
%                0      accelY     0
%                0        0      accelZ];
% 
%     % ������������� ������� ���� ����������� ��������, ������������ �����������, 
%     % ���������� �������� accel:         
%     Q=[accel*U(1,1)    accel*U(1,2);...
%        accel*U(2,1)    accel*U(2,2)];       
% 
%     if KalmanType ==1
%         % � �������� ��������� ������������ ������ �������� �������
%          R=RMS_KALMAN .^ 2*diag(ones(1,NumPosts));                
%     else 
%         % � �������� ��������� ������������ �������� ������� � ������
%         % TLE/SGP4
%          Q=NumPosts + numel(S); % ����� ��������� � ������� ���������
%          R=diag(ones(1,Q));
%          for j=1:NumPosts
%              R(j,:)=RMS_KALMAN .^ 2*R(j,:);
%          end
% 
%          dxyz=10^2;
%          for j=NumPosts+1:NumPosts+1+2
%              R(j,:)=dxyz^2*R(j,:);
%          end
% 
%          dVxyz=10^1;
%          for j=NumPosts+1+2+1:Q 
%              R(j,:)=dVxyz^2*R(j,:);
%          end         
%     end               
%      
%      NumVitkov=length(Vitok); % ����� ������� ������
%      tru_est=sgp4;
%      tru_V_est=sgp4_V;
%      b=nan(1,NumVitkov);
%      q=nan(1,NumVitkov);
%      
%    for v=1:NumVitkov
%        
%      ind=Vitok{v};
%      M=length(ind); % ����� ��������� �� �����      
%           
%      % ��������� �������� ������� ��������� ������� =
%      S=[sgp4(ind(1),:) sgp4_V(ind(1),:)];       
%      
%     % ��������� �������� �������������� ������� 
%     % ����������� ����������  ��������� ������� ���������
%      P=[D^2  0   0   0   0   0
%         0  D^2  0   0   0   0
%         0   0  D^2  0   0   0 
%         0   0   0  V^2  0   0
%         0   0   0   0  V^2  0
%         0   0   0   0   0  V^2];         
% 
%     for t = 1:M
%          % ���������� ������� ��������� ���������� �������� ������� (UTF)
%          % � �������� ��������� ������������ �� ������������ ������� f, � 
%          % ���������� �������� v (tru_V_PROECT), ����� � �������� �� ������� 
%          % �������� ������� ���� freqDownlink. ��� ��������, ��������� 
%          % f=(v/c)*freqDownlink, ��� c - �������� �����.
% 
%          if KalmanType==1         
%             u=meas_V_PROECT(ind(t),:);     
%          else
%             u=[meas_V_PROECT(ind(t),:) sgp4(ind(t),:) sgp4_V(ind(t),:)];
%          end
% 
%           [S, P] = UKFupdate(u,xyzOP,S,P, A, Q, R, KalmanType);
% 
%           tru_est(ind(t),:)=S(1:3);
%           tru_V_est(ind(t),:)=S(4:6);      
%     
%     end     
%     
%       err=tru_est(ind,:)-tru(ind,:);
%       err=err';
%       err_=sum(err.^2); % ������������ �� ��������, � ������� ���������� x y z
%       b(v)=sqrt(mean(err_));
%       
%       err=tru_V_est(ind,:)-tru_V(ind,:);
%       err=err';
%       err_=sum(err.^2); % ������������ �� ��������, � ������� ���������� x y z
%       q(v)=sqrt(mean(err_)); 
%       
%    end
% 
%   pos_err_Kalman=mean(b);
%   vel_err_Kalman=mean(q);
% % 
%  end

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
     % S=[X Y Z]
     S=zeros(1,6);
     
     D=10^-4; % ����������� ���������� ���������� �������� ������������� ��������� (� ������)      
     V=10^-4; % �/� 
     
%      accelX=10^5; % ��������� ��������� �� ���������� X
%      accelY=10^5; % ��������� ��������� �� ���������� Y
%      accelZ=10^5; % ��������� ��������� �� ���������� Z    
     accelX=10^1; % ��������� ��������� �� ���������� X
     accelY=10^1; % ��������� ��������� �� ���������� Y
     accelZ=10^1; % ��������� ��������� �� ���������� Z   
 
     dT=DT_; % ���
     % ������� �������� � ����� ��������� (���������� �������� ���� dx=x+v*dt):
%      A = [1   0  0   dT  0   0
%           0   1  0   0   dT  0 
%           0   0  1   0   0   dT
%           0   0  0   1   0   0
%           0   0  0   0   1   0
%           0   0  0   0   0   1];

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
    else % KalmanType ==2
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

         dVxyz=10^-1;
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
%             u=[meas_V_PROECT(ind(t),:) sgp4(ind(t),:)];
         end

          [S, P] = UKFupdate(u,xyzOP,S,P, A, Q, R, KalmanType);
          
          S(4:6)=sgp4_V(ind(t),:)';

          tru_est(ind(t),:)=S(1:3);
          tru_V_est(ind(t),:)=S(4:6);      
    
    end     

     ind=ind(end);

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
     % ������������: 
%      Dop=tru_DOPPLER(ind);
    
%           [opt,~, exitflag]=SOMP_Doppler(SAT,...
%                                     SAT_V,...                          
%                                     Dop,...
%                                     opt,...
%                                     h0,...
%                                     tolx,...
%                                     tolf,...
%                                     maxiter,...
%                                     maxfun,...
%                                     freqDownLink); % freqDownLink ���
%                                     ���������� ��.��������� ������� SOMP_Doppler

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

% %% �������� C��� ����� �� ������� ����� ��S �� ������ ����������� ������������� ����������� �������
% 
% %  opt=[latPost, lonPost] ; % ��������� ������ ��� ������ �������� ������������� ����������� 
%  opt=[latOP(PostNumber), lonOP(PostNumber)] ; % ��������� ������ ��� ������ �������� ������������� ����������� 
%  h0=0; % ������������� ������ ��� ������ 
%  
%  NumVitkov=length(Vitok); % ����� ������� ������
%  % ��������� ���:
%  latPost_Doppler=nan(NumVitkov,1);
%  lonPost_Doppler=nan(NumVitkov,1);
% 
%  for v=1:NumVitkov
% 
%      ind=Vitok{v};
%      ind=ind(1); % ��������� ���������
%      % ������� ���������� sgp4: 
% %      SAT=sgp4(ind,:); 
% %      SAT_V=sgp4_V(ind,:); 
% 
%      % ������ ���������� ����� ��S: 
%      SAT=tru_est_OKS(ind,:);
%      SAT_V=tru_V_est_OKS(ind,:);
% 
%      % ��������� ������� �� �������� ����������:
%      %Dop=tru_V_PROECT(ind,PostNumber); % �������� ������ 
%      Dop=meas_V_PROECT(ind,PostNumber); % ��������� (� �������)
%      % � �������� ��������� ������������ �� 
%      % ������������ ������� f (tru_DOPPLER), �  ���������� �������� v, 
%      % ����� � �������� SOMP_Doppler �� ������� ������� ������� ����
%      % freqDownLink. ��� ��������, ��������� f=(v/c)*freqDownLink, ��� c - 
%      % �������� �����
%           [opt,~, exitflag]=SOMP_Doppler(SAT,...
%                                          SAT_V,...                          
%                                          Dop,...
%                                          opt,...
%                                          h0,...
%                                          tolx,...
%                                          tolf,...
%                                          maxiter,...
%                                          maxfun);  
%      % ������������: 
% %      Dop=tru_DOPPLER(ind);
%     
% %           [opt,~, exitflag]=SOMP_Doppler(SAT,...
% %                                     SAT_V,...                          
% %                                     Dop,...
% %                                     opt,...
% %                                     h0,...
% %                                     tolx,...
% %                                     tolf,...
% %                                     maxiter,...
% %                                     maxfun,...
% %                                     freqDownLink); % freqDownLink ���
% %                                     ���������� ��.��������� ������� SOMP_Doppler
% 
%      if ~ exitflag
%           MinFlagStr='���� �� ������� ����� ��S: min ������������� ����������� ������� �� ������ !!!';
%           disp(MinFlagStr)
%      end
% 
%     latPost_Doppler(v)=opt(1); % ������ ������ �����
%     lonPost_Doppler(v)=opt(2); % ������ ������� �����    
%  
% end
% 
% %%  ������ ������ C��� ����� �� ������� ����� ��S
% 
% errPost_Doppler=nan(NumVitkov,1);
% 
% if ~ isempty(Vitok{1})
%     
%     for v=1:NumVitkov
% 
%         % ����� ��������� ����������:
%         h0=0;
%         [xIRIest, yIRIest, zIRIest]=llh2xyz( latPost_Doppler(v), lonPost_Doppler(v), h0 );             
%         errPost_Doppler(v)=sqrt((xIRIest-xyzOP(1,PostNumber)).^2 + (yIRIest-xyzOP(2,PostNumber)).^2+(zIRIest-xyzOP(3,PostNumber)).^2); 
% 
%         if boolean(errPost_Doppler(v))
%             % ����� ������������� ���������� ����������:                               
%              errPost_Doppler(v)= vdist(latPost_Doppler(v),lonPost_Doppler(v),latOP(PostNumber),lonOP(PostNumber));
%         end
%     end
% end
% 
% meanErr_Doppler_afterOKS=nanmean(errPost_Doppler);
% disp(['������� �������� ������ ���� ����� �� ������� ����� ��S = ' num2str(meanErr_Doppler_afterOKS/UnitDist)  ' �� ' ])
% disp(' ')
    
%% ������� ����� ���������

plotDoppler=boolean(1); 
plotSGP4=boolean(0);     % ��������� ��������������� ������ �� ��������� SGP4 (������, �������, ������ �� �������)  
plotSGP4_TRU=boolean(0); % ��������� ���������� SGP4 � TRU

incloKA1=satrecKA1.inclo*180/pi; % ���������� (���������) ������ ��� � ��������
eccoKA1=satrecKA1.ecco;   % ���������������� ������ ���

SAT=Norad2name( NORAD_ID_KA1);
% titl=['�������� ������������� SGP4' sprintf('\n')...
%       '���: ' SAT ', ���������� = ' num2str(incloKA1) ' ����, ���������������� =  '...
%        num2str(eccoKA1), sprintf('\n')...
%       '���� ' titlePost ' , ��� ' titleIRI  sprintf('\n')...
%       '��. �� ������ ������ ��� = ' num2str(meanErr/UnitDist) ' ��' sprintf('\n')...
%       '����� UTC �� TLE : ' utc0_KA1  sprintf('\n')...
%       '���� ������������� : ' datestr(now)];

titl=['�������� ������������� SGP4' sprintf('\n')...
      '���: ' SAT ', ���������� = ' num2str(incloKA1) ' ����, ���������������� =  '...
       num2str(eccoKA1), sprintf('\n')...
      '���� ' titlePost  sprintf('\n')...
      '����� UTC �� TLE : ' utc0_KA1  sprintf('\n')...
      '���� ������������� : ' datestr(now)];
  
% ����� ������������ �����, ��� � �������������� ���������

if NumPosts==1
        
        figure(1)
        hold on   
        plot(lonPost, latPost,'.b','MarkerSize',20)         
        lon=lonSAT1;
        lat=latSAT1;       
        plot(lon,lat,'.r','MarkerSize',10) 
%         plot(lonIRI,latIRI,'Marker','o', 'Color','g','MarkerSize',10,'LineStyle','non','LineWidth',2)
%         plot(lonIRI_DDM,latIRI_DDM,'.k','MarkerSize',10)       
        if ~ isempty(Vitok{1})
%             legend('����',['��������� ������� ������ � ���' sprintf('\n')...
%                                         '�������������� ����������'], '���','������ ��� ���')  
            legend('����',['��������� ������� ������' sprintf('\n')...
                                        '�������������� ����������']) 
        else
%             legend('����', '���')  
          legend('����')  
        end
%         plot_google_map
        xlabel('�������, ����')
        ylabel('������, ����')
        xlim([-180 +180]) 
        grid on        
        title(titl)
        
else % NumPosts > 1
    
        figure(1)
        hold on   
        plot(lonOP_0, latOP_0,'.b','MarkerSize',20)         
        plot(lonOP,latOP,'xb','MarkerSize',10)        
        lon=lonSAT1;
        lat=latSAT1;       
        plot(lon,lat,'.r','MarkerSize',10) 
%         plot(lonIRI,latIRI,'Marker','o', 'Color','g','MarkerSize',10,'LineStyle','non','LineWidth',2)
%         plot(lonIRI_DDM,latIRI_DDM,'.k','MarkerSize',10)         
        if ~ isempty(Vitok{1})
%             legend('����� ������','�����',['��������� ������� ������� � ���' sprintf('\n')...
%                                         '�������������� ����������'], '���','������ ��� ���') 
            legend('����� ������','�����',['��������� ������� �������' sprintf('\n')...
                                        '�������������� ����������']) 
        else
%             legend('����� ������', '�����', '���') 
            legend('����� ������', '�����')  
        end
%         plot_google_map
        xlabel('�������, ����')
        ylabel('������, ����')
        xlim([-180 +180]) 
        grid on        
        title(titl)            
        
end
        
        figure
        plot(dt,losCUBIK,'Marker','+', 'Color','r')
        grid on
        xlabel('�����, ���')
%         title('���������� ��������� ��� ������ � ���')
        title('���������� ��������� ��� �������')
        
%         figure
%         stem(err_DDM/UnitDist,'Marker','o', 'Color','r')
%         grid on
%         xlabel('����� �������� �����')
%         ylabel('��')        
%         title(['������ ��� ������� ���' sprintf('\n')...
%                '�������� ��������� = ' num2str(DT_) 'c'])         
%         
%         figure
%         stem(NumdiffD,'Marker','o', 'Color','r')
%         grid on
%         xlabel('����� �������� �����')
%         title(['����� ���������' sprintf('\n')...
%                '�������� ��������� = ' num2str(DT_) 'c'])
        
       
% ������ ������� ��� ����� ���, ����������� ������ �������
     if plotDoppler   
    
        titl=['�������: ' SAT sprintf('\n')...
              ' ���� ������������� ' datestr(now)];
          
%         % 5  ���������� ������� ���
%         V=sqrt(sgp4_V(:,1).^2+sgp4_V(:,2).^2+sgp4_V(:,2).^2);            
%         figure
%         plot(dt,V/UnitDist,'Marker','.','LineStyle','none')
%         xlabel('�����, ���')
%         ylabel('��/c')        
%         title('���������� ������� ���')
%         grid on        
        
        % 6 ������� ������� ����
        figure
        plot(dt, sgp4_DOPPLER(:,1)/UnitFreq,'Marker','o','MarkerSize',4,'LineStyle','none')
        hold on        
        plot(dt, tru_DOPPLER(:,1)/UnitFreq,'Marker','s','MarkerSize',4,'Color','r','LineStyle','none')      
        plot(dt, meas_DOPPLER(:,1)/UnitFreq,'Marker','x','MarkerSize',4,'Color','k','LineStyle','none')         
        legend('SGP4','TRUE','���������')
        xlabel('�����, ���')
        ylabel('���')
        title([titl sprintf('\n') '������� ������� ����'])
        grid on         
        
        
%         % 7 ����.������� ������� ����
%         figure
%         plot(dt, sgp4_DOPPLER(:,1)/freqDownLink,'Marker','o','MarkerSize',4,'LineStyle','none')
%         hold on        
%         plot(dt, tru_DOPPLER(:,1)/freqDownLink,'Marker','s','MarkerSize',4,'Color','r','LineStyle','none')         
%         legend('SGP4','TRUE')
%         xlabel('�����, ���')   
%         title([titl sprintf('\n') '����.������� ������� ���� (v/c)'])
%         grid on
        
        % 8 TRU-SGP4 ������� ������� ����
        figure
        plot(dt, (tru_DOPPLER(:,1)-sgp4_DOPPLER(:,1))/UnitFreq,'Marker','*','MarkerSize',4,'LineStyle','none','Color','k')
        hold on            
        xlabel('�����, ���')
        ylabel('���')
        title([titl sprintf('\n') '�������=TRU-SGP4 ������� ������� ����'])
        grid on         
%         
%         % 9 TRU-SGP4 ����.������� ������� ����
%         figure
%         plot(dt, (tru_DOPPLER(:,1)-sgp4_DOPPLER(:,1))/freqDownLink,'Marker','*','MarkerSize',4,'LineStyle','none','Color','k')
%         hold on            
%         xlabel('�����, ���')
%         title([titl sprintf('\n') 'TRU-SGP4 ����.������� ������� ����'])
%         grid on          
  
        % 10 ����������� ������� ������� ����
         figure        
         d1=diff(sgp4_DOPPLER(:,1));      
         td=diff(dt);
         dd1=d1./td;     
         hold on
         plot(dt(2:end), (dd1/UnitFreq)/sec_per_mim,'Marker','o','MarkerSize',4,'LineStyle','none')          
         d1=diff(tru_DOPPLER(:,1));      
         td=diff(dt);
         dd1=d1./td;     
         hold on
         plot(dt(2:end), (dd1/UnitFreq)/sec_per_mim,'Marker','s','MarkerSize',4,'Color','r','LineStyle','none')                  
         legend('SGP4','TRUE')
         title([titl sprintf('\n') '����������� ������� ������� ����'])
         xlabel('�����, ���')
         ylabel('��� / c')    
         grid on           
       
%          % 11 ���������� ��������
%          figure         
%          dt=dt/1; % ���
%          hold on
%          plot(dt, sgp4_V_PROECT(:,1)/UnitDist,'Marker','o','MarkerSize',4,'LineStyle','none')
%          plot(dt, tru_V_PROECT(:,1)/UnitDist,'Marker','s','MarkerSize',4,'Color','r','LineStyle','none')          
%          legend('SGP4','TRUE')
%          title([titl sprintf('\n') '���������� �������� ���'])
%          xlabel('�����, ���')
%          ylabel('��/�')    
%          grid on
% 
%          % 12
%          % ����������� ���������� ��������
%          figure         
%          d=diff(sgp4_V_PROECT(:,1));
%          td=diff(dt);
%          dd=d./td;
%          plot(dt(2:end), dd/UnitDist,'Marker','o','MarkerSize',4,'LineStyle','none')
%          hold on
%          d=diff(tru_V_PROECT(:,1));
%          td=diff(dt);
%          dd=d./td;
%          plot(dt(2:end), dd/UnitDist,'Marker','s','MarkerSize',4,'Color','r','LineStyle','none')                
%          legend('SGP4','TRUE')
%          title([titl sprintf('\n') '����������� ���������� �������� ���'])
%          xlabel('�����, ���')
%          ylabel('(��/c) / ���')    
%          grid on          
     end
    
        if Kalman   
            titl=['������ ������� = '  num2str(pos_err_Kalman) ' �' sprintf('\n')...
                  '������ c������� = ' num2str(vel_err_Kalman) ' �/c' sprintf('\n')];
            
            figure 
            plot(tru_est(:,1),'r','Marker','o','MarkerSize',4,'LineStyle','none')
            hold on;
            plot(tru(:,1),'Marker','s','MarkerSize',4,'LineStyle','none')          
            title([titl 'X'])
            grid on
            plot(sgp4(:,1),'k','Marker','+','MarkerSize',4,'LineStyle','none')             
            legend('������ ��� �������','�������� ����������','TLE/SGP4')                     

            figure 
            plot(tru_est(:,2),'r','Marker','o','MarkerSize',4,'LineStyle','none')
            hold on;
            plot(tru(:,2),'Marker','s','MarkerSize',4,'LineStyle','none')
            title([titl 'Y'])
            grid on
            plot(sgp4(:,2),'k','Marker','+','MarkerSize',4,'LineStyle','none')             
            legend('������ ��� �������','�������� ����������','TLE/SGP4')         

            figure 
            plot(tru_est(:,3),'r','Marker','o','MarkerSize',4,'LineStyle','none')
            hold on;
            plot(tru(:,3),'Marker','s','MarkerSize',4,'LineStyle','none')        
            title([titl 'Z'])
            grid on            
            plot(sgp4(:,3),'k','Marker','+','MarkerSize',4,'LineStyle','none')             
            legend('������ ��� �������','�������� ����������','TLE/SGP4')                     

            figure 
            plot(tru_V_est(:,1),'r','Marker','o','MarkerSize',4,'LineStyle','none')
            hold on;
            plot(tru_V(:,1),'Marker','s','MarkerSize',4,'LineStyle','none')
            title([titl 'Vx'])
            grid on
            plot(sgp4_V(:,1),'k','Marker','+','MarkerSize',4,'LineStyle','none')             
            legend('������ ��� �������','�������� ����������','TLE/SGP4')            

            figure 
            plot(tru_V_est(:,2),'r','Marker','o','MarkerSize',4,'LineStyle','none')
            hold on;
            plot(tru_V(:,2),'Marker','s','MarkerSize',4,'LineStyle','none')
            title([titl 'Vy'])
            grid on
            plot(sgp4_V(:,2),'k','Marker','+','MarkerSize',4,'LineStyle','none')             
            legend('������ ��� �������','�������� ����������','TLE/SGP4')          

            figure 
            plot(tru_V_est(:,3),'r','Marker','o','MarkerSize',4,'LineStyle','none')
            hold on;
            plot(tru_V(:,3),'Marker','s','MarkerSize',4,'LineStyle','none')
            title([titl 'Vz'])
            grid on
            plot(sgp4_V(:,3),'k','Marker','+','MarkerSize',4,'LineStyle','none')             
            legend('������ ��� �������','�������� ����������','TLE/SGP4')         
        
        end
        
        if OKS   
            titl=['������ ��S ������� = '  num2str(pos_err_OKS) ' �' sprintf('\n')...
                  '������ ��S c������� = ' num2str(vel_err_OKS) ' �/c' sprintf('\n')];
            
            figure 
            plot(tru_est_OKS(:,1),'r','Marker','o','MarkerSize',4,'LineStyle','none')
            hold on;
            plot(tru(:,1),'Marker','s','MarkerSize',4,'LineStyle','none')          
            title([titl 'X'])
            grid on
            plot(sgp4(:,1),'k','Marker','+','MarkerSize',4,'LineStyle','none')             
            legend('������ ��S','�������� ����������','TLE/SGP4')                     

            figure 
            plot(tru_est_OKS(:,2),'r','Marker','o','MarkerSize',4,'LineStyle','none')
            hold on;
            plot(tru(:,2),'Marker','s','MarkerSize',4,'LineStyle','none')
            title([titl 'Y'])
            grid on
            plot(sgp4(:,2),'k','Marker','+','MarkerSize',4,'LineStyle','none')             
            legend('������ ��S','�������� ����������','TLE/SGP4')         

            figure 
            plot(tru_est_OKS(:,3),'r','Marker','o','MarkerSize',4,'LineStyle','none')
            hold on;
            plot(tru(:,3),'Marker','s','MarkerSize',4,'LineStyle','none')        
            title([titl 'Z'])
            grid on            
            plot(sgp4(:,3),'k','Marker','+','MarkerSize',4,'LineStyle','none')             
            legend('������ ��S','�������� ����������','TLE/SGP4')                     

            figure 
            plot(tru_V_est_OKS(:,1),'r','Marker','o','MarkerSize',4,'LineStyle','none')
            hold on;
            plot(tru_V(:,1),'Marker','s','MarkerSize',4,'LineStyle','none')
            title([titl 'Vx'])
            grid on
            plot(sgp4_V(:,1),'k','Marker','+','MarkerSize',4,'LineStyle','none')             
            legend('������ ��S','�������� ����������','TLE/SGP4')            

            figure 
            plot(tru_V_est_OKS(:,2),'r','Marker','o','MarkerSize',4,'LineStyle','none')
            hold on;
            plot(tru_V(:,2),'Marker','s','MarkerSize',4,'LineStyle','none')
            title([titl 'Vy'])
            grid on
            plot(sgp4_V(:,2),'k','Marker','+','MarkerSize',4,'LineStyle','none')             
            legend('������ ��S','�������� ����������','TLE/SGP4')          

            figure 
            plot(tru_V_est_OKS(:,3),'r','Marker','o','MarkerSize',4,'LineStyle','none')
            hold on;
            plot(tru_V(:,3),'Marker','s','MarkerSize',4,'LineStyle','none')
            title([titl 'Vz'])
            grid on
            plot(sgp4_V(:,3),'k','Marker','+','MarkerSize',4,'LineStyle','none')             
            legend('������ ��S','�������� ����������','TLE/SGP4')         
        
        end        
        
    if plotSGP4_TRU
        
        ind=1;
        tl='X';
        figure
        hold on
        plot(dt, SGP4(:,ind))
        plot(dt, TRU(:,ind),'r')
        legend('SGP4','TRU')
        grid on
        xlabel('�����, ���')    
        title(tl)
        
        ind=2;
        tl='Y';
        figure
        hold on
        plot(dt, SGP4(:,ind))
        plot(dt, TRU(:,ind),'r')
        legend('SGP4','TRU')
        grid on
        xlabel('�����, ���')    
        title(tl)  
        
        ind=3;
        tl='Z';
        figure
        hold on
        plot(dt, SGP4(:,ind))
        plot(dt, TRU(:,ind),'r')
        legend('SGP4','TRU')
        grid on
        xlabel('�����, ���')    
        title(tl)          
        
    end

   
    if plotSGP4 % ������� ��������������� ���������� ���
          
       % ��������������� SPG4 �� ��������� + - = ���� �� ����� TLE � ����� 1 ���
       % ������ ���:
       
        [predictT_KA1,predictLAT_KA1,predictLON_KA1,predictH_KA1] = SPG4twoLLH(lonstr1KA1, lonstr2KA1);       
        SAT=Norad2name( NORAD_ID_KA1);  

           titl_SGP4=['�������� ������������� SGP4' sprintf('\n')...
                 '���: ' SAT ', ���������� = ' num2str(incloKA1) ' ����, ���������������� =  '...
                  num2str(eccoKA1), sprintf('\n')...                 
                 '����� UTC �� TLE : ' utc0_KA1  sprintf('\n')...
                 '���� ������������� : ' datestr(now)];
             
            a=min(adrKA1);
            b=max(adrKA1);  

            figure
            plot(predictT_KA1{1},predictLON_KA1{1}) 
            grid on
            xlabel('�����, ���')
            ylabel('�������, ����')    
            title(titl_SGP4)
            hold on

            c=min(predictLON_KA1{1});
            d=max(predictLON_KA1{1});
            e=linspace(c,d);
            plot(predictT_KA1{1}(a),e)
            plot(predictT_KA1{1}(b),e)

            figure
            plot(predictT_KA1{1},predictLAT_KA1{1}) 
            grid on
            xlabel('�����, ���')
            ylabel('������, ����')    
            title(titl_SGP4)
            hold on            
            
            c=min(predictLAT_KA1{1});
            d=max(predictLAT_KA1{1});
            e=linspace(c,d);
            plot(predictT_KA1{1}(a),e)
            plot(predictT_KA1{1}(b),e)
            
            figure
            plot(predictT_KA1{1},predictH_KA1{1}/UnitDist) 
            grid on
            xlabel('�����, ���')
            ylabel('������, ��')    
            title(titl_SGP4)
            hold on            
            
            c=min(predictH_KA1{1}/UnitDist);
            d=max(predictH_KA1{1}/UnitDist);
            e=linspace(c,d);
            plot(predictT_KA1{1}(a),e)
            plot(predictT_KA1{1}(b),e)                                       

    end
    
    % ��� ������� � ����� ����:
    f=findobj(0,'IntegerHandle','on');
    set(f,'WindowStyle', 'docked')

end % GEOLOCATION_LEO_DOPPLER_N_POST
