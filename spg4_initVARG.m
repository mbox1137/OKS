function [satrec, jd0, utc0, NORAD_ID, r, v ] = spg4_initVARG(lonstr1, lonstr2)
%
% ������������� ��������� ��������������� �� ��������� SGP4
%
% var@mail.spbstu.ru, ������� 2016

global UnitDist

whichconst=84;     % ��������� ������ ����� � ���� ���������� WGS-84, ������� ������������ � GPS 
% whichconst= 72;
global  opsmode                
opsmode= 'a'; % % ��������: a, e   

%����� ���������������:
typerun='c';   % c = prpagates at 20 min timesteps frm
               %     one day before epoch to one day after
               %     ��������: c, v, m 
               % (�� ������ ������ ��������� �������� ������ � ������ c)
               
%����� ������������� �������:
typeinput = 'e'; % e - ����� � ���� �����
                 %    ��������: e, m, d

satrec = twoline2rv_VARG( whichconst, ...
                       lonstr1, lonstr2, typerun, typeinput);
                   
% satrec - ����������� ��������� (� ���� ��������� ������) � �����, ��������� � TLE
% startmfe - ��������� ����� ������������ �����, ��������� � TLE, � �������
% stopmfe - �������� ����� ������������ �����, ��������� � TLE, � �������
% deltami - ��������, � �������                      
                   
disp('����� �������� (NORAD CAT ID): ')
NORAD_ID=satrec.satnum;
fprintf(1,' %d\n', NORAD_ID); 

jd = satrec.jdsatepoch; % ����� ���� � ��������� ����
jd0=jd;
[year,mon,day,hr,minute,sec] = invjday ( jd ); % ���� � ������� ������� UTC
disp('UTC ����� �� TLE : ')
dtm = datenum(year,mon,day,hr,minute,sec);
utc0=datestr(dtm);
disp(utc0)
disp(' ') 
        
% ��������� �������� ������� ��������� rv � �����, ��������� � TLE:
                   
[satrec,r,v] = sgp4 (satrec,  0.0);

r=r*UnitDist; %  �������������� r � �����
v=v*UnitDist; %  �������������� v � �/�     

end