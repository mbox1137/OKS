function [azIRI,elevIRI,indexLos,distIRI]=peleng(refLat,refLong,refH,Lat,Long,H,min_elev_deg)
%
% var@mail.spbstu.ru, V.A. Vargauzin, ������ 2020, ������ 2021
%**************************************************************************

if nargin < 7 
    min_elev_deg =0;
end
%% ������������� �������� ��������       

% ENU:
% refLat=latPost;
% refLong=lonPost;
% refH=hPost;

[xIRI, yIRI, zIRI]=llh2xyz( Lat, Long, H );

% ��������� ��������� ���������� ��� � ENU 
[xEastIRI,yNorthIRI,zUpIRI] = xyz2enu(refLat,refLong,refH, xIRI, yIRI, zIRI);

% % ��������e ����-���
distIRI=sqrt(xEastIRI.^2 + yNorthIRI.^2 + zUpIRI.^2);

% �������� ��������� �� ��������� EN
distIRIprojectENU=sqrt(xEastIRI.^2+yNorthIRI.^2);

% ���� �����
elevIRI=atand(zUpIRI./distIRIprojectENU); % �������

% �������� ������ ���������
indexLos=elevIRI > min_elev_deg;

% ������
azIRI=asind(xEastIRI./distIRIprojectENU); % �������

for j=1:length(azIRI)
    
if xEastIRI(j) > 0 && yNorthIRI(j) < 0
     azIRI(j)=180-azIRI(j);
elseif xEastIRI(j) < 0 && yNorthIRI(j) < 0
     azIRI(j)=abs(azIRI(j))-180;
end

if azIRI(j) < 0
    azIRI(j)=azIRI(j)+360;
end

end

end