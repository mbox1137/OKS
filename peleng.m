function [azIRI,elevIRI,indexLos,distIRI]=peleng(refLat,refLong,refH,Lat,Long,H,min_elev_deg)
%
% var@mail.spbstu.ru, V.A. Vargauzin, январь 2020, январь 2021
%**************************************************************************

if nargin < 7 
    min_elev_deg =0;
end
%% Моделирование истинных пеленгов       

% ENU:
% refLat=latPost;
% refLong=lonPost;
% refH=hPost;

[xIRI, yIRI, zIRI]=llh2xyz( Lat, Long, H );

% локальные декартовы координаты ИРИ в ENU 
[xEastIRI,yNorthIRI,zUpIRI] = xyz2enu(refLat,refLong,refH, xIRI, yIRI, zIRI);

% % растояниe Пост-ИРИ
distIRI=sqrt(xEastIRI.^2 + yNorthIRI.^2 + zUpIRI.^2);

% Проекция растояния на плоскость EN
distIRIprojectENU=sqrt(xEastIRI.^2+yNorthIRI.^2);

% Уогл места
elevIRI=atand(zUpIRI./distIRIprojectENU); % градусы

% Контроль прямой видимосим
indexLos=elevIRI > min_elev_deg;

% Азимут
azIRI=asind(xEastIRI./distIRIprojectENU); % градусы

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