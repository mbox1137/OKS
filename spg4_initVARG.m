function [satrec, jd0, utc0, NORAD_ID, r, v ] = spg4_initVARG(lonstr1, lonstr2)
%
% Инициализация процедуры прогнозирования по алгоритму SGP4
%
% var@mail.spbstu.ru, декабрь 2016

global UnitDist

whichconst=84;     % параметры модели Земли в виде эллипсоида WGS-84, который используется в GPS 
% whichconst= 72;
global  opsmode                
opsmode= 'a'; % % возможно: a, e   

%Режим прогнозирования:
typerun='c';   % c = prpagates at 20 min timesteps frm
               %     one day before epoch to one day after
               %     возможно: c, v, m 
               % (на данный момент программа работает только в режиме c)
               
%Режим интерпретации времени:
typeinput = 'e'; % e - время в виде эпохи
                 %    возможно: e, m, d

satrec = twoline2rv_VARG( whichconst, ...
                       lonstr1, lonstr2, typerun, typeinput);
                   
% satrec - орбитальные параметры (в виде структуры данных) в эпоху, указанную в TLE
% startmfe - начальное время относительно эпохи, указанной в TLE, в минутах
% stopmfe - конечное время относительно эпохи, указанной в TLE, в минутах
% deltami - интервал, в минутах                      
                   
disp('номер спутника (NORAD CAT ID): ')
NORAD_ID=satrec.satnum;
fprintf(1,' %d\n', NORAD_ID); 

jd = satrec.jdsatepoch; % число дней в юлианской дате
jd0=jd;
[year,mon,day,hr,minute,sec] = invjday ( jd ); % дата в системе времени UTC
disp('UTC эпохи из TLE : ')
dtm = datenum(year,mon,day,hr,minute,sec);
utc0=datestr(dtm);
disp(utc0)
disp(' ') 
        
% Начальные значения вектора состояния rv в эпоху, указанную в TLE:
                   
[satrec,r,v] = sgp4 (satrec,  0.0);

r=r*UnitDist; %  преобразование r в метры
v=v*UnitDist; %  преобразование v в м/с     

end