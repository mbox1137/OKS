function [loshadka] =prm(params)
loshadka = GEOLOCATION_LEO_DOPPLER_N_POST_2023
function []=GEOLOCATION_LEO_DOPPLER_N_POST_2023(varargin)
   
%**************************************************************************
%
% Функция GEOLOCATION_LEO_DOPPLER_N_POST_2023 моделирует алгоритм определения 
% координат спутника (OKS) по результат измерений частоты Доплера N постами
% (N_POST). Число постов в модели может быть любым.
%
% В качестве спутника используется низкоорбитальный космический аппарат (НКА, LEO). 
% Каждый пост(приемное устройство), в принципе, может сделать K измерений частоты Доплера "вниз" от маяка на
% борту НКА. Измерения производся на интервале видимости НКА постами. 
% За сутки НКА совершает несколько витков вокруг Земли, и, соответсвенно
% имеется несколько интервалов видимости.
%
% В соотвествии с техническим заданием в модели используется лишь ОДНО 
% (первое) ИЗМЕРЕНИЕ на каждом витке.
%
% Для моделирования истинной (TRU) траектории движения спутника используется 
% модель со смещением от "грубой" траектории TLE/SGP4. При этом сама 
% "грубая" траектория использутся в качестве начальных условий для
% поиска минимума функционала невязки. При этом модель демонстрирует 
% отличие в Доплере при смешении истинной траектории от грубой 
% (за счёт различных проекций вектора скорости НКА на направление на Пост).
%
% Модель использует 4 вспомогательные функции: 
% - dopplerSAT1_downLink 
% - OKS_Doppler
% - SOMP_Doppler
% - UKFupdate
%
% Фукция dopplerSAT1_downLink моделирует измерения частоты Доплера вниз (и
% несколько расширяет возможности ранее разработанной функции dopplerSAT1)
%
% Функция OKS_Doppler - произодит оценку истинной траектории спутника. Для
% и поиска минимума функционала невязки используется алгоритм Нелдера-Мида.
%
% Функция SOMP_Doppler решает задачу само ОМП (СОМП) Поста, т.е. по K измерениям
% частоты Доплера от маяка НКА и сведений о грубой траектурии (TLE/SGP4) делает 
% оценку координат Поста (по умочанию -первого). Для поиска минимума функционала 
% невязки используется алгоритм Нелдера-Мида. Функция, в принципе, позволяет 
% судить степени отличия траекторий (истинной и грубой). Среднее по видимым 
% интервалам значение ошибки СОМП предоставляет переменная meanErr_Doppler.
%
% В качестве альтернативы функции OKS_Doppler модель позволяет использовать
% алгоритм Калмана. Для этого используется функция UKFupdate, которая обновленет
% вектор состояния нелинейного сигма-точечного фильтра (нелинейная версия 
% фильтра Калмана). В качестве измерений для фильтра используются измерения 
% частоты Доплера вниз, а вектор состояния имеет 6 координат: 3 координаты 
% местоположения и 3 координаты скорости НКА. 
%
% На данный момент нелинейный фильтр Калмана для выделения траектории по  измерениям 
% Доплера не даёт положительного результата. По-видимому, причина неудачи в
% том, что вектор измерений выливается в скаляр (в один момент времени лишь
% одно измерение Доплера), а вектор состояния имеет размер 6. (Ранее я 
% подобный фильтр делал для задачи слежения за траекторией  гипераппарата
% (задача для точек А и Б) и там фильтр работал. Однако там вектор измерения
% (в один момент времени) имел две координаты: азимут  и угол места.) 
%
% ЗАМЕЧАНИЕ
% В этой модели ИРИ не нужен. Однако для совместимости с предыдущими
% версиями в модель введён ИРИ, местоположение которого, совпадает с 
% Местоположением первого поста.
%
% При этом (для развития модели в будущем) 
% в качестве опции используется режим ОМП ИРИ на основе диффиринциально-дальномерного 
% метода (ДДМ, DDM, см. переменную DDM_OMP), анологичный режиму в функции GEOLOCATION_LEO_DDM. 
%
%%*************************************************************************
% Последняя версия библиотеки программ по геолокации расположена по адресу:
% https://yadi.sk/d/dNrJubwO3JFoVB/
%**************************************************************************
%
% var@mail.spbstu.ru, V.A. Vargauzin, декабрь-январь 2020, 
%
%**************************************************************************

if isempty(varargin)
    help(mfilename) % печать help
    close all       % закрыть все графики
    clc
end

%% GPU

n = gpuDeviceCount;
for ii = 1:n
    g = gpuDevice(ii);
    fprintf(1, 'Device %i has ComputeCapability %s \n', ...
            g.Index, g.ComputeCapability);
end

%% Глабальные переменные

global min_per_day sec_per_day UnitDist...
       Xp Yp TAI_UTC LOD TT_TAI  GPS_UTC UT1_UTC...
       omega_Earth JD1858 JD1950 JD2000 MJD2000

%% Управление программой

GUIinput=boolean(1); % использовать GUI (иначе - ввод парааметров в командной строке)
Geoid=boolean(0); % использовать или не использовать геоид
DDM_OMP=boolean(0); % исключить алгоритм DDM для ОМП ИРИ
rng('default')

%% Константы для пересчёта в другие единицы измерений

UnitDist=10^3; % 1 км  = 1000 м
Arcs      = 3600*180/pi;         % Arcseconds per radian
Rad       = pi/180;              % Radians per degree
UnitFreq=10^3; % 1 кГц =1000 Гц
OneGHz=UnitFreq^3; % 1 ГГц

%% Физические константы

% Earth rotation (derivative of GMST at J2000; differs from inertial period by precession)
omega_Earth = 15.04106717866910/3600*Rad;  % [rad/s]; WGS-84
% omega_Earth=7.29211514670698e-05
c = 299792458; % скорость света, м/c

%% Константы времени

min_per_day=60*24;    % число минут в день = 1440
sec_per_day=60*24*60; % число секунд в день = 8640
sec_per_mim=60;       % число секунд в минуте

[year,mon,day,hr,minute,sec] = datevec('18-Nov-1858 00:00');
JD1858=jday(year,mon,day,hr,minute,sec); % % юлианская дата (число дней от 1 января 4713 г.д.н.э.)
% J1858=2400000.5; % смещение для получения из юлианской даты модифицированной илианской даты 

[year,mon,day,hr,minute,sec] = datevec('1-Jan-2000 12:00');
JD2000=jday(year,mon,day,hr,minute,sec);
% JD2000=2451545.0

[year,mon,day,hr,minute,sec] = datevec('1-Jan-1950 00:00');
JD1950=jday(year,mon,day,hr,minute,sec);
% JD1950=2433281.5

MJD2000=JD2000- JD1858;  % модифицированная юлианская дата JD2000, к которой привязана инерциальная звёздная система координат IERC
% MJD_J2000 = 51544.5; 

TAI_UTC= +37; % сек, опережение атомного времени от UTC на 2017 год (Beginning 1 January 2017:)

TT_TAI  = +32.184; % сек, разница между земным и атомным временем https://ru.wikipedia.org/wiki/%D0%AD%D1%84%D0%B5%D0%BC%D0%B5%D1%80%D0%B8%D0%B4%D0%BD%D0%BE%D0%B5_%D0%B2%D1%80%D0%B5%D0%BC%D1%8F

GPS_UTC= +8*0; % сек, метки GPS по определению опережат метки UTC на 16 сек !!! 

% Март 2017:     
UT1_UTC= 0.51363; % сек
% коэффициенты полярного сжатия:
Xp= 0.0804/Arcs; % рад 
Yp= 0.2640/Arcs; % рад 
LOD= 0; %                                  
     
%% Геодезические координаты центра постов
% широта lat, град (-90...+90), долгота lon, град (-180...+180)

hPostSee=0; % геодезическая высота измерительного поста P


               latPost=60.065843; 
               lonPost=30.246163; 
            
          
      
if Geoid       
% Преобразование высоты над уровнем моря (геоидом 'egm96') в высоту над эллипсоидом
% для этого используется карта высот Земли,записанная в файл geoidegm96grid.mat        

       hPostGeoid = geoidheight( latPost, lonPost); % высота геоида над эллипсоидом в месте расположения поста 
       hPost=hPostSee+hPostGeoid; % геодезическая высота измерительного поста     
       
end
 
%% Выбор числа постов
 
NumPosts= 3; % число постов                          


 %% Формирование области постов (OP) вокруг центра

    % Далее формируется матрица 3 x NumPosts (координаты x посты) - декартовы
    % координаты измерительных постов P в (правой) геоцентрической 
    % системе координат (ГЦСК)
    % Центр ГЦСК совпадает с центром масс Земли
    % Ось OZ - лежит на малой полуоси эллипсоида и направлена на северный полюс
    % Ось OX - направлена в точку пересечения экватора (окружности) с меридианом Гринвича.     
    
    % Центр OP:
    latOP_0=latPost;
    lonOP_0=lonPost;
%     hOP_0=hPost;
    distOP=0; % случай одного поста

if NumPosts > 1
    
    distOP=350*UnitDist; % растояние от центра Постов, м
    
    
 end    
    
    xyzOP=nan(3,NumPosts);
    az=linspace(-180,+180,NumPosts+1);
    az=az(2:end);
    latOP=nan(size(az));
    lonOP=nan(size(az));                  
         
    hOP=zeros(size(az)); % высота постов над геоидом
    hSeeOP=zeros(size(az)); % высота над уровнем моря постов
    hAntennaOP=zeros(size(az)); % высота антенн постов          
         
    for j=1:length(az)         
        [latOP(j), lonOP(j)]=vreckon(latOP_0,lonOP_0, distOP, az(j)); 
        [x, y, z] = xyz( latOP(j), lonOP(j), hOP(j), Geoid, hSeeOP(j), hAntennaOP(j) ); 
        xyzOP(:,j)=[x, y, z]';              
    end

% для совместимости с предыдущей версией:
latPost=latOP(1);
lonPost=lonPost(1);
  
%% Выбор НКА           

    
      %TLE:
      % https://www.space-track.org/#/tle
      %
      % Здесь можно получить TLE для эпохи, наиболее близкой ко времени измерения.
      %
      % Последние TLE, однако, всегда можно свободно (без регистрации, в отличие от того,
      % как требуется на сайте space-track) получить по адресу:
      %
      % https://www.n2yo.com
      %
      % При этом, на этом же сайте организован и СЕРВИС СЛЕЖЕНИЯ ЗА СПУТНИКАМИ В РЕАЛЬНОМ ВРЕМЕНИ, 
      % где широта, долгота, высота спутника обновляются почти ежесекундно (это - просто видно). 
      % В частности, есть слежение и за UFO-10 (TLE которого, как выяснилось, обновляются ОЧЕНЬ редко),
      % см.:
      %
      % https://www.n2yo.com/satellite/?s=25967    
    
           
               % IRIDIUM 151               
                lonstr1KA1='1 43074U 17083E   20352.49900893  .00000110  00000-0  32273-4 0  9994';
                lonstr2KA1='2 43074  86.3979  81.0759 0001886 103.9895 256.1511 14.34217107156373';               
                            

           
%% ФОРМИРОВАНИЕ КООРДИНАТ ИРИ 
% В этой модели ИРИ не нужен. Однако для совместимости с предыдущими
% версиями считаем Местоположение ИРИ, совпадающим с Местоположением первого
% поста

hIRI=0;    % геодезическая высота ИРИ
hIRIsee=0; % высота над уровнем моря

% Menu= { 'совпадает с первым Постом'...
%         'совпадает с центром Постов'...
%         'отстоит от центра Постов на север'...
%         'отстоит от центра Постов на юг'...
%         'отстоит от центра Постов на восток'...
%         'отстоит от центра Постов на запад'};
% 
%       
%       titl='Местоположение ИРИ';
%       [select,OK]=listdlg('ListString',Menu,...
%                             'Name',titl,...
%                             'SelectionMode','single',... 
%                             'ListSize',[350 100]); % [width height]
%        if ~OK
%            disp('Дата измерения не выбрана')           
%            return                     
%        end       

     
               
           % совпадает с 1-м Постом
               latIRI=latOP(1); 
               lonIRI=lonOP(1); 
               hIRIsee=0;
               az=0;           
              
                            
    

    % Выбор расстония до Поста
    distOA=350*UnitDist; % растояние области анализа (OA) от Поста, м
   

    [latIRI, lonIRI]=vreckon(latOP_0,lonOP_0,distOA,az);

    % проверка растояния ПОСТ-ИРИ через геодезические координаты: 
    d = vdist(latIRI,lonIRI,latOP_0,lonOP_0);
    disp(['растояние центр Постов-ИРИ recon = ' num2str(d/UnitDist) ' км'])
    disp(' ')
       
     
% Преобразование высоты над уровнем моря (геоидом 'egm96') в высоту над эллипсоидом
% для этого используется карта высот геоида на эллипсоидом,записанная в файл geoidegm96grid.mat        

       hIRIGeoid = geoidheight( latIRI, lonIRI); % высота геоида над эллипсоидом в месте расположения поста 
       hIRI=hIRIsee+hIRIGeoid; % геодезическая высота измерительного поста             


% % геоцентрические координаты xyzIRI: 
 [xIRI, yIRI, zIRI]=llh2xyz( latIRI, lonIRI, hIRI );
 xyzIRI=[xIRI, yIRI, zIRI];

    
%% Инициализация моделирования (прогнозирования) траекторий движения спутников по алгоритму SPG4

    [satrecKA1,~, utc0_KA1, NORAD_ID_KA1] = spg4_initVARG(lonstr1KA1, lonstr2KA1);   
    
     % Структура satrecKA1 (параметров орбиты):
     
%          inclo: 1.7027 - inclination - наклонение орбиты в радианах, (180/pi) * angleInRadians= градусы
%           ecco: 0.0059 - eccentricity - эксентртрисетет орбиты
%          nodeo: 0.5558 - right ascension of ascending node - долгота                        
%                          восходящего узла,  (180/pi) * angleInRadians= градусы 
%          argpo: 0.5347 - argument of perigee (output if ds) - аргумент
%                          перигея
%             mo: 5.7566 - mean anomaly (output if ds)
%             no: 0.0647 - mean motion - среднее движение
%              a: 1.0977 - большая полуось эллипса ? - непохоже
%           alta: 0.1041 - высота в аппогее ?
%           altp: 0.0912 - высота в перигее-?

% Метод преобразования ( -> ) инерциальной системы координат (ECI), 
% используемой в алгоритме SGP4 (TEME), в фиксированную 
% земную геоцентрическую систему (ECEF)
    TEME2ECEF='ECI SGP4 (TEME, True Equator Mean Equinox) -> PEF (Pseudo Earth Fixed) + TOD (True of Date)  = ECEF'; 
    % возможен вариант TEME2ECEF='ECI TEME -> PEF' ; 

%%  Выбор параметров алгоритма поиска минимума функционала невязки методом Нелдера-Мида
%
% математика метода поиска минимума функции нескольких переменных 
% симплексным методом Нелдера-Мида (Nelder-Mead) хороша описана здесь:
%
% http://www.jasoncantarella.com/downloads/NelderMeadProof.pdf
%
% и здесь:
%
% http://www.mathworks.com/help/matlab/math/optimizing-nonlinear-functions.html#bsgpq6p-11
%
% Для двух переменных (широта и долгота) симплекс - равносторонний
% треугольник.
% 
% Для трёх переменных (широта, долгота и высота) симплекс - равносторонняя
% пирамида.
%
% Алгоритм ищет оптимальный симплекс. Оптимальным считается симплекс, в 
% котором расстояния между вершинами не превышают заданную величину tolx_TDOA,
% и, при этом различие между минимальным значением функционала в одной из 
% вершин симплекса (её номер всегда 1) между значениями функционала в остальных 
% вершинах симплекса не превышала величину tolf_TDOA.
%
% Параметрами алгоритма также являются: 
% - максимальное число итераций  maxiter для поиска минимума
% - максимальное число вычислений значения функционала maxfun для 
%   его поиска минимума
%
% Достоинством алгоритма является отсутствие в необходимости вычислять
% производные, градиенты и т.п.
%
% Алгоритм Нелдера-Мида на С++:
% http://www.codeguru.com/cpp/article.php/c17505/Simplex-Optimization-Algorithm-and-Implemetation-in-C-Programming.htm   
       
tolx = 1e-4;
tolf=1e-3;    
maxfun = 200*3; % максимальное число вычислений значения функционала
                % 3 - это размер числа переменных для оптимизации
maxiter= 200*3; % максимальное число итераций                          



%% Параметры модели

timeEst=cell(1);
% Дата начала измерений:
timeEst{1}=utc0_KA1; % эпоха TLE

% DT_=30; % интервал измерений, сек
% DT_ALL= 1*min_per_day/2; % интервал наблюдения, мин (1440 мин = min_per_day = сутки)  
% dX=1;  % км
% dY=-1; % км
% dZ=1;   % км
% DT_=60*10; % интервал измерений, сек



DT_=abs(30*1); % сек
DT_ALL=abs(1*min_per_day/1); % мин
dXYZ=UnitDist*[1 1 -1]; % метры 

freqDownLink=7*OneGHz; % Гц
RMS_DOPPLER=1; % Гц
OKSType=0;
min_elev_deg=10;

if OKSType > 3
    OKSType=2;
end   
    
%% Моделирование возможных моментов времени возможных измерений 
% (реальные измерений могут быть произведены при условии прямой видимости
%  ИРИ-НКА-Пост)
        
        [year,mon,day,hr,minute,sec] = datevec(timeEst{1}); % дата прогноза
        jdEst = jday(year,mon,day,hr,minute,sec); % юлианская дата прогноза, дни    
        jd1 = jday(year,mon,day,hr,minute,sec); % юлианская дата первого измерения прогноза, дни                                            
        DT=0; % смещение времени прогноза относительно первого измерения, мин           
        i=1;
        DT_min=DT_/60; % интервал измерений, мин
        jdEst=(jdEst*min_per_day +0)/min_per_day;

        while DT < DT_ALL
          
          i=i+1;           
          jdEst=(jdEst*min_per_day +DT_min)/min_per_day;         
          [year,mon,day,hr,minute,sec] = invjday(jdEst);
          dtm = datenum(year,mon,day,hr,minute,sec);
          timeEst{i}=datestr(dtm);            
          DT=(jdEst-jd1)*min_per_day; % смещение времени прогноза относительно первого измерения, мин         

        end        
                
%% Моделирование траектории и скоростей НКА

        Nexp=numel(timeEst); % Число измерений      
        adrKA1=zeros(1,Nexp);
        
        % траектория SGP4:        
        SGP4=zeros(Nexp,3);  % "грубая" траектория (X,Y,Z), полученная моделью прогнозирования TLE/SGP4            
        SGP4_V=zeros(Nexp,3);%  траектория вектора скорости (Vx,Vy,Vz)
        
        % Доплер для области постов (OP):
        SGP4_DOPPLER = zeros(Nexp,NumPosts); % траектория доплеровского смещения
        SGP4_V_PROECT=zeros(size(SGP4_DOPPLER)); % траектория проекции вектора скорости в направлении на Пост                       
        
        % истинная траектория TRU(E):
        TRU=zeros(Nexp,3);  % истинная траектория (X,Y,Z), смещённая от TLE/SGP4            
        TRU_V=zeros(Nexp,3);%  траектория вектора скорости (Vx,Vy,Vz)
        
         % Доплер для области постов (OP):
        TRU_DOPPLER = zeros(Nexp,NumPosts); % траектория доплеровского смещения
        TRU_V_PROECT=zeros(size(SGP4_DOPPLER)); % траектория проекции вектора скорости в направлении на Пост 
        MEAS_DOPPLER=zeros(size(SGP4_DOPPLER)); % измерения доплеровского смещения
        MEAS_V_PROECT=zeros(size(SGP4_DOPPLER)); % измерения проекции вектора скорости в направлении на Пост         
                
       % вычисление dt - смещение времени измерения (в мин)
       % относительно первого измерения        
        dt=zeros(Nexp,1); 
        [year,mon,day,hr,minute,sec] = datevec(timeEst{1}); % дата измерения
        jd1 = jday(year,mon,day,hr,minute,sec); % юлианская дата первого измерения прогноза, дни 
         
        % цикл по кол-ву измерений для вычисления даты измерения:       
        if Nexp >1
            for iexp = 2:Nexp
                % вычисление момента времени прогнозирования
                [year,mon,day,hr,minute,sec] = datevec(timeEst{iexp}); % дата прогноза
                jdEst = jday(year,mon,day,hr,minute,sec); % юлианская дата прогноза, дни                                  
                dt(iexp)=(jdEst-jd1)*min_per_day; % смещение времени прогноза относительно первого изимерения, мин           
            end         
        end 
        
     % цикл по кол-ву измерений для вычисления траектории НКА:       
     for iexp = 1:Nexp             

                [ r_ecf, v_ecf, satrecKA1, adrKA1(iexp) ] = SPG4_ecef(satrecKA1, timeEst{iexp},TEME2ECEF );
                sat_TLE_xyz(1) = r_ecf(1); % м   
                sat_TLE_xyz(2) = r_ecf(2); % м
                sat_TLE_xyz(3) = r_ecf(3); % м                
                SGP4(iexp,:)=sat_TLE_xyz;                
                              
                sat_TLE_Vxyz(1) = v_ecf(1); %  м/с    
                sat_TLE_Vxyz(2) = v_ecf(2); %  м/с   
                sat_TLE_Vxyz(3) = v_ecf(3); %  м/с                
                
                SGP4_V(iexp,:)=sat_TLE_Vxyz;
                
               [SGP4_DOPPLER(iexp,:), SGP4_V_PROECT(iexp,:)] =...
                   dopplerSAT1_downLink(sat_TLE_xyz, sat_TLE_Vxyz, xyzOP, freqDownLink); 
               
                % Истинная траектория:
                sat_TRUE_xyz=sat_TLE_xyz+dXYZ; % смещение траектории TRUE относительно SGP4
                TRU(iexp,:)=sat_TRUE_xyz;                
                               
                % считаем, что скорости у траектории TRU те же векторы
                % скорости, что и у траектории SGP4:                
                TRU_V(iexp,:)=SGP4_V(iexp,:);
               
                % Доплер для истинной траектории:
                [TRU_DOPPLER(iexp,:), TRU_V_PROECT(iexp,:)] =...
                    dopplerSAT1_downLink(sat_TRUE_xyz, sat_TLE_Vxyz, xyzOP, freqDownLink);
                
                MEAS_DOPPLER(iexp,:)=TRU_DOPPLER(iexp,:)+ RMS_DOPPLER*randn(1,NumPosts);
                MEAS_V_PROECT(iexp,:)=c*MEAS_DOPPLER(iexp,:)/freqDownLink;                
              
    end       
    
%% Анализ совместной прямой видимости (los) НКА Постом и ИРИ

% min_elev_deg=10; % minimal elevation, deg

[lat,lon,h]=xyz2llh(TRU(:,1), TRU(:,2), TRU(:,3) ); % lat,lon,h - векторы-столбцы 
losCUBIK=ones(size(lat));

losAnaliz=boolean(1); % если losAnaliz=0, то далее анализируется вся 
% траектория спутника независимо от его видимости постами и ИРИ (такая
% траектория формально считается "одним витком")

if losAnaliz

    % для круга постов из NumPosts штук:
    losPost=nan(length(lat), NumPosts);
    losCUBIK_Post=ones(size(lat)); % совместная видимость всеми постами

    for j=1:NumPosts

        % ENU на Посту:
        refLat=latOP(j);
        refLong=lonOP(j);
        refH=hOP(j);
        [~,~,losPost(:,j)]=peleng(refLat,refLong,refH,lat,lon,h,min_elev_deg);
        losCUBIK_Post=losCUBIK_Post & losPost(:,j);
    end

    % анализ прямой видимости ИРИ - НКА 

    % ENU на ИРИ:
    refLat=latIRI;
    refLong=lonIRI;
    refH=hIRI;
    [~,~,losCUBIK_IRI]=peleng(refLat,refLong,refH,lat,lon,h,min_elev_deg);

    % совместная видимость:
    losCUBIK=losCUBIK_Post & losCUBIK_IRI;
    % оставляем координаты НКА при совместной прямой видимости Постом и ИРИ:

end

% Траектория SGP4 - > sgp4:
sgp4=SGP4;
sgp4(losCUBIK==0,:)=nan; % видимые декартовые координаты НКА

sgp4_V=SGP4_V;
sgp4_V(losCUBIK==0,:)=nan; % видимые декартовые координаты СКОРОСТИ НКА

sgp4_DOPPLER=SGP4_DOPPLER;
sgp4_DOPPLER(losCUBIK==0,:)=nan; % видимые декартовые координаты ДОПЛЕРА НКА

sgp4_V_PROECT=SGP4_V_PROECT;
sgp4_V_PROECT(losCUBIK==0,:)=nan; % видимые декартовые координаты проекции СКОРОСТИ НКА

% Траектория TRU - > tru:
tru=TRU;
tru(losCUBIK==0,:)=nan; % видимые декартовые координаты НКА

tru_V=TRU_V;
tru_V(losCUBIK==0,:)=nan; % видимые декартовые координаты вектора СКОРОСТИ НКА

tru_DOPPLER=TRU_DOPPLER;
tru_DOPPLER(losCUBIK==0,:)=nan; % видимые декартовые частоты Доплера

tru_V_PROECT=TRU_V_PROECT;
tru_V_PROECT(losCUBIK==0,:)=nan; % видимые декартовые координаты проекции СКОРОСТИ НКА

% Измерения:
meas_DOPPLER=MEAS_DOPPLER;
meas_DOPPLER(losCUBIK==0,:)=nan; % видимые декартовые частоты Доплера

meas_V_PROECT=MEAS_V_PROECT;
meas_V_PROECT(losCUBIK==0,:)=nan; % видимые декартовые координаты проекции СКОРОСТИ НКА

latSAT1=lat;
latSAT1(losCUBIK==0,:)=nan;  % видимая широта НКА

lonSAT1=lon;
lonSAT1(losCUBIK==0,:)=nan; % видимая долгота НКА

%% Разделение моментов времени на видимые витки

Vitok=cell(1);
ind=find(losCUBIK);

if ~isempty(ind) && length(ind) > 1
    
    v=1; % первый видимый виток
    Vitok{v}=ind(1);

    for i=2: length(ind)
        if (ind(i)-ind(i-1))==1
           Vitok{v}=[Vitok{v} ind(i)] ; % пополнение индексов текущего видимого витка
        else
            v=v+1; % следующий видимый виток
            Vitok{v}=ind(i);
        end
    end  
    
end
% 
%% Алгоритм Калмана оценки траектории на витке по измерениям смещения Доплера

Kalman=boolean(0);
if OKSType==1 || OKSType==2
  Kalman=boolean(1);  
  KalmanType=OKSType;
end

 if Kalman
     
     if RMS_DOPPLER ~=0
        RMS_KALMAN=RMS_DOPPLER*c/freqDownLink; % Погрешность измерений частоты Доплера     
     else         
        RMS_KALMAN=10; % Гц        
     end
     
     % вектор состояния фильтра =
     % S=[X Y Z Vx Vy Vz]
     S=zeros(1,6);
     
     D=10^4; % стандартное отклонение априорного гауссова распределения координат (в метрах)      
     V=10^4; % м/с 
     
     accelX=10^5; % дисперсия ускорения по координате X
     accelY=10^5; % дисперсия ускорения по координате Y
     accelZ=10^5; % дисперсия ускорения по координате Z     
 
     dT=DT_; % сек
     % матрица перехода в новое состояние (кинематика движения вида dx=x+v*dt):
     A = [1   0  0   dT  0   0
          0   1  0   0   dT  0 
          0   0  1   0   0   dT
          0   0  0   1   0   0
          0   0  0   0   1   0
          0   0  0   0   0   1];

    % ковариционная матрица двумерного шума процесса с единичными дисперсиями
    % ускорения по координатам X и Y:  
      U=[dT^4/4   dT^3/2
         dT^3/2   dT^2];

    % матрица ускорения для трехмерного шума процесса:   
      accel=[accelX     0        0
               0      accelY     0
               0        0      accelZ];

    % ковариционная матрица шума трехмерного процесса, порожденного ускорениями, 
    % задавемыми матрицей accel:         
    Q=[accel*U(1,1)    accel*U(1,2);...
       accel*U(2,1)    accel*U(2,2)];       

    if KalmanType ==1
        % в качестве измерений используются только смещение Доплера
         R=RMS_KALMAN .^ 2*diag(ones(1,NumPosts));                
    else 
        % в качестве измерений используются смещение Доплера и данные
        % TLE/SGP4
         Q=NumPosts + numel(S); % число эдементов в векторе измерений
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
     
     NumVitkov=length(Vitok); % число видимых витков
     tru_est=sgp4;
     tru_V_est=sgp4_V;
     b=nan(1,NumVitkov);
     q=nan(1,NumVitkov);
     
   for v=1:NumVitkov
       
     ind=Vitok{v};
     M=length(ind); % число измерений на витке      
          
     % начальное значение вектора состояния фильтра =
     S=[sgp4(ind(1),:) sgp4_V(ind(1),:)];       
     
    % начальное значение ковариационной матрицы 
    % погрешности оценивания  координат вектора состояния
     P=[D^2  0   0   0   0   0
        0  D^2  0   0   0   0
        0   0  D^2  0   0   0 
        0   0   0  V^2  0   0
        0   0   0   0  V^2  0
        0   0   0   0   0  V^2];         

    for t = 1:M
         % Обновление вектора состояния нелинейным фильтром Калмана (UTF)
         % В качестве измерений используется не доплеровская частота f, а 
         % радиальная скорость v (tru_V_PROECT), чтобы в алгоритм не вводить 
         % нсесущую частоту вниз freqDownlink. Это возможно, поскольку 
         % f=(v/c)*freqDownlink, где c - скорость света.

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
      err_=sum(err.^2); % суммирование по столбцам, в которых координаты x y z
      b(v)=sqrt(mean(err_));
      
      err=tru_V_est(ind,:)-tru_V(ind,:);
      err=err';
      err_=sum(err.^2); % суммирование по столбцам, в которых координаты x y z
      q(v)=sqrt(mean(err_)); 
      
   end

  pos_err_Kalman=mean(b);
  vel_err_Kalman=mean(q);
% 
 end

%% Алгоритм оценки кординат Спутника (OKS) по Доплеру от N постов по Нелдеру

OKS=boolean(0);
if OKSType==0
    OKS=boolean(1);
end

if OKS
 NumVitkov=length(Vitok); % число видимых витков
 % Результат ОKS:
 tru_est_OKS=nan(size(sgp4));
 tru_V_est_OKS=nan(size(sgp4));
 b=nan(1,NumVitkov);
 q=nan(1,NumVitkov); 

 for v=1:NumVitkov

     ind=Vitok{v};
     ind=ind(1); % ОДИНОЧНОЕ ИЗМЕРЕНИЕ
     % известа траектория sgp4: 
     SAT=sgp4(ind,:); 
     SAT_V=sgp4_V(ind,:); 
     % измерения сделаны по истинной траектории:
     Dop=meas_V_PROECT(ind,:); % В качестве измерений используется не      
    %  Dop=tru_V_PROECT(ind,:); % В качестве измерений используется не 
     % доплеровская частота f (tru_DOPPLER), а  радиальная скорость v, 
     % чтобы не вводить несущую частоту вниз freqDownLink. Это возможно,
     % поскольку f=(v/c)*freqDownLink, где c - скорость света
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
      err_=sum(err.^2); % суммирование по столбцам, в которых координаты x y z
      b(v)=sqrt(mean(err_));
      
      err=tru_V_est_OKS(ind,:)-tru_V(ind,:);
      err=err';
      err_=sum(err.^2); % суммирование по столбцам, в которых координаты x y z
      q(v)=sqrt(mean(err_));         
 end
 
   pos_err_OKS=mean(b);
   vel_err_OKS=mean(q);

   disp(['среднее значение ошибки ОКS = ' num2str(pos_err_OKS/UnitDist)  ' км ' ])
   disp(' ')
 
end

%% номер поста для ОМП
PostNumber=1;

if DDM_OMP


  %% Моделировние измерений для ОМП ИРИ методом ДДМ
  
  diffD=DDM_measure(tru,xyzOP(:,PostNumber),xyzIRI); % вектор-столбец идеальных измерений (в моменты видимости)
  noise_Diff_D=0; % шум
  diffD_measure=diffD +noise_Diff_D; % измерения с шумом 
  
%% Алгоритм ОМП ИРИ методом ДДМ на основе минимазации квадратичного функционала невязки
           
 opt=[latOP(PostNumber), lonOP(PostNumber)] ; % априорные данные для поиска минимума квадратичного функционала                           
 h0=0; % геодезическая высота для поиска 
 
 NumVitkov=length(Vitok); % число видимых витков
 latIRI_DDM=nan(NumVitkov,1);
 lonIRI_DDM=nan(NumVitkov,1);
 NumdiffD=nan(NumVitkov,1); % число измерений на видимом витке 

 for v=1:NumVitkov

     ind=Vitok{v};
     % для ОМП используется не истинная траектория tru, а грубая sgp4
     % !!!!!!!!!!!!! После оценки траектории надо бы заменить sgp4 на tru_est:
%          SAT=tru_est(ind,:); 
     SAT=sgp4(ind,:); % матрица размера N x 3
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
          MinFlagStr='ОМП ИРИ методом ДДМ: min квадратичного функционала невязки не найден !!!';
          disp(MinFlagStr)
     end

    latIRI_DDM(v)=opt(1); % оценка широты ИРИ
    lonIRI_DDM(v)=opt(2); % оценка долготы ИРИ    
 
end

%%  расчёт ошибки ОМП ИРИ методом ДДМ

err_DDM=nan(NumVitkov,1);

if ~ isempty(Vitok{1})
    
    for v=1:NumVitkov

        % через декартовы координаты:
        h0=0;
        [xIRIest, yIRIest, zIRIest]=llh2xyz( latIRI_DDM(v), lonIRI_DDM(v), h0 );             
        err_DDM(v)=sqrt((xIRIest-xIRI).^2 + yIRIest-yIRI).^2+(zIRIest-zIRI).^2; 

%         if boolean(err_DDM(v))
%             % через геодезические координаты:                               
%              err_DDM(v)= vdist(latIRI_DDM(v),lonIRI_DDM(v),latIRI,lonIRI);
%         end
    end
end

meanErr=nanmean(err_DDM);
disp(['среднее значение ОМП ИРИ методом ДДМ, км,  = ' num2str(meanErr/UnitDist)])
disp(' ')

end % if DDM_OMP

  
%% Алгоритм CОМП Поста по Доплеру до ОКS на основе минимизации квадратичного функционала невязки

%  opt=[latPost, lonPost] ; % априорные данные для поиска минимума квадратичного функционала 
 opt=[latOP(PostNumber), lonOP(PostNumber)] ; % априорные данные для поиска минимума квадратичного функционала 
 h0=0; % геодезическая высота для поиска 
 
 NumVitkov=length(Vitok); % число видимых витков
 % Результат ОМП:
 latPost_Doppler=nan(NumVitkov,1);
 lonPost_Doppler=nan(NumVitkov,1);

 for v=1:NumVitkov

     ind=Vitok{v};
     ind=ind(1); % ОДИНОЧНОЕ ИЗМЕРЕНИЕ
     % известа траектория sgp4: 
     SAT=sgp4(ind,:); 
     SAT_V=sgp4_V(ind,:); 

%      SAT=tru_est_OKS(ind,:);
%      SAT_V=tru_V_est_OKS(ind,:);

     % измерения сделаны по истинной траектории:
     %Dop=tru_V_PROECT(ind,PostNumber); % истинный Доплер 
     Dop=meas_V_PROECT(ind,PostNumber); % измерения (с ошибкой)
     % В качестве измерений используется не 
     % доплеровская частота f (tru_DOPPLER), а  радиальная скорость v, 
     % чтобы в алгоритм SOMP_Doppler не вводить несущую частоту вниз
     % freqDownLink. Это возможно, поскольку f=(v/c)*freqDownLink, где c - 
     % скорость света
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
          MinFlagStr='СОМП по Доплеру: min квадратичного функционала невязки не найден !!!';
          disp(MinFlagStr)
     end

    latPost_Doppler(v)=opt(1); % оценка широты Поста
    lonPost_Doppler(v)=opt(2); % оценка долготы ПОСТА    
 
end

%%  расчёт ошибки ОМП Поста по Доплеру до ОКS

errPost_Doppler=nan(NumVitkov,1);

if ~ isempty(Vitok{1})
    
    for v=1:NumVitkov

        % через декартовы координаты:
        h0=0;
        [xIRIest, yIRIest, zIRIest]=llh2xyz( latPost_Doppler(v), lonPost_Doppler(v), h0 );             
        errPost_Doppler(v)=sqrt((xIRIest-xyzOP(1,PostNumber)).^2 + (yIRIest-xyzOP(2,PostNumber)).^2+(zIRIest-xyzOP(3,PostNumber)).^2); 

        if boolean(errPost_Doppler(v))
            % через геодезические координаты координаты:                               
             errPost_Doppler(v)= vdist(latPost_Doppler(v),lonPost_Doppler(v),latOP(PostNumber),lonOP(PostNumber));
        end
    end
end

meanErr_Doppler=nanmean(errPost_Doppler);
disp(['среднее значение ошибки СОМП Поста по Доплеру = ' num2str(meanErr_Doppler/UnitDist)  ' км ' ])
disp(' ')


    
%% ГРАФИКА ПОСЛЕ ОБРАБОТКИ



end % GEOLOCATION_LEO_DOPPLER_N_POST
end