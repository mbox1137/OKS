function [ r_ecef, v_ecef, satrec, adrKA ] = SPG4_ecef( satrec, timeEst, toECEF)

global min_per_day sec_per_day UnitDist...
       TAI_UTC  TT_TAI GPS_UTC UT1_UTC JD1858

jd0 = satrec.jdsatepoch;
[year,mon,day,hr,minute,sec] = datevec(timeEst); % дата измерения
jdEst = jday(year,mon,day,hr,minute,sec); % юлианская дата измерения, дни

dt=(jdEst-jd0)*min_per_day ; % смещение времени измерения относительно эпохи TLE, мин                   
               
% [satrec, r, v] = sgp4 (satrec,  dt);   % dt - мин,  r - км, v - км/с 

switch toECEF
    
    case 'ECI SGP4 (TEME, True Equator Mean Equinox) -> PEF (Pseudo Earth Fixed) + TOD (True of Date)  = ECEF' % из ВАЛЛАДО
     
        utc = hms2sec( hr, minute,sec ); % сек
        ut1= utc + UT1_UTC - GPS_UTC;
        [hrtemp,mintemp,sectemp] = sec2hms(  ut1 );
        JDut1 = jday( year,mon,day, hrtemp, mintemp, sectemp );
        
        TT_UTC= TAI_UTC + TT_TAI ;   % sec
        [hrtemp,mintemp,sectemp] = sec2hms( TT_UTC );
        JDtt = jday( year,mon,day, hrtemp, mintemp, sectemp);
        
        dt=(JDut1-jd0)*min_per_day;  % смещение времени измерения относительно эпохи TLE, мин                    
               
        [satrec, r, v] = sgp4 (satrec,  dt);   % dt - мин,  r - км, v - км/с 
        
         MJDtt= (JDtt - JD1858 )/ 36525.0; 
        [r_ecef,v_ecef] = teme2ecefVARG(r',v',[],MJDtt,JDut1); % ECI (TEME) -> ECEF     
    
    case 'ECI TEME -> PEF'
        
         [satrec, r, v] = sgp4 (satrec,  dt);   % dt - мин,  r - км, v - км/с         
         jd1 = jd0 + dt/min_per_day + UT1_UTC/sec_per_day - GPS_UTC/sec_per_day ;  % число дней в юлианской дате измерения
         gst=gstimeVARG(jd1); % gst - гринвичское седерическое время (в радианах) относительно эпохи J2000 
         % угол считается относительно J2000 = 2451545.0 юлианских дней
         [r_ecef, v_ecef]=teme2pef(r, v, gst); % преобразование из ECI в ECEF     
         
    case ''
        
         [satrec, r, v] = sgp4 (satrec,  dt);   % dt - мин,  r - км, v - км/с         
         jd1 = jd0 + dt/min_per_day + UT1_UTC/sec_per_day;  % число дней в юлианской дате измерения
         gst=gstimeVARG(jd1); % gst - гринвичское седерическое время (в радианах) относительно эпохи J2000   
         [r_ecef, v_ecef]=teme2pef(r, v, gst); % преобразование из ECI в ECEF
         
    otherwise 
        
         [satrec, r, v] = sgp4 (satrec,  dt);   % dt - мин,  r - км, v - км/с                                   
         jd1 = jd0 + dt/min_per_day + UT1_UTC/sec_per_day;  % число дней в юлианской дате измерения
         gst=gstimeVARG(jd1);  % gst - гринвичское седерическое время (в радианах) относительно эпохи J2000 
         % угол считается относительно J2000 = 2451545.0 юлианских дней         
         [r_ecef, v_ecef]=eci2ecfVARG(r, v, gst); % преобразование из ECI в ECEF
end

        adrKA=round(min_per_day + dt +1);% адрес позиции KA в прогнозе с точностью до минуты   
        r_ecef=r_ecef*UnitDist;
        v_ecef=v_ecef*UnitDist;

end

