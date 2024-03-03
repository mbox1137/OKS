function [ r_ecef, v_ecef, satrec, adrKA ] = SPG4_ecef( satrec, timeEst, toECEF)

global min_per_day sec_per_day UnitDist...
       TAI_UTC  TT_TAI GPS_UTC UT1_UTC JD1858

jd0 = satrec.jdsatepoch;
[year,mon,day,hr,minute,sec] = datevec(timeEst); % ���� ���������
jdEst = jday(year,mon,day,hr,minute,sec); % ��������� ���� ���������, ���

dt=(jdEst-jd0)*min_per_day ; % �������� ������� ��������� ������������ ����� TLE, ���                   
               
% [satrec, r, v] = sgp4 (satrec,  dt);   % dt - ���,  r - ��, v - ��/� 

switch toECEF
    
    case 'ECI SGP4 (TEME, True Equator Mean Equinox) -> PEF (Pseudo Earth Fixed) + TOD (True of Date)  = ECEF' % �� �������
     
        utc = hms2sec( hr, minute,sec ); % ���
        ut1= utc + UT1_UTC - GPS_UTC;
        [hrtemp,mintemp,sectemp] = sec2hms(  ut1 );
        JDut1 = jday( year,mon,day, hrtemp, mintemp, sectemp );
        
        TT_UTC= TAI_UTC + TT_TAI ;   % sec
        [hrtemp,mintemp,sectemp] = sec2hms( TT_UTC );
        JDtt = jday( year,mon,day, hrtemp, mintemp, sectemp);
        
        dt=(JDut1-jd0)*min_per_day;  % �������� ������� ��������� ������������ ����� TLE, ���                    
               
        [satrec, r, v] = sgp4 (satrec,  dt);   % dt - ���,  r - ��, v - ��/� 
        
         MJDtt= (JDtt - JD1858 )/ 36525.0; 
        [r_ecef,v_ecef] = teme2ecefVARG(r',v',[],MJDtt,JDut1); % ECI (TEME) -> ECEF     
    
    case 'ECI TEME -> PEF'
        
         [satrec, r, v] = sgp4 (satrec,  dt);   % dt - ���,  r - ��, v - ��/�         
         jd1 = jd0 + dt/min_per_day + UT1_UTC/sec_per_day - GPS_UTC/sec_per_day ;  % ����� ���� � ��������� ���� ���������
         gst=gstimeVARG(jd1); % gst - ����������� ������������ ����� (� ��������) ������������ ����� J2000 
         % ���� ��������� ������������ J2000 = 2451545.0 ��������� ����
         [r_ecef, v_ecef]=teme2pef(r, v, gst); % �������������� �� ECI � ECEF     
         
    case ''
        
         [satrec, r, v] = sgp4 (satrec,  dt);   % dt - ���,  r - ��, v - ��/�         
         jd1 = jd0 + dt/min_per_day + UT1_UTC/sec_per_day;  % ����� ���� � ��������� ���� ���������
         gst=gstimeVARG(jd1); % gst - ����������� ������������ ����� (� ��������) ������������ ����� J2000   
         [r_ecef, v_ecef]=teme2pef(r, v, gst); % �������������� �� ECI � ECEF
         
    otherwise 
        
         [satrec, r, v] = sgp4 (satrec,  dt);   % dt - ���,  r - ��, v - ��/�                                   
         jd1 = jd0 + dt/min_per_day + UT1_UTC/sec_per_day;  % ����� ���� � ��������� ���� ���������
         gst=gstimeVARG(jd1);  % gst - ����������� ������������ ����� (� ��������) ������������ ����� J2000 
         % ���� ��������� ������������ J2000 = 2451545.0 ��������� ����         
         [r_ecef, v_ecef]=eci2ecfVARG(r, v, gst); % �������������� �� ECI � ECEF
end

        adrKA=round(min_per_day + dt +1);% ����� ������� KA � �������� � ��������� �� ������   
        r_ecef=r_ecef*UnitDist;
        v_ecef=v_ecef*UnitDist;

end

