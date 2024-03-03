% ----------------------------------------------------------------------------
%
%                           function teme2ecef
%
%  this function trsnforms a vector from the true equator mean equniox frame
%    (teme), to an earth fixed (ITRF) frame.  the results take into account
%    the effects of sidereal time, and polar motion.
%
%  author        : david vallado                  719-573-2600   10 dec 2007
%
%  revisions
%
%  inputs          description                    range / units
%    rteme       - position vector teme           km
%    vteme       - velocity vector teme           km/s
%    ateme       - acceleration vector teme       km/s2
%    LOD         - excess length of day           sec
%    MJDtt         - julian centuries of tt         centuries
%    JDut1       - julian date of ut1             days from 4713 bc
%    Xp          - polar motion coefficient       rad
%    Yp          - polar motion coefficient       rad
%
%  outputs       :
%    recef       - position vector earth fixed    km
%    vecef       - velocity vector earth fixed    km/s
%    aecef       - acceleration vector earth fixedkm/s2
%
%  locals        :
%    st          - matrix for pef - tod 
%    pm          - matrix for ecef - pef 
%
%  coupling      :
%   gstime       - greenwich mean sidereal time   rad
%   polarm       - rotation for polar motion      pef - ecef
%
%  references    :
%    vallado       2007, 219-228
%
% [recef,vecef,aecef] = teme2ecef  ( rteme,vteme,ateme,MJDtt,JDut1,LOD,Xp,Yp );
% ----------------------------------------------------------------------------

function [recef,vecef,aecef] = teme2ecefVARG  ( rteme,vteme,ateme,MJDtt,JDut1)

global omega_Earth Xp Yp LOD
% omega_Earth=7.29211514670698e-05 

        %% ------------------------ седерическое время gmst --------------------------
        gmst= gstimeVARG( JDut1 );
        
        %% TEME (True Equator Mean Equinox) to PEF (Pseudo Earth Fixed), без учёта движения полюсов:

%         thetasa    = 7.29211514670698e-05 * (1.0  - LOD/86400.0 );
        thetasa    = omega_Earth * (1.0  - LOD/86400.0 );

        st(1,1) =  cos(gmst);
        st(1,2) = -sin(gmst);
        st(1,3) =  0.0;
        st(2,1) =  sin(gmst);
        st(2,2) =  cos(gmst);
        st(2,3) =  0.0;
        st(3,1) =  0.0;
        st(3,2) =  0.0;
        st(3,3) =  1.0;
        
%         st=[cos(gmst)  -sin(gmst)  0.0
%             sin(gmst)   cos(gmst)  0.0
%             0.0         0.0        1.0];

        rpef  = st'*rteme;
%         thetasa= 7.29211514670698e-05 * (1.0  - LOD/86400.0 );
        omegaearth = [0; 0; thetasa;];     
        vpef  = st'*vteme - cross( omegaearth,rpef );

        %%  PEF to TOD (True of Date)  c учётом движения полюсов:

        [pm] = polarm( Xp, Yp, MJDtt,'80');
        recef = pm'*rpef;
        vecef = pm'*vpef;
        
        if ~isempty(ateme)

            temp  = cross(omegaearth,rpef);

            aecef = pm'*(st'*ateme - cross(omegaearth,temp) ...
                    - 2.0*cross(omegaearth,vpef));
        else
            aecef=[];
        end

