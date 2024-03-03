function [ xIRI, yIRI, zIRI] = xyz( latIRI, lonIRI, hIRI, Geoid, hIRIsee, hIRI_Antenna )

% ��� ������� ��������� ������� xyzIRI � ���������� �� �� ���� ������
%
% var@mail.spbstu.ru, ��� 2019, ������ 2021

                if Geoid && ~isnan(latIRI)       
                % �������������� ������ ��� ������� ���� (������� 'egm96') � ������ ��� �����������
                % ��� ����� ������������ ����� ����� ������ �� �����������,���������� � ���� geoidegm96grid.mat
 
                       hIRIGeoid = geoidheight( latIRI, lonIRI); % ������ ������ ��� ����������� � ����� ������������ ����� 
                       hIRI=hIRIsee+hIRIGeoid; % ������������� ������ �������������� ����� P

                end 

                hIRI=hIRI+hIRI_Antenna;
               % ��������������� ��������� ���������� IRI:               
                [xIRI, yIRI, zIRI]=llh2xyz( latIRI, lonIRI, hIRI );

end % xyz

