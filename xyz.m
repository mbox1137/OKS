function [ xIRI, yIRI, zIRI] = xyz( latIRI, lonIRI, hIRI, Geoid, hIRIsee, hIRI_Antenna )

% эта функция идентична функции xyzIRI и отличается от неё лишь именем
%
% var@mail.spbstu.ru, май 2019, январь 2021

                if Geoid && ~isnan(latIRI)       
                % Преобразование высоты над уровнем моря (геоидом 'egm96') в высоту над эллипсоидом
                % для этого используется карта высот геоида на эллипсоидом,записанная в файл geoidegm96grid.mat
 
                       hIRIGeoid = geoidheight( latIRI, lonIRI); % высота геоида над эллипсоидом в месте расположения поста 
                       hIRI=hIRIsee+hIRIGeoid; % геодезическая высота измерительного поста P

                end 

                hIRI=hIRI+hIRI_Antenna;
               % геоцентрические декартовы координаты IRI:               
                [xIRI, yIRI, zIRI]=llh2xyz( latIRI, lonIRI, hIRI );

end % xyz

