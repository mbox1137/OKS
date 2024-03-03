function SAT = Norad2name( NORAD_ID_KA)
%
% var@mail.spbstu.ru, Vargauzin V.A., июль, октябрь 2018, декабрь 2020,
% декабрь 2021

           switch NORAD_ID_KA        
                case 39022
                     SAT='YAMAL 402';
                case 40277
                     SAT='EXPRESS AM-6'; 
                case 37749 
                     SAT='KAZSAT 2';
                case 28089
                     SAT='YAMAL 202';
                case 40505
                     SAT='EXPRESS AM-7';
                case 34710
                     SAT='Eutelsat 10A';
                case 26243
                     SAT='Eutelsat 16C';
                case 38992
                     SAT='Eutelsat 21B';
                case 27948
                     SAT='Eutelsat 31A';
                case 37836
                     SAT='Eutelsat 16A';
               case 36101
                    SAT='Eutelsat 36B';   
               case 33056
                    SAT='Turksat 3A';                     
                case 39773
                     SAT='Eutelsat 3B';                    
                case 39163
                     SAT='Eutelsat 7B';                                          
                case 41310
                     SAT='Eutelsat 9B';                                       
                case 33459
                    SAT='Eutelsat Hot Bird 13C'; 
                case 39233 
                    SAT='Eutelsat 25B';
                case 39617
                    SAT='Astra 5B';
                case 36745
                    SAT='Arabsat 5A';
                case 38741
                    SAT='Hylas 2';
                case 40613 % Западное полушарие
                    SAT='Thor 7 (1 West)';                      
                case 37816 
                    SAT='Eutelsat 7 West B';                   
                case 40875
                    SAT='Eutelsat 8 West B'; 
                case 42950
                    SAT='Intelsat 37e (18 West)';
                case 43039
                    SAT='AlComSat 1 (25 West)'; 
                case 37264
                    SAT='Hispasat 30W-5 (30 West)';
               case 25967 
                    SAT='UFO 10';                      
               case 35943                   
                    SAT='COMSATBW-1';
               case 40374
                     SAT='MUOS 15.5 West';
               case 39444
                   SAT='FUNcube-1 (AO-73)';
               case 32875 % падает
                   SAT='COSMOS 2421 DEB';
               case 32785
                   SAT='Cute-1.7+APD II';
               case 27844
                   SAT='Cute-1';
               case 25544
                   SAT='МКС';
               case 38708
                   SAT='BKA 2';
               case 29449
                   SAT='HINODE (SOLAR-B)';
               case 40966
                   SAT='AEROCUBE 7';
               case 40967
                   SAT='FOX-1A (AO-85)';
               case 40660
                   SAT='AEROCUBE 8B';
               case 39403 
                   SAT='HO OPONOPONO 2';
               case 41168
                   SAT='ATHENOXAT 1';
               case 44854
                   SAT='OBJECT C';
               case 41340
                   SAT='HORYU 4';
               case 44371
                   SAT='SPACEBEE-8';
               case 40936
                   SAT='EXACTVIEW 9';
               case 43049
                   SAT='ASGARDIA 1';
               case 39268
                   SAT='POPACS';
               case 27607 
                   SAT='SO-50';
               case 37190
                   SAT='GLOBALSTAR M076';
               case 37191
                   SAT='GLOBALSTAR M077';
               case 37742
                   SAT='GLOBALSTAR M085';
               case 43074
                   SAT='IRIDIUM 151';
               case 43929
                   SAT='IRIDIUM 171';                
               case 43925 
                   SAT='IRIDIUM 173';
               case 25163
                   SAT ='GLOBALSTAR M004';
               case 35008
                   SAT ='MERIDIAN 2';
               case 37398
                   SAT ='MERIDIAN 4';                   
               case 38995
                   SAT ='MERIDIAN 6';                   
               case 40296
                   SAT ='MERIDIAN 7';                   
               case 44453
                   SAT ='MERIDIAN 8';                   
               case 45254
                   SAT ='MERIDIAN 9';
               case 28885
                   SAT ='SYRACUSE 3A'; % 2005
               case 29273
                   SAT ='SYRACUSE 3B'; % 2006
               case 49333
                   SAT ='SYRACUSE 4A'; % 2021
               otherwise           
                    SAT= ['номер (NORAD CAT ID) ' int2str(NORAD_ID_KA)];
           end            
end

