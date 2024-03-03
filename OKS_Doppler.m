%% ����������� ������ ������� � �������� �������� � ������� ����������:                                         
%  function [SAT, SAT_V]= OKS_Doppler(SAT,...
%                                     SAT_V,...                          
%                                     Doppler_measure,...
%                                     xyzOP,...                                   
%                                     tolx,...
%                                     tolf,...
%                                     maxiter,...
%                                     maxfun,...
%                                     freqDownLink) 
%                                     
%  
% % var@mail.spbstu.ru, V.A. Vargauzin, ������ 2020
% 
%  if  isempty(SAT) || isempty(Doppler_measure)
%      return
%  end
%  
%  if nargin < 9
%      freqDownLink=0;
%  end
%  
%  M=size(Doppler_measure,1); % ����� ����������
%  for j=1:M 
%      
%      opt=[SAT(j,:) SAT_V(j,:)];
%      DOP=Doppler_measure(j,:); 
%   
%     % ���������� ������� �� ����� ��������� ������ ��������:
%      printtype='none'; % 'none' 'iter' 'final' 'simplex'             
%                                           
%     % ���������� min ������������� ����������� �������
%      opt=fminsearchVARG(@fun, opt,...
%                         optimset('TolFun',tolf,...
%                                  'TolX',  tolx,...
%                                  'MaxFunEvals',maxfun,...
%                                  'MaxIter',maxiter,...
%                                  'Display', printtype ));
%        SAT(j,:)=opt(1:3);
%        SAT_V(j,:)=opt(4:6);
%  end
%                        function dR=fun(opt)
%                              
%                              if freqDownLink==0
% %                              freqDownLink=0; % ������� ����� ����, �������� � �������� ��������� 
% %                               % ������������ �� ������� �������, � ������� � ��� ��������� ����������
% %                               % ��������
%                              [ ~,Doppler_hipoteza]=dopplerSAT1_downLink(opt(1:3),...
%                                                                         opt(4:6),...
%                                                                         xyzOP,...
%                                                                         freqDownLink);  
%                              else
%                              [ Doppler_hipoteza]=dopplerSAT1_downLink(opt(1:3),...
%                                                                         opt(4:6),...
%                                                                         xyzOP,...
%                                                                         freqDownLink);
%                              end
%                              d=Doppler_hipoteza-DOP; % ������-������� �������                             
%                              dR=nansum(d.^2);
%                              %dR=(nansum(d))^2;
% 
%                        end   
%                      
% end % OKS_Doppler


%% ����������� ������ ������� ��������, ������ ������ �������� ���������:
%                                         
%  function [SAT, SAT_V]= OKS_Doppler(SAT,...
%                                     SAT_V,...                          
%                                     Doppler_measure,...
%                                     xyzOP,...                                   
%                                     tolx,...
%                                     tolf,...
%                                     maxiter,...
%                                     maxfun,...
%                                     freqDownLink) 
%                                                                            
% % 
% % var@mail.spbstu.ru, V.A. Vargauzin, ������ 2020
% 
%  if  isempty(SAT) || isempty(Doppler_measure)
%      return
%  end
%  
%  if nargin < 9
%      freqDownLink=0;
%  end
%  
%  M=size(Doppler_measure,1); % ����� ����������
%  for j=1:M 
%      
%      opt=[SAT(j,:)];
%      DOP=Doppler_measure(j,:); 
%   
%     % ���������� ������� �� ����� ��������� ������ ��������:
%      printtype='none'; % 'none' 'iter' 'final' 'simplex'             
%                                           
%     % ���������� min ������������� ����������� �������
%      opt=fminsearchVARG(@fun, opt,...
%                         optimset('TolFun',tolf,...
%                                  'TolX',  tolx,...
%                                  'MaxFunEvals',maxfun,...
%                                  'MaxIter',maxiter,...
%                                  'Display', printtype ));
%        SAT(j,:)=opt(1:3);
% %        SAT_V(j,:)=opt(4:6);
%  end
%                        function dR=fun(opt)
%                              
%                              if freqDownLink==0
% %                              freqDownLink=0; % ������� ����� ����, �������� � �������� ��������� 
% %                               % ������������ �� ������� �������, � ������� � ��� ��������� ����������
% %                               % ��������
%                              [ ~,Doppler_hipoteza]=dopplerSAT1_downLink(opt(1:3),...
%                                                                         SAT_V(j,:),...
%                                                                         xyzOP,...
%                                                                         freqDownLink);  
%                              else
%                              [ Doppler_hipoteza]=dopplerSAT1_downLink(opt(1:3),...
%                                                                         SAT_V(j,:),...
%                                                                         xyzOP,...
%                                                                         freqDownLink);
%                              end
%                              d=Doppler_hipoteza-DOP; % ������-������� �������                             
%                              dR=nansum(d.^2);
% %                              %dR=(nansum(d))^2;
% 
%                        end   
%                      
% end % OKS_Doppler

%% 2020: ����������� ������ ������� � �������� �������� � ������� ���������� �� �Ѩ� �����:
%  function [SAT, SAT_V]= OKS_Doppler(SAT,...
%                                     SAT_V,...                          
%                                     Doppler_measure,...
%                                     xyzOP,...                                   
%                                     tolx,...
%                                     tolf,...
%                                     maxiter,...
%                                     maxfun,...
%                                     freqDownLink) 
% %                                   
% %
% % var@mail.spbstu.ru, V.A. Vargauzin, ������ 2020
% 
%  if  isempty(SAT) || isempty(Doppler_measure)
%      return
%  end
%  
%  if nargin < 9
%      freqDownLink=0;
%  end
%       
%  opt=[SAT SAT_V];
%  DOP=Doppler_measure; 
%   
% % ���������� ������� �� ����� ��������� ������ ��������:
% printtype='none'; % 'none' 'iter' 'final' 'simplex'             
%                                           
% % ���������� min ������������� ����������� �������
% % opt=fminsearchVARG(@fun, opt,...
% %                         optimset('TolFun',tolf,...
% %                                  'TolX',  tolx,...
% %                                  'MaxFunEvals',maxfun,...
% %                                  'MaxIter',maxiter,...
% %                                  'Display', printtype ));
% opt=fminsearch(@fun, opt,...
%                         optimset('TolFun',tolf,...
%                                  'TolX',  tolx,...
%                                  'MaxFunEvals',maxfun,...
%                                  'MaxIter',maxiter,...
%                                  'Display', printtype ));
% SAT=opt(:,1:3);
% SAT_V=opt(:,4:6);
%        
% 
%                        function dR=fun(opt)
%                              
%                              if freqDownLink==0
% %                              freqDownLink=0; % ������� ����� ����, �������� � �������� ��������� 
% %                               % ������������ �� ������� �������, � ������� � ��� ��������� ����������
% %                               % ��������
%                              [ ~,Doppler_hipoteza]=dopplerSAT1_downLink(opt(:,1:3),...
%                                                                         opt(:,4:6),...
%                                                                         xyzOP,...
%                                                                         freqDownLink);
% 
%                              else
%                              [ Doppler_hipoteza]=dopplerSAT1_downLink(opt(:,1:3),...
%                                                                         opt(:,4:6),...
%                                                                         xyzOP,...
%                                                                         freqDownLink);
%                              end
%                              d=Doppler_hipoteza-DOP; % ������-������� �������                             
%                              dR=nansum(nansum(d.^2));
%                              %dR=(nansum(d))^2;
% 
%                        end 

%% 2023: ����������� ������ ������� � �������� �������� � ������� ���������� �� �Ѩ� �����: 
 function [SAT, SAT_V]= OKS_Doppler(SAT,...
                                    SAT_V,...                          
                                    Doppler_measure,...
                                    xyzOP,...                                   
                                    tolx,...
                                    tolf,...
                                    maxiter,...
                                    maxfun,...
                                    freqDownLink) 
%                                    
% 
% var@mail.spbstu.ru, V.A. Vargauzin, ������ 2020, ������ 2023

 if  isempty(SAT) || isempty(Doppler_measure)
     return
 end
 
 if nargin < 9
     freqDownLink=0;
 end
      
%  opt=[SAT SAT_V];
opt=[SAT]; % 2023
 DOP=Doppler_measure; 
  
% ���������� ������� �� ����� ��������� ������ ��������:
printtype='none'; % 'none' 'iter' 'final' 'simplex'             
                                          
% ���������� min ������������� ����������� �������
% opt=fminsearchVARG(@fun, opt,...
%                         optimset('TolFun',tolf,...
%                                  'TolX',  tolx,...
%                                  'MaxFunEvals',maxfun,...
%                                  'MaxIter',maxiter,...
%                                  'Display', printtype ));

opt=fminsearch(@fun, opt,...
                        optimset('TolFun',tolf,...
                                 'TolX',  tolx,...
                                 'MaxFunEvals',maxfun,...
                                 'MaxIter',maxiter,...
                                 'Display', printtype ));

SAT=opt(:,1:3);
% SAT_V=opt(:,4:6); % 2023
       

                       function dR=fun(opt)
                             
                             if freqDownLink==0
%                              freqDownLink=0; % ������� ����� ����, �������� � �������� ��������� 
%                               % ������������ �� ������� �������, � ������� � ��� ��������� ����������
%                               % ��������
%                              [ ~,Doppler_hipoteza]=dopplerSAT1_downLink(opt(:,1:3),...
%                                                                         opt(:,4:6),...
%                                                                         xyzOP,...
%                                                                         freqDownLink);
                             [~, Doppler_hipoteza]=dopplerSAT1_downLink(opt(:,1:3),...
                                                                      SAT_V,...
                                                                       xyzOP,...
                                                                       freqDownLink);
                             else
%                              [ Doppler_hipoteza]=dopplerSAT1_downLink(opt(:,1:3),...
%                                                                         opt(:,4:6),...
%                                                                         xyzOP,...
%                                                                         freqDownLink);
                             [ Doppler_hipoteza]=dopplerSAT1_downLink(opt(:,1:3),...
                                                                      SAT_V,...
                                                                       xyzOP,...
                                                                       freqDownLink);
                             end
                             d=Doppler_hipoteza-DOP; % ������-������� �������                             
                             dR=nansum(nansum(d.^2));
                             dR=sqrt(nansum(nansum(d.^2)));
%                              dR=(nansum(d))^2;
%                              dR=abs(nansum(d));

                       end   
                     
end % OKS_Doppler