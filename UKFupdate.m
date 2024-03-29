function [S, P]= UKFupdate(u, xyzPost, S, P, A, Q, R, typ)                                       
                                         
% ������� ���������� ������ ��������� ��� �����-�������� �������� (���,
% Unscented Kalman Filter,)����������� ������� ������� ����������� 
% ������������� ������
%
% u - ������ ��������� 
% S - ������ ���������
% P - �������������� ������� ����������� ���������� 
% A - ������� �������� ���������� ������ �������� �������� ���
% Q - �������������� ������� ���� �������� ��������
% R - �������������� ������� ���� ��������� 
%
% var@mail.spbstu.ru, V.A. Vargauzin, ������� 2015, ������� 2020

xhat=S;
if size(S,2) > 1
    xhat=xhat';
end
% xhat - ������-�������

if size(u,2) > 1
    u=u';
end
% u - ������-�������

    n=6; % = numel(xhat)
    n2=2*n; 
    SIGMAm = zeros(n,n2); % ����� ��� �������� �������������
    % ����� ����� ��������� (������ ����� - ������ ���������)
%     m=1; % ������ ������� ��������� = numel(u)
    m=numel(u);
    Z = zeros(m,n2); % ����� ��� �������� ������������� �����
    % ����� ��������� (������ ����� - ������ ���������) 

    % ������-������� ���������:
    y=u; 

    % ���������� 2n ����� ����� ��������� (xhat,P -> SIGMA):
    SIGMA = ukf_sigmapoints(xhat,P);

    % ���������� ����������� ��� ������ �� ����� ����� SIGMA
    % (SIGMA -> SIGMAm):
    for j = 1:n2
        SIGMAm(:,j) = updatex(SIGMA(:,j),A);
    end

    % ���������� ������ ������������ (�������� xm � �������������� 
    % ������� Pm) �� ���� 2n ������  SIGMAm (SIGMAm -> xm, Pm):
    [xm, Pm] = ukf_mean_covariance(SIGMAm);
    Pm = Pm + Q;

    % ���������� 2n ����� ����� ��������� (������������ ��������� Z ��� 
    % ���� SIGMAm) (SIGMAm -> Z): 
    for j = 1:n2
        Z(:,j) = h(SIGMAm(:,j));                   
    end

    % ���������� ������ ������������ ��������� 
    % (�������� z � �������������� ������� Pz) ��
    % ���� 2n ������������� ���������� Z (Z -> z, Pz):
    [z, Pz] = ukf_mean_covariance(Z);
    Pz = Pz + R;
    
    % ���������� �������-�������������� ������� ����� ����� ������� 
    % SIGMAm � Z: 
    Pxz = ukf_covariance(SIGMAm, Z); % ������-������ ������� n

    K = Pxz/Pz;            % ������� �������� �������
    xhat = xm + K*(y - z); % ���������� ������� ���������
    P = Pm - K*Pz*K';      % ���������� �������������� �������
    P = 0.5*(P+P');        % ����������� �������� ������������

    S=xhat;
    if size(S,2) > 1
        S=S';
    end
     
%--------------------------------------------------

% ������� ����� ����� �������� ��������� � ������������ ������� ���������:
function y = h(xhat) 

    sat=[xhat(1) xhat(2) xhat(3)];
    sat_v=[xhat(4) xhat(5) xhat(6)];
    freqDownlink=0; % ������� ����� ����, �������� � �������� ��������� 
    % ������������ �� ������� �������, � ������� � ��� ��������� ��������
    % ������� �������� 
    [~,a]=dopplerSAT1_downLink(sat,...
                               sat_v,...
                               xyzPost,...
                               freqDownlink);
    if typ==1    
         y=a';
    else
        y=[a sat sat_v]';                         
    end
   
end

%--------------------------------------------------
% ������� ���������� ������� ���������:
function xo = updatex(x,A)
    xo = A*x;
end

%--------------------------------------------------
% ������� ���������� ����� �����:
function X = ukf_sigmapoints(xbar,Pi) 
    nn = size(xbar,1);
    B = chol(nn*Pi)'; % R = chol(B) where R'*R = B
    X = [B -B] + repmat(xbar,1,nn*2);
end
 
%--------------------------------------------------
% ������� ���������� �������� (� ����������) �� ��������� �����-�����:
function [m, P] = ukf_mean_covariance(sigmaX) 
    n_2 = size(sigmaX,2);
    m = sum(sigmaX,2)/n_2;
    P = zeros(size(sigmaX,1)); % size n
    for k = 1:n_2
        P = P + ((sigmaX(:,k)-m)*(sigmaX(:,k)-m)');
    end
    P = P/n_2;
end
 
%--------------------------------------------------
% ������� ���������� �������� ���������� ����� ����� ����������� 
% ����� ����� (��������� � ���������):
function P = ukf_covariance(sigmaX, sigmaY) 
    n_2 = size(sigmaX,2); 
    mx = sum(sigmaX,2)/n_2;
    my = sum(sigmaY,2)/n_2;
    P = zeros(size(sigmaX,1),size(sigmaY,1)); % size n
    for k = 1:n_2
        P = P + ((sigmaX(:,k)-mx)*(sigmaY(:,k)-my)');
    end
    P = P/n_2;
end 
            
end
       
