

% Generate trajectories (CML trajectories) with singular destination distribution

clear all
close all
clc

%--------------------

N = 100;

T = 15;

F = [1 T 0 0;...
    0 1 0 0;...
    0 0 1 T;...
    0 0 0 1];


q = 0.01;

Q = [q*T^3/3 q*T^2/2 0 0;...
    q*T^2/2 q*T 0 0;...
    0 0 q*T^3/3 q*T^2/2;...
    0 0 q*T^2/2 q*T];


D = chol(Q,'lower');

%------------------------------

X0_m_m = [2000;5;2000;20];

X0_Cov_m = [100000 40 0 0;
            40 10 0 0;
            0 0 100000 40;
            0 0 40 10];
        
%---------------------------

XN_m_m = [15000;5;2000;-20];

XN_Cov_m = [0 0 0 0;
            0 1 0 0;
            0 0 0 0;
            0 0 0 1];
        
C_N0 = [0 0 0 0;
        0 2 0 0;
        0 0 0 0;
        0 0 0 2];


%---------------------------

iteration = 50;

Xr = zeros(4,N+1);


for iter=1:iteration
    
    %       iter
    
    
%     Dxn = chol(X0_Cov_m,'lower');
%     
%     Xr(:,1) = X0_m_m + Dxn*randn(4,1);
%     
%     Dxn = [0 0 0 0;
%         0 1 0 0;
%         0 0 0 0;
%         0 0 0 1];
%     
%     Xr(:,N+1) = XN_m_m + Dxn*randn(4,1);
    
    %----------------------------------
    
%     Dxn = chol(XN_Cov_m,'lower');
      Dxn = [0 0 0 0;
             0 1 0 0;
             0 0 0 0;
             0 0 0 1];

      Xr(:,N+1) = XN_m_m + Dxn*randn(4,1);
      
      
      C0 = X0_Cov_m - C_N0*pinv(XN_Cov_m)*C_N0';
     Dxn = chol(C0,'lower');

      Xr(:,1) = X0_m_m + C_N0*pinv(XN_Cov_m)*(Xr(:,N+1) - XN_m_m) + Dxn*randn(4,1);
      
      %-------------------------------
    
    for k=1:N-1
        
        
        %----------CN|k
        
        CNk = zeros(4,4);
        
        for ii=0:N-k-1
            
            CNk = CNk + F^ii*Q*(F^ii)';
            
        end
        
        %.....
        Gk = Q - Q*F^(N-k)'/(CNk + F^(N-k)*Q*F^(N-k)')*F^(N-k)*Q;
        DG = chol(Gk,'lower');
        
        Gk_N = Gk*F^(N-k)'/(CNk);
        Gk_km1 = F - Gk_N*F^(N-k+1);
        
        
        Xr(:,k+1) = Gk_km1*Xr(:,k) + Gk_N*Xr(:,N+1) + DG*randn(4,1);
        
        
    end
    
    
    
    figure(1)
    hold on
    plot(Xr(1,:),Xr(3,:),'b')
    grid
    
    
    kk=0:1:100;
    
    figure(2)
    hold on
    plot(kk,Xr(2,:),'b')
    grid
    
    figure(3)
    hold on
    plot(kk,Xr(4,:),'b')
    grid
    
end



