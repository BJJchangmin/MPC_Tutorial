
% A = [1 1 0; 0 1 0; 1 1 1];
% B = [0.5; 1; 0.5];
% C = [0 0 1];

A_e = [1 1 0; 0 1 0; 1 1 1];
B_e = [0.5; 1; 0.5];
C_e = [0 0 1];
Nc =5;
Np=30;

R_wei_var = 0.1;
R_wei = R_wei_var*eye(Nc);
Xk = [0.0;0.2;0.0];


h(1,:)=C_e;
F(1,:)=C_e*A_e;

for kk=2:Np
    h(kk,:)=h(kk-1,:)*A_e;
    F(kk,:)= F(kk-1,:)*A_e;
end

v=h*B_e;
Phi=zeros(Np,Nc); %declare the dimension of Phi
Phi(:,1)=v; % first column of Phi

for i=2:Nc
    Phi(:,i)=[zeros(i-1,1);v(1:Np-i+1,1)]; %Toeplitz matrix
end

BarRs=ones(Np,1);
Phi_Phi= Phi'*Phi;
Phi_F= Phi'*F;
Phi_R=Phi'*BarRs;

dleta_U = inv(Phi_Phi+R_wei)*(Phi_R-Phi_F*Xk);

Ky_I = zeros(1,Nc);
Ky_I(1) = 1;
Ky_MPC = zeros(1,Nc);
Ky_MPC(1) = 1;

Ky = Ky_I*inv(Phi_Phi+R_wei)*Phi_R;
Kmpc = Ky_MPC*inv(Phi_Phi+R_wei)*Phi_F;

N_sim = 100;



Pole = [0.01 0.0105 0.011];
K_ob = place(Ac', Cc', Pole)';
% Kmpc = [1.9685 0.9688 2.9685];

