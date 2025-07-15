clear all;

a=0.8;
b = 0.1;
Np=10;
Nc=4;

Ae = [a 0;a 1];
Be = [b;b];
Ce = [0 1];

Xk = [0.1; 0.2];
Rs = ones(Np,1);
R_wei = eye(Nc);

F(1,:) = Ce*Ae;
for i = 2:Np
    F(i,:) = F(i-1,:)*Ae;
end

v = zeros(Np,1);
A_phi = eye(size(Ae));
for i = 1:Np
    v(i) = Ce * A_phi * Be; %Ce*A^k-1*Be
    A_phi = A_phi * Ae;
end

Phi = zeros(Np, Nc);
for i = 1:Nc
    Phi(i:end, i) = v(1:Np-i+1);
end

pTp = transpose(Phi)*Phi;
pTF = transpose(Phi)*F;
Rs_F=Rs-F*Xk;

deltaU = inv(pTp+10*R_wei)*transpose(Phi)*Rs_F;

%%%%%  Closed Loop Control System %%%%%
Ky_I = zeros(1,Nc);
Ky_I(1) = 1;
Ky_MPC = zeros(1,Nc);
Ky_MPC(1) = 1;

Ky = Ky_I*inv(pTp+10*R_wei)*transpose(Phi)*Rs;
Kmpc = Ky_MPC*inv(pTp+10*R_wei)*transpose(Phi)*F;

pole_mpc = eig(Ae-Be*Kmpc);
