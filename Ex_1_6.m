omega = 10;
dt = 0.01;
Nc =3;
Np=200;
R_wei_var = 0.5;
R_wei = R_wei_var*eye(Nc);
Xk = [0.1;0.2;0.3];

s = tf('s');
numc = [0 0 omega^2];
denc = [1 0.1*omega omega^2];
[Ac,Bc,Cc,Dc] = tf2ss(numc,denc);

[Ap,Bp,Cp,Dp]=c2dm(Ac,Bc,Cc,Dc,dt);
[Phi_Phi,Phi_F,Phi_R,A_e, B_e,C_e] =mpcgain(Ap,Bp,Cp,Nc,Np);

dleta_U = inv(Phi_Phi+R_wei)*(Phi_R-Phi_F*Xk);

%%%%%  Closed Loop Control System %%%%%
Ky_I = zeros(1,Nc);
Ky_I(1) = 1;
Ky_MPC = zeros(1,Nc);
Ky_MPC(1) = 1;

Ky = Ky_I*inv(Phi_Phi+R_wei)*Phi_R;
Kmpc = Ky_MPC*inv(Phi_Phi+R_wei)*Phi_F;

pole_mpc = eig(A_e-B_e*Kmpc);

H = Phi_Phi+R_wei;
con_num = cond(H);
