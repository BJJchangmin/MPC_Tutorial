Jm = 0.01;
Bm = 0.001;
Jl = 1;
Bl = 2;
ks = 7000;
N = 100;

% Ac = [0 1 0 0;-ks/(N^2*Jm) -Bm/Jm ks/(Jm*N) 0; 0 0 0 1; ks/(Jl*N) 0 -ks/Jl -Bl/Jl];
Ac = [0 1 0 0;-ks/(N^2*Jm) -Bm/Jm ks/(Jm*N) 0; 0 0 0 1;ks/(Jl*N) 0 -ks/Jl -Bl/Jl;];
Bc = [0; 0; 0; 1;];
Cc = [1 0 0 0];
Dc = 0;
dt = 0.01;
Nc = 3;
Np= 60;
R_wei_var = 0.1;
R_wei = R_wei_var*eye(Nc);

[Ap,Bp,Cp,Dp]=c2dm(Ac,Bc,Cc,Dc,dt);
[Phi_Phi,Phi_F,Phi_R,A_e, B_e,C_e] =mpcgain(Ap,Bp,Cp,Nc,Np);

%%%%% Observable Check %%%%%
Num_s = 4;
Ob = Cp;
for i = 1:Num_s-1
    Ob = [Ob; Cp * Ap^i];
end
rank_Ob = rank(Ob);

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


