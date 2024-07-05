

function t = planerThreeLinksDynamics()
% t = planerThreeLinksDynamics();
% tau = round(double(subs(t)));

%% calculate P
syms q1 a1 l1 q2 a2 l2 q3 a3 l3

P0 = [0, 0, 0];

A02_ = trotz(q1)*transl(l1, 0, 0);
Pl1 = A02_(1:3,4);

A01 = trotz(q1)*transl(a1, 0, 0);
P1 = A01(1:3,4);

A12_ = trotz(q2)*transl(l2, 0, 0);
A02_ = A01 * A12_;
Pl2 = A02_(1:3,4);

A12 = trotz(q2)*transl(a2, 0, 0);
A02 = A01 * A12;
P2 = A02(1:3,4);


A23_ = trotz(q3)*transl(l3, 0, 0);
A23 = trotz(q3)*transl(a3, 0, 0);
A03_ = A02 * A23_;
Pl3 = A03_(1:3,4);
A03 = A02 * A23;
P3 = A03(1:3,4);

Pl1 = simplify(Pl1);
P1 = simplify(P1);

Pl2 = simplify(Pl2);
P2 = simplify(P2);

Pl3 = simplify(Pl3);
% P3 = simplify(P3);

%% calculate JlP

z = [0, 0, 1];
Jl1P1 = func(z, (Pl1 - P0));

Jl2P1 = func(z, (Pl2 - P0));
Jl2P2 = func(z, (Pl2 - P1));

Jl3P1 = func(z, (Pl3 - P0));
Jl3P2 = func(z, (Pl3 - P1));
Jl3P3 = func(z, (Pl3 - P2));

Jl1P = [Jl1P1, zeros(3,1), zeros(3,1)];
Jl2P = [Jl2P1, Jl2P2, zeros(3,1)];
Jl3P = [Jl3P1, Jl3P2, Jl3P3];


%% calculate JlO

Jl1O1 = [0; 0; 1];
Jl1O = [Jl1O1, zeros(3,1), zeros(3,1)];

Jl2O1 = [0; 0; 1];
Jl2O2 = [0; 0; 1];
Jl2O = [Jl2O1, Jl2O2, zeros(3,1)];

Jl3O1 = [0; 0; 1];
Jl3O2 = [0; 0; 1];
Jl3O3 = [0; 0; 1];
Jl3O = [Jl3O1, Jl3O2, Jl3O3];

%% calculate JmP
Pm1 = [0, 0, 0];
Jm1P1 = func(z, (Pm1-P0));

Pm2 = P1;
Jm2P1 = func(z, (Pm2-P0));
Jm2P2 = func(z, (Pm2-P1));

Pm3 = P2;
Jm3P1 = func(z, (Pm3-P0));
Jm3P2 = func(z, (Pm3-P1));
Jm3P3 = func(z, (Pm3-P2));

Jm1P = [Jm1P1, zeros(3,1), zeros(3,1)];
Jm2P = [Jm2P1, Jm2P2, zeros(3,1)];
Jm3P = [Jm3P1, Jm3P2, Jm3P3];

%% calculate JmO

syms k1 k2 k3

Jm1O1 = k1 * z';

Jm2O1 = Jl2O1;
Jm2O2 = k2 * z';

Jm3O1 = Jl3O1;
Jm3O2 = Jl3O2;
Jm3O3 = k3 * z';

Jm1O = [Jm1O1, zeros(3,1), zeros(3,1)];
Jm2O = [Jm2O1, Jm2O2, zeros(3,1)];
Jm3O = [Jm3O1, Jm3O2, Jm3O3];


%% calaulate B

syms ml1 ml2 ml3 mm1 mm2 mm3 Il1 Im1 Il2 Im2 Il3 Im3
R1 = A01(1:3,1:3);
Rm1 = eye(3);
R2 = A02(1:3,1:3);
% Rm7 = A56(1:3,1:3);
Rm1 = eye(3);
R3 = A03(1:3,1:3);
% Rm8 = A57(1:3,1:3);
Rm3 = eye(3);

% for linear velocity

B = ml1*(Jl1P.'*Jl1P) + Jl1O.'*R1*Il1*R1.'*Jl1O + mm1*(Jm1P.'*Jm1P) +  Jm1O.'*Rm1*Im1*Rm1.'*Jm1O + ...
    ml2*(Jl2P.'*Jl2P) + Jl2O.'*R2*Il2*R2.'*Jl2O + mm2*(Jm2P.'*Jm2P) +  Jm2O.'*Rm1*Im2*Rm1.'*Jm2O + ...
    ml3*(Jl3P.'*Jl3P) + Jl3O.'*R3*Il3*R3.'*Jl3O + mm3*(Jm3P.'*Jm3P) +  Jm3O.'*Rm3*Im3*Rm3.'*Jm3O;


b11 = collect(simplify(B(1,1)));
b12 = collect(simplify(B(1,2)));
b13 = collect(simplify(B(1,3)));
b21 = collect(simplify(B(2,1)));
b22 = collect(simplify(B(2,2)));
b23 = collect(simplify(B(2,3)));
b31 = collect(simplify(B(3,1)));
b32 = collect(simplify(B(3,2)));
b33 = collect(simplify(B(3,3)));

B = [b11, b12, b13;
    b21, b22, b23;
    b31, b32, b33];

%% calculate C

syms dq1 dq2 dq3

c111 = 1/2*(diff(b11,q1)+diff(b11,q1)-diff(b11,q1));
c112 = 1/2*(diff(b11,q2)+diff(b12,q1)-diff(b12,q1));
c113 = 1/2*(diff(b11,q3)+diff(b13,q1)-diff(b13,q1));

c11 = c111*dq1 + c112*dq2 + c113*dq3;
c11 = collect(simplify(c11));

c121 = 1/2*(diff(b12,q1)+diff(b11,q2)-diff(b21,q1));
c122 = 1/2*(diff(b12,q2)+diff(b12,q2)-diff(b22,q1));
c123 = 1/2*(diff(b12,q3)+diff(b13,q2)-diff(b23,q1));

c12 = c121*dq1 + c122*dq2 + c123*dq3;
c12 = collect(simplify(c12));

c131 = 1/2*(diff(b13,q1)+diff(b11,q3)-diff(b31,q1));
c132 = 1/2*(diff(b13,q2)+diff(b12,q3)-diff(b32,q1));
c133 = 1/2*(diff(b13,q3)+diff(b13,q3)-diff(b33,q1));

c13 = c131*dq1 + c132*dq2 + c133*dq3;
c13 = collect(simplify(c13));

c211 = 1/2*(diff(b21,q1)+diff(b21,q1)-diff(b11,q2));
c212 = 1/2*(diff(b21,q2)+diff(b22,q1)-diff(b12,q2));
c213 = 1/2*(diff(b21,q3)+diff(b23,q1)-diff(b13,q2));

c21 = c211*dq1 + c212*dq2 + c213*dq3;
c21 = collect(simplify(c21));

c221 = 1/2*(diff(b22,q1)+diff(b21,q2)-diff(b21,q2));
c222 = 1/2*(diff(b22,q2)+diff(b22,q2)-diff(b22,q2));
c223 = 1/2*(diff(b22,q3)+diff(b23,q2)-diff(b22,q2));

c22 = c221*dq1 + c222*dq2 + c223*dq3;
c22 = collect(simplify(c22));

c231 = 1/2*(diff(b23,q1)+diff(b21,q3)-diff(b31,q2));
c232 = 1/2*(diff(b23,q2)+diff(b22,q3)-diff(b32,q2));
c233 = 1/2*(diff(b23,q3)+diff(b23,q3)-diff(b33,q2));

c23 = c231*dq1 + c232*dq2 + c233*dq3;
c23 = collect(simplify(c23));

c311 = 1/2*(diff(b31,q1)+diff(b31,q1)-diff(b11,q3));
c312 = 1/2*(diff(b31,q2)+diff(b32,q1)-diff(b12,q3));
c313 = 1/2*(diff(b31,q3)+diff(b33,q1)-diff(b13,q3));

c31 = c311*dq1 + c312*dq2 + c313*dq3;
c31 = collect(simplify(c31));

c321 = 1/2*(diff(b32,q1)+diff(b31,q2)-diff(b21,q3));
c322 = 1/2*(diff(b32,q2)+diff(b32,q2)-diff(b22,q3));
c323 = 1/2*(diff(b32,q3)+diff(b33,q2)-diff(b23,q3));

c32 = c321*dq1 + c322*dq2 + c323*dq3;
c32 = collect(simplify(c32));

c331 = 1/2*(diff(b33,q1)+diff(b31,q3)-diff(b31,q3));
c332 = 1/2*(diff(b33,q2)+diff(b32,q3)-diff(b32,q3));
c333 = 1/2*(diff(b33,q3)+diff(b33,q3)-diff(b33,q3));

c33 = c331*dq1 + c332*dq2 + c333*dq3;
c33 = collect(simplify(c33));

C = [c11, c12, c13;
    c21, c22, c23;
    c31, c32, c33];

%% dB
dB = B;

for i = 1:3
    for j = 1:3
        dB(i,j) = diff(B(i,j), q1) * dq1 + diff(B(i,j), q2) * dq2 + diff(B(i,j), q3) * dq3;
    end
end

%% calculate N
N = collect(simplify(dB - 2*C));

A = collect(simplify(N+N.')); % A = 0 -> N is skew-symmetric, the development of B and C are right

%% calculate potential energy
syms g0
g = [0; -g0; 0];

G1 = ml1*g.'*Jl1P1 + mm1*g.'*Jm1P1 + ...
    ml2*g.'*Jl2P1 + mm2*g.'*Jm2P1 + ...
    ml3*g.'*Jl3P1 + mm3*g.'*Jm3P1;

G2 = ml2*g.'*Jl2P2 + mm2*g.'*Jm2P2 + ...
    ml3*g.'*Jl3P2 + mm3*g.'*Jm3P2;

G3 = ml3*g.'*Jl3P3 + mm3*g.'*Jm3P3;

G = [G1; G2; G3];

%% calculate torque

syms ddq1 ddq2 ddq3

ddq = [ddq1; ddq2; ddq3];
dq = [dq1; dq2; dq3];

t = collect(B*ddq + C*dq + G);




















