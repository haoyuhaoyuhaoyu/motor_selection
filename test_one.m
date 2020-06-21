%% init
clear; close all; clc;
restoredefaultpath;
gait_path = 'gait';
addpath(gait_path);
casia_gait = load('gait/casia_gait_0.4.mat');

%% param init
transmission_range = 100;
motor_type = 1;
switch motor_type
    case 1 %TQ 8526
        Kt = 0.127;
        Km = 0.254;
        Im = 0.000115;
        bm = 0;
        %Tm_max = 8.3;
        Tm_max = 8.3;
        joint_index = 10;
        u_index = 4;
        cassie_trans = 37;
    case 2 %7615
        Kt = 0.131;
        Tm_max = 2.86;
        joint_index = [7,8,14,15];
        u_index = [1,2,6,7];
    case 3 %6013
        Kt = 0.11;
        Tm_max = 1.37;
        joint_index = [13,20];
        u_index = [5,10];
end 

%% compute load rms
u_all = [casia_gait.gait(1).inputs.u, casia_gait.gait(3).inputs.u];
dx_all = [casia_gait.gait(1).states.dx, casia_gait.gait(3).states.dx];
ddx_all = [casia_gait.gait(1).states.ddx, casia_gait.gait(3).states.ddx];

u = cassie_trans*u_all(u_index, :);
dx = dx_all(joint_index, :);
ddx = ddx_all(joint_index, :);
%%
[TL_max, TL_max_dx_index] = max(abs(u));
dx_max = max(abs(dx));
TL_rms = rms(u);
ddx_rms = rms(ddx);

%% compute total E
C1 = sum(u*dx', 2);
C2 = sum(dx*ddx', 2);
C3 = sum(dx*dx', 2);
C4 = sum(ddx*ddx', 2);
C5 = sum(u*u', 2);
C6 = sum(u*ddx', 2);
C7 = sum(u*dx', 2);


E_tot_all = [];
for belta = ceil(TL_max/Tm_max):0.2:transmission_range
    E_tot = C1 + Im*(1+2*bm/Km^2)*belta^2*C2 +...
        bm*(1+bm/Km^2)*belta^2*C3 + ...
        Im^2/Km^2*belta^2*C4 + 1/(Km^2*belta^2)*C5 +...
        2*Im/Km^2*C6 + 2*bm/Km^2*C7;    
    E_tot_all = [E_tot_all, E_tot];
end

%% compute bandwidth
w_all = [];
for belta = ceil(TL_max/Tm_max):0.2:transmission_range
    w = sqrt((Tm_max/belta-TL_max/belta^2)/(dx_max^2)-bm^2)/Im;
    w_all = [w_all, w];
end

plot(w_all(1,:),E_tot_all(1,:))

belta1 = sqrt(TL_rms/(Im*ddx_rms));
belta2 = 2*TL_max/Tm_max;

