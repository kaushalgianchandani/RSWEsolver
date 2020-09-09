% rectangular domain. Wind-driven circulation, no diffusion, only
% bototm-friction
% f-plane, constant easterly
% % % % Creates 'wind_driven_t2.gif' as gif file and 'wind_driven_t2.mat' as data file.

close all;
clear
clc

%% Physical constants

Omega = 2 * pi / 86400;
R_earth = 6371*1000;
g = 9.81;
plotting = 1; %1 for not displaying output with every 'plot_rate' interation, 0 for no display

%% Solver parameters
deg2rad = pi/180.0;   % index to convert degree measure to radian
deg2meter = (pi/180.0)*R_earth ; % equivalent distance for 1 degree
lat = 30;                          % Mid-latitude (deg)

% A = 1;                             % Wave amplitude
nl = 0;                            % Linear or non-linear
H0 = 200;                          %  depth (m)
F = 0.2;                         % maximum wind stress (Pa) (1 Pa = 10 dyne/cm2)
rho0 = 1027;
Lx = 10000*1000;                          % domain length (m)
dx = 20*1000;
Ly = 2*pi*1000*1000;                   % domain width (m)
dy = 50*1000;
ah = 1e4;                           % Diffusion constant
bf = 0;

%To be used for Stommel's model
% ah = 0;
% bf = 1/(10*24*60*60);                         % friction

dt = 60;                         % Time step
N = 2000000;                     % Number of iterations
% disp(N*dt/86400);
% N = 1;

plot_rate = 100;                   % sampling every plot_rate iterations
save_data = 0;                     % Save data or not

x = 0:dx:Lx;
y = 0:dy:Ly;

% Grid size
nx = length(x);
ny = length(y);

% Relevant matrix indices.
i = 2 : (nx - 1);
im = i - 1;
ip = i + 1;

j = 2 : (ny - 1);
jm = j - 1;
jp = j + 1;

% Coriolis parameters
% f0 = 2 * Omega * sind(lat);
f0 = 0;
beta = 2 * Omega * cosd(lat) / R_earth;
% beta = 0. ; % comment if want a beta plane
fU = repmat(f0 + beta * (y - 0.5 * dy), nx, 1);
fV = repmat(f0 + beta * y, nx, 1);

% Rossby radius of deformation
R = abs(sqrt(g * H0) / f0);

% Other variables
d2x = 2 * dx;
d4x = 4 * dx;
dx2 = dx^2;
d2y = 2 * dy;
d4y = 4 * dy;
dy2 = dy^2;
tm = 1; tn = 2; tp = 3; % Time indices (tm - previous (minus)step, tn - present (now), tp - future (plus) step
delt = dt;
U = zeros(nx, ny, 3);
V = zeros(nx, ny, 3);
H = zeros(nx, ny, 3);
grav = zeros(nx, ny);
vel = zeros(nx, ny);
adv = zeros(nx, ny);
xadv1 = zeros(nx, ny);
xadv2 = zeros(nx, ny);
xadv = zeros(nx, ny);
diff_x = zeros(nx, ny);
diff_y = zeros(nx, ny);
DH1 = zeros(nx, ny);
DH2 = zeros(nx, ny);

% Stability checks
staba = 1 / dx2 + 1 / dy2;
stabc = bf / H0^2;
stabd = (stabc^2 + f0^2) / stabc;
stab1 = (sqrt(stabc^2 + 16 * g * H0 * staba) - stabc) / (4 * g * H0 * staba);
stab2 = (sqrt(stabd^2 + 8 * g * H0 * staba) - stabd) / (2 * g * H0 * staba);
stab3 = 1 / stabd;
stable = dt < min([stab1 stab2 stab3]);

% For data logging
data_idx = 1;
data_logs = 1 + floor(N / plot_rate) + 1 * (mod(N, plot_rate) ~= 0);
if (save_data)
    dataf = matfile(dataname, 'writable', true);
    dataf.deformation = R;
    dataf.friction = bf;
    dataf.windstress = F;
    dataf.time = zeros(1, data_logs);
    dataf.mean_err = zeros(1, data_logs);
    dataf.x = x(2: (nx - 1));
    dataf.y = y(2: (ny - 1));
    dataf.H = zeros(nx - 2, ny - 2, data_logs);
    dataf.U = zeros(nx - 2, ny - 2, data_logs);
    dataf.V = zeros(nx - 2, ny - 2, data_logs);
    
    
end

if ah == 0
    alpha = bf*(Lx/beta)*(1/Ly^2);
else
    alpha = -ah*(Lx/beta)*(1/Ly^4);
end
disp(alpha)

%% Initial conditions - motionless ocean
% to restart from a previous state
% IniState = load('fileName.mat','H','U','V'); 
% H(:, :, :) = IniState.H;
% U(:,:,:) =  IniState.U;
% V(:,:,:) =  IniState.V;

H(:, :, :) = H0;

normf = 1;

if ~stable
    disp('Vlues do not satisfy stability conditions');
end
hold off;

% added by KG
% wind-stress
taux=-ones(nx,ny)*F.*cos(pi*y/Ly);
tauy(i,j)=0.;

Taux = (1/rho0)*taux;
Tauy = (1/rho0)*tauy;


%% Time stepping loop

for t = 0 : N
    disp(100*t/N);
    % Plotting
    if plotting == 1
        if (mod(t, plot_rate) == 0) || t == N
            
            subplot(121)
            cla;
            contourf(x(i) / 1000, y(j) / 1000, (H(i, j, tn) - H0)',10);
            xlim([0 Lx/1000])
            ylim([0 Ly/1000])
            colorbar;
            set(gca,'FontSize',16)
            
            drawnow;
            
            hold on;
            xlabel('x (km)');
            ylabel('y (km)');
            title('Sea surface height - \eta^* (m)','FontWeight','normal')
                       

            u1 = (U(i,j,tn)./H(i, j, tn));
            v1 = (V(i,j,tn)./H(i, j, tn));
            
            subplot(122)
            cla;
            contourf(x(i)/1000,y(j)/1000,-dx*cumsum(v1,1)',15)
            shading flat
            xlabel('x (km)');
            ylabel('y (km)');
            title("Stream function - \psi^* (m^2 s^{-1})",'FontWeight','normal')
            colorbar
            set(gca,'FontSize',16)
            xlim([0 Lx/1000])
           
            
            if (save_data)
                %             % GIF creation
                %             frame = getframe(1);
                %             img = frame2im(frame);
                %             [imind, cm] = rgb2ind(img, 256);
                %             if (t == 0)
                %                 imwrite(imind, cm, gifname, 'gif', 'Loopcount', inf, 'DelayTime', 0);
                %             else
                %                 imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime', 0);
                %             end
                
                % Data logging
                dataf.time(1, data_idx) = t * dt;
                %            dataf.mean_err(1, data_idx) = diff;
                dataf.H(i, j, data_idx) = H(i, j, tn);
                dataf.U(i, j, data_idx) = U(i, j, tn);
                dataf.V(i, j, data_idx) = V(i, j, tn);
                data_idx = data_idx + 1;
                %             data.x(i) = x(i);
                %             data.y(j) = y(j);
                t*dt/86400
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Advancing H.        psi_x(i,:) = psi_x(i-1,:) - V1.*(x(i) - x(i-1));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    H(i, j, tp) = H(i, j, tm) - delt * ((U(i, j, tn) - U(im, j, tn)) / dx + (V(i, j, tn) - V(i, jm, tn)) / dy);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Advancing U.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    grav(i, j) = (H(ip, j, tn) - H(i, j, tn)) .* (H(ip, j, tn) + H(i, j, tn)) / d2x; % Gravitation
    vel(i, j) = (V(i, j, tn) + V(i, jm, tn) + V(ip, j, tn) + V(ip, jm, tn)) / 4; % Average velocity
    adv(i, j) = ((U(ip, j, tn) + U(i, j, tn)).^2 ./ H(ip, j, tn) - (U(i, j, tn) + U(im, j, tn)).^2 ./ H(i, j, tn)) / d4x; % Advection
    
    % Cross advection
    DH1(i, j) = H(i, j, tn) + H(i, jp, tn) + H(ip, j, tn) + H(ip, jp, tn);
    DH2(i, j) = H(i, j, tn) + H(i, jm, tn) + H(ip, j, tn) + H(ip, jm, tn);
    xadv1(i, j) = (U(i, jp, tn) + U(i, j, tn)) .* (V(i, j, tn) + V(ip, j, tn)) ./ DH1(i, j);
    xadv1(:, ny - 1) = 0;
    xadv2(i, j) = (U(i, j, tn) + U(i, jm, tn)) .* (V(ip, jm, tn) + V(i, jm, tn)) ./ DH2(i, j);
    xadv2(:, 2) = 0;
    xadv(i, j) = (xadv1(i, j) - xadv2(i, j)) / dy;
    % Diffusion
    diff_x(i, j) = (U(im, j, tn) - 2 * U(i, j, tn) + U(ip, j, tn)) / dx2;
    diff_y(i, j) = (U(i, jm, tn) - 2 * U(i, j, tn) + U(i, jp, tn)) / dy2;
    diff_y(:, 2) = 0; diff_y(:, ny - 1) = 0;
    
    % Finite difference equation
    U(i, j, tp) = U(i, j, tm) - delt *(g * grav(i, j) - fU(i, j) .* vel(i, j) + bf * U(i, j, tn) - ah * (diff_x(i, j) + diff_y(i, j)) - Taux(i,j)) - ...
        nl * delt * (adv(i, j) + xadv(i, j));
    
    U(1, :, tp) = 0; U(nx-1, :, tp) = 0; % Rigid boundaries in zonal direction
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Advancing V.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    grav(i, j) = (H(i, jp, tn) - H(i, j, tn)) .* (H(i, jp, tn) + H(i, j, tn)) / d2y; % Gravitation
    vel(i, j) = (U(i, j, tn) + U(im, j, tn) + U(i, jp, tn) + U(im, jp, tn)) / 4; % Average velocity
    adv(i, j) = ((V(i, jp, tn) + V(i, j, tn)).^2 ./ H(i, jp, tn) - (V(i, j, tn) + V(i, jm, tn)).^2 ./ H(i, j, tn)) / d4y; % Advection
    
    % Cross advection
    DH1(i, j) = H(i, j, tn) + H(i, jp, tn) + H(ip, j, tn) + H(ip, jp, tn);
    DH2(i, j) = H(i, j, tn) + H(im, j, tn) + H(i, jp, tn) + H(im, jp, tn);
    xadv1(i, j) = (U(i, jp, tn) + U(i, j, tn)) .* (V(i, j, tn) + V(ip, j, tn)) ./ DH1(i, j);
    xadv1(nx-1, :) = 0;
    xadv2(i, j) = (U(im, j, tn) + U(im, jp, tn)) .* (V(i, j, tn) + V(im, j, tn)) ./ DH2(i, j);
    xadv2(2, :) = 0;
    
    xadv(i, j) = (xadv1(i, j) - xadv2(i, j)) / dx;
    
    % Diffusion
    diff_x(i, j) = (V(im, j, tn) - 2 * V(i, j, tn) + V(ip, j, tn)) / dx2;
    diff_y(i, j) = (V(i, jm, tn) - 2 * V(i, j, tn) + V(i, jp, tn)) / dy2;
    diff_y(:, 2) = 0; diff_y(:, ny - 2) = 0;
    
    % Finite difference equation
    V(i, j, tp) = V(i, j, tm) - delt *(g * grav(i, j) + fV(i, j) .* vel(i, j) + bf * V(i, j, tn) - ah * (diff_x(i, j) + diff_y(i, j)) - Tauy(i,j)) - ...
        nl * delt * (adv(i, j) + xadv(i, j));
    
    V(:, 1, tp) = 0; V(:, ny - 1, tp) = 0; % Rigid boundaries in
    % meridional direction
    
    
    % Robert-Asselin Filter.
    %   if (t > 5)
    U(:, :, tn) = U(:, :, tn) + 0.1 * (U(:, :, tp) - 2 * U(:, :, tn) + U(:, :, tm));
    V(:, :, tn) = V(:, :, tn) + 0.1 * (V(:, :, tp) - 2 * V(:, :, tn) + V(:, :, tm));
    H(:, :, tn) = H(:, :, tn) + 0.1 * (H(:, :, tp) - 2 * H(:, :, tn) + H(:, :, tm));
    %    end
    
    %sum(sum(H(:,:,tp)))
    
    
    % Moving on in time.
    delt = 2 * dt;
    tm = mod(tm, 3) + 1;
    tn = mod(tn, 3) + 1;
    tp = mod(tp, 3) + 1;
    
end
