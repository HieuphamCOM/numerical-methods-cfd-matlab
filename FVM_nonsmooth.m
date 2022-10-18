%% FVM for 1D linear scalar conservation laws : NON-SMOOTH cases
%% q_t + u*q_x = 0, x \in (-1,1)
%% IC: q(x,0) = q0(x) = sin(pi*x)
%% BC: periodic boundary condition
%%

clear all
close all
clc

% problem setting
xa = -1; xb = 1;    % computational domain
ubar = 1;
f = @(q) ubar * q;    % flux function 

% spartial grid data
nx = 100;
dx = (xb - xa) / nx;
xnode = (xa:dx:xb);
for j = 1:nx
    x(j).left = xnode(j);
    x(j).right = xnode(j+1);
    x(j).center = 0.5 * ( x(j).left + x(j).right );
    x(j).dx = x(j).right - x(j).left;
end
xcenter = zeros(nx,1);
for j = 1:nx
   xcenter(j) =  x(j).center;
end

% Riemann initial condition
ql = 1; qr = -0.5;
q0 = zeros(nx,1);
qexact = zeros(nx,1);
for j = 1:nx
    if xcenter(j) < 0
        q0(j) = ql;
    else
        q0(j) = qr;
    end
end
dq0 = @(x) q0 * x; 

% temporal grid data
tfinal = 0.7;
cfl = 0.5;
dt = cfl * dx / ubar;
t = (0:dt:tfinal);
if t(end) < tfinal
    t = [t, tfinal];
end
nt = length(t);
dt = dt * ones(nt,1);
dt(end) = t(end) - t(end-1);

% FVM setting
% upwind scheme
QcurrentUP = zeros(nx,1);   % FV solution at current time step
QnextUP = zeros(nx,1);      % FV solution at next time step
% Lax-Friedrichs scheme
QcurrentLxF = zeros(nx,1);   % FV solution at current time step
QnextLxF = zeros(nx,1);      % FV solution at next time step
% Lax-Wendroff scheme
QcurrentLW = zeros(nx,1);   % FV solution at current time step
QnextLW = zeros(nx,1);      % FV solution at next time step

% high-resolution flux-limiter TVD schemes
% minmod
QcurrentMM = zeros(nx,1);   % FV solution at current time step
QnextMM = zeros(nx,1);      % FV solution at next time step
% superbee
QcurrentBEE = zeros(nx,1);   % FV solution at current time step
QnextBEE = zeros(nx,1);      % FV solution at next time step
% MC
QcurrentMC = zeros(nx,1);   % FV solution at current time step
QnextMC = zeros(nx,1);      % FV solution at next time step

% initial condition
for j = 1:nx
    QcurrentUP(j) = q0(j);
    QcurrentLxF(j) = q0(j);
    QcurrentLW(j) = q0(j);
    QcurrentMM(j) = q0(j);
    QcurrentBEE(j) = q0(j);
    QcurrentMC(j) = q0(j);
end

% time looping
figure;
for n = 2:nt
    % fixed BC
    Qnp1UP = qr;
    Q0UP = ql;
    Qnp1LxF = qr;
    Q0LxF = ql;
    Qnp1LW = qr;
    Q0LW = ql;
    Qnp1MM = qr;
    Qnp2MM = qr;
    Q0MM = ql;
    Qm1MM = ql;
    Qnp1BEE = qr;
    Qnp2BEE = qr;
    Q0BEE = ql;
    Qm1BEE = ql;
    Qnp1MC = qr;
    Qnp2MC = qr;
    Q0MC = ql;
    Qm1MC = ql;
    
    %=== 1st-order upwind FV ===
    % Qnext(j) = Qcurrent(j) - dt/dx*( Fcurrent(j+1/2) - Fcurrent(j-1/2) )
    % Fcurrent(j+1/2) = Qcurrent(j)   
    Fright = f(QcurrentUP(1));
    Fleft = f(Q0UP);
    QnextUP(1) = QcurrentUP(1) - dt(n)/x(1).dx * ( Fright - Fleft );
    for j = 2:nx
        % numerical flux functions
        Fright = f(QcurrentUP(j));
        Fleft = f(QcurrentUP(j-1));
        
        % FVM:
        QnextUP(j) = QcurrentUP(j) - dt(n)/x(j).dx * ( Fright - Fleft );
    end    
    
    %=== 1st-order Lax-Friedrichs FV ===
    % Qnext(j) = Qcurrent(j) - dt/dx*( Fcurrent(j+1/2) - Fcurrent(j-1/2) )
    % Fcurrent(j+1/2) = 0.5*(f(Q(j)) + f(Q(j+1))) - 0.5*dx/dt*(Q(j+1) - Q(j))
    Fright = 0.5*( f(QcurrentLxF(1)) + f(QcurrentLxF(2)) ) - 0.5*x(1).dx/dt(n)*( QcurrentLxF(2) - QcurrentLxF(1) );
    Fleft = 0.5*( f(Q0LxF) + f(QcurrentLxF(1)) ) - 0.5*x(1).dx/dt(n)*( QcurrentLxF(1) - Q0LxF );
    QnextLxF(1) = QcurrentLxF(1) - dt(n)/x(1).dx * ( Fright - Fleft );
    for j = 2:nx-1
        % numerical flux functions
        Fright = 0.5*( f(QcurrentLxF(j)) + f(QcurrentLxF(j+1)) ) - 0.5*x(j).dx/dt(n)*( QcurrentLxF(j+1) - QcurrentLxF(j) );
        Fleft = 0.5*( f(QcurrentLxF(j-1)) + f(QcurrentLxF(j)) ) - 0.5*x(j).dx/dt(n)*( QcurrentLxF(j) - QcurrentLxF(j-1) );
        
        % FVM:
        QnextLxF(j) = QcurrentLxF(j) - dt(n)/x(j).dx * ( Fright - Fleft );
    end    
    Fright = 0.5*( f(QcurrentLxF(nx)) + f(Qnp1LxF) ) - 0.5*x(nx).dx/dt(n)*( Qnp1LxF - QcurrentLxF(nx) );
    Fleft = 0.5*( f(QcurrentLxF(nx-1)) + f(QcurrentLxF(nx)) ) - 0.5*x(nx).dx/dt(n)*( QcurrentLxF(nx) - QcurrentLxF(nx-1) );
    QnextLxF(nx) = QcurrentLxF(nx) - dt(n)/x(nx).dx * ( Fright - Fleft );
    
    %=== 2nd-order Lax-Wendroff ===
    % Qnext(j) = Qcurrent(j) - dt/dx*( Fcurrent(j+1/2) - Fcurrent(j-1/2) )
    % Fcurrent(j+1/2) = 0.5*(f(Q(j)) + f(Q(j+1))) - 0.5*ubar*dt/dx*(f(Q(j+1)) - f(Q(j)))
    Fright = 0.5*( f(QcurrentLW(1)) + f(QcurrentLW(2)) ) - 0.5*ubar*dt(n)/x(1).dx*( f(QcurrentLW(2)) - f(QcurrentLW(1)) );
    Fleft = 0.5*( f(Q0LW) + f(QcurrentLW(1)) ) - 0.5*ubar*dt(n)/x(j).dx*( f(QcurrentLW(1)) - f(Q0LW) );
    QnextLW(1) = QcurrentLW(1) - dt(n)/x(1).dx * ( Fright - Fleft );
    for j = 2:nx-1
        % numerical flux functions
        Fright = 0.5*( f(QcurrentLW(j)) + f(QcurrentLW(j+1)) ) - 0.5*ubar*dt(n)/x(j).dx*( f(QcurrentLW(j+1)) - f(QcurrentLW(j)) );
        Fleft = 0.5*( f(QcurrentLW(j-1)) + f(QcurrentLW(j)) ) - 0.5*ubar*dt(n)/x(j).dx*( f(QcurrentLW(j)) - f(QcurrentLW(j-1)) );
        
        % FVM:
        QnextLW(j) = QcurrentLW(j) - dt(n)/x(j).dx * ( Fright - Fleft );
    end 
    Fright = 0.5*( f(QcurrentLW(nx)) + f(Qnp1LW) ) - 0.5*ubar*dt(n)/x(nx).dx*( f(Qnp1LW) - f(QcurrentLW(nx)) );
    Fleft = 0.5*( f(QcurrentLW(nx-1)) + f(QcurrentLW(nx)) ) - 0.5*ubar*dt(n)/x(nx).dx*( f(QcurrentLW(nx)) - f(QcurrentLW(nx-1)) );
    QnextLW(nx) = QcurrentLW(nx) - dt(n)/x(nx).dx * ( Fright - Fleft );
    
    
    %=== High-resolution TVD schemes ===
    %% minmod:
    % smoothness indicators
    dQright = QcurrentMM(1) - QcurrentMM(1);
    dQleft = QcurrentMM(1) - Q0MM;
    theta = dQleft / dQright;
    % flux limiter
    phi = max( 0, min(1, theta) );
    delta = phi * dQright;
    Fright = f(QcurrentMM(1)) + 0.5*ubar*(1 - ubar*dt(n)/x(1).dx)*delta;

    % smoothness indicators
    dQright = QcurrentMM(1) - Q0MM;
    dQleft = Q0MM - Qm1MM;
    theta = dQleft / dQright;
    % flux limiter
    phi = max( 0, min(1, theta) );
    delta = phi * dQright;
    Fleft = f(Q0MM) + 0.5*ubar*(1 - ubar*dt(n)/x(1).dx)*delta;
    % FVM:
    QnextMM(1) = QcurrentMM(1) - dt(n)/x(1).dx * ( Fright - Fleft );
    
    % smoothness indicators
    dQright = QcurrentMM(3) - QcurrentMM(2);
    dQleft = QcurrentMM(2) - QcurrentMM(1);
    theta = dQleft / dQright;
    % flux limiter
    phi = max( 0, min(1, theta) );
    delta = phi * dQright;
    Fright = f(QcurrentMM(2)) + 0.5*ubar*(1 - ubar*dt(n)/x(2).dx)*delta;

    % smoothness indicators
    dQright = QcurrentMM(2) - QcurrentMM(1);
    dQleft = QcurrentMM(1) - Q0MM;
    theta = dQleft / dQright;
    % flux limiter
    phi = max( 0, min(1, theta) );
    delta = phi * dQright;
    Fleft = f(QcurrentMM(1)) + 0.5*ubar*(1 - ubar*dt(n)/x(2).dx)*delta;
    % FVM:
    QnextMM(2) = QcurrentMM(2) - dt(n)/x(2).dx * ( Fright - Fleft );
    
    for j = 3:nx-1
        % smoothness indicators
        dQright = QcurrentMM(j+1) - QcurrentMM(j);
        dQleft = QcurrentMM(j) - QcurrentMM(j-1);
        theta = dQleft / dQright;
        % flux limiter
        phi = max( 0, min(1, theta) );  % minmod limiter
        delta = phi * dQright;
        Fright = f(QcurrentMM(j)) + 0.5*ubar*(1 - ubar*dt(n)/x(j).dx)*delta;
        
        % smoothness indicators
        dQright = QcurrentMM(j) - QcurrentMM(j-1);
        dQleft = QcurrentMM(j-1) - QcurrentMM(j-2);
        theta = dQleft / dQright;
        % flux limiter
        phi = max( 0, min(1, theta) );
        delta = phi * dQright;
        Fleft = f(QcurrentMM(j-1)) + 0.5*ubar*(1 - ubar*dt(n)/x(j).dx)*delta;
        
        % FVM:
        QnextMM(j) = QcurrentMM(j) - dt(n)/x(j).dx * ( Fright - Fleft );
    end
    % smoothness indicators
    dQright = Qnp1MM - QcurrentMM(nx);
    dQleft = QcurrentMM(nx) - QcurrentMM(nx-1);
    theta = dQleft / dQright;
    % flux limiter
    phi = max( 0, min(1, theta) );
    delta = phi * dQright;
    Fright = f(QcurrentMM(nx)) + 0.5*ubar*(1 - ubar*dt(n)/x(nx).dx)*delta;

    % smoothness indicators
    dQright = QcurrentMM(nx) - QcurrentMM(nx-1);
    dQleft = QcurrentMM(nx-1) - QcurrentMM(nx-2);
    theta = dQleft / dQright;
    % flux limiter
    phi = max( 0, min(1, theta) );
    delta = phi * dQright;
    Fleft = f(QcurrentMM(nx-1)) + 0.5*ubar*(1 - ubar*dt(n)/x(nx).dx)*delta;
    % FVM:
    QnextMM(nx) = QcurrentMM(nx) - dt(n)/x(nx).dx * ( Fright - Fleft );
    
    %% superbee
    % smoothness indicators
    dQright = QcurrentBEE(1) - QcurrentBEE(1);
    dQleft = QcurrentBEE(1) - Q0BEE;
    theta = dQleft / dQright;
    % flux limiter
    phi = max( 0, max( min(2*theta,1), min(theta,2) ) );   % superbee limiter
    delta = phi * dQright;
    Fright = f(QcurrentBEE(1)) + 0.5*ubar*(1 - ubar*dt(n)/x(1).dx)*delta;

    % smoothness indicators
    dQright = QcurrentBEE(1) - Q0BEE;
    dQleft = Q0BEE - Qm1BEE;
    theta = dQleft / dQright;
    % flux limiter
    phi = max( 0, max( min(2*theta,1), min(theta,2) ) );   % superbee limiter
    delta = phi * dQright;
    Fleft = f(Q0BEE) + 0.5*ubar*(1 - ubar*dt(n)/x(1).dx)*delta;
    % FVM:
    QnextBEE(1) = QcurrentBEE(1) - dt(n)/x(1).dx * ( Fright - Fleft );
    
    % smoothness indicators
    dQright = QcurrentBEE(3) - QcurrentBEE(2);
    dQleft = QcurrentBEE(2) - QcurrentBEE(1);
    theta = dQleft / dQright;
    % flux limiter
    phi = max( 0, max( min(2*theta,1), min(theta,2) ) );   % superbee limiter
    delta = phi * dQright;
    Fright = f(QcurrentBEE(2)) + 0.5*ubar*(1 - ubar*dt(n)/x(2).dx)*delta;

    % smoothness indicators
    dQright = QcurrentBEE(2) - QcurrentBEE(1);
    dQleft = QcurrentBEE(1) - Q0BEE;
    theta = dQleft / dQright;
    % flux limiter
    phi = max( 0, max( min(2*theta,1), min(theta,2) ) );   % superbee limiter
    delta = phi * dQright;
    Fleft = f(QcurrentBEE(1)) + 0.5*ubar*(1 - ubar*dt(n)/x(2).dx)*delta;
    % FVM:
    QnextBEE(2) = QcurrentBEE(2) - dt(n)/x(2).dx * ( Fright - Fleft );
    
    for j = 3:nx-1
        % smoothness indicators
        dQright = QcurrentBEE(j+1) - QcurrentBEE(j);
        dQleft = QcurrentBEE(j) - QcurrentBEE(j-1);
        theta = dQleft / dQright;
        % flux limiter
        phi = max( 0, max( min(2*theta,1), min(theta,2) ) );   % superbee limiter
        delta = phi * dQright;
        Fright = f(QcurrentBEE(j)) + 0.5*ubar*(1 - ubar*dt(n)/x(j).dx)*delta;
        
        % smoothness indicators
        dQright = QcurrentBEE(j) - QcurrentBEE(j-1);
        dQleft = QcurrentBEE(j-1) - QcurrentBEE(j-2);
        theta = dQleft / dQright;
        % flux limiter
        phi = max( 0, max( min(2*theta,1), min(theta,2) ) );   % superbee limiter
        delta = phi * dQright;
        Fleft = f(QcurrentBEE(j-1)) + 0.5*ubar*(1 - ubar*dt(n)/x(j).dx)*delta;
        
        % FVM:
        QnextBEE(j) = QcurrentBEE(j) - dt(n)/x(j).dx * ( Fright - Fleft );
    end
    % smoothness indicators
    dQright = Qnp1BEE - QcurrentBEE(nx);
    dQleft = QcurrentBEE(nx) - QcurrentBEE(nx-1);
    theta = dQleft / dQright;
    % flux limiter
    phi = max( 0, max( min(2*theta,1), min(theta,2) ) );   % superbee limiter
    delta = phi * dQright;
    Fright = f(QcurrentBEE(nx)) + 0.5*ubar*(1 - ubar*dt(n)/x(nx).dx)*delta;

    % smoothness indicators
    dQright = QcurrentBEE(nx) - QcurrentBEE(nx-1);
    dQleft = QcurrentBEE(nx-1) - QcurrentBEE(nx-2);
    theta = dQleft / dQright;
    % flux limiter
    phi = max( 0, max( min(2*theta,1), min(theta,2) ) );   % superbee limiter
    delta = phi * dQright;
    Fleft = f(QcurrentBEE(nx-1)) + 0.5*ubar*(1 - ubar*dt(n)/x(nx).dx)*delta;
    % FVM:
    QnextBEE(nx) = QcurrentBEE(nx) - dt(n)/x(nx).dx * ( Fright - Fleft );
    
    %=== exact solution ===
    for j = 1:nx
        if ( xcenter(j) - ubar * t(n) ) < 0
            qexact(j) = ql;
        else
            qexact(j) = qr;
        end
    end
    
    
    %=== plotting current solution ===
%     plot(xcenter, QcurrentUP,'-db',...
%          xcenter, QcurrentLxF,'-vr',...
%          xcenter, QcurrentLW,'->m',...
%          xcenter, qexact,'-k');
%     %axis([xa xb, -1, 1]);
%     xlabel('x'); ylabel('q(x,t)');
%     title(sprintf('Time = %f',t(n)));
%     legend('Upwind','LxF','Lax-Wendroff','Exact','Location','SouthWest');
%     pause(0.5);

    plot(xcenter, QcurrentUP,'-db',...
         xcenter, QcurrentMM,'-vr',...
         xcenter, QcurrentBEE,'->m',...
         xcenter, qexact,'-k');
    %axis([xa xb, -1, 1]);
    xlabel('x'); ylabel('q(x,t)');
    title(sprintf('Time = %f',t(n)));
    legend('upwind','minmod','superbee','Exact','Location','SouthWest');
    pause(0.5);
    
    % update next time step
    QcurrentUP = QnextUP;
    QcurrentLxF = QnextLxF;
    QcurrentLW = QnextLW;
    QcurrentMM = QnextMM;
    QcurrentBEE = QnextBEE;
end

% plotting at final time
% plot(xcenter, QcurrentUP,'-db',...
%      xcenter, QcurrentLxF,'-vr',...
%      xcenter, QcurrentLW,'->g',...
%      xcenter, qexact,'-k');
% %axis([xa xb, -1, 1]);
% xlabel('x'); ylabel('q(x,t)');
% title(sprintf('Time = %f',t(n)));
% legend('Upwind','LxF','Lax-Wendroff','Exact','Location','SouthWest');

% plot(xcenter, QcurrentUP,'-db',...
%      xcenter, QcurrentMM,'-vr',...
%      xcenter, QcurrentBEE,'->g',...
%      xcenter, qexact,'-k');
% %axis([xa xb, -1, 1]);
% xlabel('x'); ylabel('q(x,t)');
% title(sprintf('Time = %f',t(n)));
% legend('upwind','minmod','superbee','Exact','Location','SouthWest');