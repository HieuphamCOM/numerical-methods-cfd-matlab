%% FVM for 1D linear scalar conservation laws : smooth cases
%% q_t + u*q_x = 0, x \in (-1,1)
%% IC: q(x,0) = q0(x) = sin(pi*x)
%% BC: periodic boundary condition
%%

clear all
close all
clc

% problem setting
xa = -pi; xb = pi;    % computational domain
ubar = 1;
f = @(q) ubar * q;    % flux function 
q0 = @(x) sin( x );   % initial condition
dq0 = @(x) -cos( x ); % antiderivative of IC
qexact = @(x,t) q0( x - ubar * t );


% spartial grid data
nx = 50;
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

% temporal grid data
tfinal = 2;
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

% initial condition
for j = 1:nx
    QcurrentUP(j) = ( dq0(x(j).right) - dq0(x(j).left) ) / x(j).dx;
    QcurrentLxF(j) = ( dq0(x(j).right) - dq0(x(j).left) ) / x(j).dx;
    QcurrentLW(j) = ( dq0(x(j).right) - dq0(x(j).left) ) / x(j).dx;
end

% time looping
figure;
for n = 2:nt
    % periodic BC
    Qnp1UP = QcurrentUP(1);      % Q(nx+1) = Q(1)    
    Q0UP = QcurrentUP(end);      % Q(0) = Q(nx)
    Qnp1LxF = QcurrentLxF(1);    % Q(nx+1) = Q(1)
    Q0LxF = QcurrentLxF(end);    % Q(0) = Q(nx)
    Qnp1LW = QcurrentLW(1);      % Q(nx+1) = Q(1)  
    %Qnp2LW = QcurrentLW(2);      % Q(nx+2) = Q(2)  
    Q0LW = QcurrentLW(end);      % Q(0) = Q(nx)
    %Qm1LW = QcurrentLW(end-1);   % Q(-1) = Q(nx-1)    
    
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
    Fright = 0.5*( f(QcurrentLxF(1)) + f(QcurrentLxF(2)) )...
        - 0.5*x(1).dx/dt(n)*( QcurrentLxF(2) - QcurrentLxF(1) );
    Fleft = 0.5*( f(Q0LxF) + f(QcurrentLxF(1)) )...
        - 0.5*x(1).dx/dt(n)*( QcurrentLxF(1) - Q0LxF );
    QnextLxF(1) = QcurrentLxF(1) - dt(n)/x(1).dx * ( Fright - Fleft );
    for j = 2:nx-1
        % numerical flux functions
        Fright = 0.5*( f(QcurrentLxF(j)) + f(QcurrentLxF(j+1)) )...
            - 0.5*x(j).dx/dt(n)*( QcurrentLxF(j+1) - QcurrentLxF(j) );
        Fleft = 0.5*( f(QcurrentLxF(j-1)) + f(QcurrentLxF(j)) )...
            - 0.5*x(j).dx/dt(n)*( QcurrentLxF(j) - QcurrentLxF(j-1) );
        
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
    
    
    % plotting current solution
    plot(xcenter, QcurrentUP,'-db',...
         xcenter, QcurrentLxF,'-vr',...
         xcenter, QcurrentLW,'->g',...
         xcenter, qexact(xcenter,t(n)),'-k');
    axis([xa xb, -1, 1]);
    xlabel('x'); ylabel('q(x,t)');
    title(sprintf('Time = %f',t(n)));
    legend('Upwind','LxF','Lax-Wendroff','Exact');
    pause(0.5);
    
    % update next time step
    QcurrentUP = QnextUP;
    QcurrentLxF = QnextLxF;
    QcurrentLW = QnextLW;
end

% plotting at final time
plot(xcenter, QcurrentUP,'-db',...
     xcenter, QcurrentLxF,'-vr',...
     xcenter, QcurrentLW,'->g',...
     xcenter, qexact(xcenter,t(n)),'-k');
axis([xa xb, -1, 1]);
xlabel('x'); ylabel('q(x,t)');
title(sprintf('Time = %f',t(n)));
legend('Upwind','LxF','Lax-Wendroff','Exact');


























































