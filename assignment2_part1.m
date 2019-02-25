%ELEC4700 Assignment 2 Part 1
%Xiaochen Xin 100989338
global C
    C.q_0 = 1.60217653e-19;             % electron charge
    C.hb = 1.054571596e-34;             % Dirac constant
    C.h = C.hb * 2 * pi;                % Planck constant
    C.m_0 = 9.10938215e-31;             % electron mass
    C.kb = 1.3806504e-23;               % Boltzmann constant
    C.eps_0 = 8.854187817e-12;          % vacuum permittivity
    C.mu_0 = 1.2566370614e-6;           % vacuum permeability
    C.c = 299792458;                    % speed of light
    C.g = 9.80665;                      % metres (32.1740 ft) per s²

nx = 75;
ny = 50;

dx = 1;
dy = 1;
G = sparse(nx*ny, ny*nx);
alpha = (C.hb^2) / (2 * C.m_0);

map = @(i,j) j + (i - 1)*ny;

%part a
% Setting up boundary conditions
% Go through Xs
for i=1:nx
    % Go through Ys
    for j=1:ny
        n = map(i,j);
        nxm = map(i-1,j);
        nxp = map(i+1,j);
        nym = map(i,j-1);
        nyp = map(i,j+1);
        
        % left column
        if i == 1            
            V(n) = 1;
            G(n,:) = 0;
            G(n,n) = 1;
        % right column
        elseif i == nx
            V(n) = 0;
            G(n,:) = 0;
            G(n,n) = 1;            
        % bottom row
        elseif j == 1
            G(n,:) = 0;
            G(n,n) = -3;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nyp) = 1;
        % top row
        elseif j == ny
            G(n,:) = 0;
            G(n,n) = -3;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nyp) = 1;
        else
            G(n,:) = 0;
            G(n,n) = -4;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nym) = 1;
            G(n,nyp) = 1;
        end
    end
end

F = G\V;

% Surf Plot
surfa = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = map(i,j);
        surfa(i,j) = F(n);
    end
end

%Plot
figure(1)
surf(surfa)
title('Electrostatic Potential in Rectangle')
xlabel('Width')
ylabel('Length')
zlabel('Voltage')


%Part b --------------------------------------------------------

% Start setting up boundary conditions
% Go through Xs
for i=1:nx
    % Go through Ys
    for j=1:ny
        n = map(i,j);
        nxm = map(i-1,j);
        nxp = map(i+1,j);
        nym = map(i,j-1);
        nyp = map(i,j+1);
        
        % left column
        if i == 1
            V(n) = 1;
            G(n,:) = 0;
            G(n,n) = 1;            
        % right column
        elseif i == nx
            V(n) = 1;
            G(n,:) = 0;
            G(n,n) = 1;            
        % bottom row
        elseif j == 1
            % More BCs related to V are defined for part 2
            V(n) = 0;
            G(n,:) = 0;
            G(n,n) = 1;            
        % top row
        elseif j == ny
            V(n) = 0;
            G(n,:) = 0;
            G(n,n) = 1;           
        else
            G(n,:) = 0;
            G(n,n) = -4;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nym) = 1;
            G(n,nyp) = 1;
        end
    end
end

F = G\V;

% Surf Plot
surfb = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = map(i,j);
        surfb(i,j) = F(n);
    end
end

%Plot
figure(2)
surf(surfb)
title('Electrostatic Potential in Rectangle')
xlabel('Width')
ylabel('Length')
zlabel('Voltage')
