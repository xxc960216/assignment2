%ELEC4700 Assignment 2 Part 2
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
    
% Part a

%working area same as part 1
W = 20;
L = 30;

halfX = L/2;
halfY = W/2;

G = zeros(L*W,L*W);
B = zeros(L*W,1);

%conductivity
s1 = 1;
s2 = 0.01;

%resistive regions size(randomly assigned)
rL = L*1/3;
rW = W*4/7;

map = @(i,j) j + (i - 1)*W;

%map containing conductivity
Smap = zeros(L,W);
for i =1:1:L
    for j =1:1:W
        n = map(i,j);
        nxm = map(i-1, j);
        nxp = map(i+1, j);
        nyp = map(i, j+1);
        nym = map(i, j-1);
        
        %left column
        if(i==1)
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = 1;
            Smap(i,j) = s1;
        %right column
        elseif(i==L)
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = 0;
            Smap(i,j) = s1;
        %bottom row
        elseif(j==1)
            G(n,:) = 0;
            %bc- within the bottle-neck
            if(i > halfX-(rL/2) && i < halfX+(rL/2))
                G(n,nxm) = s2;
                G(n,nxp) = s2;
                G(n,nyp) = s2;
                G(n,n) = -3*s2;
                Smap(i,j) = s2;
            else
                G(n,nxm) = s1;
                G(n,nxp) = s1;
                G(n,nyp) = s1;
                G(n,n) = -3;
                Smap(i,j) = s1;
            end
        %top row
        elseif(j==W)
            G(n,:) = 0;
            % boundry condition-- within the bottle-neck
            if(i > halfX-(rL/2) && i < halfX+(rL/2))
                G(n,nxm) = s2;
                G(n,nxp) = s2;
                G(n,nym) = s2;
                G(n,n) = -3*s2;
                Smap(i,j) = s2;
            else
                G(n,nxm) = s1;
                G(n,nxp) = s1;
                G(n,nym) = s1;
                G(n,n) = -3;
                Smap(i,j) = s1;
            end
        % in the middle    
        else
            G(n,:) = 0;
            G(n,n) = -4;
            %  boundry condition-- within the bottle-neck for both X and Y
            if((i > halfX-(rL/2) && i < halfX+(rL/2)) && (j > halfY+(rW/2) || j < halfY-(rW/2)))
           
                G(n,nxp) = s2;
                G(n,nxm) = s2;
                G(n,nyp) = s2;
                G(n,nym) = s2;
                Smap(i,j) = s2;
            % not within bottle-neck region
            else
                G(n,nxp) = s1;
                G(n,nxm) = s1;
                G(n,nyp) = s1;
                G(n,nym) = s1;
                Smap(i,j) = s1;
            end     
        end
    end
end

V = G\B;

%map
Vmap = zeros(L,W);
for i =1:1:L
    for j =1:1:W
        n=map(i,j);
        Vmap(i,j) =V(n);
    end
end

[Ey,Ex] = gradient(Vmap);

E = gradient(Vmap);

J = Smap.*E;

figure(3)
surf(Smap)
colorbar
title('Sigma map'),xlabel('X'),ylabel('Y'),zlabel('Sigma');

figure(4)
surf(Vmap)
colorbar
title('Voltage map'),xlabel('X'),ylabel('Y'),zlabel('Voltage')

figure(5)
surf(Ex)
colorbar
title('Ex'),xlabel('X'),ylabel('Y'),zlabel('E Field');

figure(6)
surf(Ey)
colorbar
title('Ey'),xlabel('X'),ylabel('Y'),zlabel('E Field');

figure(7)
surf(J)
colorbar
title('Current Density'),xlabel('X'),ylabel('Y'),zlabel('Current/m^2');


%Part 2-b
loopb = 8;

currentb = zeros(1,loopb);
for i =1:1:loopb
    currentb(i) = part2b(i);
end

figure(8)
plot(1:1:loopb,currentb);
title('Mesh Density')
xlabel('Number of Points Per Unit Spacing')
ylabel('Current')
ylim([-0.5e-6 0])


%part 2-c
loopc = 10;

currentc = zeros(1,loopc);
for i =2:1:loopc
    currentc(i) = part2c(1/i);
end

figure(9)
plot(1:1:loopc,currentc);
title('Bottleneck Narrowing')
xlabel('Spacing of L*1/x')
ylabel('Current')
xlim([2 loopc])

%part 2-d
loopd = 100;

currentd = zeros(1,loopd);
for i =1:1:loopd
    sig = i*1e-2;
    currentd(i) = part2d(sig);
end

figure(10)
plot(1e-2:1e-2:1,currentd);
title('Varying Sigma')
xlabel('Sigma')
ylabel('Current')

function I = part2b(mesh_density)

%working area
W = 20;
L = 30;

G = sparse(L*W,L*W);
B = zeros(L*W,1);

%conductivity
s1 = 1;
s2 = 0.01;

%for mesh density
dr = 1/mesh_density;

pointsNumberL = mesh_density*L;
pointsNumberW = mesh_density*W;

%resistive regions size
rL = L*1/3;
rW = W*4/7;

halfX = L/2;
halfY = W/2;

map = @(i,j) j + (i - 1)*pointsNumberW;

Smap = zeros(L,W);
for i =1:1:pointsL
    for j =1:1:pointsNumberW
        n = map(i,j);
        nxm = map(i-1, j);
        nxp = map(i+1, j);
        nyp = map(i, j+1);
        nym = map(i, j-1);
        
        if(i==1)
            G(n,:) = 0;
            G(n,n) = 1/(dr^2);
            B(n) = 1;
            Smap(i,j) = s1;
        elseif(i==pointsNumberL)
            G(n,:) = 0;
            G(n,n) = 1/(dr^2);
            B(n) = 1;
            Smap(i,j) = s1;
        elseif(j==1)
            G(n,:) = 0;
            G(n,n) = -3/(dr^2);
            if(i/mesh_density > halfX-(rL/2) && i/mesh_density < halfX+(rL/2))
                G(n,nxm) = s2/(dr^2);
                G(n,nxp) = s2/(dr^2);
                G(n,nyp) = s2/(dr^2);
                Smap(i,j) = s2;
            else
                G(n,nxm) = s1/(dr^2);
                G(n,nxp) = s1/(dr^2);
                G(n,nyp) = s1/(dr^2);
                Smap(i,j) = s1;
            end
        elseif(j==pointsNumberW)
            G(n,:) = 0;
            G(n,n) = -3/(dr^2);
            if(i/mesh_density > halfX-(rL/2) && i/mesh_density < halfX+(rL/2))
                G(n,nxm) = s2/(dr^2);
                G(n,nxp) = s2/(dr^2);
                G(n,nym) = s2/(dr^2);
                Smap(i,j) = s2;
            else
                G(n,nxm) = s1/(dr^2);
                G(n,nxp) = s1/(dr^2);
                G(n,nym) = s1/(dr^2);
                Smap(i,j) = s1;
            end
        else          
            G(n,:) = 0;
            G(n,n) = -4/(dr^2);
        
            if((i/mesh_density > halfX-(rL/2) && i/mesh_density < halfX+(rL/2)) && (j/mesh_density > halfY+(rW/2) || j/mesh_density < halfY-(rW/2)))
          
                G(n,nxp) = s2/(dr^2);
                G(n,nxm) = s2/(dr^2);
                G(n,nyp) = s2/(dr^2);
                G(n,nym) = s2/(dr^2);
                Smap(i,j) = s2;
            else
                G(n,nxp) = s1/(dr^2);
                G(n,nxm) = s1/(dr^2);
                G(n,nyp) = s1/(dr^2);
                G(n,nym) = s1/(dr^2);
                Smap(i,j) = s1;
            end     
        end
    end
end

V = G\B;

%map
Vmap = zeros(pointsL,pointsW);
for i =1:1:pointsL
    for j =1:1:pointsW
        n=j+(i-1)*pointsW;
        Vmap(i,j) =V(n);
    end
end

E = gradient(Vmap);

%a matrix of current density
J = Smap.*E;

%size of the area
area = L*W;
%average current density
Javg = sum(sum(J))/(pointsL*pointsW);
I = Javg/area;

end

function I = part2c(bottle_neck_narrow_ratio)

%size of working area
W = 20;
L = 30;

G = zeros(L*W,L*W);
B = zeros(L*W,1);

%conductivity
s1 = 1;
s2 = 0.01;

%centre point of the working area
halfX = L/2;
halfY = W/2;

%resistive regions size
rL = L*1/3;
rW = W*bottle_neck_narrow_ratio;
map = @(i,j) j + (i - 1)*W;
Smap = zeros(L,W);
for i =1:1:L
    for j =1:1:W
        n = map(i,j);
        nxm = map(i-1, j);
        nxp = map(i+1, j);
        nyp = map(i, j+1);
        nym = map(i, j-1);
        
        %left column
        if(i==1)
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = 1;
            Smap(i,j) = s1;
        %right column
        elseif(i==L)
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = 0;
            Smap(i,j) = s1;
        %bottom row
        elseif(j==1)
            G(n,:) = 0;
            G(n,n) = -3;
            if(i > halfX-(rL/2) && i < halfX+(rL/2))
                G(n,nxm) = s2;
                G(n,nxp) = s2;
                G(n,nyp) = s2;
                Smap(i,j) = s2;
            else
                G(n,nxm) = s1;
                G(n,nxp) = s1;
                G(n,nyp) = s1;
                Smap(i,j) = s1;
            end
        %top row 
        elseif(j==W)
            G(n,:) = 0;
            G(n,n) = -3;
            if(i > halfX-(rL/2) && i < halfX+(rL/2))
                G(n,nxm) = s2;
                G(n,nxp) = s2;
                G(n,nym) = s2;
                Smap(i,j) = s2;
            else
                G(n,nxm) = s1;
                G(n,nxp) = s1;
                G(n,nym) = s1;
                Smap(i,j) = s1;
            end
        else          
            G(n,:) = 0;
            G(n,n) = -4;
           
            if((i > halfX-(rL/2) && i < halfX+(rL/2)) && (j > halfY+(rW/2) || j < halfY-(rW/2)))         
                G(n,nxp) = s2;
                G(n,nxm) = s2;
                G(n,nyp) = s2;
                G(n,nym) = s2;
                Smap(i,j) = s2;
            else
                G(n,nxp) = s1;
                G(n,nxm) = s1;
                G(n,nyp) = s1;
                G(n,nym) = s1;
                Smap(i,j) = s1;
            end     
        end
    end
end

V = G\B;

%map
Vmap = zeros(L,W);
for i =1:1:L
    for j =1:1:W
        n=j+(i-1)*W;
        Vmap(i,j) =V(n);
    end
end

E = gradient(Vmap);

J = Smap.*E;

area = L*W;
Javg = sum(sum(J))/(L*W);
I = Javg/area;

end
function I = part2d(varying_sigma)

%size of working area
W = 20;
L = 30;

halfX = L/2;
halfY = W/2;

G = zeros(L*W,L*W);
B = zeros(L*W,1);

%conductivity
s1 = 1;
s2 = varying_sigma;

%resistive regions size
rL = L*1/3;
rW = W*4/7;

map = @(i,j) j + (i - 1)*W;

Smap = zeros(L,W);
for i =1:1:L
    for j =1:1:W
        n = map(i,j);
        nxm = map(i-1, j);
        nxp = map(i+1, j);
        nyp = map(i, j+1);
        nym = map(i, j-1);
        %left column
        if(i==1)
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = 1;
            Smap(i,j) = s1;
        %right column
        elseif(i==L)
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = 0;
            Smap(i,j) = s1;
        %bottom row
        elseif(j==1)
            G(n,:) = 0;
            G(n,n) = -3;
            if(i > halfX-(rL/2) && i < halfX+(rL/2))
                G(n,nxm) = s2;
                G(n,nxp) = s2;
                G(n,nyp) = s2;
                Smap(i,j) = s2;
            else
                G(n,nxm) = s1;
                G(n,nxp) = s1;
                G(n,nyp) = s1;
                Smap(i,j) = s1;
            end
        %top row    
        elseif(j==W)
            G(n,:) = 0;
            G(n,n) = -3;
            if(i > halfX-(rL/2) && i < halfX+(rL/2))
                G(n,nxm) = s2;
                G(n,nxp) = s2;
                G(n,nym) = s2;
                Smap(i,j) = s2;
            else
                G(n,nxm) = s1;
                G(n,nxp) = s1;
                G(n,nym) = s1;
                Smap(i,j) = s1;
            end
        else          
            G(n,:) = 0;
            G(n,n) = -4;
            if((i > halfX-(rL/2) && i < halfX+(rL/2)) && (j > halfY+(rW/2) || j < halfY-(rW/2)))
                G(n,nxp) = s2;
                G(n,nxm) = s2;
                G(n,nyp) = s2;
                G(n,nym) = s2;
                Smap(i,j) = s2;
            else
                G(n,nxp) = s1;
                G(n,nxm) = s1;
                G(n,nyp) = s1;
                G(n,nym) = s1;
                Smap(i,j) = s1;
            end     
        end
    end
end

V = G\B;

%map
Vmap = zeros(L,W);
for i =1:1:L
    for j =1:1:W
        n=j+(i-1)*W;
        Vmap(i,j) =V(n);
    end
end


E = gradient(Vmap);

J = Smap.*E;

area = L*W;
Javg = sum(sum(J))/(L*W);
I = Javg/area;

end
