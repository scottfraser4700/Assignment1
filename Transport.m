%constants
C.q_0 = 1.60217653e-19;
C.m_0 = 9.10938215e-31;
C.eps_0 = 8.854187817e-12;
C.kb = 1.3806504e-23;
C.T = 300;
C.am = 1.66053892e-27;
frameWidth = 200e-9;
frameHeight = 100e-9;
Vth = sqrt(C.kb * C.T / (0.26*C.m_0));
Tstop = 1e-10;
t = 0;
dt = 1e-12;
%setting up frame
figure
xlim([0 frameWidth])
ylim([0 frameHeight])
%setting up global variables, shame me later
global nAtoms;      %number of particles used
global X;           %x positions of all particles
global Y;           %y positions of all particles
global VX;          %X component of velocities
global VY;          %Y component of velocities
global Xnext;       %Next x coordinate of particles
global Ynext;       %Next y coordinate of particles
global Xp;          %previous x coord
global Yp;          %previous y coord


nAtoms = 100;
Xp = zeros(1,nAtoms);
Yp = zeros(1,nAtoms);
Xnext = zeros(1,nAtoms);
Ynext = zeros(1,nAtoms);
direction = 2*pi*rand(1, nAtoms);
VX = Vth .* cos(direction);
VY = Vth .* sin(direction);
X = frameWidth .* rand(1, nAtoms);
Y = frameHeight .* rand(1, nAtoms);
History

%electrons placed randomly, with uniform velocity and random directions

while t < Tstop
    %for each cycle in time, move electrons in accordance with their
    %velocities
    Xp = X;
    Yp = Y;
    Xnext = X + VX*dt;
    Ynext = Y + VY*dt;
    %X boundary conditions set, and x positions updated
    right = Xnext>frameWidth;
    left = Xnext<0;
    Xnext(right) = Xnext(right)-frameWidth;
    Xnext(left) = Xnext(left) + frameWidth;
    X = Xnext;
    Y = Ynext;
    %Y boundary conditions set
    top = Ynext > frameHeight;
    bottom = Ynext < 0;
    VY(top | bottom) = VY(top | bottom) * -1;
    
    
    t = t+dt;
    pause(0.1);
    plot([Xp(1:3) X(1:3)], [Yp(1:3) Y(1:3)], 'b')
    %scatter(X(1,:), Y(1,:))
    xlim([0 frameWidth])
    ylim([0 frameHeight])
end
%TODO add calculations for temperature, and plot temperature
%TODO reimplement plotting setup
