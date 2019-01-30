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
dt = 1e-14;
%setting up global variables, shame me later
global nAtoms;      %number of particles used
global X;           %x positions of all particles
global Y;           %y positions of all particles
global VX;          %X component of velocities
global VY;          %Y component of velocities
global Xnext;       %Next x coordinate of particles
global Ynext;       %Next y coordinate of particles
global EK;
global Temperature;
%initialization of vectors
nAtoms = 100;
Xnext = zeros(1,nAtoms);
Ynext = zeros(1,nAtoms);
direction = 2*pi*rand(1, nAtoms);
VX = Vth .* cos(direction);
VY = Vth .* sin(direction);
X = frameWidth .* rand(1, nAtoms);
Y = frameHeight .* rand(1, nAtoms);
EK = 0;
Temperature = 0;
%electrons placed randomly, with uniform velocity and random directions

while t < Tstop
    %for each cycle in time, move electrons in accordance with their
    %velocities
   
    Xnext = X + VX*dt;
    Ynext = Y + VY*dt;
    %X boundary conditions set
    right = Xnext>frameWidth;
    left = Xnext<0;
    Xnext(right) = Xnext(right)-frameWidth;
    Xnext(left) = Xnext(left) + frameWidth;
    %Y boundary conditions set
    top = Ynext > frameHeight;
    bottom = Ynext < 0;
    VY(top | bottom) = VY(top | bottom) * -1;
    %EK = mean(VX.*VX+VY.*VY);
    %Temperature = EK / nAtoms / C.kb;
    %plotting, but trying to avoid plotting the full horizontal jump
    
   %subplot(2,1,2); 
   %figure(2)
  %plot(t, Temperature);
   %hold on
   %xlim([0 1e-10]) 
  % ylim([0 1e32])
    
    
    if abs((Xnext(1) - X(1)) < 2*abs(VX(1))*dt)
        %subplot(2,1,1);
        figure(1)
        plot([Xnext(1) X(1)], [Ynext(1) Y(1)], 'blue')
        xlim([0 frameWidth])
        ylim([0 frameHeight])
        hold on
    end
    if abs((Xnext(2) - X(2)) < 2*abs(VX(2))*dt)
        %subplot(2,1,1);
        figure(1)
        plot([Xnext(2) X(2)], [Ynext(2) Y(2)], 'red')
        xlim([0 frameWidth])
        ylim([0 frameHeight])
        hold on
    end
    if abs((Xnext(3) - X(3)) < 2*abs(VX(3))*dt)
        %subplot(2,1,1);
        figure(1)
        plot([Xnext(3) X(3)], [Ynext(3) Y(3)], 'green')
        xlim([0 frameWidth])
        ylim([0 frameHeight])
        hold on
    end
    
    
    %updating positions, and advancing time a step forward so the while
    %loop works
    X = Xnext;
    Y = Ynext;
    t = t+dt;
    
 
    pause(0.01);
    
    
end
%TODO add calculations for temperature, and plot temperature
