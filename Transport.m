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
nAtoms = 100;
figure
xlim([0 frameWidth])
ylim([0 frameHeight])

atomsX = frameWidth .* rand(1, nAtoms);
atomsY = frameHeight .* rand(1, nAtoms);
direction = 2*pi .* rand(1, nAtoms);
atomsVX = Vth .* cos(direction);
atomsVY = Vth .* sin(direction);
scatter(atomsX(1,:), atomsY(1,:))

%electrons placed randomly, with uniform velocity and random directions


%TODO - make function to calculate forces
%TODO - make function to initialize atom arrangement
%TODO - make time-loop to update velocities and positions
%TODO - add plotting to time-loop
%TODO - add boundary conditions