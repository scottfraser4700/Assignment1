function = GetForces(PhiCutoff,Epsilon, Sigma)
%GETFORCES assigns values to the global vectors FX and FY
%

for i = 1:nAtoms
   FX(i) = 0;
   FY(i) = 0;
   phi(i) = 0;
   
   for j = 1:nAtoms
       %skip for cases where it checks itself
       if i == j, continue, end
       %checks whether the particles are neighbours
       dx = x(i) - x(j);
       dy = y(i) - y(j);
       r = sqrt(dx^2+dy^2);
       if r > PhiCutoff, continue, end
       
       [aphi, dPhidr] = LJ(r, Epsilon, sigma);
       
       angle = atan2(dy, dx);
       dFX = dPhidr*cos(angle);
       dFY = dPhidr*sin(angle);
       phi(i) = phi(i) + aphi;
       FX(i) = FX(i) + dFx;
       FY(i) = FY(i) + dFy;
   end
end



end

