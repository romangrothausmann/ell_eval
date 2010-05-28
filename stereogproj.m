function r = stereogproj(theta, phi, radius, phi0, lambda0)

for n=1:1:size(theta,2)

  #stereographic projection
  #http://mathworld.wolfram.com/StereographicProjection.html
  k= 2 * radius / (1 + sin(phi0) * sin(phi(n)) + cos(phi0) * cos(phi(n)) * cos(theta(n) - lambda0));
  x= k * cos(phi(n)) * sin(theta(n) - lambda0);
  y= k *              (cos(phi0) * sin(phi(n)) - sin(phi0) * cos(phi(n)) * cos(theta(n) - lambda0));
  r(:,n)= [x, y];#stereographic projected point list
end;
