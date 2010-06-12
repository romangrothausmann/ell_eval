clear all

t= load octave_test02.txt;

N= size(t, 1)
m= zeros(N,12)
#e= zeros(N,3)
num= 101
phi0=pi(1)/4 
lambda0=pi(1)/4

for n=1:1:N

  i= [t(n,1),t(n,4),t(n,5)
      t(n,4),t(n,2),t(n,6)
      t(n,5),t(n,6),t(n,3)]

  [v, l] = eig (i)
  [theta, phi, r]= cart2sph(l(1,1), l(2,2), l(3,3))
  e1(n,:)= [theta, phi, r]
  #[x,y,z]= sph2cart (theta, phi, 1)
  #e2(n,:)=  [x,y,z]
  #p1(n,:)= [x / (1 - z), y / (1 - z)]

  r= horzcat (reshape (v, 1, 9), l(1,1), l(2,2), l(3,3))
  m(n,:)=r
  #fprintf(result.txt, "

end

a=[linspace(0,pi(1)/2,num); zeros(1, num)]
b=[linspace(0,pi(1)/2,num); ones(1, num) * pi(1)/2]
c=[zeros(1, num); linspace(0,pi(1)/2,num)]
d= horzcat(a, b, c)
d= [d(1,:); d(2,:)]
d= horzcat(d, e1(:,1:2)')
#d= e1(:,1:2)'

for n=1:1:size(d,2)
  k= 2 / (1 + sin(phi0) * sin(d(1,n)) + cos(phi0) * cos(d(1,n)) * cos(d(2,n) - lambda0))
  x= k * cos(d(1,n)) * sin(d(2,n) - lambda0)
  y= k *     (cos(phi0) * sin(d(1,n)) - sin(phi0) * cos(d(1,n)) * cos(d(2,n) - lambda0))
  p2(n,:)= [x, y]
  [x,y,z]= sph2cart (d(2,n), d(1,n), 1)
  e3(n,:)=  [x,y,z]
end

#e4= vertcat(e2,e3)
#e4=e3

#N= size(e4, 1)
#for n=1:1:N
#  p2(n,:)= [e4(n,1) / (1 - e4(n,3) ),  e4(n,2) / (1 -  e4(n,3))]
#  k= 2 / ( 1 + sin(phi1)sin(phi
#  p2(n,:)= [
#end

#scatter3 (e3(:,1),e3(:,2),e3(:,3), [], e3(:,1))
scatter (p2(:,1), p2(:,2), [], 1)# p2(:,1))

#save -ascii m.txt m





