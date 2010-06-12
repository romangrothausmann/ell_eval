#!/net/home/ftd/localusr/bin/octave -qf

#skript to visualize ellipsoid axis ratios
#this is useful to identify amount of "smarties" and "cigars"
######
#REALLY MAKE SURE TO NOT MIX UP THETA AND PHI!!!!!!!


clear all;

#awk 'NR>1 {print $6, $7, $8, $9, $10, $11, $3, $4, $5}' pudel/az239b_080725a_bak02_3D-ana.csv > octave/az239b_080725a_bak02_3D-ana_I.txt
#t=rand(200,6);

lbb= 0;

arg_list = argv ();
if nargin
  if (nargin > 1)
    lbb= 1;
    printf("Doing BBox fitting with %s...\n", arg_list{2});
    bb= load(arg_list{2}); #ellipsoid_bbox.txt;
  else
    printf("Doing I fitting with %s...\n", arg_list{1});
    t= load(arg_list{1}); #octave_test02.txt;
  endif
else
printf("Usage: %s <I-fit-file> [BBox-fit-file]\n", program_name);
exit(1)
endif



#awk 'NR>1 {print $12, $13, $14}' az239b_080725a_bak02_3D-ana.csv > octave/az239b_080725a_bak02_3D-ana_BBox.txt   

#bb= sort(bb')'; #don't sort here!!! it's done below

N= size(t, 1);
#m= zeros(N,12);
#e= zeros(N,3);
num= 101;
radius= 1;
phi0=pi(1)/4;
lambda0=pi(1)/4;
mscale= 20/9;
mass= 1;

if (lbb)
  [fid, msg] = fopen ("BBox-fit.txt", "w");
else
  [fid, msg] = fopen ("I-fit.txt", "w");
endif

for n=1:1:N;

  I= [t(n,1),t(n,4),t(n,5);
      t(n,4),t(n,2),t(n,6);
      t(n,5),t(n,6),t(n,3)];
  #I= [t(n,2)+t(n,3),-t(n,4),-t(n,5);  #avizo calculates the math. 2order moments not the pyhsical!
  #    -t(n,4),t(n,1)+t(n,3),-t(n,6); #they are also normalized by their volumes
  #    -t(n,5),-t(n,6),t(n,1)+t(n,2)] * t(n,10);
    #i=rand(3,3);
  [v, l]= eig(I);
    #l=rand(3,3);

  #test if eigenvalues are correlated to their eigenvectors
  if (I*v != [l(1,1)*v(:,1)'; l(2,2)*v(:,2)'; l(3,3)*v(:,3)']')
    n
    I*v
    [l(1,1)*v(:,1)'; l(2,2)*v(:,2)'; l(3,3)*v(:,3)']'
    fflush(stdout);
  end

  if (lbb)
    ax= [bb(n,1), bb(n,2), bb(n,3)]; #axis fitted to BBox
  else
    #ax= [l(1,1), l(2,2), l(3,3)]; #half axis fitted to I
    #ll= 2*sqrt(inv(l));
    #ax= mscale * 2*sum(sqrt(l));
    ax= sum(sqrt(l));
    #ax= ceil(sum(l)*1000)/1000;
    #ax= abs(sqrt(5/2/mass*sum(l*(ones(3,3)-2*eye(3))))); #http://en.wikipedia.org/wiki/Ellipsoid#Mass_properties
    #ax= [ll(1,1), ll(2,2), ll(3,3)];#axis fitted to I
  endif
  
  vol= 4/3*pi*ax(1)*ax(2)*ax(3);
  if(vol>0)
    ax= 2*(t(n,10)/vol)^(1/3)*ax;
  else
    n
    t(n,:)
    I
    l
    ax
    vol
    printf("Volum <= 0! Aborting\n")
    break
  endif

  [axs, axi]= sort (ax); #making a < b < c if axis-names are not specially assigned!
  #eu= euler_angles(v(:,axi(1)), v(:,axi(2)), v(:,axi(3))); #index ordered v

  [theta, phi, r]= cart2sph(axs(1), axs(2), axs(3)); 
  e1(n,:)= [theta, phi, r];

 # if (n == 37)
 #   i
 #   l
 #   ll
 #   e1(n,:)
 #   fflush(stdout);
 # end

 #fprintf(fid, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", ax(1), ax(2), ax(3), eu(1), eu(2), eu(3), t(n,7), t(n,8), t(n,9));
  #fprintf(fid, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", ax, t(n,7:9), v);
  fprintf(fid, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", axs, t(n,7:9), v(:,axi(1)), v(:,axi(2)), v(:,axi(3)));
end;

fclose(fid);

w =acos(dot([1,1,1], [1,1,0])/(norm([1,1,1]) * norm([1,1,0])));
#a0=[linspace(0,pi/2,num); zeros(1, num)];# 1/4 circle in ab-plane
#a0=[pi/4, pi/2]';#[1,0,0] point?
#b0=[ones(1, num) * pi/2; linspace(0,pi/2,num)];#arc from [0, 1, 0] to a0
b0=[ones(1, num) * pi/2; linspace(pi/4,pi/2,num)];#arc from b0 to a0 
#b0=[pi/2, pi/4]';#[0, 1/sqrt(2), 1/sqrt(2)] point
#c0=[zeros(1, num); linspace(0,pi/2,num)];#arc from [1, 0, 0] to a0
c0=[pi/4, w]';#[1/sqrt(3), 1/sqrt(3), 1/sqrt(3)] point
#a1=[linspace(0,pi/2,num); ones(1, num) * w];#arc from [1/sqrt(3), 0, sqrt(2/3)] to  [0, 1/sqrt(3), sqrt(2/3)]
#c1=[ones(1, num) * pi/4; linspace(0,pi/2,num)]; #arc from [1/sqrt(2), 1/sqrt(2), 0] to a0; only good for 2D view!
#c1=
#hl= horzcat(a0, b0, c0, a1, c1);
hl= horzcat(b0, c0);
#d= [d(1,:); d(2,:)];
d= horzcat(hl, e1(:,1:2)');
#d= e1(:,1:2)';

for n=1:1:size(d,2);
  #stereographic projection
  #http://mathworld.wolfram.com/StereographicProjection.html
  k= 2 * radius / (1 + sin(phi0) * sin(d(2,n)) + cos(phi0) * cos(d(2,n)) * cos(d(1,n) - lambda0));
  x= k * cos(d(2,n)) * sin(d(1,n) - lambda0);
  y= k *              (cos(phi0) * sin(d(2,n)) - sin(phi0) * cos(d(2,n)) * cos(d(1,n) - lambda0));
  p2d(n,:)= [x, y];#stereographic projected point list

  [x1,y1,z1]= sph2cart (d(1,n), d(2,n), 1);#projection of 3D points onto unit sphere
  p3d(n,:)= [x1,y1,z1];
end;

gp3d=  p3d(1:size(hl,2),:);  #guide points in 3D
gp2d=  p2d(1:size(hl,2),:);  #guide points in 2D

dp3d=  p3d(size(hl,2)+1:size(d,2),:);  #data points in 3D
dp2d=  p2d(size(hl,2)+1:size(d,2),:);  #data points in 2D

[Theta, R]= cart2pol(p2d(:,1), p2d(:,2));
ra= atan((gp2d(1,2) - gp2d(size(gp2d,1),2)) / (gp2d(1,1) - gp2d(size(gp2d,1),1))); #ra= 31.325°; why not 30°, stereographic projection error?
[p2dr(:,1), p2dr(:,2)]= pol2cart(Theta - ra, R); #rotation for second 2D-hist
#[p2dr(:,1), p2dr(:,2)]= pol2cart(Theta - 30 / 180 * pi, R); #rotation of 30° for second 2D-hist

gp2dr= p2dr(1:size(hl,2),:); #guide points in 2D rotated
dp2dr= p2dr(size(hl,2)+1:size(d,2),:); #data points in 2D rotated

#scatter3 (e1(:,1),e1(:,2),e1(:,3), [], 1);
scatter3 (dp3d(:,1),dp3d(:,2),dp3d(:,3), [], 1); #"markersize", 3, 1);
w3= vertcat(gp3d, gp3d(1,:)); #create guide line matrix
hold on
plot3 (w3(:,1), w3(:,2), w3(:,3))        #plot guide line matrix
hold off
#axis ("square");
axis ([0,1,0,1,0,1],"square");

azimuth= 135;
#azimuth= 315;
elevation= acosd(dot([1,1,1], [1,1,0])/(norm([1,1,1]) * norm([1,1,0])));
#elevation= elevation + 90;
view(azimuth, elevation);
#break
####printing now...

#paper_size = [640, 480];
#set (gcf, "paperunits", "inches")
#set (gcf, "papertype", "<custom>")
#set (gcf, "papersize", paper_size)
#set (gcf, "paperposition", [0, 0, paper_size])

#figure('Position', [0, 0, 600, 400]); 

print('ellipsoid_oct02_01.png', '-dpng');#, '-r100');
print('ellipsoid_oct02_01.svg', '-dsvg');

####printing end

#scatter (p2d(:,1), p2d(:,2), [], 1)# p2d(:,1));
scatter (dp2d(:,1), dp2d(:,2), [], 1)
w2= vertcat(gp2d, gp2d(1,:)); #create guide line matrix
hold on
plot (w2(:,1), w2(:,2))        #plot guide line matrix
hold off
#axis ("square");
axis ([gp2d(size(gp2d,1),1),gp2d(1,1),gp2d(size(gp2d,1)-1,2),gp2d(size(gp2d,1),2)],"square");
#break

####printing again...

print('ellipsoid_oct02_02.png', '-dpng');#, '-r100');
print('ellipsoid_oct02_02.svg', '-dsvg');

####printing end

#save -ascii m.txt m;

bin= 50;
dgp2d= vertcat(dp2d, gp2d(1,:), gp2d(size(gp2d,1)-1,:), gp2d(size(gp2d,1),:));

xbin=bin;
ybin=round( (gp2d(size(gp2d,1)-1,2) - gp2d(size(gp2d,1),2)) / (gp2d(1,1) - gp2d(size(gp2d,1),1)) * xbin);

vXEdge = linspace(gp2d(size(gp2d,1),1)-0/xbin,gp2d(1,1)+1/xbin,xbin); #dgp2d(3,:)<>c0; dgp2d(2,:)<>b0; dgp2d(1,:)<>a0
vYEdge = linspace(gp2d(size(gp2d,1),2)-1/ybin,gp2d(size(gp2d,1)-1,2)+2/ybin,ybin);
#mHist2d = hist2d([dgp2d(:,2),dgp2d(:,1)],vYEdge,vXEdge); #2D-hist with guide points
mHist2d = hist2d([dp2d(:,2),dp2d(:,1)],vYEdge,vXEdge); #2D-hist without guide points


nXBins = length(vXEdge);
nYBins = length(vYEdge);
vXLabel = 0.5*(vXEdge(1:(nXBins-1))+vXEdge(2:nXBins));
vYLabel = 0.5*(vYEdge(1:(nYBins-1))+vYEdge(2:nYBins));
pcolor(vXLabel, vYLabel, mHist2d); 
#w2= vertcat(gp2d, gp2d(1,:)); #create guide line matrix
hold on
plot (w2(:,1), w2(:,2))        #plot guide line matrix
hold off
shading flat; #means no border around each hist rectangle

####create a colour map for many small values

N= 256;
j= jet(N+1);
for n=1:1:N
  i=round((n/N)^(1/2) * N + 1);
  cmap(n,:)= j(i,:);
end

cmap=vertcat([1,1,1],cmap);

####end colour map cration

colormap(cmap)
#colormap(hsv(128))
#caxis([0, 10])#ignore extreme etremes
colorbar #show colorbar
#axis ([0,gp2d(size(1,1)),gp2d(size(gp2d,1)-1,2),gp2d(size(gp2d,1),2)],"square");
axis ("square");

####printing again...

print('ellipsoid_oct02_03.png', '-dpng');#, '-r100');
print('ellipsoid_oct02_03.svg', '-dsvg');

####printing end

#break

clear bin xbin ybin vXEdge vYEdge mHist2d nXBins nYBins vXLabel vYLabel

bin= 50;
dgp2dr= vertcat(dp2dr, gp2dr(1,:), gp2dr(size(gp2d,1)-1,:), gp2dr(size(gp2d,1),:));

xbin=bin;
ybin=round( (gp2dr(size(gp2dr,1)-1,2) - gp2dr(size(gp2dr,1),2)) / (gp2dr(1,1) - gp2dr(size(gp2dr,1),1)) * xbin);

vXEdge = linspace(gp2dr(size(gp2dr,1),1)-0/xbin,gp2dr(1,1)+1/xbin,xbin); #dgp2dr(3,:)<>c0; dgp2dr(2,:)<>b0; dgp2dr(1,:)<>a0
#vYEdge = linspace(gp2dr(size(gp2dr,1),2),gp2dr(size(gp2dr,1)-1,2),ybin);#here extreme points are missing
vYEdge = linspace(gp2dr(size(gp2dr,1),2)-1/ybin,gp2dr(size(gp2dr,1)-1,2)+2/ybin,ybin); #here not
#mHist2d = hist2d([dgp2dr(:,2),dgp2dr(:,1)],vYEdge,vXEdge);#2D-hist with guide points
mHist2d = hist2d([dp2dr(:,2),dp2dr(:,1)],vYEdge,vXEdge); #2D-hist without guide points


nXBins = length(vXEdge);
nYBins = length(vYEdge);
vXLabel = 0.5*(vXEdge(1:(nXBins-1))+vXEdge(2:nXBins));
vYLabel = 0.5*(vYEdge(1:(nYBins-1))+vYEdge(2:nYBins));
pcolor(vXLabel, vYLabel, mHist2d); 
shading flat; #means no border around each hist rectangle
w2r= vertcat(gp2dr, gp2dr(1,:)); #create guide line matrix
hold on
plot (w2r(:,1), w2r(:,2))        #plot guide line matrix
hold off

####create a colour map for many small values

N= 256;
j= jet(N+1);
for n=1:1:N
  i=round((n/N)^(1/2) * N + 1);
  cmap(n,:)= j(i,:);
end

cmap=vertcat([1,1,1],cmap);

####end colour map cration

colormap(cmap)
#colormap(hsv(128))
#caxis([0, 10])#ignore extreme etremes
colorbar #show colorbar
axis ("square");

####printing again...

print('ellipsoid_oct02_04.png', '-dpng');#, '-r100');
print('ellipsoid_oct02_04.svg', '-dsvg');

####printing end
