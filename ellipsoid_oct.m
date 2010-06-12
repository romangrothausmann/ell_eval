#!/net/home/ftd/localusr/bin/octave -qf

#skript to visualize ellipsoid axis ratios
#this is useful to identify amount of "smarties" and "cigars"
######
#REALLY MAKE SURE TO NOT MIX UP THETA AND PHI!!!!!!!


clear all;

#awk 'NR>1 {print $6, $7, $8, $9, $10, $11, $3, $4, $5}' pudel/az239b_080725a_bak02_3D-ana.csv > octave/az239b_080725a_bak02_3D-ana_I.txt
t= load I.txt; #octave_test02.txt;
#t=rand(200,6);

lbb= 0;

arg_list = argv ();
if nargin
  if (arg_list{1} == "-b")
    lbb= 1;
    printf("Doing BBox fitting...\n");
    endif
else
  printf("Doing I fitting...\n");
endif



#awk 'NR>1 {print $12, $13, $14}' az239b_080725a_bak02_3D-ana.csv > octave/az239b_080725a_bak02_3D-ana_BBox.txt   
if (lbb)
  bb= load BBox.txt; #ellipsoid_bbox.txt;
endif

#bb= sort(bb')'; #don't sort here!!! it's done below

N= size(t, 1);
#m= zeros(N,12);
#e= zeros(N,3);
num= 101;
phi0=pi(1)/4;
lambda0=pi(1)/4;

if (lbb)
  [fid, msg] = fopen ("BBox-fit.txt", "w");
else
  [fid, msg] = fopen ("I-fit.txt", "w");
endif

for n=1:1:N;

  I= [t(n,1),t(n,4),t(n,5);
      t(n,4),t(n,2),t(n,6);
      t(n,5),t(n,6),t(n,3)];
    #i=rand(3,3);
  [v, l]= eig(I);
    #l=rand(3,3);

  #test if eigenvalues are correlated to their eigenvectors
  if (I*v != [l(1,1)*v(:,1)'; l(2,2)*v(:,2)'; l(3,3)*v(:,3)']')
    I*v
    [l(1,1)*v(:,1); l(2,2)*v(:,2); l(3,3)*v(:,3)]
    fflush(stdout);
  end

  if (lbb)
    ax= [bb(n,1), bb(n,2), bb(n,3)]; #axis fitted to BBox
  else
    #ax= [l(1,1), l(2,2), l(3,3)]; #half axis fitted to I
    ax= [2*l(1,1), 2*l(2,2), 2*l(3,3)];#axis fitted to I
  endif
  
  [axs, axi]= sort (ax); #making a < b < c if axis-names are not specially assigned!
  #eu= euler_angles(v(:,1), v(:,2), v(:,3));
  eu= euler_angles(v(:,axi(1)), v(:,axi(2)), v(:,axi(3))); #index ordered v

  #ll= sort(ax); #making a < b < c if axis-names are not specially assigned!
  [theta, phi, r]= cart2sph(axs(1), axs(2), axs(3)); 

  e1(n,:)= [theta, phi, r];

 # if (n == 37)
 #   i
 #   l
 #   ll
 #   e1(n,:)
 #   fflush(stdout);
 # end

  #r= horzcat (reshape (v, 1, 9), l(1,1), l(2,2), l(3,3));
  #r= horzcat (l(1,1), l(2,2), l(3,3), eu, t(n,7), t(n,8), t(n,9));
  #m(n,:)=r;
  #fprintf(fid, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", ax(1), ax(2), ax(3), eu(1), eu(2), eu(3), t(n,7), t(n,8), t(n,9));
  fprintf(fid, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", axs(1), axs(2), axs(3), eu(1), eu(2), eu(3), t(n,7), t(n,8), t(n,9));
end;

fclose(fid);

w =acos(dot([1,1,1], [1,1,0])/(norm([1,1,1]) * norm([1,1,0])));
#a0=[linspace(0,pi(1)/2,num); zeros(1, num)];
a0=[pi(1)/4, pi(1)/2]';
#b0=[ ones(1, num) * pi(1)/2; linspace(0,pi(1)/2,num)];
#b0=[ ones(1, num) * pi(1)/2; linspace(pi(1)/4,pi(1)/2,num)];
b0=[pi(1)/2, pi(1)/4]';
#c0=[zeros(1, num); linspace(0,pi(1)/2,num)];
c0=[pi(1)/4, w]';
#a1=[linspace(0,pi(1)/2,num); ones(1, num) * w];
#c1=[ones(1, num) * pi(1)/4; linspace(0,pi(1)/2,num)]; #only good for 2D view!
#c1=
#d= horzcat(a0, b0, c0, a1, c1);
d= horzcat(a0, b0, c0);
#d= [d(1,:); d(2,:)];
d= horzcat(d, e1(:,1:2)');
#d= e1(:,1:2)';

for n=1:1:size(d,2);
  k= 2 / (1 + sin(phi0) * sin(d(2,n)) + cos(phi0) * cos(d(2,n)) * cos(d(1,n) - lambda0));
  x= k * cos(d(2,n)) * sin(d(1,n) - lambda0);
  y= k *     (cos(phi0) * sin(d(2,n)) - sin(phi0) * cos(d(2,n)) * cos(d(1,n) - lambda0));
  p2(n,:)= [x, y];

  [x1,y1,z1]= sph2cart (d(1,n), d(2,n), 1);
  e3(n,:)= [x1,y1,z1];
end;

[Theta, R]= cart2pol(p2(:,1), p2(:,2));
[p3(:,1), p3(:,2)]= pol2cart(Theta - 30 / 180 * pi, R);


#scatter3 (e1(:,1),e1(:,2),e1(:,3), [], 1);
scatter3 (e3(:,1),e3(:,2),e3(:,3), [], 1); #"markersize", 3, 1);
axis ("square");

azimuth= 135;
#azimuth= 315;
elevation= acosd(dot([1,1,1], [1,1,0])/(norm([1,1,1]) * norm([1,1,0])));
#elevation= elevation + 90;
view(azimuth, elevation);

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

scatter (p2(:,1), p2(:,2), [], 1)# p2(:,1));
axis ("square");

####printing again...

print('ellipsoid_oct02_02.png', '-dpng');#, '-r100');
print('ellipsoid_oct02_02.svg', '-dsvg');

####printing end

#save -ascii m.txt m;

bin= 50;

xbin=bin;
ybin=round( (p2(1,2) - p2(3,2)) / (p2(2,1) - p2(3,1)) * xbin);

vXEdge = linspace(p2(3,1)-0/xbin,p2(2,1)+1/xbin,xbin); #p2(3,:)<>c0; p2(2,:)<>b0; p2(1,:)<>a0
vYEdge = linspace(p2(3,2)-1/ybin,p2(1,2)+2/ybin,ybin);
mHist2d = hist2d([p2(:,2),p2(:,1)],vYEdge,vXEdge);


nXBins = length(vXEdge);
nYBins = length(vYEdge);
vXLabel = 0.5*(vXEdge(1:(nXBins-1))+vXEdge(2:nXBins));
vYLabel = 0.5*(vYEdge(1:(nYBins-1))+vYEdge(2:nYBins));
pcolor(vXLabel, vYLabel, mHist2d); 
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
axis ("square");

####printing again...

print('ellipsoid_oct02_03.png', '-dpng');#, '-r100');
print('ellipsoid_oct02_03.svg', '-dsvg');

####printing end


clear bin xbin ybin vXEdge vYEdge mHist2d nXBins nYBins vXLabel vYLabel

bin= 50;

xbin=bin;
ybin=round( (p3(1,2) - p3(3,2)) / (p3(2,1) - p3(3,1)) * xbin);

vXEdge = linspace(p3(3,1)-0/xbin,p3(2,1)+1/xbin,xbin); #p3(3,:)<>c0; p3(2,:)<>b0; p3(1,:)<>a0
vYEdge = linspace(p3(3,2)-1/ybin,p3(1,2)+2/ybin,ybin);
mHist2d = hist2d([p3(:,2),p3(:,1)],vYEdge,vXEdge);


nXBins = length(vXEdge);
nYBins = length(vYEdge);
vXLabel = 0.5*(vXEdge(1:(nXBins-1))+vXEdge(2:nXBins));
vYLabel = 0.5*(vYEdge(1:(nYBins-1))+vYEdge(2:nYBins));
pcolor(vXLabel, vYLabel, mHist2d); 
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
axis ("square");

####printing again...

print('ellipsoid_oct02_04.png', '-dpng');#, '-r100');
print('ellipsoid_oct02_04.svg', '-dsvg');

####printing end
