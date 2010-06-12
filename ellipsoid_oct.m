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

gnuplot_binary ("GNUTERM=wxt gnuplot -geometry 800x800"); 
set (0, 'defaulttextfontname', 'arial');
#use: cat .Xresources
#! gnuplot options
#! modify this for a convenient window size
#gnuplot*geometry: 600x600


N= size(t, 1);
#m= zeros(N,12);
#e= zeros(N,3);
num= 101;
radius= 1;
phi0=pi(1)/4;
lambda0=pi(1)/4;
mscale= 20/9;
mass= 1;

Ns= 0;
Nz= 0;
Nsz= 0;


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
#  if (I*v != [l(1,1)*v(:,1)'; l(2,2)*v(:,2)'; l(3,3)*v(:,3)']')
#    n
#    I*v
#    [l(1,1)*v(:,1)'; l(2,2)*v(:,2)'; l(3,3)*v(:,3)']'
#    fflush(stdout);
#    #break
# end

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
    fflush(stdout);
    break
  endif

  [axs, axi]= sort (ax); #making a < b < c if axis-names are not specially assigned!
  #eu= euler_angles(v(:,axi(1)), v(:,axi(2)), v(:,axi(3))); #index ordered v

  if (axs(2) / axs(3) > axs(1) / axs(2))
    Ns++;
  else 
    if (axs(2) / axs(3) < axs(1) / axs(2))
      Nz++;
    else
      Nsz++;
    endif
  endif

  [theta, phi, r]= cart2sph(axs(1), axs(2), axs(3)); 
  dp3ds(:,n)= [theta, phi, r];

  #fprintf(fid, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", ax, t(n,7:9), v);
  fprintf(fid, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", axs, t(n,7:9), v(:,axi(1)), v(:,axi(2)), v(:,axi(3)));
end;

printf("# smarty-like: %d; # cigar-like: %d; # on edge: %d\n", Ns, Nz, Nsz);

fclose(fid);
#break

ws =acos(dot([1,1,1], [1,1,0])/(norm([1,1,1]) * norm([1,1,0])));
#a0=[linspace(0,pi/2,num); zeros(1, num)];# 1/4 circle in ab-plane
#a0=[pi/4, pi/2]';#[1,0,0] point?
#b0=[ones(1, num) * pi/2; linspace(0,pi/2,num)];#arc from [0, 1, 0] to a0
b0s=[ones(1, num) * pi/2; linspace(pi/4,pi/2,num)];#arc from b0 to a0 
#b0=[pi/2, pi/4]';#[0, 1/sqrt(2), 1/sqrt(2)] point
c00= [1/sqrt(3); 1/sqrt(3); 1/sqrt(3)];
c00s=[pi/4; ws]';#[1/sqrt(3), 1/sqrt(3), 1/sqrt(3)] point
c0s=[ones(1, num) * pi/4; linspace(pi/2,ws,num)];#arc from w to a0
#c1=[ones(1, num) * pi/4; linspace(0,pi/2,num)]; #arc from [1/sqrt(2), 1/sqrt(2), 0] to a0; only good for 2D view!
#a1=[linspace(0,pi/2,num); ones(1, num) * w];#arc from [1/sqrt(3), 0, sqrt(2/3)] to  [0, 1/sqrt(3), sqrt(2/3)]


#### guide points for b/c > a/b
xsn= 100;
#xs= linspace(sqrt((sqrt(33)-1)/2),100,xsn); #start point corr. to c0
#xs= linspace(1/xsn,1,xsn); 
xs= linspace(1,0.001,xsn); 
xs= xs.*xs; #make'm more evenly spaced; element by element multiplication

for n=1:1:xsn
  #den= (1 + xs(n)^2 + xs(n)^4)^(1/4);
  den= sqrt(1/xs(n)^2+1/xs(n)+1);
  sx= 1 / den;
  sy= 1 / den / sqrt(xs(n));
  sz= 1 / den / xs(n);
  s0(:,n)= [sx, sy, sz]; #separation points for b/c > a/b 
end
s0= horzcat(s0, [0,0,1]');#add c0 at end

#### guide points for b== 1/sqrt(3)
clear xs
xs= linspace(1/sqrt(3),0,num);# b== 1/sqrt(3) line
for n=1:1:num
  #s1x= xs;
  #s1y= 1/sqrt(3);
  s1z= sqrt(2/3 - xs(n)^2);
  s1(:,n)= [xs(n), 1/sqrt(3), s1z];
end

#### guide points for b== c
clear xs
xs= linspace(1/sqrt(3),0,num);# b== 1/sqrt(3) line
for n=1:1:num
  a0y= sqrt((1-xs(n)^2)/2);
  a0(:,n)= [xs(n), a0y, a0y];
end

#d0= horzcat(c0, b0, c0); #start from c0 over b0s to c0 over s0
#d0s= horzcat(c00s, b0s, c1s); #start from c0 over b0s to c1 over s0
#gp3d0s= vertcat (d0s, ones(1,size(d0s,2))); #combine sph. guide points

#cart. for plotting
clear xt yt zt;
[xt, yt, zt]= sph2cart (b0s(1,:), b0s(2,:), ones(1,size(b0s,2)));#projection of 3D points onto unit sphere
b0= vertcat (xt, yt, zt);

clear xt yt zt;
[xt, yt, zt]= sph2cart (c0s(1,:), c0s(2,:), ones(1,size(c0s,2)));#projection of 3D points onto unit sphere
c0= vertcat (xt, yt, zt);

clear xt yt zt;
[xt, yt, zt]= sph2cart (dp3ds(1,:), dp3ds(2,:), ones(1,size(dp3ds,2)));#projection of 3D data points onto unit sphere
dp3d= vertcat (xt, yt, zt);

#sph. for stereog. proj.
clear xt yt zt;
[xt, yt, zt]= cart2sph (s0(1,:), s0(2,:), s0(3,:));#projection of 3D guide points onto unit sphere
s0s= vertcat (xt, yt, zt);

clear xt yt zt;
[xt, yt, zt]= cart2sph (s1(1,:), s1(2,:), s1(3,:));#projection of 3D points onto unit sphere
s1s= vertcat (xt, yt, zt);

clear xt yt zt;
[xt, yt, zt]= cart2sph (a0(1,:), a0(2,:), a0(3,:));#projection of 3D points onto unit sphere
a0s= vertcat (xt, yt, zt);

#stereographic projection on unit sphere
c00p= stereogproj(c00s(1), c00s(2), 1, phi0, lambda0);
c0p= stereogproj(c0s(1,:), c0s(2,:), 1, phi0, lambda0);
b0p= stereogproj(b0s(1,:), b0s(2,:), 1, phi0, lambda0);
a0p= stereogproj(a0s(1,:), a0s(2,:), 1, phi0, lambda0);
s0p= stereogproj(s0s(1,:), s0s(2,:), 1, phi0, lambda0);
dp2d= stereogproj(dp3ds(1,:), dp3ds(2,:), 1, phi0, lambda0);


###plotting 3D
scatter3 (dp3d(1,:),dp3d(2,:),dp3d(3,:), [], "b"); #"markersize", 3, 1);
hold on
plot3 (a0(1,:), a0(2,:), a0(3,:), "k")
plot3 (b0(1,:), b0(2,:), b0(3,:), "k")
plot3 (c0(1,:), c0(2,:), c0(3,:), "k")
plot3 (s0(1,:), s0(2,:), s0(3,:), "k")
#plot3 (s1(1,:), s1(2,:), s1(3,:), "k")
hold off

#axis ("square");
axis ([0,1,0,1,0,1],"square");

azimuth= 135;
#azimuth= 315;
elevation= acosd(dot([1,1,1], [1,1,0])/(norm([1,1,1]) * norm([1,1,0])));
#elevation= elevation + 90;
view(azimuth, elevation);

xlabel("a");
ylabel("b");
zlabel("c");
text (c00(1,1), c00(2,1), c00(3,1), "sphere\npoint", "horizontalalignment", "right"); #looks nicer
text (b0(1,1), b0(2,1), b0(3,1), "circle\npoint");
text (c0(1,1), c0(2,1), c0(3,1), "line\npoint", "horizontalalignment", "right");
text (a0(1,floor(num/4)), a0(2,floor(num/4)) + .05, a0(3,floor(num/4)), "oblate arc", "rotation", 30);
text (b0(1,floor(num/2)) - .05, b0(2,floor(num/2)), b0(3,floor(num/2)), "ellipse arc", "rotation", -50);
text (c0(1,floor(num/4*3)) + .05, c0(2,floor(num/4*3)), c0(3,floor(num/4*3)), "prolate arc", "rotation", 90);
text (s0(1,size(s0,2)-30) - .03, s0(2,size(s0,2)-30), s0(3,size(s0,2)-30), "separation curve", "rotation", -75);

#break
####printing now...

#paper_size = [640, 480];
#set (gcf, "paperunits", "inches")
#set (gcf, "papertype", "<custom>")
#set (gcf, "papersize", paper_size)
#set (gcf, "paperposition", [0, 0, paper_size])

#figure('Position', [0, 0, 600, 400]); 

print('ellipsoid_oct02_01.png', '-dpng', '-S800,800');#, '-F/usr/X11R6/lib/X11/fonts/msttf/arial.ttf');#, '-r100');
print('ellipsoid_oct02_01.svg', '-dsvg', '-S800,800');


view(110, 10);

print('ellipsoid_oct02_01b.png', '-dpng', '-S800,800');#, '-F/usr/X11R6/lib/X11/fonts/msttf/arial.ttf');#, '-r100');
print('ellipsoid_oct02_01b.svg', '-dsvg', '-S800,800');

####printing end

scatter (dp2d(1,:), dp2d(2,:), [], 1)
hold on
plot (a0p(1,:), a0p(2,:), "k")
plot (b0p(1,:), b0p(2,:), "k")
plot (c0p(1,:), c0p(2,:), "k")
plot (s0p(1,:), s0p(2,:), "k")
#plot3 (s1(1,:), s1(2,:), "k")
hold off
axis ([c00p(1,1), b0p(1,1), c00p(2,1), c0p(2,1), ],"square");

text (c00p(1,1) - .02, c00p(2,1) + .02, "sphere\npoint", "horizontalalignment", "right"); #looks nicer
text (b0p(1,1) + .02, b0p(2,1), "circle\npoint");
text (c0p(1,1) - .02, c0p(2,1), "line\npoint", "horizontalalignment", "right");
text (a0p(1,floor(num/4)) + .02, a0p(2,floor(num/4)), "oblate line", "rotation", 30);
text (b0p(1,floor(num/2)) + .02, b0p(2,floor(num/2)), "ellipse arc", "rotation", -50);
text (c0p(1,floor(num/4*3)) - .02, c0p(2,floor(num/4*3)), "prolate line", "rotation", 90);
text (s0p(1,size(s0p,2)-20) + .02, s0p(2,size(s0p,2)-20), "separation curve", "rotation", -75);

#xlabel("");
#ylabel("");
set (gca, 'xtick', "");#the ticks aren't correct!
set (gca, 'ytick', "");


#break

####printing again...

print('ellipsoid_oct02_02.png', '-dpng');#, '-r100');
print('ellipsoid_oct02_02.svg', '-dsvg');

####printing end

#save -ascii m.txt m;

bin= 30;
#dgp2d= vertcat(dp2d, gp2d(1,:), gp2d(size(gp2d,1)-1,:), gp2d(size(gp2d,1),:));

dx= (b0p(1,1) - c00p(1,1));
dy= (c0p(2,1) - c00p(2,1));

xbin=bin;
ybin=round( dy / dx * xbin); #make'm squares

vXEdge = linspace(c00p(1,1)-0/xbin, b0p(1,1)+1/xbin, xbin); #dgp2d(3,:)<>c0; dgp2d(2,:)<>b0; dgp2d(1,:)<>a0
vYEdge = linspace(c00p(2,1)-1/ybin, c0p(2,1)+2/ybin, ybin);
#mHist2d = hist2d([dgp2d(:,2),dgp2d(:,1)],vYEdge,vXEdge); #2D-hist with guide points
mHist2d = hist2d([dp2d(2,:)',dp2d(1,:)'],vYEdge,vXEdge); #2D-hist without guide points


nXBins = length(vXEdge);
nYBins = length(vYEdge);
vXLabel = 0.5*(vXEdge(1:(nXBins-1))+vXEdge(2:nXBins));
vYLabel = 0.5*(vYEdge(1:(nYBins-1))+vYEdge(2:nYBins));

set (gca, 'xtick', "");#the ticks aren't correct!
set (gca, 'ytick', "");

pcolor(vXLabel, vYLabel, mHist2d); #mHist2D acts as color value
hold on
plot (a0p(1,:), a0p(2,:), "k")
plot (b0p(1,:), b0p(2,:), "k")
plot (c0p(1,:), c0p(2,:), "k")
plot (s0p(1,:), s0p(2,:), "k")
#plot3 (s1(1,:), s1(2,:), "k")
hold off

shading flat; #means no border around each hist rectangle

#axis ([c00p(1,1), b0p(1,1), c00p(2,1), c0p(2,1), ],"square");#setting axis range here can be bad!

text (c00p(1,1) - .02, c00p(2,1) + .02, "sphere\npoint", "horizontalalignment", "right"); #looks nicer
text (b0p(1,1) + .02, b0p(2,1), "circle\npoint");
text (c0p(1,1) - .02, c0p(2,1), "line\npoint", "horizontalalignment", "right");
text (a0p(1,floor(num/4)) + .02, a0p(2,floor(num/4)), "oblate line", "rotation", 30);
text (b0p(1,floor(num/2)) + .02, b0p(2,floor(num/2)), "ellipse arc", "rotation", -50);
text (c0p(1,floor(num/4*3)) - .02, c0p(2,floor(num/4*3)), "prolate line", "rotation", 90);
text (s0p(1,size(s0p,2)-20) + .02, s0p(2,size(s0p,2)-20), "separation curve", "rotation", -75);

#xlabel("");
#ylabel("");
set (gca, 'xtick', "");#the ticks aren't correct!
set (gca, 'ytick', "");


####create a colour map for many small values

#N= 256;
N=max (max (mHist2d));
cmap= jet(N + 1);
cmap(1,:)=[1,1,1];
#cmap=vertcat([1,1,1],cmap);

#j= jet(N+1);
#for n=1:1:N
#  i=round((n/N)^(1/2) * N + 1); #sqrt gradient
#  cmap(n,:)= j(i,:);
#end
#cmap=vertcat([1,1,1],cmap);

####end colour map cration

colormap(cmap)
#colormap(hsv(128))
#caxis([0, 10])#ignore extreme etremes
colorbar #show colorbar
axis ("square");#setting axis range here can be bad!

####printing again...

print('ellipsoid_oct02_03.png', '-dpng');#, '-r100');
print('ellipsoid_oct02_03.svg', '-dsvg');

####printing end

#break


####rot. for 2nd 2D-hist

function res= rotate(data, angle);

[Theta, R]= cart2pol(data(1,:), data(2,:));
[datar(1,:), datar(2,:)]= pol2cart(Theta + angle, R);
res= vertcat (datar(1,:), datar(2,:));
endfunction

ra= -atan((b0p(2,1) - c00p(2,1)) / (b0p(1,1) - c00p(1,1))); #ra= 31.325°; why not 30°, stereographic projection error?
#ra= -30 / 180 * pi;#Here the hist squares form a straight

c00r= rotate(c00p, ra);
c0r= rotate(c0p, ra);
b0r= rotate(b0p, ra);
a0r= rotate(a0p, ra);
s0r= rotate(s0p, ra);
dp2dr= rotate(dp2d, ra);
#dp2dr= rotate(horzcat(dp2d, a0p), ra);#for testing hist squares straight


clear xbin ybin vXEdge vYEdge mHist2d nXBins nYBins vXLabel vYLabel

dxr= (b0r(1,1) - a0r(1,1));
dyr= (c0r(2,1) - a0r(2,1));

xbin= round(bin / dx * dxr); #make the squares the same size as before
ybin= round(dyr / dxr * xbin); #make'm squares

vXEdge = linspace(a0r(1,1)-0/xbin, b0r(1,1)+1/xbin, xbin); #dgp2d(3,:)<>c0; dgp2d(2,:)<>b0; dgp2d(1,:)<>a0
vYEdge = linspace(a0r(2,1)-1/ybin, c0r(2,1)+2/ybin, ybin);
#mHist2d = hist2d([dgp2d(:,2),dgp2d(:,1)],vYEdge,vXEdge); #2D-hist with guide points
mHist2d = hist2d([dp2dr(2,:)',dp2dr(1,:)'],vYEdge,vXEdge); #2D-hist without guide points


nXBins = length(vXEdge);
nYBins = length(vYEdge);
vXLabel = 0.5*(vXEdge(1:(nXBins-1))+vXEdge(2:nXBins));
vYLabel = 0.5*(vYEdge(1:(nYBins-1))+vYEdge(2:nYBins));

set (gca, 'xtick', "");#the ticks aren't correct!
set (gca, 'ytick', "");

pcolor(vXLabel, vYLabel, mHist2d); #mHist2D acts as color value
hold on
plot (a0r(1,:), a0r(2,:), "k")
plot (b0r(1,:), b0r(2,:), "k")
plot (c0r(1,:), c0r(2,:), "k")
plot (s0r(1,:), s0r(2,:), "k")
#plot3 (s1(1,:), s1(2,:), "k")
hold off

shading flat; #means no border around each hist rectangle

text (c00r(1,1) - .02, c00r(2,1) + .02, "sphere\npoint", "horizontalalignment", "right"); #looks nicer
text (b0r(1,1) + .02, b0r(2,1), "circle\npoint");
text (c0r(1,1) - .02, c0r(2,1), "line\npoint", "horizontalalignment", "right");
text (a0r(1,floor(num/4)), a0r(2,floor(num/4)) - .01, "oblate line");
text (b0r(1,floor(num/4*3)) + .02, b0r(2,floor(num/4*3)), "ellipse arc", "rotation", -80);
text (c0r(1,floor(num/4*3)) - .02, c0r(2,floor(num/4*3)), "prolate line", "rotation", 60);
text (s0r(1,size(s0r,2)-20) + .02, s0r(2,size(s0r,2)-20), "separation curve", "rotation", -105);

#xlabel("");
#ylabel("");
set (gca, 'xtick', "");#the ticks aren't correct!
set (gca, 'ytick', "");


N=max (max (mHist2d)); #max can be different
cmap= jet(N + 1);
cmap(1,:)=[1,1,1];
#cmap=vertcat([1,1,1],cmap);

colormap(cmap)
#colormap(hsv(128))
#caxis([0, 10])#ignore extreme etremes
colorbar #show colorbar
axis ("square");#setting axis range here can be bad!

####printing again...

print('ellipsoid_oct02_04.png', '-dpng');#, '-r100');
print('ellipsoid_oct02_04.svg', '-dsvg');

####printing end


printf("# smarty-like: %d; # cigar-like: %d; # on edge: %d\n", Ns, Nz, Nsz);
