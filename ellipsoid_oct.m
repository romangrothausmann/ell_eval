#!/net/home/ftd/localusr/bin/octave -qf
#this script needs the output from the ITK-program: analyse02 or label_ana01 

#script to visualize ellipsoid axis ratios
#this is useful to identify amount of "smarties" and "cigars"
######
#REALLY MAKE SURE TO NOT MIX UP THETA AND PHI!!!!!!!

#from ellipsoid_oct08.m
#adding global axes orientation evaluation
#2D-hist colormap according to max of all three axes
#corrected 2D-hist range for all 2D-hists!!!
#added two further separation lines
#no sep. lines according abs. error in a,b,c possible (unsolv. polynomial or 4th degree)
#doing it numerically (no error propagation taken into account) and with colouring
#outputing particle surface for blender
#taking erode/dilate error into account
#introduced fuzzy logic for ell-type evaluation
#using only delated error because labels are missing in eroded ws 

clear all;



ee_min= 0.00000000000001;

arg_list = argv ();
if nargin != 3
  printf("Usage: %s <analysis.txt> <erode_analysis.txt> <dilate_analysis.txt>\n", program_name);
  exit(1)
else
  printf("Evaluating ellipsoids from %s...\n", arg_list{1});
  t= load(arg_list{1}); #octave_test02.txt;
  printf("Evaluating eroded analysis from %s...\n", arg_list{2});
  te= load(arg_list{2}); #octave_test02.txt;
  printf("Evaluating dilated analysis from %s...\n", arg_list{3});
  td= load(arg_list{3}); #octave_test02.txt;
endif



addpath("~/octave")
gnuplot_binary ("gnuplot -geometry 800x800"); 
#gnuplot_binary ("GNUTERM=wxt gnuplot -geometry 800x800"); 
#gnuplot_binary ("sed 's/ pt 6 / pt 5 /g' | gnuplot -geometry 800x800"); 
#gnuplot_binary ("tee octave.gp | gnuplot -V");
#gnuplot_binary ('tee octave.gp');#this is not possible, octave checks for gnuplot version!
#gset terminal dump
set (0, 'defaulttextfontname', 'arial');
#use: cat .Xresources
#! gnuplot options
#! modify this for a convenient window size
#gnuplot*geometry: 600x600

#da= 1*.26 #const. abs. error in a; could be read from file for each a individually
#daxs= [da,da,da] #same abs. error for all axes
sm_min=6 #has to be > length(es)/2
ci_min=sm_min

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
tsum= 0;
#es=[ 1, 1, 1;
#     1,-1, 1;
#     1,-1, 1;]

es=[0,0,0];
for x=0:1:1
  for y=0:1:1
    for z=0:1:1
      es=vertcat([x*2-1,y*2-1,z*2-1],es);
    end
  end
end

es=es(1:length(es)-1,:);

if sm_min <= length(es)/2 #has to be > length(es)/2
  exit(1)
endif

[fid, msg] = fopen ("I-fit.txt", "w");
fprintf(fid, \
        "#ell_a\tell_b\tell_c\tell_x\tell_y\tell_z\ta_x\ta_y\ta_z\tb_x\tb_y\tb_z\tc_x\tc_y\tc_z\tell_t\tindex\tp_surf\n");

for n=1:1:N;

  if (t(n,1) != td(n,1))
    printf("Indexs don't match! Aborting!\n")
    exit(1)
  endif
  if (t(n,1) != te(n,1))
    printf("Indexs don't match! Aborting!\n")
    exit(1)
  endif

  p_index= t(n,1);
  p_pos=  [t(n,2),t(n,3),t(n,4)];
  ax= [t(n,5),t(n,6),t(n,7)];
  v=  [t(n,8),t(n,9),t(n,10);
       t(n,11),t(n,12),t(n,13);
       t(n,14),t(n,15),t(n,16),]';
  p_V= t(n,17);
  
  vol= 4/3*pi*ax(1)*ax(2)*ax(3);
  if(vol>0)
    ax= 2*(p_V/vol)^(1/3)*ax;
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

  tdax= [td(n,5),td(n,6),td(n,7)];
  teax= [te(n,5),te(n,6),te(n,7)];
  [tdaxs, tdaxi]= sort (tdax); #making a < b < c if axis-names are not specially assigned!
  [teaxs, teaxi]= sort (teax); #making a < b < c if axis-names are not specially assigned!


  is_sm= 0;
  is_ci= 0;

  #this error evaluation can be done with fuzzy logic!
  #this makes sence if the error criterion is very hard and often only a single exception causes uncertainty
  #with fuzzy logic a minimum of certain evaluations have to be met
  #(fuzzy logic with reals is possible but very difficult: the ratio of the two volumes created by the cut of the separation surface through the error box)
  esum= 0;
  for i=1:1:length(es) 
    for j=1:1:size(es,2)
      daxs(j)= abs(axs(j) - tdaxs(j)); #since axes may change abs necessary to have a positive error
      printf("t(n,1): %f; n: %d; e: %f; axs: %f; tdaxs: %f\n", t(n,1),n, daxs(j), axs(j), tdaxs(j))
      esum= esum + daxs(j);
    endfor

    ee(1)= (axs(1) + es(i,1)*daxs(1));
    ee(2)= (axs(2) + es(i,2)*daxs(2));
    ee(3)= (axs(3) + es(i,3)*daxs(3));

    for ii=1:1:length(ee)
      if ee(ii) < 0
        ee(ii)= ee_min;
      endif
    endfor

    #keine Projektion vorher nötig, da jeder Faktor sich rauskürzt!
    is_sm+= ((ee(1) / ee(2)) < (ee(2) / ee(3)));
    is_ci+= ((ee(1) / ee(2)) > (ee(2) / ee(3)));
    printf("t(n,1): %f; i: %d; a/b: %f; b/c: %f; sm: %d; ci: %d\n", t(n,1),i,(ee(1) / ee(2)), (ee(2) / ee(3)),is_sm,is_ci)
  endfor
  mee= esum / size(es,1) / size(es,2)
  tsum= tsum + mee;
  
  if ((is_sm >= sm_min) && (is_ci >= ci_min))
    printf("is_sm is_ci true! This shlould never happen! t(n,1): %f\n", t(n,1))
    return #stop for --persist
  endif

  if is_sm >= sm_min
    if (axs(1) / axs(2) > axs(2) / axs(3))
      printf("is_sm wrong! This can happen for very big errors. \
          Couting as uncertain! t(n,1): %f\n", t(n,1))
      Nsz++;
      ce(:,n)= [0,0,1];
      et= 0;
      #return #stop for --persist
      continue
    endif
    Ns++;
    ce(:,n)= [0,1,0];
    et= 1;
  else 
    if is_ci >= ci_min
      if (axs(1) / axs(2) < axs(2) / axs(3))
        printf("is_ci wrong! t(n,1): %f\n", t(n,1))
        return #stop for --persist
      endif
      Nz++;
      ce(:,n)= [1,0,0];
      et= -1;
    else
      Nsz++;
      ce(:,n)= [0,0,1];
      et= 0;
    endif
  endif

  #if (axs(2) / axs(3) > axs(1) / axs(2))
  #  u(1,:,n)= [0, 0, 0];
  #  u(2,:,n)= [0, 0, 0];
  #  u(3,:,n)= [0, 0, 0];
  #else 
  #  u(1,:,n)= v(:,1);
  #  u(2,:,n)= v(:,2);
  #  u(3,:,n)= v(:,3);
  #endif

  u(1,:,n)= v(:,axi(1));
  u(2,:,n)= v(:,axi(2));
  u(3,:,n)= v(:,axi(3));
  #u(1,:,n)= v(:,1);
  #u(2,:,n)= v(:,2);
  #u(3,:,n)= v(:,3);


  [theta, phi, r]= cart2sph(axs(1), axs(2), axs(3)); 
  dp3ds(:,n)= [theta, phi, r];


  #fprintf(fid, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", ax, t(n,7:9), v);
  fprintf(fid, \
          "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\n", \
           axs, p_pos, v(:,axi(1)), v(:,axi(2)), v(:,axi(3)), et, p_index, t(n,18));
end;

fclose(fid);
printf("# smarty-like: %d; # cigar-like: %d; # uncertain: %d; ratio: %.2f\n", Ns, Nz, Nsz, Ns/Nz);
tsum / N

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
#scatter3 (dp3d(1,:),dp3d(2,:),dp3d(3,:), [], "b"); #"markersize", 3, 1);
color=vertcat(sqrt(2)*dp3d(2,:),sqrt(3)*dp3d(1,:),1/(1-sqrt(1/3))*(dp3d(3,:)-sqrt(1/3)));
#color=vertcat(sqrt(3)*dp3d(1,:),sqrt(2)*dp3d(2,:),1/(1-sqrt(1/3))*(dp3d(3,:)-sqrt(1/3)));



#z= ones(1,size(x,1));
#scatter3(ua, ub, uc, 50, color, 's')
#axis ([-1,1,-1,1,-1,1],"square");

#break

#scatter3 (dp3d(1,:),dp3d(2,:),dp3d(3,:), 50, color, 's'); #'d'
scatter3 (dp3d(1,:),dp3d(2,:),dp3d(3,:), 50, ce, 's'); #'d'
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

#break
print('ellipsoid_oct02_01.svg', '-dsvg', '-S800,800');
print('ellipsoid_oct02_01.png', '-dpng', '-S800,800');#, '-F/usr/X11R6/lib/X11/fonts/msttf/arial.ttf');#, '-r100');


view(110, 10);

print('ellipsoid_oct02_01b.png', '-dpng', '-S800,800');#, '-F/usr/X11R6/lib/X11/fonts/msttf/arial.ttf');#, '-r100');
print('ellipsoid_oct02_01b.svg', '-dsvg', '-S800,800');

####printing end
#break




###global axes orientations 
##with sinusoidal projection: http://mathworld.wolfram.com/SinusoidalProjection.html

clear phi0 the0
clear xbin ybin vXEdge vYEdge mHist2d nXBins nYBins vXLabel vYLabel

l00= pi/2; #projection centre
bin= 21;#make sure binning fits scatter plot impression!



phi0= linspace(-pi/2,pi/2,num);
the0= pi .* ones(1,num);
l0(1,:)=   (the0 - pi*0/2) .* cos(phi0);
l0(2,:)=   phi0;
l90(1,:)=  (the0 - pi*1/2) .* cos(phi0);
l90(2,:)=  phi0;
l180(1,:)= (the0 - pi*2/2) .* cos(phi0);
l180(2,:)= phi0;
l270(1,:)= (the0 - pi*3/2) .* cos(phi0);
l270(2,:)= phi0;
l360(1,:)= (the0 - pi*4/2) .* cos(phi0);
l360(2,:)= phi0;
 
for n=1:1:3

  ua= squeeze(u(n,1,:));
  ub= squeeze(u(n,2,:));
  uc= squeeze(u(n,3,:));

  ##ellipsoids are symmetry along their axes so restrict the orientation
  ##to hemisphere
  for i=1:1:size(ua)
    if ub(i) < 0
      ua(i)= -ua(i);
      ub(i)= -ub(i);
      uc(i)= -uc(i);
    endif
  endfor

  clear theta phi r x y
  [theta, phi, r]= cart2sph(ua, ub, uc); 

  ##ellipsoids are symmetry along their axes so restrict the orientation
  ##to hemisphere
  #for i=1:1:size(theta)
  #  if theta(i) < 0
  #    theta(i)= theta(i) + pi;
  #    phi(i)= -phi(i);
  #  endif
  #endfor

  x= (theta - l00) .* cos(phi);
  y= phi;

  #scatter(x, y, 500, color, 's')
  scatter(x, y, 500, ce, 's')
  hold on
  #plot (  l0(1,:),   l0(2,:), "k")
  plot ( l90(1,:),  l90(2,:), "k")
  plot (l180(1,:), l180(2,:), "k")
  plot (l270(1,:), l270(2,:), "k")
  #plot (l360(1,:), l360(2,:), "k")
  hold off

  axis ("square");

  ####printing now...

  fn= sprintf("%s_%d.svg","ellipsoid_oct02_05" , n);
  print(fn, '-dsvg', '-S800,800');

  fn= sprintf("%s_%d.png","ellipsoid_oct02_05" , n);
  print(fn, '-dpng', '-S800,800');

  ####printing end

  ###2D-hist now

  #lef= min(x);
  #rig= max(x);
  #bot= min(y);
  #top= max(y);
  lef= min(l270(1,:));
  rig= max( l90(1,:));
  bot= min(l270(2,:));
  top= max(l270(2,:));
  dx= (rig - lef);
  dy= (top - bot);

  xbin=bin;
  ybin=round( dy / dx * xbin); #make'm squares
  #xbin= round(bin / dx * dxr); #make the squares the same size as before
  #ybin= round(dyr / dxr * xbin); #make'm squares

  #x= vertcat(x, l0(1,:)');#for testing guide line pos
  #y= vertcat(y, l0(2,:)');

  vXEdge = linspace(lef, rig+(rig-lef)/(xbin-1), xbin+1);#add one entry at end
  vYEdge = linspace(bot, top+(top-bot)/(ybin-1), ybin+1);
  mHist2d = hist2d([y, x],vYEdge,vXEdge); #2D-hist without guide points

  vXEdget(n,:) = vXEdge;
  vYEdget(n,:) = vYEdge;
  mHist2dt(n,:,:) = mHist2d;

  
  Nt(n)= max(max(mHist2d));

endfor

N= max(Nt);
cmapg= jet(N + 1);
cmapg(1,:)=[1,1,1];

for n=1:1:3

  vXEdge = vXEdget(n,:);
  vYEdge = vYEdget(n,:);
  mHist2d = squeeze(mHist2dt(n,:,:));

  nXBins = length(vXEdge);
  nYBins = length(vYEdge);
  #vXLabel = 0.5*(vXEdge(1:(nXBins-1))+vXEdge(2:nXBins));
  #vYLabel = 0.5*(vYEdge(1:(nYBins-1))+vYEdge(2:nYBins));
  vXLabel = vXEdge(1:(nXBins-1));
  vYLabel = vYEdge(1:(nYBins-1));
  
  set (gca, 'xtick', "");#the ticks aren't correct!
  set (gca, 'ytick', "");
  
  pcolor(vXLabel, vYLabel, mHist2d); #mHist2D acts as color value
  #imagesc(mHist2d);
  hold on
  #plot (  l0(1,:),   l0(2,:), "k")
  plot ( l90(1,:),  l90(2,:), "k")
  plot (l180(1,:), l180(2,:), "k")
  plot (l270(1,:), l270(2,:), "k")
  #plot (l360(1,:), l360(2,:), "k")
  hold off

  shading flat; #means no border around each hist rectangle
  #shading faceted

  #caxis([0, Nt(n)]);
  cmap= cmapg(1:Nt(n)+1,:);
  colormap(cmap);
  colorbar #show colorbar
  axis ("square");#setting axis range here can be bad!

  ####printing now...

  fn= sprintf("%s_%d.svg","ellipsoid_oct02_05h" , n);
  print(fn, '-dsvg', '-S800,800');

  fn= sprintf("%s_%d.png","ellipsoid_oct02_05h" , n);
  print(fn, '-dpng', '-S800,800');

  ####printing end


end

#break


#scatter (dp2d(1,:), dp2d(2,:), [], 1)
#scatter (dp2d(1,:), dp2d(2,:), 500, color, 's')
scatter (dp2d(1,:), dp2d(2,:), 500, ce, 's')
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
lef= c00p(1,1);
rig= b0p(1,1);
bot= c00p(2,1);
top= c0p(2,1);
dx= (rig - lef);
dy= (top - bot);

xbin=bin;
ybin=round( dy / dx * xbin); #make'm squares

vXEdge = linspace(lef, rig+(rig-lef)/(xbin-1), xbin+1);
vYEdge = linspace(bot, top+(top-bot)/(ybin-1), ybin+1);
mHist2d = hist2d([dp2d(2,:)',dp2d(1,:)'],vYEdge,vXEdge); #2D-hist without guide points


nXBins = length(vXEdge);
nYBins = length(vYEdge);
vXLabel = vXEdge(1:(nXBins-1));
vYLabel = vYEdge(1:(nYBins-1));

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

function res= shear_y(data, m);
res= vertcat (data(1,:), data(2,:) - m * data(1,:));
endfunction

m= (b0p(2,1) - c00p(2,1)) / (b0p(1,1) - c00p(1,1));
#ra= -atan((b0p(2,1) - c00p(2,1)) / (b0p(1,1) - c00p(1,1))); #ra= 31.325°; why not 30°, stereographic projection error?

c00r= shear_y(c00p, m);
c0r= shear_y(c0p, m);
b0r= shear_y(b0p, m);
a0r= shear_y(a0p, m);
s0r= shear_y(s0p, m);
dp2dr= shear_y(dp2d, m);
#dp2dr= shear_y(horzcat(dp2d, a0p), m);#for testing hist squares straight


clear xbin ybin vXEdge vYEdge mHist2d nXBins nYBins vXLabel vYLabel

lef= a0r(1,1);
rig= b0r(1,1);
bot= a0r(2,1);
top= c0r(2,1);
dxr= (rig - lef);
dyr= (top - bot);

xbin= round(bin / dx * dxr); #make the squares the same size as before
ybin= round(dyr / dxr * xbin); #make'm squares

vXEdge = linspace(lef, rig+(rig-lef)/(xbin-1), xbin+1);
vYEdge = linspace(bot, top+(top-bot)/(ybin-1), ybin+1);

#vXEdge = linspace(a0r(1,1)-1/xbin, b0r(1,1)+1/xbin, xbin); 
#vYEdge = linspace(a0r(2,1)-1/ybin, c0r(2,1)+2/ybin, ybin);
mHist2d = hist2d([dp2dr(2,:)',dp2dr(1,:)'],vYEdge,vXEdge); 


nXBins = length(vXEdge);
nYBins = length(vYEdge);
vXLabel = vXEdge(1:(nXBins-1));
vYLabel = vYEdge(1:(nYBins-1));

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

text (c00p(1,1) - .02, c00p(2,1) + .02, "sphere\npoint", "horizontalalignment", "right"); #looks nicer
text (b0p(1,1) + .02, b0p(2,1), "circle\npoint");
text (c0p(1,1) - .02, c0p(2,1), "line\npoint", "horizontalalignment", "right");
text (a0p(1,floor(num/4)) + .02, a0p(2,floor(num/4)), "oblate line", "rotation", 30);
text (b0p(1,floor(num/2)) + .02, b0p(2,floor(num/2)), "ellipse arc", "rotation", -50);
text (c0p(1,floor(num/4*3)) - .02, c0p(2,floor(num/4*3)), "prolate line", "rotation", 90);
text (s0p(1,size(s0p,2)-20) + .02, s0p(2,size(s0p,2)-20), "separation curve", "rotation", -75);
text (c00p(1,1), c00p(2,1) - .1, sprintf("# oblate-like: %d; # prolate-like: %d; # uncertain: %d; ratio: %.2f\n", Ns, Nz, Nsz, Ns/Nz));

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



