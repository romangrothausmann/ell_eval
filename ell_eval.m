#!/usr/bin/octave -q
##this script needs the output from ITK-CLIs analyse_binary or analyse_labels

##script to visualize ellipsoid axis ratios
##this is useful to identify amount of "smarties" and "cigars"
#######
##REALLY MAKE SURE TO NOT MIX UP THETA AND PHI!!!!!!!

##from ellipsoid_oct08.m
##adding global axes orientation evaluation
##2D-hist colormap according to max of all three axes
##corrected 2D-hist range for all 2D-hists!!!
##added two further separation lines
##no sep. lines according abs. error in a,b,c possible (unsolv. polynomial or 4th degree)

#doing it numerically (no error propagation taken into account) and with colouring
#outputing particle surface for blender
#introduced fuzzy logic for ell-type evaluation
#rgb for global alignment and stereographic projection 
#intruducing sphere category blue; uncertain yellow
#draw special lines longer in 3D
#draw also lines for abr an bcr
#using local sphere sampling for global ell ori
#using c/sqrt(ab) as weights
#normalizing weights
#normalizing individual weights

#ellipsoid_oct29.m:
#modified to work @mhh on pc3G20768 with octave 3.6.2
#expects axis lenths (a,b,c) > 0; can often indirectly be achieved by vol > 64:
#e.g. awk 'NR>1 {if ( 64 < $17) {print $0}}' az239b_080725a_sirt_th-12707286.31_lr+255+0+1.ana > az239b_080725a_sirt_th-12707286.31_lr+255+0+1_v64.ana
#_30: fixed limits from dra for nw to [0,0.1]
#_31: creating a point-mesh for spherical coords as vis for publication using local_hemisphere_sampling_hist2d
#_32: added separate drawing of 3D hemispheres for pub vis

clear all;


ee_min= 0.00000000000001;

quiet= 0;
arg_list = argv ();
if nargin != 1
  printf("Usage: %s <analysis.txt> [-q]\n", program_name);
  exit(1)
else
  printf("Evaluating ellipsoids from %s...\n", arg_list{1});
  t= load(arg_list{1}); #octave_test02.txt;
  if nargin == 2
    if (arg_list{2} == "-q")
      quiet= 1;
    endif
  endif 
endif



addpath("~/octave/functions")
#graphics_toolkit fltk; #for octave 3.6.2, nice for viewing but not saving
graphics_toolkit gnuplot;
#gnuplot_binary ("/home/grothama/skripte/octave_gnuplot_exec"); #for octave 3.6.2
#set size square
#set (gca, "PlotBoxAspectRatioMode", "manual", "PlotBoxAspectRatio", [1 1 1]);

#figure, set(gca,'Size',"square"); #for octave 3.6.2, has no the ef
#figure, set(gcf,'position',[100 100 800 800]); #for octave 3.6.2, has no the effect of gnuplot> set size square
#gnuplot_binary ("gnuplot", "set size square"); #for octave 3.6.2
#gnuplot_binary ("gnuplot", "geometry 800x800"); #for octave 3.6.2
#gnuplot_binary ("/usr/bin/gnuplot -geometry 800x800 "); 
#gnuplot_binary ("GNUTERM=wxt gnuplot -geometry 800x800"); 
#gnuplot_binary ("sed 's/ pt 6 / pt 5 /g' | gnuplot -geometry 800x800"); 
#gnuplot_binary ("tee octave.gp | gnuplot -V");
#gnuplot_binary ('tee octave.gp');#this is not possible, octave checks for gnuplot version!
#gset terminal dump

#figure (1)
#clf ()

set (0, 'defaulttextfontname', 'arial');
#set (gca (), "plotboxaspectratio", [1 1 1.45])#empirical ratio, needs to be set after axis ([0,1,0,1,0,1],"square");!!!, to make output look as if axis ([0,1,0,1,0,1], "equal"), axis ([0,1,0,1,0,1], "equal"); #resets plotboxaspectratio to 1 1 1!!! axis ([0,1,0,1,0,1], "square"); does not reset plotboxaspectratio!!!!!
#use: cat .Xresources
#! gnuplot options
#! modify this for a convenient window size
#gnuplot*geometry: 600x600


#ps3d=  40;
#ps2d= 400;
ps3d= 3;#for octave 3.6.2
ps2d= 3;#for octave 3.6.2
out3D=sprintf("%s_AR", arg_list{1});
outGO=sprintf("%s_GO", arg_list{1});
nplot= 0;


da= 2 #const. abs. error in a; could be read from file for each a individually
daxs= [da,da,da] #same abs. error for all axes
sm_min=6 #has to be > length(es)/2, can only reach up to 8!!!
ci_min=sm_min
sp_min=6 #can only reach up to 6!!!
spe=.2 #half sphere cube width
ev= ones(1,3)/norm(ones(1,3)); #unit vector in [1,1,1] direction

N= size(t, 1);
#m= zeros(N,12);
#e= zeros(N,3);
num= 101;
radius= 1;
phi0=acos(dot(ev,[1,1,0]/norm([1,1,0]))) #acos(dot(ev,[0,0,1])) #0 #pi(1)/4;
lambda0=pi(1)/4
abr3=1/3 #a/b ratio at which to draw a line for vis
bcr3=1/3 #b/c ratio at which to draw a line for vis
abr2=1/2 #a/b ratio at which to draw a line for vis
bcr2=1/2 #b/c ratio at which to draw a line for vis
abr1=1/1.5 #a/b ratio at which to draw a line for vis
bcr1=1/1.5 #b/c ratio at which to draw a line for vis

mscale= 20/9;
mass= 1;

Nss= 0;
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
es

if sm_min <= length(es)/2 #has to be > length(es)/2
  exit(1)
endif

[fid, msg] = fopen ("I-fit.txt", "w");
fprintf(fid, \
        "#ell_a\tell_b\tell_c\tell_x\tell_y\tell_z\ta_x\ta_y\ta_z\tb_x\tb_y\tb_z\tc_x\tc_y\tc_z\tell_t\tindex\n");

for n=1:1:N;

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
#    l
    ax
    vol
    printf("Volum <= 0! Aborting\n")
    fflush(stdout);
    #exit
    return #break is not enough any more !!!!!
  endif

  [axs, axi]= sort (ax); #making a < b < c if axis-names are not specially assigned!
  #eu= euler_angles(v(:,axi(1)), v(:,axi(2)), v(:,axi(3))); #index ordered v

  is_sp= 0;
  is_sm= 0;
  is_ci= 0;
  esum= 0;
  
  eaxs= ev*norm(axs); #point within individual sphere limit?

  #this error evaluation can be done with fuzzy logic!
  #this makes sence if the error criterion is very hard and often only a single exception causes uncertainty
  #with fuzzy logic a minimum of certain evaluations have to be met
  #(fuzzy logic with reals is possible but very difficult: the ratio of the two volumes created by the cut of the separation surface through the error box)
  for i=1:1:length(es) 
    for j=1:1:size(es,2)
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

## error box containing sphere point?
    if (es(i,:)==[-1, -1, -1]) 
      for jj=1:1:size(es,2)
        is_sp+= (ee(jj) < eaxs(jj));
      endfor
    else if (es(i,:)==[1, 1, 1])
        for jj=1:1:size(es,2)
          is_sp+= (ee(jj) > eaxs(jj));
        endfor
      endif
    endif
   
    is_sm+= ((ee(1) / ee(2)) < (ee(2) / ee(3)));
    is_ci+= ((ee(1) / ee(2)) > (ee(2) / ee(3)));
    ## if !quiet
    ##   printf("t(n,1): %f; n: %d; i: %d; a/b: %f; b/c: %f; sm: %d; ci: %d; sp: %d\n", t(n,1),n,i,(ee(1) / ee(2)), (ee(2) / ee(3)),is_sm,is_ci,is_sp)
    ## endif
  endfor
  mee= esum / size(es,1) / size(es,2);
  tsum= tsum + mee;

  #is_sp= norm(axs/norm(axs) - ev/norm(ev)) < spe; #point within const. sphere limit?

  #if ==-case can be neglected this can be shortened:
  #is_ci= length(es)-is_sm
  #because of this is_sm > sm_min && is_ci > ci_min can never happen

  if is_sp >= sp_min
    Nss++;
    ce(:,n)= [0,0,1];
    et= 0;
  else if is_sm >= sm_min
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
        ce(:,n)= [1,1,0];
        et= 2;
      endif
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

  u_axes(:,n)= [axs(1), axs(2), axs(3)];

  [theta, phi, r]= cart2sph(axs(1), axs(2), axs(3)); 
  dp3ds(:,n)= [theta, phi, r];


  #fprintf(fid, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", ax, t(n,7:9), v);
  fprintf(fid, \
          "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\n", \
           axs, p_pos, v(:,axi(1)), v(:,axi(2)), v(:,axi(3)), et, p_index);
end;

fclose(fid);
printf("# sphere-like: %d; # smarty-like: %d; # cigar-like: %d; # uncertain: %d; ratio: %.2f\n", Nss,Ns, Nz, Nsz, Ns/Nz);
Nss+Ns+Nz+Nsz
tsum / N

#break

ws =acos(dot([1,1,1], [1,1,0])/(norm([1,1,1]) * norm([1,1,0])));
#a0=[linspace(0,pi/2,num); zeros(1, num)];# 1/4 circle in ab-plane
#a0=[pi/4, pi/2]';#[1,0,0] point?
#b0=[ones(1, num) * pi/2; linspace(0,pi/2,num)];#arc from [0, 1, 0] to a0
b0s=[ones(1, num) * pi/2; linspace(pi/4,pi/2,num)];#arc from b0 to a0 
m0s=[ones(1, num) * pi/2; linspace(0,pi/4,num)];# rest of 1/4 circle in bc-plane 
#b0=[pi/2, pi/4]';#[0, 1/sqrt(2), 1/sqrt(2)] point
k0s=[zeros(1, num); linspace(0,pi/2,num)];# 1/4 circle in ac-plane 
l0s=[linspace(0,pi/2,num); zeros(1, num)];# 1/4 circle in ab-plane
c00=ev'; #def. at beginning for sphere eval.
#c00= [1/sqrt(3); 1/sqrt(3); 1/sqrt(3)]; 
c00s=[pi/4; ws]';#[1/sqrt(3), 1/sqrt(3), 1/sqrt(3)] point
c0s=[ones(1, num) * pi/4; linspace(pi/2,ws,num)];#arc from w to a0
n0s=[ones(1, num) * pi/4; linspace(ws,0,num)];#rest of arc from ab-plane to a0
cstart= 1/sqrt(2+abr3^2);
vstart= [abr3*cstart, cstart, cstart];
phistart= pi/2-acos(dot(vstart, [0,0,1])/(norm(vstart) * norm([0,0,1])));
theta= acot(abr3);
u0s=[ones(1, num) * theta; linspace(pi/2,phistart,num)];#arc for a/b==.5
cstart= 1/sqrt(2+abr2^2);
vstart= [abr2*cstart, cstart, cstart];
phistart= pi/2-acos(dot(vstart, [0,0,1])/(norm(vstart) * norm([0,0,1])));
theta= acot(abr2);
o0s=[ones(1, num) * theta; linspace(pi/2,phistart,num)];#arc for a/b==.5
cstart= 1/sqrt(2+abr1^2);
vstart= [abr1*cstart, cstart, cstart];
phistart= pi/2-acos(dot(vstart, [0,0,1])/(norm(vstart) * norm([0,0,1])));
theta= acot(abr1);
q0s=[ones(1, num) * theta; linspace(pi/2,phistart,num)];#arc for a/b==.5
#c1=[ones(1, num) * pi/4; linspace(0,pi/2,num)]; #arc from [1/sqrt(2), 1/sqrt(2), 0] to a0; only good for 2D view!
#a1=[linspace(0,pi/2,num); ones(1, num) * w];#arc from [1/sqrt(3), 0, sqrt(2/3)] to  [0, 1/sqrt(3), sqrt(2/3)]

####creating a point-mesh for spherical coords
pbin= 31;
tbin= 2*pbin-1;
phis=[linspace(-pi/2,pi/2,pbin)];# circle in ac-plane 
thes=[linspace(-pi,pi,tbin)];# circle in ab-plane

rphis= repmat(phis, tbin, 1);
rthes= repmat(thes, pbin, 1);
#bm=cat (3,rthes, rphis'); #tbin x pbin x 2 matrix!
bl= [reshape(rthes,1,[])', reshape(rphis',1,[])']'; #point-coords of bm as vertcat
####done

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
#xs= linspace(1,0,num);# b== 1/sqrt(3) line
for n=1:1:num
  #s1x= xs;
  #s1y= 1/sqrt(3);
  s1z= sqrt(2/3 - xs(n)^2);
  s1(:,n)= [xs(n), 1/sqrt(3), s1z];
end

#### guide points for b== c
clear xs
xs= linspace(1/sqrt(3),0,num);#a variable
for n=1:1:num
  a0y= sqrt((1-xs(n)^2)/2); #b==c
  a0(:,n)= [xs(n), a0y, a0y];
end

#### guide points for b== c: rest of 1/4 arc
clear xs
#xs= linspace(1/sqrt(3),0,num);
xs= linspace(1,1/sqrt(3),num);
for n=1:1:num
  ja0y= sqrt((1-xs(n)^2)/2);
  ja0(:,n)= [xs(n), ja0y, ja0y];
end

#### guide points for b/c== bcr
clear xs
astart= 1/sqrt(2+1/bcr2^2);
xs= linspace(astart,0,num);#a variable
for n=1:1:num
  p0z= sqrt((1-xs(n)^2)/(bcr2^2+1));
  p0(:,n)= [xs(n), bcr2*p0z, p0z];
end

#### guide points for b/c== bcr
clear xs
astart= 1/sqrt(2+1/bcr1^2);
xs= linspace(astart,0,num);#a variable
for n=1:1:num
  r0z= sqrt((1-xs(n)^2)/(bcr1^2+1));
  r0(:,n)= [xs(n), bcr1*r0z, r0z];
end

#### guide points for b/c== bcr
clear xs
astart= 1/sqrt(2+1/bcr3^2);
xs= linspace(astart,0,num);#a variable
for n=1:1:num
  v0z= sqrt((1-xs(n)^2)/(bcr3^2+1));
  v0(:,n)= [xs(n), bcr3*v0z, v0z];
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
[xt, yt, zt]= sph2cart (k0s(1,:), k0s(2,:), ones(1,size(k0s,2)));#projection of 3D points onto unit sphere
k0= vertcat (xt, yt, zt);

clear xt yt zt;
[xt, yt, zt]= sph2cart (l0s(1,:), l0s(2,:), ones(1,size(l0s,2)));#projection of 3D points onto unit sphere
l0= vertcat (xt, yt, zt);

clear xt yt zt;
[xt, yt, zt]= sph2cart (m0s(1,:), m0s(2,:), ones(1,size(m0s,2)));#projection of 3D points onto unit sphere
m0= vertcat (xt, yt, zt);

clear xt yt zt;
[xt, yt, zt]= sph2cart (n0s(1,:), n0s(2,:), ones(1,size(n0s,2)));#projection of 3D points onto unit sphere
n0= vertcat (xt, yt, zt);

clear xt yt zt;
[xt, yt, zt]= sph2cart (o0s(1,:), o0s(2,:), ones(1,size(o0s,2)));#projection of 3D points onto unit sphere
o0= vertcat (xt, yt, zt);

clear xt yt zt;
[xt, yt, zt]= sph2cart (q0s(1,:), q0s(2,:), ones(1,size(q0s,2)));#projection of 3D points onto unit sphere
q0= vertcat (xt, yt, zt);

clear xt yt zt;
[xt, yt, zt]= sph2cart (u0s(1,:), u0s(2,:), ones(1,size(u0s,2)));#projection of 3D points onto unit sphere
u0= vertcat (xt, yt, zt);

clear xt yt zt;
[xt, yt, zt]= sph2cart (dp3ds(1,:), dp3ds(2,:), ones(1,size(dp3ds,2)));#projection of 3D data points onto unit sphere
dp3d= vertcat (xt, yt, zt);

clear xt yt zt;
[xt, yt, zt]= sph2cart (bl(1,:), bl(2,:), ones(1,size(bl,2)));#projection of 3D points onto unit sphere
bl0= vertcat (xt, yt, zt);


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

clear xt yt zt;
[xt, yt, zt]= cart2sph (p0(1,:), p0(2,:), p0(3,:));#projection of 3D points onto unit sphere
p0s= vertcat (xt, yt, zt);

clear xt yt zt;
[xt, yt, zt]= cart2sph (r0(1,:), r0(2,:), r0(3,:));#projection of 3D points onto unit sphere
r0s= vertcat (xt, yt, zt);

clear xt yt zt;
[xt, yt, zt]= cart2sph (v0(1,:), v0(2,:), v0(3,:));#projection of 3D points onto unit sphere
v0s= vertcat (xt, yt, zt);

#stereographic projection on unit sphere
c00p= stereogproj(c00s(1), c00s(2), 1, phi0, lambda0);
c0p= stereogproj(c0s(1,:), c0s(2,:), 1, phi0, lambda0);
b0p= stereogproj(b0s(1,:), b0s(2,:), 1, phi0, lambda0);
a0p= stereogproj(a0s(1,:), a0s(2,:), 1, phi0, lambda0);
s0p= stereogproj(s0s(1,:), s0s(2,:), 1, phi0, lambda0);
o0p= stereogproj(o0s(1,:), o0s(2,:), 1, phi0, lambda0);
p0p= stereogproj(p0s(1,:), p0s(2,:), 1, phi0, lambda0);
q0p= stereogproj(q0s(1,:), q0s(2,:), 1, phi0, lambda0);
r0p= stereogproj(r0s(1,:), r0s(2,:), 1, phi0, lambda0);
u0p= stereogproj(u0s(1,:), u0s(2,:), 1, phi0, lambda0);
v0p= stereogproj(v0s(1,:), v0s(2,:), 1, phi0, lambda0);
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
#scatter3 (dp3d(1,:),dp3d(2,:),dp3d(3,:), 50, ce, 's'); #'d'
scatter3 (dp3d(1,:),dp3d(2,:),dp3d(3,:), ps3d, ce', 's', 'filled')#for octave 3.6.2
hold on
plot3 (a0(1,:), a0(2,:), a0(3,:), "k")
plot3 (b0(1,:), b0(2,:), b0(3,:), "k")
plot3 (c0(1,:), c0(2,:), c0(3,:), "k")
plot3 (s0(1,:), s0(2,:), s0(3,:), "k")
#plot3 (s1(1,:), s1(2,:), s1(3,:), "k") #arc for b== ws
plot3 (k0(1,:), k0(2,:), k0(3,:), "k")
plot3 (l0(1,:), l0(2,:), l0(3,:), "k")
plot3 (m0(1,:), m0(2,:), m0(3,:), "k")
plot3 (n0(1,:), n0(2,:), n0(3,:), "k") 
plot3 (o0(1,:), o0(2,:), o0(3,:), "k") #a/b==abr2
plot3 (p0(1,:), p0(2,:), p0(3,:), "k") #b/c==bcr2
plot3 (q0(1,:), q0(2,:), q0(3,:), "k") #a/b==abr1
plot3 (r0(1,:), r0(2,:), r0(3,:), "k") #b/c==bcr1
plot3 (u0(1,:), u0(2,:), u0(3,:), "k") #a/b==abr3
plot3 (v0(1,:), v0(2,:), v0(3,:), "k") #b/c==bcr3
plot3 (ja0(1,:), ja0(2,:), ja0(3,:), "k")
hold off

#set (gca (), "plotboxaspectratio", [1 1 1.45])#empirical ratio

#axis ("square");
axis ([0,1,0,1,0,1],"square");
set (gca (), "plotboxaspectratio", [1 1 1.45])#empirical ratio, needs to be set after axis ([0,1,0,1,0,1],"square");!!!

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
#text (a0(1,floor(num/3*2)), a0(2,floor(num/3*2)) + .05, a0(3,floor(num/3*2)), "oblate arc", "rotation", 30);
text (b0(1,floor(num/2)) - .05, b0(2,floor(num/2)), b0(3,floor(num/2)), "ellipse arc", "rotation", -50);
#text (b0(1,floor(num/4*3)) - .05, b0(2,floor(num/4*3)), b0(3,floor(num/4*3)), "ellipse arc", "rotation", -50);
text (c0(1,floor(num/4*3)) + .05, c0(2,floor(num/4*3)), c0(3,floor(num/4*3)), "prolate arc", "rotation", 90);
#text (c0(1,floor(num/3)) + .05, c0(2,floor(num/3)), c0(3,floor(num/3)), "prolate arc", "rotation", 90);
text (s0(1,size(s0,2)-30) - .03, s0(2,size(s0,2)-30), s0(3,size(s0,2)-30), "separation curve", "rotation", -75);


#break
####printing now...

#paper_size = [640, 480];
#set (gcf, "paperunits", "inches")
#set (gcf, "papertype", "<custom>")
#set (gcf, "papersize", paper_size)
#set (gcf, "paperposition", [0, 0, paper_size])

#figure('Position', [0, 0, 600, 400]); 


####printing now...

nplot= nplot + 1;
if !quiet
  printf("Printing plot # %d", nplot)
endif
print(sprintf("%s_%.2d.png", out3D, nplot), '-dpng', '-S800,800');#, '-F/usr/X11R6/lib/X11/fonts/msttf/arial.ttf');#, '-r100');
print(sprintf("%s_%.2d.svg", out3D, nplot), '-dsvg', '-S800,800');#has to be there for axis ("square") to work even with svg (-S not possible any more with gnuplot > 4.3.0 ???)
if !quiet
  printf(" done.\n", nplot)
endif

####printing end

#return

view(110, 10);

####printing now...

nplot= nplot + 1;
if !quiet
  printf("Printing plot # %d", nplot)
endif
print(sprintf("%s_%.2d.png", out3D, nplot), '-dpng', '-S800,800');#, '-F/usr/X11R6/lib/X11/fonts/msttf/arial.ttf');#, '-r100');
print(sprintf("%s_%.2d.svg", out3D, nplot), '-dsvg', '-S800,800');#has to be there for axis ("square") to work even with svg (-S not possible any more with gnuplot > 4.3.0 ???)
if !quiet
  printf(" done.\n", nplot)
endif

####printing end


# box("off");



# ####printing now...

# nplot= nplot + 1;
# if !quiet
#   printf("Printing plot # %d", nplot)
# endif
# print(sprintf("%s_%.2d.png", out3D, nplot), '-dpng', '-S800,800');#, '-F/usr/X11R6/lib/X11/fonts/msttf/arial.ttf');#, '-r100');
# print(sprintf("%s_%.2d.svg", out3D, nplot), '-dsvg', '-S800,800');#has to be there for axis ("square") to work even with svg (-S not possible any more with gnuplot > 4.3.0 ???)
# if !quiet
#   printf(" done.\n", nplot)
# endif

# ####printing end













#######################global axes orientations BEGIN
##with local sphere sampling

clear phi0 the0
clear xbin ybin vXEdge vYEdge mHist2d nXBins nYBins vXLabel vYLabel

l00= pi/2; #projection centre
#bin= round(pi*10);#make sure binning fits scatter plot impression!
bin= 31


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

  uao(n,:)= squeeze(u(n,1,:));
  ubo(n,:)= squeeze(u(n,2,:));
  uco(n,:)= squeeze(u(n,3,:));

  ##ellipsoids are symmetry along their axes so mirror each point on opposite side of the unit sphere;-), see Diss for explanation
  uam(n,:)= -uao(n,:);
  ubm(n,:)= -ubo(n,:);
  ucm(n,:)= -uco(n,:);
  
  ua(n,:)= horzcat(uao(n,:), uam(n,:));
  ub(n,:)= horzcat(ubo(n,:), ubm(n,:));
  uc(n,:)= horzcat(uco(n,:), ucm(n,:));

  n_1=mod(n-1,3)+1;
  n_2=mod(n+0,3)+1;
  n_3=mod(n+1,3)+1;

  #uW(n,:)= [[u_axes(n_1,:)./sqrt(u_axes(n_2,:).*u_axes(n_3,:))]';[u_axes(n_1,:)./sqrt(u_axes(n_2,:).*u_axes(n_3,:))]'];#for weighted density

  #uW= ones(size(ua(n,:),2),1)#for number density
  uW= ones(3, size(ua(n,:),2));#for number density

  #uW= uW./sum(uW);
  uW(n,:)= uW(n,:)./size(uao(n,:),2);

  u_lss(n,:)= local_sphere_sampling([ua(n,:);ub(n,:);uc(n,:)]', uW(n,:), 10, 1);#'
  Nt(n)= max(u_lss(n,:))*size(uao(n,:),2); #save max value for unified colour gradient of all 3 plots

  max(u_lss(n,:))
endfor
#break

#Nmax= max(Nt);
Nmax= 0.3*size(uao(n,:),2);
#midcmap= jet(Nmax + 0.2 * Nmax);#dark red region is about 1/5th of jet length
midcmap= jet(1.15 * Nmax + 1);#just to be on the save side that midcmap is larger than Nt(n)+1

##some octave magick to delete dark red colours;-)
del_p = midcmap(:,1) < 1 & midcmap(:,2) == 0  & midcmap(:,3) == 0;
size(midcmap)
midcmap(del_p,:)= [];
size(midcmap)
##magick done

#cmapg= midcmap;
Nc= 10000;
cmapg= vertcat(midcmap, ones(Nc - Nmax,1) * midcmap(end,:));#extend with last colour of midcmap for values mapped above Nmax
#cmapg= jet();
#cmapg(1,:)=[1,1,1];

for n=1:1:3

  N= size(ua(n,:),2);

  printf("Doing %d run...\n", n)
  #cmap= cmapg(1:Nt(n)+1,:);
  cmap=jet(128); #artificial, to gen. qualtitative images

  c_i= round(u_lss(n,:) ./ max(u_lss(n,:)) .* (length(cmap) - 1) + 1); #for abs values: hight has to be normalized for this
  #[c_r, c_g, c_b]= ind2rgb(gray2ind(wM), cmapg);
  m_c= zeros(N,3);
  for nnn=1:1:N;
    p_c= cmap(c_i(nnn),:);
    m_c(nnn,:)= p_c;
    #if (p_c == [1,1,1]) #remove white for the mesh colouring
      #p_c= [0,0,0] #black
    #  p_c= cmapg(2,:); #colour after white mapped twice!
    #endif
  endfor

  figure (1)
  clf ()

  posu= ua(n,:) > 0;
  negu= !posu;
  posb= bl0(1,:) > 0;
  negb= !posb;

  #scatter3(ua(n,:), ub(n,:), uc(n,:), ps3d, m_c , 's', 'filled');
  #hold on
  #scatter3(bl0(1,:), bl0(2,:), bl0(3,:), 1, 'black', 's', 'filled')
  scatter3(ua(n,negu), ub(n,negu), uc(n,negu), ps3d, m_c(negu,:) , 's', 'filled');
  hold on
  scatter3(bl0(1,negb), bl0(2,negb), bl0(3,negb), 1, 'black', 's', 'filled')
  hold off

  axis ([-1,1,-1,1,-1,1],"square"); #view(0, 90);
  #axis ([-1,1,-1,0,-1,1],"square");#since only hemisphere matters ##doesn't look nice since square seems not implemented for 3D:-(
  set (gca (), "plotboxaspectratio", [1 1 1.45])#empirical ratio, needs to be set after axis ([0,1,0,1,0,1],"square");!!!

  #azimuth= 135;
  #elevation= acosd(dot([1,1,1], [1,1,0])/(norm([1,1,1]) * norm([1,1,0])));
  #view(azimuth, elevation);
  view(110, 10);

  xlabel("x");
  ylabel("y");
  zlabel("z");


####printing now...

nplot= nplot + 1;
if !quiet
  printf("Printing plot # %d", nplot)
endif

filename = sprintf ("%s_%.2d_pbaspect=[%f,%f,%f].pdf", outGO, nplot, get (gca (), "plotboxaspectratio"));
filename(filename=="0") = [];
print ("-dpdfwrite", filename)
filename = sprintf ("%s_%.2d_pbaspect=[%f,%f,%f].png", outGO, nplot, get (gca (), "plotboxaspectratio"));
filename(filename=="0") = [];
print ("-dpng", filename)
filename = sprintf ("%s_%.2d_pbaspect=[%f,%f,%f].svg", outGO, nplot, get (gca (), "plotboxaspectratio"));
filename(filename=="0") = [];
print ("-dsvg", filename)

print(sprintf("%s_%.2d_%d.png", outGO, nplot, n), '-dpng', '-S800,800');#, '-F/usr/X11R6/lib/X11/fonts/msttf/arial.ttf');#, '-r100');
print(sprintf("%s_%.2d_%d.svg", outGO, nplot, n), '-dsvg', '-S800,800');#has to be there for axis ("square") to work even with svg (-S not possible any more with gnuplot > 4.3.0 ???)
if !quiet
  printf(" done.\n", nplot)
endif

####printing end


  figure (1)
  clf ()


  scatter3(ua(n,posu), ub(n,posu), uc(n,posu), ps3d, m_c(posu,:) , 's', 'filled');
  hold on
  scatter3(bl0(1,posb), bl0(2,posb), bl0(3,posb), 1, 'black', 's', 'filled')
  hold off

  axis ([-1,1,-1,1,-1,1],"square"); #view(0, 90);
  #axis ([-1,1,-1,0,-1,1],"square");#since only hemisphere matters ##doesn't look nice since square seems not implemented for 3D:-(
  set (gca (), "plotboxaspectratio", [1 1 1.45])#empirical ratio, needs to be set after axis ([0,1,0,1,0,1],"square");!!!

  #azimuth= 135;
  #elevation= acosd(dot([1,1,1], [1,1,0])/(norm([1,1,1]) * norm([1,1,0])));
  #view(azimuth, elevation);
  view(110, 10);

  xlabel("x");
  ylabel("y");
  zlabel("z");


####printing now...

nplot= nplot + 1;
if !quiet
  printf("Printing plot # %d", nplot)
endif

filename = sprintf ("%s_%.2d_pbaspect=[%f,%f,%f].pdf", outGO, nplot, get (gca (), "plotboxaspectratio"));
filename(filename=="0") = [];
print ("-dpdfwrite", filename)
filename = sprintf ("%s_%.2d_pbaspect=[%f,%f,%f].png", outGO, nplot, get (gca (), "plotboxaspectratio"));
filename(filename=="0") = [];
print ("-dpng", filename)
filename = sprintf ("%s_%.2d_pbaspect=[%f,%f,%f].svg", outGO, nplot, get (gca (), "plotboxaspectratio"));
filename(filename=="0") = [];
print ("-dsvg", filename)

print(sprintf("%s_%.2d_%d.png", outGO, nplot, n), '-dpng', '-S800,800');#, '-F/usr/X11R6/lib/X11/fonts/msttf/arial.ttf');#, '-r100');
print(sprintf("%s_%.2d_%d.svg", outGO, nplot, n), '-dsvg', '-S800,800');#has to be there for axis ("square") to work even with svg (-S not possible any more with gnuplot > 4.3.0 ???)
if !quiet
  printf(" done.\n", nplot)
endif

####printing end



  #if(n==3)
  #   return
  #endif



  #return

  clear theta phi r x y
  [theta, phi, r]= cart2sph(ua(n,:), ub(n,:), uc(n,:)); 
  theta=theta.*180./pi;
  phi= phi.*180./pi;

  scatter3(theta, phi, u_lss(n,:), ps3d, m_c , 's', 'filled');
  #axis ([-pi,pi,-pi/2,pi/2],"square");
  #axis ([-pi/2,pi/2,-pi/2,pi/2],"square");#since only hemisphere matters
  axis ([-90,90,-90,90],"square");#since only hemisphere matters
  set (gca, 'xtick', [-90:10:90]);
  set (gca, 'ytick', [-90:10:90]);


####printing now...

nplot= nplot + 1;
if !quiet
  printf("Printing plot # %d", nplot)
endif
print(sprintf("%s_%.2d_%d.png", outGO, nplot, n), '-dpng', '-S800,800');#, '-F/usr/X11R6/lib/X11/fonts/msttf/arial.ttf');#, '-r100');
print(sprintf("%s_%.2d_%d.svg", outGO, nplot, n), '-dsvg', '-S800,800');#has to be there for axis ("square") to work even with svg (-S not possible any more with gnuplot > 4.3.0 ???)
if !quiet
  printf(" done.\n", nplot)
endif

####printing end




  scatter(theta, phi, ps2d, m_c , 's', 'filled');
  #text (-1.5, -2, sprintf("# cells: %d", N));
  text (-90, -100, sprintf("# cells: %d", N));

  cmap= cmapg(1:Nt(n)+1,:);
  #cmap(1,:)=[1,1,1];#only for field depiction
  colormap(cmap);
  colorbar #show colorbar

  #axis ([-pi,pi,-pi/2,pi/2],"square");
  #axis ([-pi/2,pi/2,-pi/2,pi/2],"square");#since only hemisphere matters
  axis ([-90,90,-90,90],"square");#since only hemisphere matters
  set (gca, 'xtick', [-90:10:90]);
  set (gca, 'ytick', [-90:10:90]);


####printing now...

nplot= nplot + 1;
if !quiet
  printf("Printing plot # %d", nplot)
endif
print(sprintf("%s_%.2d_%d.png", outGO, nplot, n), '-dpng', '-S800,800');#, '-F/usr/X11R6/lib/X11/fonts/msttf/arial.ttf');#, '-r100');
print(sprintf("%s_%.2d_%d.svg", outGO, nplot, n), '-dsvg', '-S800,800');#has to be there for axis ("square") to work even with svg (-S not possible any more with gnuplot > 4.3.0 ???)
if !quiet
  printf(" done.\n", nplot)
endif

####printing end



###2D-hist now

# lef= min(theta);
# rig= max(theta);
# bot= min(phi);
# top= max(phi);
#lef= -pi;
#rig=  pi;
#lef= -pi/2; #since only hemisphere matters
#rig=  pi/2; #since only hemisphere matters
#bot= -pi/2;
#top=  pi/2;
lef= -90; #since only hemisphere matters
rig=  90; #since only hemisphere matters
bot= -90;
top=  90;
#lef= min(l270(1,:));
#rig= max( l90(1,:));
#lef= min(l360(1,:));
#rig= max(  l0(1,:));
#lef= -repr;
#rig=  repr;
#bot= min(l270(2,:));
#top= max(l270(2,:));
dx= (rig - lef);
dy= (top - bot);

xbin=bin;
ybin=round( dy / dx * xbin); #make'm squares
#ybin= round( dy / dx * xbin) + 1; #make'm odd!
#xbin= round(bin / dx * dxr); #make the squares the same size as before
#ybin= round(dyr / dxr * xbin); #make'm squares

#x= vertcat(x, l0(1,:)');#for testing guide line pos
#y= vertcat(y, l0(2,:)');



vXEdge = linspace(lef, rig+(rig-lef)/(xbin-1), xbin+1);#add one entry at end
vYEdge = linspace(bot, top+(top-bot)/(ybin-1), ybin+1);
#[mHist2d, wM] = hist2d_wm([y, x],vYEdge,vXEdge,rel_c_area); #2D-hist without guide points
#mHist2d = avg_hist2d_w([phi, theta],vYEdge,vXEdge,u_lss(n,:));###why does this not work any more???
mHist2d = avg_hist2d_w([phi; theta]',vYEdge,vXEdge,u_lss(n,:));

#Hmax= max(max(mHist2d))
#Hmin= min(min(mHist2d)) #== 0 because areas are positve!
#Hsum= sum(sum(mHist2d)) #!= 1 becauese it's not normalized any more!


nXBins = length(vXEdge);
nYBins = length(vYEdge);
  #vXLabel = 0.5*(vXEdge(1:(nXBins-1))+vXEdge(2:nXBins));
  #vYLabel = 0.5*(vYEdge(1:(nYBins-1))+vYEdge(2:nYBins));
vXLabel = vXEdge(1:(nXBins-1));
vYLabel = vYEdge(1:(nYBins-1));
  
#set (gca, 'xtick', "");#the ticks aren't correct!
#set (gca, 'ytick', "");
#set (gca, 'xtick', []);#3.6.2 #the ticks aren't correct!
#set (gca, 'ytick', []);

pcolor(vXLabel, vYLabel, mHist2d); #mHist2D acts as color value
#imagesc(mHist2d); #no x and y scale! would need meshgrid first!
#surf(vXLabel, vYLabel, mHist2d); #3D but no bars
#surfc(vXLabel, vYLabel, mHist2d); #surf + contours
#mesh(vXLabel, vYLabel, mHist2d); #similar to surf but wireframe
#plot3(vXLabel, vYLabel, mHist2d); #plot points connected


shading flat; #means no border around each hist rectangle
  #shading faceted

#text (-3, -1.5, sprintf("# cells: %d", N));

#caxis([Nmin,Nmax]); #color between icosaedron and tetraedron
#caxis([0,1]); #color between icosaedron and tetraedron
#caxis('auto')
caxis([0,0.1]) #fixed limits from dra for nw

  #cmap= cmapg(1:Nt(n)+1,:);
  cmap= midcmap; #for fixed limits from dra for nw 
  cmap(1,:)=[1,1,1];
  colormap(cmap);
  colorbar #show colorbar
  #axis ("square");#setting axis range here can be bad!
  axis ("image");#square and tight ;-)
  #set (gca, 'xtick', [-90:30:90]);
  #set (gca, 'ytick', [-90:30:90]);
  ##no range can be specified for minor tick (not as in gnuplot mxticks)
  #set (gca, 'xminortick', [-90:10:90]);
  #set (gca, 'yminortick', [-90:10:90]);
  #set (gca, "xminortick", "on", "yminortick", "on")
  ##using main ticks and remove them in inkscape...
  set (gca, 'xtick', [-90:10:90]);
  set (gca, 'ytick', [-90:10:90]);

# printf("Before break")
# break
# return
# printf("After break")

####printing now...

nplot= nplot + 1;
if !quiet
  printf("Printing plot # %d", nplot)
endif
print(sprintf("%s_%.2d_%d.png", outGO, nplot, n), '-dpng', '-S800,800');#, '-F/usr/X11R6/lib/X11/fonts/msttf/arial.ttf');#, '-r100');
print(sprintf("%s_%.2d_%d.svg", outGO, nplot, n), '-dsvg', '-S800,800');#has to be there for axis ("square") to work even with svg (-S not possible any more with gnuplot > 4.3.0 ???)
if !quiet
  printf(" done.\n", nplot)
endif

###pure image for reprojection with e.g. G.projector;-)
imwrite (uint8(flipud(mHist2d(1:end-1,1:end-1))*length(cmap)/0.1), cmap, sprintf("%s_%.2d_%d.png", outGO, nplot, n)); # removing last row (empty) and last column (point mirror of first column) https://www.gnu.org/software/octave/doc/v4.0.1/Index-Expressions.html # upper cmap fixed to 0.1 (as caxis([0,0.1]) above), still needs extending cmap to all values to avoid "filling with black"-warning # indexed image needs to be uint indexing a color [0;1] of the cmap: https://www.gnu.org/software/octave/doc/v4.0.1/Representing-Images.html  # saveimage replaced by imwrite: https://www.gnu.org/software/octave/doc/v4.0.1/Loading-and-Saving-Images.html https://www.gnu.org/software/octave/doc/v4.0.1/Obsolete-Functions.html  # if image type is not known use im2uint8: http://octave.sourceforge.net/image/function/im2uint8.html  http://stackoverflow.com/questions/30188519/octave-how-to-load-grayscale-images-in-double-format#30189527

####printing end


endfor

#######################global axes orientations END


##### BEGIN 2D point cloud with RGB, Y point colouring

scatter (dp2d(1,:), dp2d(2,:), ps2d, ce', 's', 'filled')#3.6.2
hold on
plot (a0p(1,:), a0p(2,:), "k")
plot (b0p(1,:), b0p(2,:), "k")
plot (c0p(1,:), c0p(2,:), "k")
plot (s0p(1,:), s0p(2,:), "k")
#plot3 (s1(1,:), s1(2,:), "k")
plot (o0p(1,:), o0p(2,:), "k") #a/b==abr2
plot (p0p(1,:), p0p(2,:), "k") #b/c==abr2
plot (q0p(1,:), q0p(2,:), "k") #a/b==abr1
plot (r0p(1,:), r0p(2,:), "k") #b/c==abr1
plot (u0p(1,:), u0p(2,:), "k") #a/b==abr3
plot (v0p(1,:), v0p(2,:), "k") #b/c==abr3
hold off
axis ([c00p(1,1), b0p(1,1), c00p(2,1), c0p(2,1), ], "equal");

text (c00p(1,1) - .02, c00p(2,1) + .02, "sphere\npoint", "horizontalalignment", "right"); #looks nicer
text (b0p(1,1) + .02, b0p(2,1), "circle\npoint");
text (c0p(1,1) - .02, c0p(2,1), "line\npoint", "horizontalalignment", "right");
text (a0p(1,floor(num/4)) + .02, a0p(2,floor(num/4)), "oblate line", "rotation", 30);
text (b0p(1,floor(num/2)) + .02, b0p(2,floor(num/2)), "ellipse arc", "rotation", -50);
text (c0p(1,floor(num/4*3)) - .02, c0p(2,floor(num/4*3)), "prolate line", "rotation", 90);
text (s0p(1,size(s0p,2)-20) + .02, s0p(2,size(s0p,2)-20), "separation curve", "rotation", -75);
text (c00p(1,1), c00p(2,1) - .1, sprintf("# oblate-like: %d; # prolate-like: %d; ratio: %.2f\n# spere-like: %d; # uncertain: %d", Nz, Ns, Ns/Nz,Ns,Nsz));

#xlabel("");
#ylabel("");
#set (gca, 'xtick', "");#the ticks aren't correct!
#set (gca, 'ytick', "");
set (gca, 'xtick', []);#3.6.2 #the ticks aren't correct!
set (gca, 'ytick', []);


####printing now...

nplot= nplot + 1;
if !quiet
  printf("Printing plot # %d", nplot)
endif
print(sprintf("%s_%.2d.png", out3D, nplot), '-dpng', '-S800,800');#, '-F/usr/X11R6/lib/X11/fonts/msttf/arial.ttf');#, '-r100');
print(sprintf("%s_%.2d.svg", out3D, nplot), '-dsvg');#has to be there for axis ("square") to work even with svg (-S not possible any more with gnuplot > 4.3.0 ???)
if !quiet
  printf(" done.\n", nplot)
endif

####printing end

##### END 2D point cloud with RGB, Y point colouring


##### BEGIN 2D-hist of point cloud with jet cmap

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

#set (gca, 'xtick', "");#the ticks aren't correct!
#set (gca, 'ytick', "");
set (gca, 'xtick', []);#3.6.2 #the ticks aren't correct!
set (gca, 'ytick', []);

pcolor(vXLabel, vYLabel, mHist2d); #mHist2D acts as color value
hold on
plot (a0p(1,:), a0p(2,:), "k")
plot (b0p(1,:), b0p(2,:), "k")
plot (c0p(1,:), c0p(2,:), "k")
plot (s0p(1,:), s0p(2,:), "k")
#plot3 (s1(1,:), s1(2,:), "k")
plot (o0p(1,:), o0p(2,:), "k") #a/b==abr2
plot (p0p(1,:), p0p(2,:), "k") #b/c==abr2
plot (q0p(1,:), q0p(2,:), "k") #a/b==abr1
plot (r0p(1,:), r0p(2,:), "k") #b/c==abr1
plot (u0p(1,:), u0p(2,:), "k") #a/b==abr3
plot (v0p(1,:), v0p(2,:), "k") #b/c==abr3
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
#set (gca, 'xtick', "");#the ticks aren't correct!
#set (gca, 'ytick', "");
set (gca, 'xtick', []);#3.6.2 #the ticks aren't correct!
set (gca, 'ytick', []);


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
axis ("equal");#setting axis range here can be bad!



####printing now...

nplot= nplot + 1;
if !quiet
  printf("Printing plot # %d", nplot)
endif
print(sprintf("%s_%.2d.png", out3D, nplot), '-dpng', '-S800,800');#, '-F/usr/X11R6/lib/X11/fonts/msttf/arial.ttf');#, '-r100');
print(sprintf("%s_%.2d.svg", out3D, nplot), '-dsvg');#has to be there for axis ("square") to work even with svg (-S not possible any more with gnuplot > 4.3.0 ???)
if !quiet
  printf(" done.\n", nplot)
endif

####printing end

##### END 2D-hist of point cloud with jet cmap


##### BEGIN shearY(30°) 2D-hist of point cloud with jet cmap

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
o0r= shear_y(o0p, m);
p0r= shear_y(p0p, m);
q0r= shear_y(q0p, m);
r0r= shear_y(r0p, m);
u0r= shear_y(u0p, m);
v0r= shear_y(v0p, m);
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

#set (gca, 'xtick', "");#the ticks aren't correct!
#set (gca, 'ytick', "");
set (gca, 'xtick', []);#3.6.2 #the ticks aren't correct!
set (gca, 'ytick', []);

pcolor(vXLabel, vYLabel, mHist2d); #mHist2D acts as color value
hold on
plot (a0r(1,:), a0r(2,:), "k")
plot (b0r(1,:), b0r(2,:), "k")
plot (c0r(1,:), c0r(2,:), "k")
plot (s0r(1,:), s0r(2,:), "k")
#plot3 (s1(1,:), s1(2,:), "k")
plot (o0r(1,:), o0r(2,:), "k") #a/b==abr2
plot (p0r(1,:), p0r(2,:), "k") #b/c==abr2
plot (q0r(1,:), q0r(2,:), "k") #a/b==abr1
plot (r0r(1,:), r0r(2,:), "k") #b/c==abr1
plot (u0r(1,:), u0r(2,:), "k") #a/b==abr3
plot (v0r(1,:), v0r(2,:), "k") #b/c==abr3
hold off
# not working with 3.8.2 or 4.0.0: set (gca, 'children', flipud (get (gca, 'children'))) # Octave has no uistack(hP,'top'): http://octave.1599824.n4.nabble.com/moving-graphic-objecs-td2997314.html http://stackoverflow.com/questions/7674700/how-to-change-the-order-of-lines-in-a-matlab-figure

shading flat; #means no border around each hist rectangle

text (c00p(1,1) - .02, c00p(2,1) + .02, "sphere\npoint", "horizontalalignment", "right"); #looks nicer
text (b0p(1,1) + .02, b0p(2,1), "circle\npoint");
text (c0p(1,1) - .02, c0p(2,1), "line\npoint", "horizontalalignment", "right");
text (a0p(1,floor(num/4)) + .02, a0p(2,floor(num/4)), "oblate line", "rotation", 30);
text (b0p(1,floor(num/2)) + .02, b0p(2,floor(num/2)), "ellipse arc", "rotation", -50);
text (c0p(1,floor(num/4*3)) - .02, c0p(2,floor(num/4*3)), "prolate line", "rotation", 90);
text (s0p(1,size(s0p,2)-20) + .02, s0p(2,size(s0p,2)-20), "separation curve", "rotation", -75);
#text (c00p(1,1), c00p(2,1) - .1, sprintf("# oblate-like: %d; # prolate-like: %d; # uncertain: %d; ratio: %.2f\n", Ns, Nz, Nsz, Ns/Nz));
text (c00p(1,1), c00p(2,1) - .1, sprintf("# oblate-like: %d; # prolate-like: %d; ratio: %.2f\n# spere-like: %d; # uncertain: %d", Ns, Nz, Ns/Nz,Nss,Nsz));

#xlabel("");
#ylabel("");
#set (gca, 'xtick', "");#the ticks aren't correct!
#set (gca, 'ytick', "");
set (gca, 'xtick', []);#3.6.2 #the ticks aren't correct!
set (gca, 'ytick', []);


N=max (max (mHist2d)); #max can be different
cmap= jet(N + 1);
cmap(1,:)=[1,1,1];
#cmap=vertcat([1,1,1],cmap);

colormap(cmap)
#colormap(hsv(128))
#caxis([0, 10])#ignore extreme etremes
#colorbar #show colorbar
axis ("equal", "off");#setting axis range here can be bad!



####printing now...

nplot= nplot + 1;
if !quiet
  printf("Printing plot # %d", nplot)
endif
print(sprintf("%s_%.2d.png", out3D, nplot), '-dpng', '-S800,800');#, '-F/usr/X11R6/lib/X11/fonts/msttf/arial.ttf');#, '-r100');
print(sprintf("%s_%.2d.svg", out3D, nplot), '-dsvg'); # sprintf('"-S%d,%d"', (b0p(1,1)-c00p(1,1))*scale, (c0p(2,1)-c00p(2,1))*scale ) does not influence BG rec size or pos # remove BG rec, not possible from within octave as gnuplot term is overwritten by print-command? (not tried) http://stackoverflow.com/questions/18169221/gnuplot-png-file-without-border-line
if !quiet
  printf(" done.\n", nplot)
endif

####printing end

##### END shearY(30°) 2D-hist of point cloud with jet cmap



##### BEGIN shearY(30°) 2D-hist of point cloud with b2w cmap

cmap= gray(N + 1);
cmap(1,:)=[1,1,1];
colormap(cmap)


####printing now...

nplot= nplot + 1;
if !quiet
  printf("Printing plot # %d", nplot)
endif
print(sprintf("%s_%.2d.png", out3D, nplot), '-dpng', '-S800,800');#, '-F/usr/X11R6/lib/X11/fonts/msttf/arial.ttf');#, '-r100');
print(sprintf("%s_%.2d.svg", out3D, nplot), '-dsvg');#has to be there for axis ("square") to work even with svg (-S not possible any more with gnuplot > 4.3.0 ???)
if !quiet
  printf(" done.\n", nplot)
endif

####printing end

##### END shearY(30°) 2D-hist of point cloud with b2w cmap

##### BEGIN shearY(30°) 2D-hist of point cloud with w2b cmap


cmap= gray(N + 1);
#cmap(1,:)=[1,1,1];
#colormap(cmap(end:-1:1))
colormap(flipud(cmap))



####printing now...

nplot= nplot + 1;
if !quiet
  printf("Printing plot # %d", nplot)
endif
print(sprintf("%s_%.2d.png", out3D, nplot), '-dpng', '-S800,800');#, '-F/usr/X11R6/lib/X11/fonts/msttf/arial.ttf');#, '-r100');
print(sprintf("%s_%.2d.svg", out3D, nplot), '-dsvg');#has to be there for axis ("square") to work even with svg (-S not possible any more with gnuplot > 4.3.0 ???)
if !quiet
  printf(" done.\n", nplot)
endif

####printing end


##### END shearY(30°) 2D-hist of point cloud with w2b cmap

