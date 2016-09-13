##!/net/home/ftd/localusr/bin/octave -qf
##this script needs the output from the ITK-program: analyse02 or label_ana01 

##script to evaluate AND compare global ellipsoid orientations from different rec.

#######################################################
##REALLY MAKE SURE TO NOT MIX UP THETA AND PHI!!!!!!!##
##AND TO USE RAD NOT DEG FOR ANY ANGLULAR FUNCTION!!!##
#######################################################

##from ellipsoid_oct29.m
##expects axis lenths (a,b,c) > 0; can often indirectly be achieved by vol > 64:
##e.g. awk 'NR>1 {if ( 64 < $17) {print $0}}' az239b_080725a_sirt_th-12707286.31_lr+255+0+1.ana > az239b_080725a_sirt_th-12707286.31_lr+255+0+1_v64.ana


##_02: added switch for weighted and unweighted lss
##_03: using just jet as cmap
##_04: full sphere sampling
##_05: hemisphere sampling
##_06: added printing of Max ans Sum of hist
##_07: changed cw scale to fit new analysis

clear all; # prevents emacs from renaming file/function

######### BEGIN inline function definitions

######octave function to evaluate point densities on a unit sphere by the local sphere sampling method
#####sampling only at 2D hist field-points
###code taken from octave/scripts/surf_orient_ng_05.m
###code regads point weights
###added test to check if all the surface of the sphere is sampled
###i.e. if sampling increment and sampling radius lead to sufficient overlap
###for sphere-sphere intersection radius see: http://mathworld.wolfram.com/Sphere-SphereIntersection.html

###### sampling only a hemisphere, i.e. theta: [- 90; 90]
###used local_sphere_sampling_hist2d_03.m



######
#REALLY MAKE SURE TO NOT MIX UP THETA AND PHI!!!!!!!

#########not yet tested!!!!


function lssHist = local_hemisphere_sampling_hist2d (mX, mW, xBin, yBin, lss_radius, verbose=1)

#clear all;
#printf("Not yet tested!!!. Aborting!\n")
#exit(1)

warning ("error", "Octave:divide-by-zero");

#addpath("~/octave")

sp_lim= lss_radius / 180 * pi; #local sphere radius; does not have to be samller than 1 if range(norm(data point-vectors)) > 1 !!!

points= mX;
weights= mW;

nCol = xBin;
nRow = yBin;

lssHist = zeros(nRow,nCol);

Nx= size(points, 1);

quiet= !verbose;


####check if sampling increment and sampling radius lead to sufficient overlap

d= 1;
R= 1;
r= sp_lim;
a= 1/2/d*sqrt((-d+r-R)*(-d-r+R)*(-d+r+R)*(d+r+R)); #http://mathworld.wolfram.com/Sphere-SphereIntersection.html

###law of cosines of angles for a rectangular spherical triangle #http://mathworld.wolfram.com/SphericalTrigonometry.html
###cos(c)= cos(a)cos(b)

#ca= cos(2*a);
#cu= cos(pi/xBin)*cos(pi/yBin);
cc= acos(cos(pi/xBin)*cos(pi/yBin));
da= 2*a;
#cc, da

#if (cu < ca) ##NOTE using cosine invertes relation!!!
if (cc > da)
  cc, da
  printf("WARNING: The chosen sample radius and binning will not lead to a complete covering of the unit sphere!")
  fflush(stdout);
  #cc, da #does not output correctly, why???? Has to be before printf! Why???
  #return
endif

###check done.


if !quiet
  printf("Checking if points lie on sphere...")
  fflush(stdout);
endif

nn= sqrt(dot(points', points')); #norm(c_norm)
if(range(nn) >= .0001) #there is a rounding error!
  max(nn)
  min(nn)
  printf("Data points vectors not normalized! Using them as they are!!!\n")
endif


if !quiet
  printf(" done.\n")
  fflush(stdout);
endif



###local sphere average evaluation


if !quiet
  printf("Calculating weighted point densities by local sphere sampling...\n")
  fflush(stdout);
endif

Np= nRow * nCol;
prn_int= round(Np / 1000); # 1000 progress reports;-)

i= 0;
for n=1:1:nRow #phi
  for m=1:1:nCol #theta
    i+= 1;
    the= (m - 1) * (180 / (nCol - 1)) -  90; #scale range to [- 90; 90]: theta
    phi= (n - 1) * (180 / (nRow - 1)) -  90; #scale range to [- 90; 90]: phi

    #clear xt yt zt;
    [xt, yt, zt]= sph2cart (the/180*pi, phi/180*pi, 1);#projection of 3D points onto unit sphere
    xyz= vertcat(xt, yt, zt)';
    
    nnn= sqrt(dot(xyz', xyz')); #norm(c_norm)
    if(range(nnn) >= .0001) #there is a rounding error!
      max(nn)
      min(nn)
      printf("Reference point not normalized! This should not happen! Aborting!\n")
      exit(1)
    endif

    shifted_n= points - (ones(Nx,1) * xyz);
    shifted_l= sqrt(dot(shifted_n',shifted_n'));#all lengths

    for k=1:1:Nx;
      if (shifted_l(k) <= sp_lim)
	lssHist(n, m)+= weights(k);# / a_tot;
      endif
    endfor
    if !quiet
      if(mod(i, prn_int) == 0) #dont report for each n!
	printf("\r%6.2f%% ", 100 * i / Np )
	fflush(stdout);
      endif
    endif

  endfor
endfor

if !quiet
  printf(" done.\n")
  fflush(stdout);
endif

#return lssHist
endfunction # local_hemisphere_sampling_hist2d

######### END inline function definitions


global quiet= 0;

global bin= 30;

## global ps3d= 3;#for octave 3.6.2
## global ps2d= 3;#for octave 3.6.2
global nplot= 0;

addpath("~/octave/functions")
#graphics_toolkit fltk; #for octave 3.6.2, nice for viewing but not saving
graphics_toolkit gnuplot;

#global arg_list Na
arg_list = argv ();
Na= nargin;

if Na <= 0
  printf("Usage: %s <analysis1.txt> [analysis2.txt ...] [-q]\n", program_name);
  exit(1)
endif


function mhist = read_and_lss(file_name, weighted)
  ##using a function renders clearing unnecessary:
  #clear t p_index p_pos ax v vol axs axi u u_axes uao

  #global arg_list ##not possible for string list???
  global bin #needs to be declared global in main AND in function
  
  printf("Evaluating ellipsoids from %s...\n", file_name);
  t= load(file_name);
  N= size(t, 1);

  for l=1:1:N;

    p_index= t(l,1);
    p_pos=  [t(l,2),t(l,3),t(l,4)];
    ax= [t(l,5),t(l,6),t(l,7)];
    v=  [t(l,8),t(l,9),t(l,10);
	 t(l,11),t(l,12),t(l,13);
	 t(l,14),t(l,15),t(l,16),]';
    p_V= t(l,17);
    
    vol= 4/3*pi*ax(1)*ax(2)*ax(3);
    if(vol>0)
      ax= 2*(p_V/vol)^(1/3)*ax;
    else
      l
      t(l,:)
      I
      ax
      vol
      printf("Volum <= 0! Aborting\n")
      fflush(stdout);
      #exit
      return #break is not enough any more !!!!!
    endif

    [axs, axi]= sort (ax); #making a < b < c if axis-names are not specially assigned!
    #eu= euler_angles(v(:,axi(1)), v(:,axi(2)), v(:,axi(3))); #index ordered v


    u(1,:,l)= v(:,axi(1)); #xyz of a-axis
    u(2,:,l)= v(:,axi(2)); #xyz of b-axis
    u(3,:,l)= v(:,axi(3)); #xyz of c-axis
    #u(1,:,n)= v(:,1);
    #u(2,:,n)= v(:,2);
    #u(3,:,n)= v(:,3);

    u_axes(:,l)= [axs(1), axs(2), axs(3)];

  endfor#l
  clear l #just to make sure further use causes an error ;-)


  xbin=bin; #lss done on hemisphere!!!
  ybin=bin; #lss done on hemisphere!!!

  ###loop over axis a, b, c
  for n=1:1:3

    uao(n,:)= squeeze(u(n,1,:)); #x-coord.
    ubo(n,:)= squeeze(u(n,2,:)); #y-coord.
    uco(n,:)= squeeze(u(n,3,:)); #z-coord.

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

    if weighted
      uW(n,:)= [[u_axes(n_1,:)./sqrt(u_axes(n_2,:).*u_axes(n_3,:))]';[u_axes(n_1,:)./sqrt(u_axes(n_2,:).*u_axes(n_3,:))]'];#for weighted density
    else
      uW= ones(3, size(ua(n,:),2));#for number density
    endif

    #uW= uW./sum(uW);
    uW(n,:)= uW(n,:)./size(uao(n,:),2);

    #u_lss(n,:)= local_sphere_sampling([ua(n,:);ub(n,:);uc(n,:)]', uW(n,:), 10, 1);#'
    #u_lss(n,:)= local_sphere_sampling_at_pos([ua(n,:);ub(n,:);uc(n,:)]', uW(n,:), 10, 1);#'
    #hist_l= local_sphere_sampling_hist2d([ua(n,:);ub(n,:);uc(n,:)]', uW(n,:), xbin, ybin, 10, 1);#'
    hist_l= local_hemisphere_sampling_hist2d([ua(n,:);ub(n,:);uc(n,:)]', uW(n,:), xbin, ybin, 10, 1);#'

    #size(hist_l)
    printf("Max of hist: %f\n", max(reshape(hist_l,1,[])))
    printf("Sum of hist: %f (greater 7.2 (0.008*900) because it dep. on the ldd and the overlap of the adjacent lss, i.e. especially high for c as ldd is high around the poles.)\n", sum(reshape(hist_l,1,[])))
    #[x, ix] = max(hist_l)

    mhist(n,:,:)= hist_l;

  endfor#n
  clear n #just to make sure further use causes an error ;-)

endfunction #read_and_lss






function plot_hists (mHist2d, fname, n, bin, cmap, abs_max)
         
  global nplot quiet


  out3D=sprintf("%s_AR", fname);
  outGO=sprintf("%s_GO", fname);


  #mHist2d=squeeze(hists(k,n,:,:));


  ###2D-hist:

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
  #lef= -180; 
  #rig=  180; 
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

  ybin=bin;
  xbin=round(ybin / dy * dx); #make'm squares
  #ybin= round( dy / dx * xbin) + 1; #make'm odd!
  #xbin= round(bin / dx * dxr); #make the squares the same size as before
  #ybin= round(dyr / dxr * xbin); #make'm squares

  #x= vertcat(x, l0(1,:)');#for testing guide line pos
  #y= vertcat(y, l0(2,:)');



  vXEdge = linspace(lef, rig+(rig-lef)/(xbin-1), xbin+1);#add one entry at end
  vYEdge = linspace(bot, top+(top-bot)/(ybin-1), ybin+1);

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

  #return

  #text (-3, -1.5, sprintf("# cells: %d", N));

  #caxis([Nmin,Nmax]); #color between icosaedron and tetraedron
  caxis([0, abs_max]); #color between icosaedron and tetraedron
  #caxis('auto')
  colormap(cmap);
  colorbar #show colorbar
  #axis ("square");#setting axis range here can be bad!
  axis("image");#square and no extra space ;-)

  #set (gca, 'xtick', [-90:30:90]);
  #set (gca, 'ytick', [-90:30:90]);
  ##no range can be specified for minor tick (not as in gnuplot mxticks)
  #set (gca, 'xminortick', [-90:10:90]);
  #set (gca, 'yminortick', [-90:10:90]);
  #set (gca, "xminortick", "on", "yminortick", "on")
  ##using main ticks and remove them in inkscape...
  set (gca, 'xtick', [-90:10:90]);
  set (gca, 'ytick', [-90:10:90]);
  ##or set them this way ;-)
  set (gca, 'xticklabel', ['-90';'  ';'  ';'-60';'  ';'  ';'-30';'  ';'  ';'0';'  ';'  ';'30';'  ';'  ';'60';'  ';'  ';'90']);
  set (gca, 'yticklabel', ['-90';'  ';'  ';'-60';'  ';'  ';'-30';'  ';'  ';'0';'  ';'  ';'30';'  ';'  ';'60';'  ';'  ';'90']);

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
  #saveimage (sprintf("%s_%.2d_%d.ppm", outGO, nplot, n), mHist2d, "ppm");#loss of colour:-( only if c in [0;255]!!!!
  #imwrite (pci, sprintf("%s_%.2d.png", out2D, nplot));#loss of colour:-(
  #img= imagesc(mHist2d); ##cannot be printed as image?!
  
  ##individual (relative) scaling 
  #i_min= min(reshape(mHist2d,1,[]));
  #i_max= max(reshape(mHist2d,1,[]));

  ##absolut scaling
  i_min= 0; 
  i_max= abs_max; #values above 256 do not matter; see /usr/share/octave/3.4.3/m/image/ind2rgb.m at line 44
  
  mHist2d_s= (mHist2d - i_min)/(i_max-i_min)*255;
  mHist2d_s= flipud(mHist2d_s); #y-axis points down in images!

  #min(reshape(mHist2d_s,1,[]))
  #max(reshape(mHist2d_s,1,[]))

  imwrite(uint8(mHist2d_s/1+1), sprintf("%s_%.2d_%d_m2g.png", outGO, nplot, n), "png");
  imwrite(uint8(mHist2d_s/1+1), colormap(cmap), sprintf("%s_%.2d_%d_m2c.png", outGO, nplot, n), "png"); #hist values must be in [1;256]; only [1;256] from cmap is used!!!
  #imwrite(ind2rgb(gray2ind(mHist2d)), sprintf("%s_%.2d_%d_m2c.png", outGO, nplot, n), "png");
  #print(sprintf("%s_%.2d_%d_m2i.png", outGO, nplot, n), '-dpng');

  ####printing end

endfunction#plot_hists


######################## main

###color-map:



Nmax= 256; #using 256 ensures cmap to be suitable for imwrite
c_extension= 1.1429;#dark red region is about 1/5th of jet length
#midcmap= jet(Nmax*c_extension);#dark red region is about 1/5th of jet length
midcmap= jet(Nmax);#dark red region is about 1/5th of jet length

# ##some octave magick to delete dark red colours;-)
# del_p = midcmap(:,1) < 1 & midcmap(:,2) == 0  & midcmap(:,3) == 0;
# before= size(midcmap)
# midcmap(del_p,:)= [];
# after= size(midcmap)
# c_extension= before/after ##== 1.1428
# ##magick done

cmap= midcmap;



######doing it number-weighted

printf("Doing number weighted (nw) evaluation...\n")

###first read in and do lss

for k=1:1:Na;

  hist_l= read_and_lss(arg_list{k}, 0);#0: nw

  #size(hist_l)
  #max(hist_l)
  #max(reshape(hist_l,1,[]))
  hists(k,:,:,:)= hist_l;
  #max(reshape(hists(k,:,:,:),1,[]))

endfor#k


###get absolut maximum of all hists

#abs_max= max(reshape(hists,1,[])) ###calc max
abs_max= 0.1*c_extension  ###set max to 0.1 as in diss

##test
#imagesc(squeeze(hists(1,1,:,:))')
#axis ("square");#setting axis range here can be bad!
#return




###plot hist now that abs_max is known

for k=1:1:Na;
  for n=1:1:3
    plot_hists(squeeze(hists(k,n,:,:)), sprintf("%s_nw", arg_list{k}), n, bin, cmap, abs_max);
  endfor#n
endfor#k





######doing it weighted

printf("Doing weighted (cw) evaluation...\n")

###first read in and do lss

for k=1:1:Na;

  hist_l= read_and_lss(arg_list{k}, 1);#1: cw

  #size(hist_l)
  #max(hist_l)
  #max(reshape(hist_l,1,[]))
  hists(k,:,:,:)= hist_l;
  #max(reshape(hists(k,:,:,:),1,[]))

endfor#k


###get absolut maximum of all hists

#abs_max= max(reshape(hists,1,[])) ###calc max
#abs_max= 0.3*c_extension  ###set max to 0.3 as in diss
abs_max= 0.23*c_extension  ###set max according to new eval 

##test
#imagesc(squeeze(hists(1,1,:,:))')
#axis ("square");#setting axis range here can be bad!
#return


###plot hist now that abs_max is known

for k=1:1:Na;
  for n=1:1:3
    plot_hists(squeeze(hists(k,n,:,:)), sprintf("%s_cw", arg_list{k}), n, bin, cmap, abs_max);
  endfor#n
endfor#k
