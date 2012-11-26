##!/net/home/ftd/localusr/bin/octave -qf
##this script needs the output from the ITK-program: analyse02 or label_ana01 

##script to evaluate AND compare global ellipsoid orientations from different rec.
#######
##REALLY MAKE SURE TO NOT MIX UP THETA AND PHI!!!!!!!

##from ellipsoid_oct29.m
##expects axis lenths (a,b,c) > 0; can often indirectly be achieved by vol > 64:
##e.g. awk 'NR>1 {if ( 64 < $17) {print $0}}' az239b_080725a_sirt_th-12707286.31_lr+255+0+1.ana > az239b_080725a_sirt_th-12707286.31_lr+255+0+1_v64.ana


##_02: added switch for weighted and unweighted lss

clear all;

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


  xbin=bin;
  ybin=bin;

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
    mhist(n,:,:)= local_sphere_sampling_hist2d([ua(n,:);ub(n,:);uc(n,:)]', uW(n,:), xbin, ybin, 10, 1);#'
  #size( hist_l)
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

# ###pure image for reprojection with e.g. G.projector;-)
# saveimage (sprintf("%s_%.2d_%d.ppm", outGO, nplot, n), mHist2d, "ppm");#loss of colour:-( only if c in [0;255]!!!!
# #imwrite (pci, sprintf("%s_%.2d.png", out2D, nplot));#loss of colour:-(

# ####printing end

endfunction#plot_hists




###color-map:


Nmax= 128;
midcmap= jet(Nmax + 0.2 * Nmax);#dark red region is about 1/5th of jet length


##some octave magick to delete dark red colours;-)
del_p = midcmap(:,1) < 1 & midcmap(:,2) == 0  & midcmap(:,3) == 0;
size(midcmap)
midcmap(del_p,:)= [];
size(midcmap)
##magick done

cmap= midcmap;



######doing it number-weighted

printf("Doing number weighted (nw) evaluation...\n")

###first read in and do lss

for k=1:1:Na;

  hist_l= read_and_lss(arg_list{k}, 0);#0: nw

  size(hist_l)
  #max(hist_l)
  #max(reshape(hist_l,1,[]))
  hists(k,:,:,:)= hist_l;
  #max(reshape(hists(k,:,:,:),1,[]))

endfor#k


###get absolut maximum of all hists

abs_max= max(reshape(hists,1,[]))

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

  size(hist_l)
  #max(hist_l)
  #max(reshape(hist_l,1,[]))
  hists(k,:,:,:)= hist_l;
  #max(reshape(hists(k,:,:,:),1,[]))

endfor#k


###get absolut maximum of all hists

abs_max= max(reshape(hists,1,[]))

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
