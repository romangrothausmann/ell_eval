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

### ToDo
## anotate guide-lines
## replace hist2d by hist3 (http://octave.sourceforge.net/statistics/function/hist3.html) # even octave forge hist2d (http://octave.sourceforge.net/plot/function/hist2d.html) is deprecated: http://stackoverflow.com/questions/35879723/what-is-the-right-package-to-use-hist2d-on-octave-cygwin#35880354

clear all; # prevents emacs from renaming file/function

######### BEGIN inline function definitions

function r = stereogproj(theta, phi, radius, phi0, lambda0)

for n=1:1:size(theta,2)

  #stereographic projection
  #http://mathworld.wolfram.com/StereographicProjection.html
  k= 2 * radius / (1 + sin(phi0) * sin(phi(n)) + cos(phi0) * cos(phi(n)) * cos(theta(n) - lambda0));
  x= k * cos(phi(n)) * sin(theta(n) - lambda0);
  y= k *              (cos(phi0) * sin(phi(n)) - sin(phi0) * cos(phi(n)) * cos(theta(n) - lambda0));
  r(:,n)= [x, y];#stereographic projected point list
endfor;
endfunction # stereogproj

%% use hist3 instead, recommended here: http://de.mathworks.com/matlabcentral/fileexchange/1487-2d-histogram-matrix
%% from http://de.mathworks.com/matlabcentral/fileexchange/1487-2d-histogram-matrix/content/hist2d.m
%function mHist = hist2d ([vY, vX], vYEdge, vXEdge)
%2 Dimensional Histogram
%Counts number of points in the bins defined by vYEdge, vXEdge.
%size(vX) == size(vY) == [n,1]
%size(mHist) == [length(vYEdge) -1, length(vXEdge) -1]
%
%EXAMPLE
%   mYX = rand(100,2);
%   vXEdge = linspace(0,1,10);
%   vYEdge = linspace(0,1,20);
%   mHist2d = hist2d(mYX,vYEdge,vXEdge);
%
%   nXBins = length(vXEdge);
%   nYBins = length(vYEdge);
%   vXLabel = 0.5*(vXEdge(1:(nXBins-1))+vXEdge(2:nXBins));
%   vYLabel = 0.5*(vYEdge(1:(nYBins-1))+vYEdge(2:nYBins));
%   pcolor(vXLabel, vYLabel,mHist2d); colorbar
function mHist = hist2d (mX, vYEdge, vXEdge)
nCol = size(mX, 2);
if nCol < 2
    error ('mX has less than two columns')
end

nRow = length (vYEdge)-1;
nCol = length (vXEdge)-1;

vRow = mX(:,1);
vCol = mX(:,2);

mHist = zeros(nRow,nCol);

for iRow = 1:nRow
    rRowLB = vYEdge(iRow);
    rRowUB = vYEdge(iRow+1);
    
    vColFound = vCol((vRow >= rRowLB) & (vRow < rRowUB)); # to be consistent with histc
    
    if (~isempty(vColFound))
        
        
        vFound = histc (vColFound, vXEdge);
        
        nFound = (length(vFound)-1);
        
        if (nFound ~= nCol)
            disp([nFound nCol])
            error ('hist2d error: Size Error')
        end
        
        [nRowFound, nColFound] = size (vFound);
        
        nRowFound = nRowFound - 1;
        nColFound = nColFound - 1;
        
        if nRowFound == nCol
            mHist(iRow, :)= vFound(1:nFound)';
        elseif nColFound == nCol
            mHist(iRow, :)= vFound(1:nFound);
        else
            error ('hist2d error: Size Error')
        end
    end
    
end
endfunction # hist2d


function annotate2D (a0p, b0p, c0p, s0p, c00p, num)
  text (c00p(1,1) - .02, c00p(2,1), "sphere\npoint", "horizontalalignment", "right"); #looks nicer
  text (b0p(1,1) + .02, b0p(2,1), "circle\npoint");
  text (c0p(1,1) - .02, c0p(2,1), "line\npoint", "horizontalalignment", "right");
  text (a0p(1,floor(num/3)) + .01, a0p(2,floor(num/3)) - .03, "oblate line", "rotation", 30);
  text (b0p(1,floor(num/1.8)) + .01, b0p(2,floor(num/1.8)), "ellipse arc", "rotation", -52);
  text (c0p(1,floor(num/4*2.5)) - .02, c0p(2,floor(num/4*2.5)), "prolate line", "rotation", 90);
  text (s0p(1,size(s0p,2)-38) + .02, s0p(2,size(s0p,2)-38), "separation curve", "rotation", -97);
endfunction # annotate2D


######### END inline function definitions


ee_min= 0.00000000000001; # why not 0?

quiet= 0;
arg_list = argv ();
if nargin != 2
  printf("Usage: %s <analysis.txt> ell-axis_error\n", program_name);
  exit(1)
else
  printf("Evaluating ellipsoids from %s...\n", arg_list{1});
  t= load(arg_list{1}); #octave_test02.txt;
  da= str2num(arg_list{2}) # const. abs. error in a; could be read from file for each a individually
endif



#graphics_toolkit fltk; #for octave 3.6.2, nice for viewing but not saving
graphics_toolkit gnuplot;


ps3d= 3; # plotting point size
ps2d= 3; # plotting point size


daxs= [da,da,da] #same abs. error for all axes # scale for "unit" error box (uEB) to become an "individual" error box (iEB)
clear da; # just to make sure it is not used later

sm_min=6 # number of EB corners that have to lie on the oblate side of the separation "surface" #has to be > length(es)/2 (== 4 for 3D, half the corners of EB), can only reach up to 8 (for each corner of EB, i.e. fully on the oblate side)!!!
ci_min=sm_min # number of EB corners that have to lie on the prolate side of the separation "surface"
sp_min=6 #can only reach up to 6 (for each coord of [-1,-1,-1] and [1,1,1] of uEB)!!! 6 means error box contains sphere-point, values below 6 would allow exception to this stringent criteria in some dimension (can't think of a situation when this would make sense)
spe=.2 #half sphere cube width
ev= ones(1,3)/norm(ones(1,3)); #unit vector in [1,1,1] direction

N= size(t, 1);
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

#### construct "unit" error-box (uEB)
##    1   1   1
##    1   1  -1
##    1  -1   1
##    1  -1  -1
##   -1   1   1
##   -1   1  -1
##   -1  -1   1
##   -1  -1  -1

es=[];
for x=0:1:1
  for y=0:1:1
    for z=0:1:1
      es=vertcat([x*2-1,y*2-1,z*2-1],es);
    end
  end
end


if sm_min <= length(es)/2 #has to be > length(es)/2
  exit(1)
endif

[fid, msg] = fopen (sprintf("%s.ells", arg_list{1}), "w");
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
    ax
    vol
    printf("Volum <= 0! Aborting\n")
    fflush(stdout);
    return #break is not enough any more !!!!!
  endif

  [axs, axi]= sort (ax); #making a < b < c if axis-names are not specially assigned!

  is_sp= 0;
  is_sm= 0;
  is_ci= 0;
  esum= 0;

  ##### instead of projecting the axes-vectors (and its EB) onto the unit sphere the error evaluation is done on spheres that contain each specific iEB corner, i.e. "individual" spheres for each EB corner, can be regarded as using a separation surface (separation line scaled to any radius) that is only evaluated on for spedific radii ("individual" spheres)
  eaxs= ev*norm(axs); #point within individual sphere limit?

  #this error evaluation can be done with fuzzy logic!
  #this makes sence if the error criterion is very hard and often only a single exception causes uncertainty
  #with fuzzy logic a minimum of certain evaluations have to be met
  #(fuzzy logic with reals is possible but very difficult: the ratio of the two volumes created by the cut of the separation surface through the error box)
  for i=1:1:length(es) # each uEB corner
    for j=1:1:size(es,2)
      esum= esum + daxs(j);
    endfor

    ### for ith corner: scale unit error-box (es) by axis-errors (daxs) and translate by axis-lengths (axs) such that the iEB corner (not center) lies on the "individual" evaluation sphere
    ee(1)= (axs(1) + es(i,1)*daxs(1));
    ee(2)= (axs(2) + es(i,2)*daxs(2));
    ee(3)= (axs(3) + es(i,3)*daxs(3));

    ### in case scaled and translated iEB corner (ee) is outside first quadrant (error > axes)
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

    ## an if-statement might speed this up (no issue so far, so left as is)
    is_sm+= ((ee(1) / ee(2)) < (ee(2) / ee(3))); # iEB corner on oblate side?
    is_ci+= ((ee(1) / ee(2)) > (ee(2) / ee(3))); # iEB corner on prolate side?
    ## if !quiet
    ##   printf("t(n,1): %f; n: %d; i: %d; a/b: %f; b/c: %f; sm: %d; ci: %d; sp: %d\n", t(n,1),n,i,(ee(1) / ee(2)), (ee(2) / ee(3)),is_sm,is_ci,is_sp)
    ## endif
  endfor # each uEB corner
  mee= esum / size(es,1) / size(es,2);
  tsum= tsum + mee;

  #is_sp= norm(axs/norm(axs) - ev/norm(ev)) < spe; #point within const. sphere limit?

  #if ==-case can be neglected this can be shortened:
  #is_ci= length(es)-is_sm
  #because of this is_sm > sm_min && is_ci > ci_min can never happen

  if is_sp >= sp_min # isSpherType? basically, only sp_min:= 6 makes sense, see above
    Nss++; # count as spherType
    ce(:,n)= [0,0,1]; # blue
    et= 0;
  else # if not isSpherType
    if is_sm >= sm_min # isOblateType if at least n (sm_min) iEB corners (is_sm) lie in oblate side of seperation surface
      if (axs(1) / axs(2) > axs(2) / axs(3))
        printf("is_sm wrong! This can happen for very big errors. \
            Couting as uncertain! t(n,1): %f\n", t(n,1))
        Nsz++;
        ce(:,n)= [0,0,0]; # black to distinguish from true uncertainType
        et= 0;
        continue
      endif
      Ns++; # count as oblateType (smartie)
      ce(:,n)= [0,1,0]; # green
      et= 1;
    else # not isSpherType nor isOblateType
      if is_ci >= ci_min # isProlateType if at least n (ci_min) iEB corners (is_ci) lie in prolate side of seperation surface
        if (axs(1) / axs(2) < axs(2) / axs(3))
          printf("is_ci wrong! t(n,1): %f\n", t(n,1))
          return #stop for --persist
        endif
        Nz++; # count as prolateType (zigarre)
        ce(:,n)= [1,0,0]; # red
        et= -1;
      else
        Nsz++; # count as uncertainType (not spher nor smartie nor cigar)
        ce(:,n)= [1,1,0]; # yellow
        et= 2;
      endif
    endif
  endif
  
  u(1,:,n)= v(:,axi(1));
  u(2,:,n)= v(:,axi(2));
  u(3,:,n)= v(:,axi(3));

  u_axes(:,n)= [axs(1), axs(2), axs(3)];

  [theta, phi, r]= cart2sph(axs(1), axs(2), axs(3)); 
  dp3ds(:,n)= [theta, phi, r];


  fprintf(fid, \
          "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\n", \
           axs, p_pos, v(:,axi(1)), v(:,axi(2)), v(:,axi(3)), et, p_index);
end;

fprintf(fid, \
	"## sphere-type: %d; oblate-type: %d; prolate-type: %d; uncertain-type: %d; oblate/prolate ratio: %.2f\n", Nss, Ns, Nz, Nsz, Ns/Nz);
fclose(fid);
Nss+Ns+Nz+Nsz
tsum / N


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
bl= [reshape(rthes,1,[])', reshape(rphis',1,[])']'; #point-coords of bm as vertcat
####done

#### guide points for b/c > a/b
xsn= 100;
xs= linspace(1,0.001,xsn); 
xs= xs.*xs; #make'm more evenly spaced; element by element multiplication

for n=1:1:xsn
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
color=vertcat(sqrt(2)*dp3d(2,:),sqrt(3)*dp3d(1,:),1/(1-sqrt(1/3))*(dp3d(3,:)-sqrt(1/3)));


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

axis ([0,1,0,1,0,1],"square");
set (gca (), "plotboxaspectratio", [1 1 1.45])#empirical ratio, needs to be set after axis ([0,1,0,1,0,1],"square");!!!

azimuth= 135;
elevation= acosd(dot([1,1,1], [1,1,0])/(norm([1,1,1]) * norm([1,1,0])));
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


print(sprintf("%s-%s.svg", arg_list{1}, "3Dsym"), '-dsvg', '-S800,800');


#### replot for different view that needs adjusted text coords

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

axis ([0,1,0,1,0,1],"square");
set (gca (), "plotboxaspectratio", [1 1 1.45])#empirical ratio, needs to be set after axis ([0,1,0,1,0,1],"square");!!!

xlabel("a");
ylabel("b");
zlabel("c");
text (c00(1,1), c00(2,1), c00(3,1), "sphere\npoint", "horizontalalignment", "right"); #looks nicer
text (b0(1,1), b0(2,1), b0(3,1), "circle\npoint");
text (c0(1,1), c0(2,1), c0(3,1), "line\npoint", "horizontalalignment", "right");
text (a0(1,floor(num/8)), a0(2,floor(num/8)) + .05, a0(3,floor(num/8)), "oblate arc", "rotation", 33);
text (b0(1,floor(num/2)) - .05, b0(2,floor(num/2)), b0(3,floor(num/2)), "ellipse arc", "rotation", -35);
text (c0(1,floor(num/4*3)) + .05, c0(2,floor(num/4*3)), c0(3,floor(num/4*3)), "prolate arc", "rotation", 118);
text (s0(1,size(s0,2)-38) - .03, s0(2,size(s0,2)-38), s0(3,size(s0,2)-38), "separation curve", "rotation", -66);

view(110, 10); # has to be set after plotting

print(sprintf("%s-%s.svg", arg_list{1}, "3Dasym"), '-dsvg', '-S800,800');



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
axis ("equal", "off"); # specifying axis extents here leads to two rect. polygons in SVG, one filled white, one as frame

annotate2D(a0p, b0p, c0p, s0p, c00p, num);

set (gca, 'xtick', []);#3.6.2 #the ticks aren't correct!
set (gca, 'ytick', []);


print(sprintf("%s-%s.svg", arg_list{1}, "2Dpdist"), '-dsvg', '-S800,800');

##### END 2D point cloud with RGB, Y point colouring


##### BEGIN 2D-hist of point cloud with jet cmap

bin= 30;
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

set (gca, 'xtick', []);#3.6.2 #the ticks aren't correct!
set (gca, 'ytick', []);

pcolor(vXLabel, vYLabel, mHist2d); #mHist2D acts as color value
shading flat; #means no border around each hist rectangle # should be directly after pcolor otherwise leads to "octave only supports 3-D filled triangular patches" if e.g. scatter is used after "hold on"
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


#axis ([c00p(1,1), b0p(1,1), c00p(2,1), c0p(2,1), ],"square");#setting axis range here can be bad!

annotate2D(a0p, b0p, c0p, s0p, c00p, num);

set (gca, 'xtick', []);#3.6.2 #the ticks aren't correct!
set (gca, 'ytick', []);


####create a colour map for many small values
N=max (max (mHist2d));
cmap= jet(N + 1);
cmap(1,:)=[1,1,1];

colormap(cmap)
colorbar #show colorbar
axis ("equal");#setting axis range here can be bad!


print(sprintf("%s-%s.svg", arg_list{1}, "2Dhist"), '-dsvg', '-S800,800');

##### END 2D-hist of point cloud with jet cmap


##### BEGIN annotations for shearY(30°) 2D-hist of point cloud with jet cmap

figure # creates empty plot canvas to "hold on"
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

annotate2D(a0p, b0p, c0p, s0p, c00p, num);

axis ("equal", "off");#setting axis range here can be bad!

print(sprintf("%s-%s.svg", arg_list{1}, "2Danno"), '-dsvg', '-S800,800');

##### END annotations for shearY(30°) 2D-hist of point cloud with jet cmap


##### BEGIN 2D point cloud with RGB, Y point colouring, for OV with 2Dshist

scatter (dp2d(1,:), dp2d(2,:), ps2d, ce', 's', 'filled')#3.6.2
hold on
plot (a0p(1,:), a0p(2,:), "k") # outer guide-lines are essential to get the same plot extent in SVG for OV
plot (b0p(1,:), b0p(2,:), "k")
plot (c0p(1,:), c0p(2,:), "k")
hold off
axis ("equal", "off");#setting axis range here can be bad!

set (gca, 'xtick', []);#3.6.2 #the ticks aren't correct!
set (gca, 'ytick', []);


print(sprintf("%s-%s.svg", arg_list{1}, "2DpdistOV"), '-dsvg', '-S800,800');

##### END 2D point cloud with RGB, Y point colouring, for OV with 2Dshist


##### BEGIN shearY(30°) 2D-hist of point cloud

function res= rotate(data, angle);

[Theta, R]= cart2pol(data(1,:), data(2,:));
[datar(1,:), datar(2,:)]= pol2cart(Theta + angle, R);
res= vertcat (datar(1,:), datar(2,:));
endfunction

function res= shear_y(data, m);
res= vertcat (data(1,:), data(2,:) - m * data(1,:));
endfunction

m= (b0p(2,1) - c00p(2,1)) / (b0p(1,1) - c00p(1,1));

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

lef= c00r(1,1);
rig= b0r(1,1);
bot= c00r(2,1);
top= c0r(2,1);
dxr= (rig - lef);
dyr= (top - bot);

xbin= round(bin / dx * dxr); #make the squares the same size as before
ybin= round(dyr / dxr * xbin); #make'm squares

vXEdge = linspace(lef, rig+(rig-lef)/(xbin-1), xbin+1);
vYEdge = linspace(bot, top+(top-bot)/(ybin-1), ybin+1);

mHist2d = hist2d([dp2dr(2,:)',dp2dr(1,:)'],vYEdge,vXEdge); 


nXBins = length(vXEdge);
nYBins = length(vYEdge);
vXLabel = vXEdge(1:(nXBins-1));
vYLabel = vYEdge(1:(nYBins-1));

pcolor(vXLabel, vYLabel, mHist2d); #mHist2D acts as color value
shading flat; #means no border around each hist rectangle

N=max (max (mHist2d)); #max can be different
## ##### jet cmap
## cmap= jet(N + 1);
## cmap(1,:)=[1,1,1];

## ##### b2w cmap
## cmap= gray(N + 1);
## cmap(1,:)=[1,1,1];

##### w2b cmap
cmap= flipud(gray(N + 1));


colormap(cmap)
axis ("equal", "off"); # setting axis range here can be bad!


print(sprintf("%s-%s.svg", arg_list{1}, "2Dshist"), '-dsvg', '-S800,800');

##### END shearY(30°) 2D-hist of point cloud
