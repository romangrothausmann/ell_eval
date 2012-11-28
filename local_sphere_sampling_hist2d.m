######octave function to evaluate point densities on a unit sphere by the local sphere sampling method
#####sampling only at 2D hist field-points
###code taken from octave/scripts/surf_orient_ng_05.m
###code regads point weights
###added test to check if all the surface of the sphere is sampled
###i.e. if sampling increment and sampling radius lead to sufficient overlap
###for sphere-sphere intersection radius see: http://mathworld.wolfram.com/Sphere-SphereIntersection.html


######
#REALLY MAKE SURE TO NOT MIX UP THETA AND PHI!!!!!!!

#########not yet tested!!!!


function lssHist = local_sphere_sampling_hist2d (mX, mW, xBin, yBin, lss_radius, verbose=1)

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
#cu= cos(2*pi/xBin)*cos(pi/yBin);
cc= acos(cos(2*pi/xBin)*cos(pi/yBin));
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
    phi= (n - 1) * (180 / (nRow - 1)) -  90; #scale range to [- 90; 90]: phi
    the= (m - 1) * (360 / (nCol - 1)) - 180; #scale range to [-180;180]: theta

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



