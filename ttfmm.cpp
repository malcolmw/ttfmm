#include <cmath>
#include <float.h>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <deque>
#include "ttfmm.h"

void init_VelocityGrid(GeoGrid *VV){
  float value;
  std::ifstream inf;
  VV->nlat=32;
  VV->nlon=32;
  VV->nz=23;
  VV->lat0=32.38;
  VV->lon0=-118.17;
  VV->z0=-2.0;
  VV->dlat=0.07;
  VV->dlon=0.09;
  VV->dz=1.0;

  VV->value=(double ***) malloc(VV->nlat*sizeof(double));
  for (int ilat=0; ilat<VV->nlat; ilat++){
    VV->value[ilat]=(double **) malloc(VV->nlon*sizeof(double));
    for (int ilon=0; ilon<VV->nlon; ilon++){
      VV->value[ilat][ilon]=(double *) malloc(VV->nz*sizeof(double));
    }
  }

  inf.open("/Users/malcolcw/Projects/Shared/Velocity/FANG2016/original/Vp.dat");
  for (int iz=0; iz<VV->nz; iz++){
    for (int ilat=0; ilat<VV->nlat; ilat++){
      for (int ilon=0; ilon<VV->nlon; ilon++){
        inf >> value;
        VV->value[ilat][ilon][iz] = (double) value;
      }
    }
  }
  inf.close();
}

void init_TTGrid(GeoGrid *TT){
  TT->nlat=32;
  TT->nlon=32;
  TT->nz=23;
  TT->lat0=32.38;
  TT->lon0=-118.17;
  TT->z0=-2.0;
  TT->dlat=0.07;
  TT->dlon=0.09;
  TT->dz=1.0;

  TT->value=(double ***) malloc(TT->nlat*sizeof(double));
  for (int ilat=0; ilat<TT->nlat; ilat++){
    TT->value[ilat]=(double **) malloc(TT->nlon*sizeof(double));

    for (int ilon=0; ilon<TT->nlon; ilon++){
      TT->value[ilat][ilon]=(double *) malloc(TT->nz*sizeof(double));

      for (int iz=0; iz<TT->nz; iz++){
        TT->value[ilat][ilon][iz]=DBL_MAX;
      }
    }
  }
}

void init_source(GeoGrid *TT, std::deque<GridNode> &active, float lat, float lon, float z){
  GridNode node;
  node.ilat=TT->nlat/2 ;
  node.ilon=TT->nlon/2 ;
  node.iz=TT->nz/2 ;
  //node.ilat=0;
  //node.ilon=0;
  //node.iz=0;
  node.value=0;
  TT->value[node.ilat][node.ilon][node.iz]=node.value;
  active.push_back(node);
}

void free_grid(GeoGrid *GG){
  for (int ilat=0; ilat<GG->nlat; ilat++){
    for (int ilon=0; ilon<GG->nlon; ilon++){
      free(GG->value[ilat][ilon]);
    }
    free(GG->value[ilat]);
  }
  free(GG->value);
}

void dump_grid(GeoGrid GG){
  printf("%d %d %d\n", GG.nlat, GG.nlon, GG.nz);
  printf("%f %f %f\n", GG.lat0, GG.lon0, GG.z0);
  printf("%f %f %f\n", GG.dlat, GG.dlon, GG.dz);
  for (int ilat=0; ilat<GG.nlat; ilat++){
    for (int ilon=0; ilon<GG.nlon; ilon++){
      for (int iz=0; iz<GG.nz; iz++){
        printf("%5.2f ", GG.value[ilat][ilon][iz]);
      }
      printf("\n");
    }
  }
}

void print_grid(GeoGrid GG){
  printf("%d %d %d\n", GG.nlat, GG.nlon, GG.nz);
  printf("%f %f %f\n", GG.lat0, GG.lon0, GG.z0);
  printf("%f %f %f\n", GG.dlat, GG.dlon, GG.dz);
  for (int iz=0; iz<GG.nz; iz++){
    for (int ilat=GG.nlat-1; ilat>=0; ilat--){
      for (int ilon=GG.nlon-1; ilon>=0; ilon--){
        if ( GG.value[ilat][ilon][iz] != DBL_MAX ){
          printf("%5.2f ", GG.value[ilat][ilon][iz]);
        }
        else{
          printf("  -   ");
        }
      }
      printf("\n");
    }
    //printf("\n");
  }
}

void print_active(std::deque<GridNode> active){
  for (int i=0; i<active.size(); i++){
    printf("INFO:: active[%d]=(%7.2f, %d, %d, %d)\n",
           i,
           active[i].value,
           active[i].ilat,
           active[i].ilon,
           active[i].iz);
  }
}

bool is_frozen(int ilat, int ilon, int iz, std::deque<GridNode> &frozen){
  for (int i=0; i<frozen.size(); i++){
    if (   (frozen[i].ilat == ilat)
        && (frozen[i].ilon == ilon)
        && (frozen[i].iz == iz) ){
      return(true);
    }
  }
  return(false);
}

void geo2sph(double *coords){
  double lat=coords[0];
  double lon=coords[1];
  double z=coords[2];
  coords[0]=EARTH_RADIUS-z;
  coords[1]=(90.0-lat)*PI/180.0;
  coords[2]=lon*PI/180.0;
}

void update_active(GridNode node, std::deque<GridNode> &active){
  for (int i=0; i<active.size(); i++){
    if (   active[i].ilat == node.ilat
        && active[i].ilon == node.ilon
        && active[i].iz   == node.iz ){
      active[i].value = node.value;
      return;
    }
  }
  if (node.value < active[0].value) active.push_front(node);
  else active.push_back(node);
}

void update(GeoGrid *TT, GeoGrid *VV,
            std::deque<GridNode> &active,
            std::deque<GridNode> &frozen)
{
  hpsort(active.size()-1, active);
  GridNode hot=active[0], nbr;
  active.pop_front();
  int near[6][3]={{hot.ilat-1, hot.ilon,   hot.iz},
                  {hot.ilat+1, hot.ilon,   hot.iz},
                  {hot.ilat,   hot.ilon-1, hot.iz},
                  {hot.ilat,   hot.ilon+1, hot.iz},
                  {hot.ilat,   hot.ilon,   hot.iz-1},
                  {hot.ilat,   hot.ilon,   hot.iz+1}};
  double ur, ut, up, ddr2, ddt2, ddp2;
  double A, B, C, det, tt;
  double drho, dtheta, dphi;
  double rho, theta, phi;
  double coords[3];
  drho=TT->dz;
  dtheta=TT->dlat*PI/180.0;
  dphi=TT->dlon*PI/180.0;
  for (int inear=0; inear<6; inear++){
    nbr.ilat = near[inear][0];
    nbr.ilon = near[inear][1];
    nbr.iz = near[inear][2];
    // Don't update neighbours if they don't exist or are frozen.
    if ((nbr.ilat < 0) || (nbr.ilat>TT->nlat-1)) continue;
    if ((nbr.ilon < 0) || (nbr.ilon>TT->nlon-1)) continue;
    if ((nbr.iz < 0)   || (nbr.iz>TT->nz-1    )) continue;
    if (is_frozen(nbr.ilat, nbr.ilon, nbr.iz, frozen)) continue;
    nbr.value = TT->value[nbr.ilat][nbr.ilon][nbr.iz];

    coords[0]=TT->lat0+nbr.ilat*TT->dlat;
    coords[1]=TT->lon0+nbr.ilon*TT->dlon;
    coords[2]=TT->z0+nbr.iz*TT->dz;
    geo2sph(coords);
    rho=coords[0];
    theta=coords[1];
    phi=coords[2];

    // Compute a travel-time update for nbr
    ur=DBL_MAX;
    if ((nbr.iz-1>=0) && (TT->value[nbr.ilat][nbr.ilon][nbr.iz-1]<ur))
      ur=TT->value[nbr.ilat][nbr.ilon][nbr.iz-1];
    if ((nbr.iz+1<TT->nz) && (TT->value[nbr.ilat][nbr.ilon][nbr.iz+1]<ur))
      ur=TT->value[nbr.ilat][nbr.ilon][nbr.iz+1];
    if ((ur==DBL_MAX) || (ur >= nbr.value)){
      ur=0;
      ddr2=0;
    }
    else
      ddr2=1.0/pow(drho, 2);

    ut=DBL_MAX;
    if ((nbr.ilat-1>=0) && (TT->value[nbr.ilat-1][nbr.ilon][nbr.iz]<ut))
      ut=TT->value[nbr.ilat-1][nbr.ilon][nbr.iz];
    if ((nbr.ilat+1<TT->nlat) && (TT->value[nbr.ilat+1][nbr.ilon][nbr.iz]<ut))
      ut=TT->value[nbr.ilat+1][nbr.ilon][nbr.iz];
    if ((ut==DBL_MAX) || (ut>=nbr.value)){
      ut=0;
      ddt2=0;
    }
    else
      ddt2=1.0/pow(rho*dtheta, 2);

    up=DBL_MAX;
    if ((nbr.ilon-1>=0) && (TT->value[nbr.ilat][nbr.ilon-1][nbr.iz]<up))
      up=TT->value[nbr.ilat][nbr.ilon-1][nbr.iz];
    if ((nbr.ilon+1<TT->nlon) && (TT->value[nbr.ilat][nbr.ilon+1][nbr.iz]<up))
      up=TT->value[nbr.ilat][nbr.ilon+1][nbr.iz];
    if ((up==DBL_MAX) || (up >= nbr.value)){
      up=0;
      ddp2=0;
    }
    else
      ddp2=1.0/pow(rho*sin(theta)*dphi, 2);

    A=ddr2+ddt2+ddp2;
    B=-2*(ddr2*ur+ddt2*ut+ddp2*up);
    C=pow(ur,2)*ddr2+pow(ut,2)*ddt2+pow(up,2)*ddp2;
    C-=1.0/pow(VV->value[nbr.ilat][nbr.ilon][nbr.iz],2);
    det=pow(B,2)-4.0*A*C;
    if (det<0){
      printf("WARNING:: determinate was negative during update - %.5f\n", det);
      A=ddt2+ddp2;
      B=-2*(ddt2*ut+ddp2*up);
      C=pow(ut,2)*ddt2+pow(up,2)*ddp2;
      C-=1.0/pow(VV->value[nbr.ilat][nbr.ilon][nbr.iz],2);
      det=pow(B,2)-4.0*A*C;
      //for (int iz=-1; iz<2; iz++){
      //  for (int ilat=1; ilat>-2; ilat--){
      //    for (int ilon=1; ilon>-2; ilon--){
      //      if (TT->value[nbr.ilat+ilat][nbr.ilon+ilon][nbr.iz+iz] == DBL_MAX)
      //        printf("  -   ");
      //      else
      //        printf("%5.2f ", TT->value[nbr.ilat+ilat][nbr.ilon+ilon][nbr.iz+iz]);
      //    }
      //    printf("\n");
      //  }
      //  printf("\n\n");
      //}
      //printf("rho      =%.5f\n", rho);
      //printf("theta    =%.5f\n", theta);
      //printf("phi      =%.5f\n", phi);
      //printf("nbr.ilat =%d\n", nbr.ilat);
      //printf("nbr.ilon =%d\n", nbr.ilon);
      //printf("nbr.iz   =%d\n", nbr.iz);
      //printf("nbr.value=%.5f\n", nbr.value);
      //printf("VV       =%.5f\n", VV->value[nbr.ilat][nbr.ilon][nbr.iz]);
      //printf("Dr(U)    =%.5f\n", (nbr.value-ur)*ddr2);
      //printf("Dt(U)    =%.5f\n", (nbr.value-ut)*ddt2);
      //printf("Dp(U)    =%.5f\n", (nbr.value-up)*ddp2);
      //printf("ur       =%.5f\n", ur);
      //printf("ut       =%.5f\n", ut);
      //printf("up       =%.5f\n", up);
      //printf("ddr2     =%.5f\n", ddr2);
      //printf("ddt2     =%.5f\n", ddt2);
      //printf("ddp2     =%.5f\n", ddp2);
      //printf("A        =%.5f\n", A);
      //printf("B        =%.5f\n", B);
      //printf("C        =%.5f\n", C);
      //printf("det      =%.5f\n", det);
      //printf("WARNING:: negative determinant on update\n");
      //exit(-1);
    }
    tt=(-B+sqrt(det))/(2.0*A);
    if (tt<nbr.value){
      nbr.value=tt;
      TT->value[nbr.ilat][nbr.ilon][nbr.iz]=nbr.value;
      update_active(nbr, active);
    }
  }
  frozen.push_back(hot);
}

void update_dep(GeoGrid *TT, GeoGrid *VV,
            std::deque<GridNode> &active,
            std::deque<GridNode> &frozen)
{
  hpsort(active.size()-1, active);
  int ilat, ilon, iz;
  GridNode hot=active[0], test;
  double dtheta=TT->dlat*PI/180.0;
  double dphi=TT->dlon*PI/180.0;
  double drho=TT->dz;
  double ur, ut, up, ddr2, ddt2, ddp2, A, B, C;
  double coords[3];
  int near[6][3]={{hot.ilat-1, hot.ilon,   hot.iz},
                  {hot.ilat+1, hot.ilon,   hot.iz},
                  {hot.ilat,   hot.ilon-1, hot.iz},
                  {hot.ilat,   hot.ilon+1, hot.iz},
                  {hot.ilat,   hot.ilon,   hot.iz-1},
                  {hot.ilat,   hot.ilon,   hot.iz+1}};
  //printf("START update()\n");
  //print_active(active);
  active.pop_front();
  printf("hot node is (%.2f, %d, %d, %d); there are currently %zd active nodes\n",
         hot.value,
         hot.ilat,
         hot.ilon,
         hot.iz,
         active.size());
  // Determine the geographic coordinates of the hot grid node
  coords[0] = TT->lat0+TT->dlat*hot.ilat;
  coords[1] = TT->lon0+TT->dlon*hot.ilon;
  coords[2] = TT->z0+TT->dz*hot.ilon;
  // Convert geographic coordinates to spherical
  geo2sph(coords);
  for (int inear=0; inear<6; inear++){
    test.ilat=near[inear][0];
    test.ilon=near[inear][1];
    test.iz=near[inear][2];
    if ((test.ilat < 0) || (test.ilat >= TT->nlat)) continue;
    if ((test.ilon < 0) || (test.ilon >= TT->nlon)) continue;
    if ((test.iz   < 0) || (test.iz   >= TT->nz  )) continue;
    if (is_frozen(test.ilat, test.ilon, test.iz, frozen)){
      printf("node (%d, %d, %d) is frozen\n", test.ilat, test.ilon, test.iz);
      continue;
    }
    ut=std::min(TT->value[std::max(test.ilat-1,0)][test.ilon][test.iz],
                TT->value[std::min(test.ilat+1,TT->nlat-1)][test.ilon][test.iz]);
    up=std::min(TT->value[test.ilat][std::max(test.ilon-1,0)][test.iz],
                TT->value[test.ilat][std::min(test.ilon+1,TT->nlon-1)][test.iz]);
    ur=std::min(TT->value[test.ilat][test.ilon][std::max(test.iz-1,0)],
                TT->value[test.ilat][test.ilon][std::min(test.iz+1,TT->nz-1)]);
    if (ur >= TT->value[test.ilat][test.ilon][test.iz]){
      ur=0;
      ddr2=0;
    }else{
      ddr2=1.0/pow(drho, 2);
    }
    if (ut >= TT->value[test.ilat][test.ilon][test.iz]){
      ut=0;
      ddt2=0;
    }else{
      ddt2=1.0/pow(coords[0]*dtheta, 2);
    }
    if (up >= TT->value[test.ilat][test.ilon][test.iz]){
      up=0;
      ddp2=0;
    }else{
      ddp2=1.0/pow(coords[0]*sin(coords[1])*dphi, 2);
    }
    A=ddr2+ddt2+ddp2;
    if (A == 0) continue;
    B=-2*(ur*ddr2+ut*ddt2+up*ddp2);
    C=pow(ur,2)*ddr2+pow(ut,2)*ddt2+pow(up,2)*ddp2;
    //C-=1/pow(VV->value[test.ilat][test.ilon][test.iz],2);
    C-=1/pow(VV->value[hot.ilat][hot.ilon][hot.iz],2);
    test.value=(-B+std::max(sqrt(pow(B,2)-4.0*A*C),0.0))/(2.0*A);
    if (std::isnan(test.value)){
      printf("NaN (%d %d %d) - A=%f, B=%f, C=%f\n",
             test.ilat,
             test.ilon,
             test.iz,
             A,
             B,
             C);
      continue;
    }
    TT->value[test.ilat][test.ilon][test.iz]=test.value;
    update_active(test, active);
    //printf("(%d %d %d) %.2f\n", test.ilat, test.ilon, test.iz, test.value);
  }
  hpsort(active.size()-1, active);
  printf("freezing node (%.2f, %d, %d, %d); there are currently %zd active nodes\n",
         hot.value,
         hot.ilat,
         hot.ilon,
         hot.iz,
         active.size());
  frozen.push_back(hot);
}

int ttfmm ()
{
  GeoGrid VV, TT;
  std::deque<GridNode> active, frozen;
  init_VelocityGrid(&VV);
  init_TTGrid(&TT);
  init_source(&TT, active, -999.99, -999.99, -999.99);
  while (active.size() > 0){
  //for(int i=0; i<20; i++){
    update(&TT, &VV, active, frozen);
    //print_grid(TT);
  }
  //print_active(active);
  //print_grid(VV);
  dump_grid(TT);
  free_grid(&VV);
  free_grid(&TT);
	return(0);
}

int ttfmm_ucalc ()
{
	return(0);
}

int main(int argc, char *argv[]){
  printf("START main()\n");
  ttfmm();
  return(0);
}
