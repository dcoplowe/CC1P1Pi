#ifndef _MINERVAUTILS_H_
#define _MINERVAUTILS_H_

#ifndef EPSILON
#define EPSILON  1e-10
#endif

class MINERVAUtils
{
 public:
  static const TVector3 * GetNuDirRec(const Double_t vertex_pos[]);

  static void PrintPDP(){cout<<"\nMINERVAUtils::PrintPDP v10r8 22/Jan/2016 PDP_XYZ=("<<GetMeanPDP_X()<<", "<<GetMeanPDP_Y()<<", "<<GetMeanPDP_Z()<<")\n"<<endl;}
  //== independent of NeutrinoUtils
  //from treePDP.C r4522 1M events from mc_minerva.txt Truth tree, weighted by mc_cvweight_total, in [m], v10r8p8: 0.2119, 45.58, -770.4, 20/Jan/2016
  //from /minerva/app/users/xlu/cmtuser/Minerva_v10r8/Ana/NukeCCQE/twoTrack/ana/include/playlists/playlist_minerva1/mc_minerva.txt 1M events, weighted, r4553, v10r8:   0.231135, 45.368069, -766.384058, 22/Jan/2016
  static Double_t GetMeanPDP_X(){return 0.231135;}
  static Double_t GetMeanPDP_Y(){return 45.368069;}
  static Double_t GetMeanPDP_Z(){return -766.384058;}

 private:
  static const TVector3 * GetNuDirSim(const Double_t nu_NuParentDecPoint[], const Double_t vertex_truepos[], const Float_t nu_truedir[], const TVector3 * savepdp=0x0);

  static void SetNeutrinoParentDecPointRec(TVector3 *vec){ vec->SetXYZ( GetMeanPDP_X(), GetMeanPDP_Y(),  GetMeanPDP_Z() ); }

  static void GlobalToLocal(const Double_t nu_NuParentDecPoint[], TVector3 * pl);

  static const TVector3 * CalcNuDir(const TVector3 * nup0Local, const TVector3 * nup1Local);
};


const TVector3 * MINERVAUtils::GetNuDirSim(const Double_t nu_NuParentDecPoint[], const Double_t vertex_truepos[], const Float_t nu_truedir[], const TVector3 * savepdp)
{
  //-------- nup0 
  //in m
  TVector3 tmp0;
  GlobalToLocal(nu_NuParentDecPoint, &tmp0);  

  if(savepdp){
    delete savepdp;
    savepdp = new TVector3(tmp0);  
  }

  //------ nup1
  //in m
  TVector3 nup1Local(vertex_truepos);
  nup1Local *= 0.001; //default mm

  const TVector3 * dirsim = CalcNuDir(&tmp0, &nup1Local);

  //test ->
  TVector3 nuDirTrue(nu_truedir);
  nuDirTrue *= 1./nuDirTrue.Mag();
  
  const Double_t ddir=1-dirsim->Dot(nuDirTrue);
  //pass 1e-6, fail at 1E-7, tested with mc_minerva_no_fsi_plastic_target.root
  //fail at 1E-6 with mc_minerva_plastic_target.root
  if(ddir>1e-5){
    printf("test bad!! %e\n", ddir);
    dirsim->Print();
    nuDirTrue.Print();
    exit(1);
  }
  //test <-
  
  return dirsim;
}

void MINERVAUtils::GlobalToLocal(const Double_t nu_NuParentDecPoint[], TVector3 * pl)
{
  //pl output in [m]

  TVector3 pg(nu_NuParentDecPoint);
  pg *= 0.001; //[mm] -> [m]

  //exactly following convert-numi-to-minerva.py -> "although it looks like someone already did the bit for changing the labeling of the axes -- Phil" -> true, because the PDP z has the exponential distribution
  const Double_t tmpx = pg.X();
  const Double_t tmpy = pg.Y();
  const Double_t tmpz = pg.Z();

  const Double_t theta=-0.0582977560;
  //in m
  const Double_t xoffset=0.2486;
  const Double_t yoffset=60.350;
  const Double_t zoffset=-1022.74;

  const Double_t xmnv=tmpx+xoffset;
  const Double_t ymnv=TMath::Cos(theta)*tmpy+TMath::Sin(theta)*tmpz+yoffset;
  const Double_t zmnv=-TMath::Sin(theta)*tmpy+TMath::Cos(theta)*tmpz+zoffset;

  pl->SetXYZ(xmnv, ymnv, zmnv);
}

const TVector3 * MINERVAUtils::CalcNuDir(const TVector3 * nup0Local, const TVector3 * nup1Local)
{
  if( nup0Local->Mag()<EPSILON || nup1Local->Mag()<EPSILON ){
    printf("MINERVAUtils::CalcNuDir bad input %f %f\n", nup0Local->Mag(), nup1Local->Mag()); exit(1);
    return 0x0;
  }

  //=====================================

  TVector3 *nuDirCalc = new TVector3( (*nup1Local) - (*nup0Local) );
  (*nuDirCalc) *= 1./nuDirCalc->Mag();

  return nuDirCalc;
}

const TVector3 * MINERVAUtils::GetNuDirRec(const Double_t vertex_pos[])
{
  //-------- nup0 
  //in m, mean of the distribution, 6B neutrino flux
  TVector3 nup0Local;
  SetNeutrinoParentDecPointRec(&nup0Local);
  
  //------ nup1
  //in m
  TVector3 nup1Local(vertex_pos);
  nup1Local *= 0.001; //default mm

  return CalcNuDir(&nup0Local, &nup1Local);
}

#endif
