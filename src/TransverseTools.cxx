#include "TransverseTools.h"

#include "TVector3.h"
#include <vector>
#include "TMath.h"
#include <iostream>

#ifndef EPSILON
#define EPSILON  1e-10
#endif

using namespace std;

TransverseTools::TransverseTools() : m_PDP_x(0.231135) , m_PDP_y(45.368069), m_PDP_z(-766.384058),
 m_Theta(-0.0582977560), m_XOffset(0.2486), m_YOffset(60.350), m_ZOffset(-1022.74) {
    //    m_PDP = new TVector3(m_PDP_x, m_PDP_y, m_PDP_z);
}

TransverseTools::~TransverseTools(){
    //    delete m_PDP;
}

double TransverseTools::GetDPTTRec(double vtx[], const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom)// const
{
    const TVector3 * nudir = GetNuDirRec(vtx);
    return GetDPTTBase(nudir, mumom, prmom,pimom);
}

double TransverseTools::GetDPTTRec(std::vector<double> vtx, const TVector3 *& mumom, const TVector3 *& prmom, 
const TVector3 *& pimom)// const
{
    double vertex[3] = { vtx[0], vtx[1], vtx[2] };
    return GetDPTTRec(vertex, mumom, prmom, pimom);
}

double TransverseTools::GetDPTTSim(double vtx[], double pdp[], const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom)// const
{
    const TVector3 * nudir = GetNuDirSim(vtx, pdp);
    return GetDPTTBase(nudir, mumom, prmom,pimom);
}

double TransverseTools::GetDPTTSim(std::vector<double> vtx, std::vector<double> pdp, const TVector3 *& mumom, const TVector3 *& prmom,
 const TVector3 *& pimom)// const
{
    double tmp_vtx[3] = { vtx[0], vtx[1], vtx[2] };
    double tmp_pdp[3] = { pdp[0], pdp[1], pdp[2] };
    return GetDPTTSim(tmp_vtx, tmp_pdp, mumom, prmom, pimom);
}

double TransverseTools::GetDPTTBase(const TVector3 *&nudir, const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom)// const
{
    TVector3 tmp1_vec = nudir->Cross(*mumom);
    tmp1_vec *= 1/tmp1_vec.Mag();   
    TVector3 sum_vec = *prmom + *pimom;
    return sum_vec.Dot(tmp1_vec);
}

TVector3 * TransverseTools::SetDPT(const TVector3 *& ptmuon, const TVector3 *& ptproton, const TVector3 *& ptpion)// const
{
    //ptmuon and ptproton already in the same plain which is perpendicular to the neutrino and already in a near back-to-back 
    //configuration
    TVector3 tmpd = (*ptmuon) + (*ptproton) + (*ptpion);
    TVector3 tmp_had = (*ptproton) + (*ptpion);
    
    double phi = TMath::ACos( ptmuon->Dot(tmp_had)*(-1)/(ptmuon->Mag()*tmp_had.Mag()) );
    
    double theta = TMath::ACos( tmpd.Dot(*ptmuon)*(-1)/(tmpd.Mag()*ptmuon->Mag())  );
    TVector3 * deltapt = new TVector3();
    deltapt->SetMagThetaPhi(tmpd.Mag(),theta, phi);
    return deltapt;
}

TVector3 * TransverseTools::GetVecT(const TVector3 *& refdir, const TVector3 *& mom)// const
{
    //w.r.t. beam direction   
    if(!refdir){
        cout << "TransverseTools::GetVecT refdir null" << endl;
        exit(1);
    }
    TVector3 vRotated(*mom);
    vRotated.Rotate(TMath::Pi(), *refdir);
    TVector3 *vt = new TVector3( (*mom - vRotated)*0.5 );
    return vt;
}

TVector3 * TransverseTools::GetPTRec(double vtx[], const TVector3 *& mom)// const
{    
    const TVector3 * nudir = GetNuDirRec(vtx);  
    return GetVecT(nudir, mom);
}

TVector3 * TransverseTools::GetPTRec(std::vector<double> vtx, const TVector3 *& mom)
{
    double tmp_vtx[3] = { vtx[0], vtx[1], vtx[2] };
    return GetPTRec(tmp_vtx, mom);
}

TVector3 * TransverseTools::GetPTSim(double vtx[], double pdp[], const TVector3 *& mom)// const
{    
    const TVector3 * nudir = GetNuDirSim(vtx, pdp);  
    return GetVecT(nudir, mom);
}
                                    
TVector3 * TransverseTools::GetPTSim(std::vector<double> vtx, std::vector<double> pdp, const TVector3 *& mom)
{
    double tmp_vtx[3] = { vtx[0], vtx[1], vtx[2] };
    double tmp_pdp[3] = { pdp[0], pdp[1], pdp[2] };
    return GetPTSim(tmp_vtx, tmp_pdp, mom);
}

TVector3 * TransverseTools::GetTransVarsRec(double vtx[], const TVector3 *& mumom, const TVector3 *& prmom,
 const TVector3 *& pimom, double &dpTT, double &dpTMag, double &dalphaT, double &dphiT)// const
{
    const TVector3 * nudir = GetNuDirRec(vtx);
    return GetTransVarsBase(nudir, mumom, prmom, pimom, dpTT, dpTMag, dalphaT, dphiT);
}

TVector3 * TransverseTools::GetTransVarsRec(std::vector<double> vtx, const TVector3 *& mumom, const TVector3 *& prmom,
 const TVector3 *& pimom, double &dpTT, double &dpTMag, double &dalphaT, double &dphiT)// const
{
    double vector[3] = { vtx[0], vtx[1], vtx[2] };
    return GetTransVarsRec(vector, mumom, prmom, pimom, dpTT, dpTMag, dalphaT, dphiT);
}

TVector3 * TransverseTools::GetNuDirBase(const TVector3 *& vtx, const TVector3 *& PDP)// const
{
    //If PDP exists make sure to rotote to minerva coords: -- May not need to do this
    TVector3 * nup1local = new TVector3(*vtx);
    (*nup1local) *= 0.001;//in meters (default mm)
    if( PDP->Mag() < EPSILON || nup1local->Mag() < EPSILON ){
        cout << "TransverseTools::CalcNuDir bad input " << PDP->Mag() << " " << nup1local->Mag() << endl;
        return 0x0;
    }
    TVector3 *nuDirCalc = new TVector3( (*nup1local) - (*PDP) );
    (*nuDirCalc) *= 1./nuDirCalc->Mag();
    return nuDirCalc;
}

TVector3 * TransverseTools::GetNuDirRec(double vtx[])// const
{
    //If PDP exists make sure to rotote to minerva coords: -- May not need to do this
    const TVector3 * VTX = new TVector3(vtx[0], vtx[1], vtx[2]);
    const TVector3 * PDP = new TVector3(m_PDP_x, m_PDP_y, m_PDP_z);
    return GetNuDirBase(VTX, PDP);
}

TVector3 * TransverseTools::GetNuDirRec(std::vector<double> vtx)
{    
    double tmp_vtx[3] = { vtx[0], vtx[1], vtx[2] };
    return GetNuDirRec(tmp_vtx);
}

TVector3 * TransverseTools::NuMiToMin(double nu_NuParentDecPoint[])// const
{
    // This converts from Numi (beam) coordinates to Minerva coordinates.
    //pl output in [m]
    TVector3 pg(nu_NuParentDecPoint);
    pg *= 0.001; //[mm] -> [m]
    const Double_t tmpx = pg.X();
    const Double_t tmpy = pg.Y();
    const Double_t tmpz = pg.Z();
    const Double_t xmnv = tmpx + m_XOffset;
    const Double_t ymnv =  TMath::Cos(m_Theta)*tmpy + TMath::Sin(m_Theta)*tmpz + m_YOffset;
    const Double_t zmnv = -TMath::Sin(m_Theta)*tmpy + TMath::Cos(m_Theta)*tmpz + m_ZOffset;
    return new TVector3(xmnv, ymnv, zmnv);
}

TVector3 * TransverseTools::GetNuDirSim(double vtx[], double pdp[])// const
{
    // May need to delete these after use? 
    const TVector3 * PDP = NuMiToMin(pdp);
    const TVector3 * VTX = new TVector3(vtx[0], vtx[1], vtx[2]);
    return GetNuDirBase(VTX, PDP);
}

TVector3 * TransverseTools::GetNuDirSim(std::vector<double> vtx, std::vector<double> pdp)// const
{   
    double tmp_vtx[3] = { vtx[0], vtx[1], vtx[2] };
    double tmp_pdp[3] = { pdp[0], pdp[1], pdp[2] };   
    return GetNuDirSim(tmp_vtx, tmp_pdp);
}

TVector3 * TransverseTools::GetTransVarsSim(double vtx[], double pdp[], const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom, double &dpTT,
                                 double &dpTMag, double &dalphaT, double &dphiT)// const;
{
    const TVector3 * nudir = GetNuDirSim(vtx, pdp);
    return GetTransVarsBase(nudir, mumom, prmom, pimom, dpTT, dpTMag, dalphaT, dphiT);
}

TVector3 * TransverseTools::GetTransVarsSim(std::vector<double> vtx, std::vector<double> pdp, const TVector3 *& mumom, const TVector3 *& prmom, 
    const TVector3 *& pimom, double &dpTT, double &dpTMag, double &dalphaT, double &dphiT)// const;
{
    double tmp_vtx[3] = { vtx[0], vtx[1], vtx[2] };
    double tmp_pdp[3] = { pdp[0], pdp[1], pdp[2] };
    return GetTransVarsSim(tmp_vtx, tmp_pdp, mumom, prmom, pimom, dpTT, dpTMag, dalphaT, dphiT);
}

TVector3 * TransverseTools::GetTransVarsBase(const TVector3 *& nudir, const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom, double &dpTT,
                                 double &dpTMag, double &dalphaT, double &dphiT)
{

    const TVector3 * mupT = GetVecT(nudir, mumom);
    const TVector3 * prpT = GetVecT(nudir, prmom);
    const TVector3 * pipT = GetVecT(nudir, pimom);
    
    TVector3 * deltapt = SetDPT(mupT, prpT, pipT);
    
    dpTMag  = deltapt->Mag();
    dalphaT = (deltapt->Theta())*TMath::RadToDeg();
    dphiT   = (deltapt->Phi())*TMath::RadToDeg();
    dpTT    = GetDPTTBase(nudir, mumom, prmom, pimom);
    
    //TODO: May cause seg fault:
    delete nudir;
    delete mupT;
    delete prpT;
    delete pipT;
    
    return deltapt;
}
