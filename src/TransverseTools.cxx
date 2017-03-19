#include "TransverseTools.h"

#include "TVector3.h"
#include <vector>
#include "TMath.h"
#include <iostream>

using namespace std;

TransverseTools::TransverseTools() : m_PDP_x(0.231135) , m_PDP_y(45.368069), m_PDP_z(-766.384058), m_Theta(-0.0582977560), m_XOffset(0.2486), m_YOffset(60.350), m_ZOffset(-1022.74) {
    //    m_PDP = new TVector3(m_PDP_x, m_PDP_y, m_PDP_z);
}

TransverseTools::~TransverseTools(){
    //    delete m_PDP;
}

double TransverseTools::GetDPTT(double vtx[], const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom, bool is_truth) const
{
    TVector3 * nudir = new TVector3();
    
    if(is_truth){
        nudir->SetXYZ(vtx[0],vtx[1],vtx[2]);
    }
    else{
        const TVector3 * tmp_vec = GetNuDirRec(vtx);
        if(tmp_vec) nudir->SetXYZ(tmp_vec->X(),tmp_vec->Y(),tmp_vec->Z());
        else {
            cout << " GetNuDirRec returned NULL " << endl;
            return -999.;
        }
    }
    
    TVector3 tmp1_vec = nudir->Cross(*mumom);
    tmp1_vec *= 1/tmp1_vec.Mag();
    
    TVector3 sum_vec = *prmom + *pimom;
    
    return sum_vec.Dot(tmp1_vec);
}

double TransverseTools::GetDPTT(std::vector<double> vtx, const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom, bool is_truth) const
{
    double vertex[3] = { vtx[0], vtx[1], vtx[2] };
    return GetDPTT(vertex, mumom, prmom, pimom, is_truth);
}

TVector3 * TransverseTools::SetDPT(const TVector3 *& ptmuon, const TVector3 *& ptproton, const TVector3 *& ptpion) const
{
    //ptmuon and ptproton already in the same plain which is perpendicular to the neutrino and already in a near back-to-back configuration
    TVector3 tmpd = (*ptmuon) + (*ptproton) + (*ptpion);
    TVector3 tmp_had = (*ptproton) + (*ptpion);
    
    double phi = TMath::ACos( ptmuon->Dot(tmp_had)*(-1)/(ptmuon->Mag()*tmp_had.Mag()) );
    
    double theta = TMath::ACos( tmpd.Dot(*ptmuon)*(-1)/(tmpd.Mag()*ptmuon->Mag())  );
    TVector3 * deltapt = new TVector3();
    deltapt->SetMagThetaPhi(tmpd.Mag(),theta, phi);
    return deltapt;
}

TVector3 * TransverseTools::GetVecT(const TVector3 *& refdir, const TVector3 *& mom) const {
    //
    //w.r.t. beam direction
    //
    if(!refdir){
        cout << "CC1P1PiAnalysis::GetVecT refdir null" << endl;
        exit(1);
    }
    
    
    TVector3 vRotated(*mom);
    vRotated.Rotate(TMath::Pi(), *refdir);
    
    TVector3 *vt = new TVector3( (*mom - vRotated)*0.5 );
    
    return vt;
}

TVector3 * TransverseTools::GetPT(double vtx[], const TVector3 *& mom, bool is_truth) const
{
    TVector3 * nudir = new TVector3();
    
    if(is_truth){
        nudir->SetXYZ(vtx[0],vtx[1],vtx[2]);
    }
    else{
        const TVector3 * tmp_vec = GetNuDirRec(vtx);
        if(tmp_vec){
            nudir->SetXYZ(tmp_vec->X(),tmp_vec->Y(),tmp_vec->Z());
            delete tmp_vec;
        }
        else{
            cout << " GetNuDirRec returned NULL " << endl;
            return 0x0;
        }
    }
    
    const TVector3 * neutrino_dir = new TVector3(nudir->X(),nudir->Y(), nudir->Z());
    
    TVector3 * pT = GetVecT(neutrino_dir, mom);
    
    delete nudir;
    delete neutrino_dir;
    
    return pT;
}

TVector3 * TransverseTools::GetTransverseVars(double vtx[], const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom, double &dpTT, double &dpTMag, double &dalphaT, double &dphiT, bool is_truth) const
{
    TVector3 * nudir = new TVector3();
    
    if(is_truth){
        nudir->SetXYZ(vtx[0],vtx[1],vtx[2]);
    }
    else{
        const TVector3 * tmp_vec = GetNuDirRec(vtx);
        if(tmp_vec) nudir->SetXYZ(tmp_vec->X(),tmp_vec->Y(),tmp_vec->Z());
        else{
            cout << " GetNuDirRec returned NULL " << endl;
            return 0x0;
        }
    }
    
    const TVector3 * neutrino_dir = new TVector3(nudir->X(),nudir->Y(), nudir->Z());
    
    const TVector3 * mupT = GetVecT(neutrino_dir, mumom);
    const TVector3 * prpT = GetVecT(neutrino_dir, prmom);
    const TVector3 * pipT = GetVecT(neutrino_dir, pimom);
    
    TVector3 * deltapt = SetDPT(mupT, prpT, pipT);
    
    dpTMag  = deltapt->Mag();
    dalphaT = (deltapt->Theta())*TMath::RadToDeg();
    dphiT   = (deltapt->Phi())*TMath::RadToDeg();
    dpTT    = GetDPTT(vtx, mumom, prmom, pimom, is_truth);
    
    //TODO: May cause seg fault:
    delete nudir;
    delete neutrino_dir;
    delete mupT;
    delete prpT;
    delete pipT;
    
    return deltapt;
}

TVector3 * TransverseTools::GetTransverseVars(std::vector<double> vtx, const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom, double &dpTT, double &dpTMag, double &dalphaT, double &dphiT, bool is_truth) const
{
    double vector[3] = { vtx[0], vtx[1], vtx[2] };
    
    return GetTransverseVars(vector, mumom, prmom, pimom, dpTT, dpTMag, dalphaT, dphiT, is_truth);
}


TVector3 * TransverseTools::GetNuDirRec(double vtx[]) const
{
    //If PDP exists make sure to rotote to minerva coords: -- May not need to do this
    
    TVector3 * PDP;
    
    if(pdp[0] != -999. || pdp[1] != -999. || pdp[2] != -999.){
        PDP->SetXYZ(pdp[0], pdp[1], pdp[2]);
    }
    else PDP->SetXYZ(m_PDP_x, m_PDP_y, m_PDP_z);
    
    TVector3 * nup1local = new TVector3(vtx[0], vtx[1], vtx[2]);
    (*nup1local) *= 0.001;//in meters (default mm)
    
    if( PDP->Mag() < EPSILON || nup1local->Mag() < EPSILON ){
        cout << "CC1P1PiAnalysis::CalcNuDir bad input " << PDP->Mag() << " " << nup1local->Mag() << endl;
        return 0x0;
    }
    
    TVector3 *nuDirCalc = new TVector3( (*nup1local) - (*PDP) );
    (*nuDirCalc) *= 1./nuDirCalc->Mag();
    
    return nuDirCalc;
}

TVector3 * TransverseTools::GlobalToLocal(const double nu_NuParentDecPoint[]) const{
    //pl output in [m]
    TVector3 pg(nu_NuParentDecPoint);
    pg *= 0.001; //[mm] -> [m]
    //exactly following convert-numi-to-minerva.py -> "although it looks like someone already did the bit for changing the labeling of the axes -- Phil" -> true, because the PDP z has the exponential distribution
    const Double_t tmpx = pg.X();
    const Double_t tmpy = pg.Y();
    const Double_t tmpz = pg.Z();
    
    const Double_t xmnv = tmpx + m_XOffset;
    const Double_t ymnv =  TMath::Cos(m_Theta)*tmpy + TMath::Sin(m_Theta)*tmpz + m_YOffset;
    const Double_t zmnv = -TMath::Sin(m_Theta)*tmpy + TMath::Cos(m_Theta)*tmpz + m_ZOffset;
    
    return new TVector3(xmnv, ymnv, zmnv);
}