#ifndef __TRANSVERSETOOLS__H
#define __TRANSVERSETOOLS__H

#include "TROOT.h"
//#include <vector>//?

class TVector3;
class vector;//Is this good?

class TransverseTools {
    
public:
    TransverseTools();
    ~TransverseTools();
    
    TVector3 * GetNuDirRec(double vtx[]);// const;
    
    TVector3 * GetNuDirSim(double vtx[], double pdp[]);// const;
    TVector3 * GetNuDirSim(std::vector<double> vtx, std::vector<double> pdp);// const;
    
    double GetDPTT(double vtx[], const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom, bool is_truth = false);// const;
    double GetDPTT(std::vector<double> vtx, const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom, bool is_truth = false);// const;
    
    TVector3 * GetVecT(const TVector3 *& refdir, const TVector3 *& mom);// const;
    
    TVector3 * GetPT(double vtx[], const TVector3 *& mom, bool is_truth = false);// const;
    TVector3 * GetPT(std::vector<double> vtx, const TVector3 *& mom, bool is_truth = false);// const;

    TVector3 * SetDPT(const TVector3 *& ptmuon, const TVector3 *& ptproton, const TVector3 *& ptpion);// const;
    
    TVector3 * GetTransverseVars(double vtx[], const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom, double &dpTT,
                                 double &dpTMag, double &dalphaT, double &dphiT, bool is_truth = false);// const;
    
    TVector3 * GetTransverseVars(std::vector<double> vtx, const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom,
                                 double &dpTT, double &dpTMag, double &dalphaT, double &dphiT, bool is_truth = false);// const;
    
private:
    const Double_t m_PDP_x;
    const Double_t m_PDP_y;
    const Double_t m_PDP_z;
    
//    TVector3 * m_PDP;
    
    const Double_t m_Theta;// = -0.0582977560;
    const Double_t m_XOffset;// = 0.2486;
    const Double_t m_YOffset;// = 60.350;
    const Double_t m_ZOffset;// = -1022.74;
    
    TVector3 * GlobalToLocal(double nu_NuParentDecPoint[]);// const;

};

#endif
