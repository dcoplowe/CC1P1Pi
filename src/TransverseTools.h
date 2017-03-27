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
    TVector3 * GetNuDirRec(std::vector<double> vtx);// const;
    
    TVector3 * GetNuDirSim(double vtx[], double pdp[]);// const;
    TVector3 * GetNuDirSim(std::vector<double> vtx, std::vector<double> pdp);// const;
    
    double GetDPTTRec(double vtx[], const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom);// const;
    double GetDPTTRec(std::vector<double> vtx, const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom);// const;

    double GetDPTTSim(double vtx[], double pdp[], const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom);
    double GetDPTTSim(std::vector<double> vtx, std::vector<double> pdp, const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom);
    
    double GetDPTTDir(double dir[], const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom);
    double GetDPTTDir(std::vector<double> dir, const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom);

    TVector3 * GetVecT(const TVector3 *& refdir, const TVector3 *& mom);// const;
    
    TVector3 * GetPTSim(double vtx[], double pdp[], const TVector3 *& mom);// const;
    TVector3 * GetPTSim(std::vector<double> vtx, std::vector<double> pdp, const TVector3 *& mom);// const;

    TVector3 * GetPTRec(double vtx[], const TVector3 *& mom);// const;
    TVector3 * GetPTRec(std::vector<double> vtx, const TVector3 *& mom);// const;

    TVector3 * GetPTDir(double dir[], const TVector3 *& mom);// const;
    TVector3 * GetPTDir(std::vector<double> dir, const TVector3 *& mom);// const;

    TVector3 * SetDPT(const TVector3 *& ptmuon, const TVector3 *& ptproton, const TVector3 *& ptpion);// const;
    
    TVector3 * GetTransVarsRec(double vtx[], const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom, double &dpTT,
                                 double &dpTMag, double &dalphaT, double &dphiT);// const;
    
    TVector3 * GetTransVarsRec(std::vector<double> vtx, const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom,
                                 double &dpTT, double &dpTMag, double &dalphaT, double &dphiT);// const;

    TVector3 * GetTransVarsDir(double dir[], const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom, double &dpTT,
                                 double &dpTMag, double &dalphaT, double &dphiT);// const;
    
    TVector3 * GetTransVarsDir(std::vector<double> dir, const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom,
                                 double &dpTT, double &dpTMag, double &dalphaT, double &dphiT);// const;

    //Use these for checking reco:
    TVector3 * GetTransVarsSim(double vtx[], double pdp[], const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom, double &dpTT,
                                 double &dpTMag, double &dalphaT, double &dphiT);// const;
    
    TVector3 * GetTransVarsSim(std::vector<double> vtx, std::vector<double> pdp, const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom,
                                 double &dpTT, double &dpTMag, double &dalphaT, double &dphiT);// const;

    void RotateToNuMi(double &y, double &z);
    
private:
    const Double_t m_PDP_x;
    const Double_t m_PDP_y;
    const Double_t m_PDP_z;
        
    const Double_t m_Theta;// = -0.0582977560;
    const Double_t m_XOffset;// = 0.2486;
    const Double_t m_YOffset;// = 60.350;
    const Double_t m_ZOffset;// = -1022.74;
    
    TVector3 * NuMiToMin(double nu_NuParentDecPoint[]);
    TVector3 * GetNuDirBase(const TVector3 *& vtx, const TVector3 *& PDP);
    TVector3 * GetTransVarsBase(const TVector3 *& nudir, const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom, double &dpTT,
                                 double &dpTMag, double &dalphaT, double &dphiT);
    double GetDPTTBase(const TVector3 *&nudir, const TVector3 *& mumom, const TVector3 *& prmom, const TVector3 *& pimom);
};

#endif
