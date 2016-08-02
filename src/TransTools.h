#ifndef TRANSTOOLS_H
#define TRANSTOOLS_H

#include <TVector3.h>

class TransTools {
public:
    TVector3 * GetVecT( TVector3 * refdir, TVector3 mom, bool is_truth = false);
    void GetTransverseVars(double vtx[], TVector3 mom1, TVector3 mom2, TVector3 mom3, bool is_truth = false);
    
    double GetdpT(double vtx[], TVector3 mom1, TVector3 mom2, TVector3 mom3, bool is_truth = false);
    double GetdpTT(double vtx[], TVector3 mom1, TVector3 mom2, TVector3 mom3, bool is_truth = false);
    double GetdalphapT(double vtx[], TVector3 mom1, TVector3 mom2, TVector3 mom3, bool is_truth = false);
    double GetdphiT(double vtx[], TVector3 mom1, TVector3 mom2, TVector3 mom3, bool is_truth = false);
    
    
private:
    
};


double vtx[3] = {evt->vtx[0], evt->vtx[1], evt->vtx[2]};
const TVector3 * fullNeutrino = MINERVAUtils::GetNuDirRec(vtx);

