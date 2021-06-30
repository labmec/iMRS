//
//  TMRSDarcyFlowWithMem.cpp
//  LinearTracer
//
//  Created by Omar Durán on 10/10/19.

#include "TMRSMemory.h"

#include "TMRSDarcyFractureFlowWithMem.h"

template <class TMEM>
TMRSDarcyFractureFlowWithMem<TMEM>::TMRSDarcyFractureFlowWithMem() : TBase() {
}

template <class TMEM>
TMRSDarcyFractureFlowWithMem<TMEM>::TMRSDarcyFractureFlowWithMem(int mat_id, int dimension) : TBase(mat_id, dimension) {
    
    
}

template <class TMEM>
TMRSDarcyFractureFlowWithMem<TMEM>::TMRSDarcyFractureFlowWithMem(const TMRSDarcyFractureFlowWithMem & other) : TMRSDarcyFlowWithMem<TMEM>(other){
}

template <class TMEM>
TMRSDarcyFractureFlowWithMem<TMEM> & TMRSDarcyFractureFlowWithMem<TMEM>::operator=(const TMRSDarcyFractureFlowWithMem & other){
    // check for self-assignment
    if(&other == this){
        return *this;
    }
    return *this;
}

template <class TMEM>
TMRSDarcyFractureFlowWithMem<TMEM>::~TMRSDarcyFractureFlowWithMem(){
    
}

template <class TMEM>
void TMRSDarcyFractureFlowWithMem<TMEM>::FillDataRequirements( TPZVec<TPZMaterialDataT<STATE>> &datavec) const{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        //datavec[idata].fNormalVec = true;
    }
}

template <class TMEM>
void TMRSDarcyFractureFlowWithMem<TMEM>::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fDeformedDirections = true;
    }
}

template <class TMEM>
void TMRSDarcyFractureFlowWithMem<TMEM>::SetDataTransfer(TMRSDataTransfer & SimData){
    this->mSimData = SimData;
}

template <class TMEM>
void TMRSDarcyFractureFlowWithMem<TMEM>::Print(std::ostream &out) const{
    TPZMaterial::Print(out);
}

template <class TMEM>
int TMRSDarcyFractureFlowWithMem<TMEM>::VariableIndex(const std::string &name) const{
    if(!strcmp("Flux",name.c_str()))            return  1;
    if(!strcmp("Pressure",name.c_str()))        return  2;
    if(!strcmp("div_q",name.c_str()))           return  3;
    if(!strcmp("kappa",name.c_str()))           return  4;
    if(!strcmp("g_average",name.c_str()))        return  5;
    if(!strcmp("p_average",name.c_str()))        return  6;
    return TPZMaterial::VariableIndex(name);
}

template <class TMEM>
int TMRSDarcyFractureFlowWithMem<TMEM>::NSolutionVariables(int var) const{
    if(var == 1) return 3;
    if(var == 2) return 1;
    if(var == 3) return 1;
    if(var == 4) return 1;
    if(var == 5) return 1;
    if(var == 6) return 1;
    return TBase::NSolutionVariables(var);
}

template <class TMEM>
void TMRSDarcyFractureFlowWithMem<TMEM>::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<REAL> &Solout) {
    TMRSDarcyFlowWithMem<TMEM>::Solution(datavec,var,Solout);
}

template <class TMEM>
void TMRSDarcyFractureFlowWithMem<TMEM>::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    int qb = 0;
    int pb = 1;
    int sb = 2;
    
    
    //    datavec[qb].Print(std::cout);
    TPZFNMatrix<100,REAL> phi_qs       = datavec[qb].phi;
    TPZFNMatrix<100,REAL> phi_ps       = datavec[pb].phi;
    TPZFNMatrix<300,REAL> dphi_qs      = datavec[qb].dphix;
    TPZFNMatrix<100,REAL> dphi_ps      = datavec[pb].dphix;
    
    
    TPZFNMatrix<40, REAL> div_phi = datavec[qb].divphi;
    REAL div_q = datavec[qb].divsol[0][0];
    
    int nphi_q       = datavec[qb].fVecShapeIndex.NElements();
    int nphi_p       = phi_ps.Rows();
    int first_q      = 0;
    int first_p      = nphi_q + first_q;
    int nvecs = datavec[qb].fDeformedDirections.Cols();
    
    // first index in fVecShapeIndex corresponding to the first transverse flux
    int first_transverse_q = 0;
    // first index in fVecShapeIndex corresponding to the second transverse flux
    int second_transverse_q = 0;
    for(int i=0; i< nphi_q; i++)
    {
        if(first_transverse_q == 0 && datavec[qb].fVecShapeIndex[i].first == nvecs-2) first_transverse_q = i;
        if(second_transverse_q == 0 && datavec[qb].fVecShapeIndex[i].first == nvecs-1) second_transverse_q = i;
    }
    if(first_transverse_q == 0 || second_transverse_q == 0 || first_transverse_q == second_transverse_q)
    {
        DebugStop();
    }
    
    
    TPZManVector<STATE,3> q  = datavec[qb].sol[0];
    STATE p                  = datavec[pb].sol[0][0];
    STATE sw                 = datavec[sb].sol[0][0];
    
    // Get the data at integrations points
    long gp_index = datavec[qb].intGlobPtIndex;
    TMEM & memory = this->GetMemory().get()->operator[](gp_index);
    
    TPZFNMatrix<3,STATE> phi_q_i(3,1,0.0), kappa_inv_phi_q_j(3,1,0.0), kappa_inv_q(3,1,0.0),kappa_inv_qFrac(3,1,0.0) ;
    
    TPZFNMatrix<9,STATE>  kappaInv(3,3,0.0);
    kappaInv=memory.m_kappa_inv;
    REAL kappaNormal = 1.0;
    REAL ad =1.0;
    REAL eps = ad;
    REAL factor = 1.0;
    
    int s_i, s_j;
    int v_i, v_j;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            kappa_inv_q(i,0) += memory.m_kappa_inv(i,j)*q[j];
        }
    }
    
   
   
    
    for (int iq = 0; iq < first_transverse_q; iq++)
    {
        //        datavec[qb].Print(std::cou        t);
        v_i = datavec[qb].fVecShapeIndex[iq].first;
        s_i = datavec[qb].fVecShapeIndex[iq].second;
        STATE kappa_inv_q_dot_phi_q_i = 0.0;
        for (int i = 0; i < 3; i++) {
            phi_q_i(i,0) = phi_qs(s_i,0) * datavec[qb].fDeformedDirections(i,v_i);
            kappa_inv_q_dot_phi_q_i        += kappa_inv_q(i,0)*phi_q_i(i,0);
        }
        
        ef(iq + first_q) += weight * ( kappa_inv_q_dot_phi_q_i - p * div_phi(iq,0));
        
        for (int jq = 0; jq < first_transverse_q; jq++)
        {
            
            v_j = datavec[qb].fVecShapeIndex[jq].first;
            s_j = datavec[qb].fVecShapeIndex[jq].second;
//            if(v_j < nvecs-2){
                kappa_inv_phi_q_j.Zero();
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        kappa_inv_phi_q_j(i,0) += (1.0/eps)* memory.m_kappa_inv(i,j) * phi_qs(s_j,0) * datavec[qb].fDeformedDirections(j,v_j);
                    }
                }
                
                STATE kappa_inv_phi_q_j_dot_phi_q_i = 0.0;
                for (int j = 0; j < 3; j++) {
                    kappa_inv_phi_q_j_dot_phi_q_i +=  kappa_inv_phi_q_j(j,0)*phi_q_i(j,0);
                }
                
                ek(iq + first_q,jq + first_q) += weight * kappa_inv_phi_q_j_dot_phi_q_i;
            }
            if(v_j >=  nvecs-2){
                DebugStop();
            }
            for (int jp = 0; jp < nphi_p; jp++)
            {
                ek(iq + first_q, jp + first_p) += weight * ( - div_phi(iq,0) ) * phi_ps(jp,0);
                ek(jp + first_p, iq + first_q) += weight * ( - div_phi(iq,0) ) * phi_ps(jp,0);
               
            }
        
    }
//    
    // compute the contribution of the hdivbound
    for (int iq = first_transverse_q; iq < second_transverse_q; iq++)
    {
        //        datavec[qb].Print(std::cou        t);
        v_i = datavec[qb].fVecShapeIndex[iq].first;
        s_i = datavec[qb].fVecShapeIndex[iq].second;
        STATE kappa_inv_q_dot_phi_q_i = 0.0;

        for (int i = 0; i < 3; i++) {
            kappa_inv_q(i,0) += (1.0/kappaNormal)*q[i];
        }
        
        for (int i = 0; i < 3; i++) {
            phi_q_i(i,0) = phi_qs(s_i,0) * datavec[qb].fDeformedDirections(i,v_i);
            kappa_inv_q_dot_phi_q_i        += kappa_inv_q(i,0)*phi_q_i(i,0);
        }

//        ef(iq + first_q) += weight * ( kappa_inv_q_dot_phi_q_i - p * div_phi(iq,0));
        for (int jq = first_transverse_q; jq < second_transverse_q; jq++)
        {

            v_j = datavec[qb].fVecShapeIndex[jq].first;
            s_j = datavec[qb].fVecShapeIndex[jq].second;
            if(v_j != v_i) DebugStop();
            kappa_inv_phi_q_j.Zero();

            // kappanormal is the orthogonal permeability
            // phi_qs is the scalar function corresponding to jq (scalar)
            // s_j is the index of the scalar function corresponding to jq
            // v_j is the index of the vector corresponding to jq
            // kappa_inv_phi_q_j = vector equal to 1/kappa phi_qs * phi_qs * vector(v_j)
            for (int j = 0; j < 3; j++) {
                REAL KappaInvVal = (1.0/kappaNormal);
                kappa_inv_phi_q_j(j,0) = (ad*KappaInvVal/2.0)* factor * phi_qs(s_j,0) * datavec[qb].fDeformedDirections(j,v_j);
            }
            

            // kappa_inv_phi_q_j_dot_phi_q_i is the dot product of both vectors
            STATE kappa_inv_phi_q_j_dot_phi_q_i = 0.0;
            for (int j = 0; j < 3; j++) {
                kappa_inv_phi_q_j_dot_phi_q_i += kappa_inv_phi_q_j(j,0)*phi_q_i(j,0);
            }
//            ek(iq + first_q,jq + first_q) += weight * kappa_inv_phi_q_j_dot_phi_q_i;
        }
        for (int jp = 0; jp < nphi_p; jp++)
        {
            ek(iq + first_q, jp + first_p) += weight * ( - div_phi(iq,0) ) * phi_ps(jp,0);
            ek(jp + first_p, iq + first_q) += weight * ( - div_phi(iq,0) ) * phi_ps(jp,0);

        }
    }
//    // compute the contribution of the hdivbound
    for (int iq = second_transverse_q; iq < nphi_q; iq++)
    {
        //        datavec[qb].Print(std::cou        t);
        v_i = datavec[qb].fVecShapeIndex[iq].first;
        s_i = datavec[qb].fVecShapeIndex[iq].second;
        STATE kappa_inv_q_dot_phi_q_i = 0.0;

        for (int i = 0; i < 3; i++) {
            kappa_inv_q(i,0) += (1.0/kappaNormal)*q[i];
        }
        
        for (int i = 0; i < 3; i++) {
            phi_q_i(i,0) = phi_qs(s_i,0) * datavec[qb].fDeformedDirections(i,v_i);
            kappa_inv_q_dot_phi_q_i        += kappa_inv_q(i,0)*phi_q_i(i,0);
        }
        
//        ef(iq + first_q) += weight * ( kappa_inv_q_dot_phi_q_i - p * div_phi(iq,0));
        for (int jq = second_transverse_q; jq < nphi_q; jq++)
        {
            v_j = datavec[qb].fVecShapeIndex[jq].first;
            s_j = datavec[qb].fVecShapeIndex[jq].second;
            if(v_j != v_i) DebugStop();
            kappa_inv_phi_q_j.Zero();

            // kappanormal is the orthogonal permeability
            // phi_qs is the scalar function corresponding to jq (scalar)
            // s_j is the index of the scalar function corresponding to jq
            // v_j is the index of the vector corresponding to jq
            // kappa_inv_phi_q_j = vector equal to 1/kappa phi_qs * phi_qs * vector(v_j)
            for (int j = 0; j < 3; j++) {
                REAL KappaInvVal = (1.0/kappaNormal);
                kappa_inv_phi_q_j(j,0) += (ad*KappaInvVal/2)* factor *phi_qs(s_j,0) * datavec[qb].fDeformedDirections(j,v_j);
            }
            
            // kappa_inv_phi_q_j_dot_phi_q_i is the dot product of both vectors
            STATE kappa_inv_phi_q_j_dot_phi_q_i = 0.0;
            for (int j = 0; j < 3; j++) {
                kappa_inv_phi_q_j_dot_phi_q_i += kappa_inv_phi_q_j(j,0)*phi_q_i(j,0);
            }
//            ek(iq + first_q,jq + first_q) += weight * kappa_inv_phi_q_j_dot_phi_q_i;

        }
        for (int jp = 0; jp < nphi_p; jp++)
        {
            ek(iq + first_q, jp + first_p) += weight * (-  div_phi(iq,0) ) * phi_ps(jp,0);
            ek(jp + first_p, iq + first_q) += weight * ( - div_phi(iq,0) ) * phi_ps(jp,0);

        }
    }
    
    
    if(this->mSimData.mTNumerics.m_four_approx_spaces_Q){
        ContributeFourSpaces(datavec,weight,ek,ef);
    }
}

template <class TMEM>
void TMRSDarcyFractureFlowWithMem<TMEM>::ContributeFourSpaces(const TPZVec<TPZMaterialDataT<STATE>> &datavec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
    
    int qb = 0;
    int pb = 1;
    int g_avgb = 2;
    int p_avgb = 3;
    
    TPZFNMatrix<100,REAL> phi_ps       = datavec[pb].phi;
    int nphi_q       = datavec[qb].fVecShapeIndex.NElements();
    int nphi_p       = phi_ps.Rows();
    
    int nphi_gb = datavec[g_avgb].phi.Rows();
    int nphi_pb = datavec[p_avgb].phi.Rows();
    if(nphi_q+nphi_p+nphi_gb+nphi_pb != ek.Rows())
    {
        DebugStop();
    }
    
    STATE p     = datavec[pb].sol[0][0];
    STATE g_avg = datavec[g_avgb].sol[0][0];
    STATE p_avg = datavec[p_avgb].sol[0][0];
    
    for(int ip=0; ip<nphi_p; ip++)
    {
        ef(nphi_q+ip,0) += weight * g_avg * phi_ps(ip,0);
        ek(nphi_q+ip,nphi_q+nphi_p) += weight * phi_ps(ip,0);
        
        ek(nphi_q+nphi_p,nphi_q+ip) += weight * phi_ps(ip,0);
    }
    
    ef(nphi_q+nphi_p+1,0) += -weight * g_avg;
    ek(nphi_q+nphi_p+1,nphi_q+nphi_p) += -weight;
    
    ef(nphi_q+nphi_p,0) += weight * (p - p_avg);
    ek(nphi_q+nphi_p,nphi_q+nphi_p+1) += -weight;
    
}

template <class TMEM>
void TMRSDarcyFractureFlowWithMem<TMEM>::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    TPZFMatrix<STATE> ekfake(ef.Rows(),ef.Rows(),0.0);
    
//    this->Contribute(datavec, weight, ekfake, ef);
    this->Contribute(datavec, weight, ef);
    
    if(TMRSDarcyFractureFlowWithMem<TMEM>::fUpdateMem){
        int qb = 0;
        int pb = 1;
        long gp_index = datavec[pb].intGlobPtIndex;
        TMEM & memory = this->GetMemory().get()->operator[](gp_index);
        TPZVec<REAL> q_n = datavec[qb].sol[0][0];
        memory.m_flux = q_n;
        
        REAL p_n = datavec[pb].sol[0][0];
        memory.m_p = p_n;
    }
    
}

template <class TMEM>
void TMRSDarcyFractureFlowWithMem<TMEM>::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc){
    TPZFMatrix<STATE> ekfake(ef.Rows(),ef.Rows(),0.0);
    this->ContributeBC(datavec, weight, ekfake, ef, bc);
}

template <class TMEM>
void TMRSDarcyFractureFlowWithMem<TMEM>::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc){
    
    
    REAL gBigNumber = 1.0e12; //TPZMaterial::gBigNumber;
    
    int qb = 0;
    TPZFNMatrix<100,REAL> phi_qs       = datavec[qb].phi;
    
    int nphi_q       = phi_qs.Rows();
    int first_q      = 0;
    
    TPZManVector<STATE,3> q  = datavec[qb].sol[0];
    
    TPZManVector<STATE,1> bc_data(1,0.0);
    bc_data[0] = bc.Val2()[0];
    
    switch (bc.Type()) {
        case 0 :    // Dirichlet BC  PD
        {
            STATE p_D = bc_data[0];
            for (int iq = 0; iq < nphi_q; iq++)
            {
                ef(iq + first_q) += weight * p_D * phi_qs(iq,0);
            }
        }
            break;
            
        case 1 :    // Neumann BC  QN
        {
            
            for (int iq = 0; iq < nphi_q; iq++)
            {
                REAL qn_N = bc_data[0];
                REAL qn = 0.0;
                qn = q[0];
                
                ef(iq + first_q) += weight * gBigNumber * (qn - qn_N) * phi_qs(iq,0);
                for (int jq = 0; jq < nphi_q; jq++)
                {
                    ek(iq + first_q,jq + first_q) += weight * gBigNumber * phi_qs(jq,0) * phi_qs(iq,0);
                }
                
            }
            
        }
            break;
            
        default: std::cout << "This BC doesn't exist." << std::endl;
        {
            
            DebugStop();
        }
            break;
    }
    
    return;
    
}

template class TMRSDarcyFractureFlowWithMem<TMRSMemory>;

























