//
//  TMRSDarcyFlowWithMem.h
//  LinearTracer
//
//  Created by Omar Durán on 10/10/19.
//

#ifndef TMRSDarcyFractureFlowWithMem_h
#define TMRSDarcyFractureFlowWithMem_h

#include <stdio.h>
#include "TPZMatWithMem.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "TMRSDataTransfer.h"
#include "TMRSDarcyFlowWithMem.h"

template <class TMEM>
class TMRSDarcyFractureFlowWithMem : public TMRSDarcyFlowWithMem<TMEM> {
    
    
public:
    
    /// Default constructor
    TMRSDarcyFractureFlowWithMem();
    
    /// Constructor based on a material id
    TMRSDarcyFractureFlowWithMem(int mat_id, int dimension);
    
    /// Constructor based on a TPBrMatMixedDarcy object
    TMRSDarcyFractureFlowWithMem(const TMRSDarcyFractureFlowWithMem & other);
    
    /// Constructor based on a TPBrMatMixedDarcy object
    TMRSDarcyFractureFlowWithMem &operator=(const TMRSDarcyFractureFlowWithMem & other);
    
    /// Default destructor
    ~TMRSDarcyFractureFlowWithMem();
    
    /// Set the required data at each integration point
    void FillDataRequirements(TPZVec<TPZMaterialData> &datavec) override;
    
    /// Set the required data at each integration point
    void FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec) override;
    
    /// Returns the name of the material
    std::string Name() override {
        return "TMRSDarcyFlowWithMem";
    }
    
    /// Returns the integrable dimension of the material */
    int Dimension() const override {return this->m_dimension;}
    
    /// Returns the number of state variables associated with the material
    int NStateVariables() const override {return 1;}
    
    virtual TPZMaterial *NewMaterial() override
    {
        return new TMRSDarcyFractureFlowWithMem<TMEM>(*this);
    }
    
    /// Set data transfer object
    void SetDataTransfer(TMRSDataTransfer & SimData);
    
    /// Print out the data associated with the material
    void Print(std::ostream &out = std::cout) override;
    
    /// Returns the variable index associated with the name
    int VariableIndex(const std::string &name) override;
    
    /// returns the number of variables associated with the variable indexed by var.
    int NSolutionVariables(int var) override;
    
    /// Returns the solution associated with the var index based on a finite element approximation (Used for TPZPostProcAnalysis)
    void Solution(TPZMaterialData &datavec, int var, TPZVec<REAL> &Solout) override {
        DebugStop();
    }
    
    /// Returns the solution associated with the var index based on a finite element approximation
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) override;
    
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
    
    void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc){
        DebugStop();
    }
    
    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }
    
    void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }
    
    void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }
    
    void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef){
        DebugStop();
    }
    
    // Contribute Methods being used
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void ContributeFourSpaces(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
};

#endif /* TMRSDarcyFlowWithMem_h */
