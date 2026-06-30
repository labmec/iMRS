//
//  TMRSTransportAnalysis.cpp
//
//  Created by Omar Durán on 10/15/19.
//

#include "TMRSTransportAnalysis.h"
#include "pzfunction.h"
#include "TPZTracerFlow.h"
#include "pzvtkmesh.h"

// Uses the new vtk function developed by Fran
#define USENEWVTK

#ifdef USENEWVTK
#include "TPZVTKGenerator.h"
#endif
#include "TPZVTKGeoMesh.h"

#ifdef USING_TBB
#include <tbb/parallel_for.h>
#endif
#include "TPZSpStructMatrix_Eigen.h"
#include "TPZSpMatrixEigen.h"
TMRSTransportAnalysis::TMRSTransportAnalysis(){
    
}

TMRSTransportAnalysis::~TMRSTransportAnalysis(){
    
}

TMRSTransportAnalysis::TMRSTransportAnalysis(TPZCompMesh * cmesh_mult,
                                             const RenumType& renumtype) : TPZLinearAnalysis(cmesh_mult, renumtype){
    
    TPZMultiphysicsCompMesh *cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(cmesh_mult);
    if(cmesh){
        m_soltransportTransfer.BuildTransferData(cmesh_mult);
    }
}

void TMRSTransportAnalysis::SetDataTransfer(TMRSDataTransfer * sim_data){
    m_sim_data = sim_data;
}

TMRSDataTransfer * TMRSTransportAnalysis::GetDataTransfer(){
    return m_sim_data;
}

int TMRSTransportAnalysis::GetNumberOfIterations(){
    return m_k_iteration;
}
auto satFunction = [](const TPZVec<REAL> &coord, TPZVec<STATE> &rhsVal, TPZFMatrix<STATE> &matVal) -> void
{
    REAL x = coord[0];
    REAL y = coord[1];
    REAL dt_istep = coord[2];
    
    REAL ax = 1.0;//Flux
    REAL ay = 1.0;//Flux

    REAL xShifted = x - (ax * dt_istep);
    REAL yShifted = y - (ay * dt_istep);

    rhsVal[0] = sin(M_PI * xShifted) * sin(M_PI * yShifted);
    rhsVal[0] = 1.0;
    
    if(rhsVal[0]<0.0){
        rhsVal[0]=0.0;
    }

    matVal.Resize(3, 1);
    matVal(0, 0) = 0.0;
    matVal(1, 0) = 0.0;
    matVal(2, 0) = 0.0;
};
void TMRSTransportAnalysis::Configure(int n_threads, bool UsePardiso_Q){
    
    if (UsePardiso_Q) {
        TPZSpStructMatrix<> matrix(Mesh());
        matrix.SetNumThreads(n_threads);
        SetStructuralMatrix(matrix);
//        TPZStepSolver<STATE> step;
//        step.SetDirect(ELU);
//        SetSolver(step);
        
//        TPZSpStructMatrixEigen matrix(Mesh());
//        matrix.SetNumThreads(n_threads);
//        SetStructuralMatrix(matrix);
//        TPZStepSolver<STATE> step;
//        step.SetDirect(ELU);
//        SetSolver(step);
        
    }else{
        TPZSkylineNSymStructMatrix matrix(Mesh());
        matrix.SetNumThreads(n_threads);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELU);
        SetSolver(step);
        SetStructuralMatrix(matrix);
    }
}

void TMRSTransportAnalysis::Assemble(){
    if (m_parallel_execution_Q) {
        fTransportSpMatrix->Assemble();
//        Assemble_parallel();
    }else{
        Assemble_serial();
    }
}

void TMRSTransportAnalysis::Assemble_serial(){
//    TPZAnalysis::Assemble();
    int ncells = fAlgebraicTransport.fCellsData.fVolume.size();
    if(!this->fSolver){
        DebugStop();
    }
    TPZMatrix<STATE> *mat = 0;
    TPZMatrixSolver<STATE> *matsol = dynamic_cast<TPZMatrixSolver<STATE> *>(fSolver);
    if(matsol && !matsol->Matrix())
    {
        TPZBaseMatrix *othermat = fStructMatrix->Create();
        mat = dynamic_cast<TPZMatrix<STATE> *>(othermat);
        matsol->SetMatrix(mat);
    }
    else if(matsol)
    {
        mat = matsol->Matrix().operator->();
    }
    else
    {
        DebugStop();
    }
//    mat->Redim(ncells, ncells);
//    Rhs().Redim(ncells,1);
    mat->Zero();
    Rhs().Zero();
    TPZFMatrix<STATE> &rhs = Rhs();
    //Volumetric Elements
    for (int ivol = 0; ivol<ncells; ivol++) {
        int eqindex = fAlgebraicTransport.fCellsData.fEqNumber[ivol];
        TPZFMatrix<double> elmat, ef;
        elmat.Resize(1, 1);
        ef.Resize(1,1);
        fAlgebraicTransport.Contribute(ivol, elmat, ef);
        mat->AddSub(eqindex, eqindex, elmat);
        rhs.AddSub(eqindex, 0,ef);
    }
   
    //Interface Elements
    int interID = m_sim_data->mTGeometry.mInterface_material_id;
    int ninterfaces = fAlgebraicTransport.fInterfaceData[interID].fFluxSign.size();
    if (ninterfaces<1) {
        DebugStop();
    }
    for (int interf = 0; interf<ninterfaces; interf++) {
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport.fInterfaceData[interID].fLeftRightVolIndex[interf];
        int left = lrindex.first;
        int right = lrindex.second;
        int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
        int righteq = fAlgebraicTransport.fCellsData.fEqNumber[right];
        TPZVec<int64_t> destinationindex(2);
        destinationindex[0]=lefteq;
        destinationindex[1]=righteq;
        TPZFMatrix<double> elmat, ef;
        elmat.Resize(2, 2);
        ef.Resize(2, 1);
        fAlgebraicTransport.ContributeInterface(interf,elmat, ef);
        mat->AddKel(elmat, destinationindex);
        rhs.AddFel(ef, destinationindex);
    }
    
    int inlet_mat_id = -2;
    //INLET
    ninterfaces = fAlgebraicTransport.fInterfaceData[inlet_mat_id].fFluxSign.size();
    for (int interf = 0; interf<ninterfaces; interf++) {
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport.fInterfaceData[inlet_mat_id].fLeftRightVolIndex[interf];
        int left = lrindex.first;
        int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
//        int right = lrindex.second;
        TPZVec<int64_t> destinationindex(1);
        destinationindex[0]=lefteq;
        TPZFMatrix<double> elmat, ef;
        ef.Resize(1, 1);
        fAlgebraicTransport.ContributeBCInletInterface(interf,ef);
        rhs.AddFel(ef, destinationindex);
    }
    
    int outlet_mat_id = -4;
    //outlet
    ninterfaces = fAlgebraicTransport.fInterfaceData[outlet_mat_id].fFluxSign.size();
    for (int interf = 0; interf<ninterfaces; interf++) {
        std::pair<int64_t, int64_t> lrindex= fAlgebraicTransport.fInterfaceData[outlet_mat_id].fLeftRightVolIndex[interf];
        int left = lrindex.first;
        int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
        TPZVec<int64_t> destinationindex(1);
        destinationindex[0]=lefteq;
        TPZFMatrix<double> elmat, ef;
        elmat.Resize(1, 1);
        ef.Resize(1, 1);
        fAlgebraicTransport.ContributeBCOutletInterface(interf,elmat,ef);
        mat->AddKel(elmat, destinationindex);
        rhs.AddFel(ef, destinationindex);
    }
}

void TMRSTransportAnalysis::AssembleResidual(){
    if (m_parallel_execution_Q) {
        AssembleResidual_Eigen();
    }else{
        AssembleResidual_serial();
    }
}


void TMRSTransportAnalysis::RunTimeStep(){
    
    TPZMultiphysicsCompMesh * cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
    TPZCompMesh * cmesh2 = dynamic_cast<TPZCompMesh *>(Mesh());
    if (!cmesh &&!cmesh2) {
        DebugStop();
    }
    
    int n = m_sim_data->mTNumerics.m_max_iter_transport;
    bool stop_criterion_Q = false;
    bool stop_criterion_corr_Q = false;
    REAL res_norm = 1.0;
    REAL corr_norm = 1.0;
    REAL res_tol = m_sim_data->mTNumerics.m_res_tol_transport;
    REAL corr_tol = m_sim_data->mTNumerics.m_corr_tol_transport;
//    TPZFMatrix<STATE> dx(Solution()),x(Solution());
    
    TPZFMatrix<STATE> dx(Solution()),x(Solution()),x2(Solution()), x3(Solution()) ;
    TPZFMatrix<STATE> correction(Solution());
    correction.Zero();
    fAlgebraicTransport.fCellsData.UpdateSaturations(x);
    std::cout<<"************"<<std::endl;
        
    if(m_sim_data->mTNumerics.m_TransportSolMethod==0){
        ComputeInitialGuess(x);
            
    }
    if(m_sim_data->mTNumerics.m_TransportSolMethod==1){
           
        RGK1(x, x2);
           // x2.Print(std::cout);
        fAlgebraicTransport.fCellsData.UpdateSaturations(x2);
    }
    if(m_sim_data->mTNumerics.m_TransportSolMethod==2){
        RGK1(x, x2);
        RGK1(x2, x3);
        auto sol = 0.5*(x2+x3);
        fAlgebraicTransport.fCellsData.UpdateSaturations(sol);
    }
    if(m_sim_data->mTNumerics.m_TransportSolMethod==3){
        UPDATERGK1(x,x2,cmesh2);
            //x2.Print(std::cout);
        fAlgebraicTransport.fCellsData.UpdateSaturations(x2);
    }
        
    
    
//    ComputeInitialGuess(x); // from the linear problem (tangent and residue)
//    bool QN_converge_Q = QuasiNewtonSteps(x,20); // assuming linear operator (tangent)
//    if(QN_converge_Q){
//        return;
//    }

    //Linear problem Benchmark
   
    if(Norm(Rhs()) < res_tol){
        std::cout << "Transport operator: Converged - (InitialGuess)" << std::endl;
        std::cout << "Number of iterations = " << 1 << std::endl;
        std::cout << "residue norm = " << Norm(Rhs()) << std::endl;
        return;
    }
    for(m_k_iteration = 1; m_k_iteration <= n; m_k_iteration++){
       
        NewtonIteration();
        dx = Solution();
//        std::cout<<"Sol Correct: "<<std::endl;
        x += dx;
        

        LoadSolution(x);
//        cmesh->LoadSolutionFromMultiPhysics();
//        PostProcessTimeStep();
        fAlgebraicTransport.fCellsData.UpdateSaturations(x);
//        fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(m_sim_data->mTPetroPhysics.mKrModel);
        NewtonIteration();
        dx = Solution();
        //       std::cout<<"Sol Correct: "<<std::endl;
       //        dx.Print(std::cout);
        x += dx;
        
        LoadSolution(x);
//        cmesh->LoadSolutionFromMultiPhysics();
//        PostProcessTimeStep();
        fAlgebraicTransport.fCellsData.UpdateSaturations(x);
        
        
        fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(m_sim_data->mTPetroPhysics.mKrModel);

//        fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(m_sim_data->mTNumerics.m_ISLinearKrModelQ);//No compila con esta linea

        
    
        AssembleResidual();
        corr_norm = Norm(dx);
        res_norm = Norm(Rhs());
        
#ifdef PZDEBUG
        {
            if(std::isnan(corr_norm) || std::isnan(res_norm))
            {
                DebugStop();
            }
        }
#endif
        std::cout << "res_norm " << res_norm << " corr_norm " << corr_norm << std::endl;
        stop_criterion_Q = (res_norm < res_tol);
        stop_criterion_corr_Q = (corr_norm < corr_tol);
        if (stop_criterion_Q || stop_criterion_corr_Q) {
            std::cout << "Transport operator: Converged" << std::endl;
            std::cout << "Number of iterations = " << m_k_iteration << std::endl;
            std::cout << "residue norm = " << res_norm << std::endl;
            break;
        }

    }
    if (!stop_criterion_Q && !stop_criterion_corr_Q) DebugStop(); //failed to converge
}


void TMRSTransportAnalysis::RGK1(TPZFMatrix<STATE> &sn, TPZFMatrix<STATE> &snp1){

    int nels = fAlgebraicTransport.fCellsData.fSaturation.size();
    
    std::cout<<"*** sn  ***"<<std::endl;
    sn.Print(std::cout);
    //x(0,0)=1.0;
    int interfaceId =100;
    int bcIletInterface= 2;
    int bcOutletInterface = 3;
    
    int ninternalInterfaces = fAlgebraicTransport.fInterfaceData[interfaceId].fIntegralFlux.size();
    int nBCInletInterfaces = fAlgebraicTransport.fInterfaceData[bcIletInterface].fIntegralFlux.size();
    int nBCOutletInterfaces = fAlgebraicTransport.fInterfaceData[bcOutletInterface].fIntegralFlux.size();
    //Contribute para las interfaces internas.
    
    
    for (int ivol =0; ivol <nels; ivol++) {
        double Vol = fAlgebraicTransport.fCellsData.fVolume[ivol];
        int volindex = fAlgebraicTransport.fCellsData.fEqNumber[ivol];
        double satlast =fAlgebraicTransport.fCellsData.fSaturationLastState[ivol];
        snp1(volindex) = satlast;
    }
    std::cout<<"*** snp1  ***"<<std::endl;
    snp1.Print(std::cout);
    
    //test gabriel 13/06
    sn=snp1;
    //test gabriel 13/06

    
    
    TPZVec<double> grad(nels);
    
    //Contribute para las interfaces INLET.
    for(int iinter=0; iinter<nBCInletInterfaces; iinter++){
        
        //SnP1 = Sn - L(sn)
        //L(sn) = DT*Sn*flux*A/vOL;
        
        
        auto LRCelIndex = fAlgebraicTransport.fInterfaceData[bcIletInterface].fLeftRightVolIndex[iinter];
        int LCelIndexFake = LRCelIndex.first;
        int LCelIndex = fAlgebraicTransport.fCellsData.fEqNumber[LCelIndexFake];
        
        //
        std::vector<REAL> centerCord = fAlgebraicTransport.fCellsData.fCenterCoordinate[LCelIndex];
        int geoIndex = fAlgebraicTransport.fCellsData.fGeoIndex[LCelIndex];;
        TPZGeoEl *gel2D = fGeoMesh->Element(geoIndex);
        
        
//        TEST SIDES
        int nsides = gel2D->NSides();
        int nnodes = gel2D->NNodes();
        for(int iside=nnodes; iside<nsides-1;iside++){
            TPZGeoElSide gelside(gel2D, iside);
            TPZVec<REAL> Xcenter(3,0.0);
            gelside.CenterX(Xcenter);
            REAL xcoord = Xcenter[0];
            REAL ycoord = Xcenter[1];
            REAL zcoord = Xcenter[2];
        }
//        TEST SIDES
        
        TPZVec<REAL> masscent(2,0.0);
        TPZVec<REAL> result(3,0.0);
        gel2D->CenterPoint(gel2D->NSides()-1, masscent);
        gel2D->X(masscent, result);
        
        if(result[0] != centerCord[0]){
            DebugStop();
        }
        auto allNeighs = FindAllElNeighs(gel2D);
        TPZGeoEl *trueBoundaryEl = nullptr;
        for (int ineig=0; ineig<allNeighs.size(); ineig++) {
            TPZGeoElSide gelNeig = allNeighs[ineig];
            int matid = gelNeig.Element()->MaterialId();
            if (matid == bcIletInterface) {
                trueBoundaryEl =gelNeig.Element();
                break;
            }
        }
        if (!trueBoundaryEl) {
            DebugStop();
        }
        
        TPZVec<REAL> masscent2(trueBoundaryEl->Dimension(),0.0);
        TPZVec<REAL> result2(3,0.0);
        trueBoundaryEl->CenterPoint(trueBoundaryEl->NSides()-1, masscent2);
        trueBoundaryEl->X(masscent2, result2);
        result2[2] = m_sim_data->mTNumerics.m_istep;
        //
       
        double fdt = m_sim_data->mTNumerics.m_dt;
        
        //double SatLnm1 = fAlgebraicTransport.fboundaryCMatVal[bcIletInterface].second;
        double SatLnm1 = 1.0 ;//InletValue
        
        //
       // const TPZVec<REAL> coord(centerCord[0], centerCord[1],centerCord[2]);
        TPZVec<STATE> center(3);
        center[0] =result2[0];
        center[1] =result2[1];
        center[2] =result2[2];
        TPZVec<STATE> rhsVal(2,0.0);
        TPZFMatrix<STATE> matVal;
         satFunction(center, rhsVal, matVal);
        REAL SatInlet = matVal(0,0); // UPDATE VALUE
        
        //
        
        double VolL = fAlgebraicTransport.fCellsData.fVolume[LCelIndex];
        double fluxInt = fAlgebraicTransport.fInterfaceData[bcIletInterface].fIntegralFlux[iinter];
       
        double SatNewL =  -1.0*SatLnm1*(1/VolL) *fluxInt* fdt;
        snp1(LCelIndex) +=SatNewL;
        
        auto LRgelIndex = fAlgebraicTransport.fInterfaceData[bcIletInterface].fLeftRightGelIndex[iinter];
        int LeftGeo = LRgelIndex.first;
        
        if (LeftGeo == 202 || LeftGeo==203) {
            int test = 1;
        }
        
       
       
    }
    
// Contribute para las interfaces INTERNAS
    
    for(int iinter=0; iinter<ninternalInterfaces; iinter++){
        
        //SnP1 = Sn - L(sn)
        //L(sn) = DT*Sn*flux*A/vOL;
        
        auto LRCelIndex = fAlgebraicTransport.fInterfaceData[interfaceId].fLeftRightVolIndex[iinter];
        int LCelIndexFake = LRCelIndex.first;
        int LCelIndex = fAlgebraicTransport.fCellsData.fEqNumber[LCelIndexFake];
        
//------------- Test ----------------
        auto geoIdInterface = fAlgebraicTransport.fInterfaceData[interfaceId].fcelindex[iinter];
        int geoIdCellLeft = fAlgebraicTransport.fCellsData.fGeoIndex[LCelIndex];
        TPZGeoEl *gel = fCompMesh->Reference()->Element(geoIdCellLeft);
        int nsides = gel->NSides();
        int firstside = gel->FirstSide(gel->Dimension()-1);
        int LeftSide = 0;
        
        for (int iside = firstside; iside<nsides-1; iside++) {
            TPZGeoElSide gelside(gel,iside);
            TPZStack<TPZGeoElSide> allneigh;
            gelside.AllNeighbours(allneigh);
            int nneighs = allneigh.size();
            for (int ineigh=0; ineigh<nneighs; ineigh++) {
                TPZGeoElSide elSideNeigh = allneigh[ineigh];
                auto elside = elSideNeigh.Element();
                int elSideMaterial = elside->MaterialId();
                int neighIndex = elSideNeigh.Id();
                if (neighIndex==geoIdInterface) {
                    LeftSide = iside;
                    break;
                }
            }
        }
        
//------------- End ----------------

    
        
        int RCelIndexFake = LRCelIndex.second;
        int RCelIndex = fAlgebraicTransport.fCellsData.fEqNumber[RCelIndexFake];
        
        double fdt = m_sim_data->mTNumerics.m_dt;
        auto LRgelIndex = fAlgebraicTransport.fInterfaceData[interfaceId].fLeftRightGelIndex[iinter];
        int LeftGeo = LRgelIndex.first;
        int RightGeo = LRgelIndex.second;
        
//------------- Test ----------------
        int geoIdCellRight = fAlgebraicTransport.fCellsData.fGeoIndex[RCelIndex];
        TPZGeoEl *gelRigh = fCompMesh->Reference()->Element(geoIdCellRight);
        int RightSide = 0;
        
        for (int iside= firstside; iside<nsides-1; iside++) {
            TPZGeoElSide gelside (gelRigh,iside);
            TPZStack<TPZGeoElSide> allneigh;
            gelside.AllNeighbours(allneigh);
            for (int ineigh = 0; ineigh<allneigh.size(); ineigh++) {
                TPZGeoElSide elSideNeigh = allneigh[ineigh];
                int neighIndex = elSideNeigh.Id();
                if (neighIndex==geoIdInterface) {
                    RightSide = iside;
                    break;
                }
            }
        }
//------------- End ----------------

        
        
        auto leftValue = fAlgebraicTransport.fCellsData.fSaturationEdges[LCelIndex][LeftSide-firstside];
        auto rightValue = fAlgebraicTransport.fCellsData.fSaturationEdges[RCelIndex][LeftSide-firstside];
        
        
    
        
        double SatLnm1 = sn(LCelIndex);
        double SatRnm1 = sn(RCelIndex);
        
        
        
        double VolL = fAlgebraicTransport.fCellsData.fVolume[LCelIndex];
        double VolR = fAlgebraicTransport.fCellsData.fVolume[RCelIndex];
        double fluxInt = fAlgebraicTransport.fInterfaceData[interfaceId].fIntegralFlux[iinter];
        double beta = 0;
        if (fluxInt>0) {
            beta =1;
        }
        double SatNewL =  -1.0*(beta*SatLnm1*(1/VolL) + (1-beta)*SatRnm1*(1/VolR))*fluxInt* fdt;
        double SatNewR =  1.0*(beta*SatLnm1*(1/VolL) + (1-beta)*SatRnm1*(1/VolR))*fluxInt* fdt;
        
        
        snp1(LCelIndex) +=SatNewL;
        snp1(RCelIndex) +=SatNewR;
        if(LCelIndex==RCelIndex){
            int stop=1;
        }
        if(SatNewL==SatNewR && (SatNewR!=0)){
            int stop=1;
        }
        if (LeftGeo == 202 || LeftGeo==203) {
            int stop=1;
        }
        if (RightGeo == 202 || RightGeo==203) {
            int stop=1;
        }
    }
    
   
    
    //Contribute para las interfaces OUTLET.
    for(int iinter=0; iinter<nBCOutletInterfaces; iinter++){
        
        //SnP1 = Sn - L(sn)
        //L(sn) = DT*Sn*flux*A/vOL;
       
        
        auto LRCelIndex = fAlgebraicTransport.fInterfaceData[bcOutletInterface].fLeftRightVolIndex[iinter];
        int LCelIndexFake = LRCelIndex.first;
        int LCelIndex = fAlgebraicTransport.fCellsData.fEqNumber[LCelIndexFake];
       // int RCelIndex = LRCelIndex.second;
        double fdt = m_sim_data->mTNumerics.m_dt;
        
        double SatLnm1 = sn(LCelIndex);
        //double SatLnm1 = 1.0 ;//InletVal
        
        
        double VolL = fAlgebraicTransport.fCellsData.fVolume[LCelIndex];
        double fluxInt = fAlgebraicTransport.fInterfaceData[bcOutletInterface].fIntegralFlux[iinter];
       
        double SatNewL =  -1.0*SatLnm1*(1/VolL) *fluxInt* fdt;
        snp1(LCelIndex) +=SatNewL;
       
    }
    std::cout<<"*** snp1 after  ***"<<std::endl;
    snp1.Print(std::cout);
    int stop = 1;
    
}



void TMRSTransportAnalysis::ComputeInitialGuess(TPZFMatrix<STATE> &x){
    
    TPZMultiphysicsCompMesh * cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
    TPZCompMesh * cmesh2 = dynamic_cast<TPZCompMesh *>(Mesh());
    if (!cmesh &&!cmesh2) {
        DebugStop();
    }
    
    fAlgebraicTransport.fCellsData.UpdateSaturations(x);
    fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(m_sim_data->mTPetroPhysics.mKrModel);
    std::cout<<"initial sats:"<<std::endl;
    int dim = fAlgebraicTransport.fCellsData.fSaturation.size();
    for (int i=0; i<dim; i++) {
        std::cout<<fAlgebraicTransport.fCellsData.fSaturation[i]<<std::endl;
    }
    
    LoadSolution(x);
    if(cmesh){
        cmesh->LoadSolutionFromMultiPhysics();
    }
    else{
        cmesh2->LoadSolution(x);
    }
 
    
    NewtonIteration();
    
    auto dx = Solution();
    std::cout<<"Sol Correct: "<<std::endl;
    int nrows = dx.Rows();
    
//    for(int i =0; i<nrows; i++){
//        std::cout<<dx. <<std::endl;
//
//    }
    
//    std::cout<<"SOLUTION: "<<std::endl;
//    std::cout<<x<<std::endl;
    x += Solution();

    LoadSolution(x);
    if(cmesh){
        cmesh->LoadSolutionFromMultiPhysics();
    }
    else{
        cmesh2->LoadSolution(x);
    }
    fAlgebraicTransport.fCellsData.UpdateSaturations(x);
    
    std::cout<<"new sats:"<<std::endl;
    for (int i=0; i<dim; i++) {
        std::cout<<fAlgebraicTransport.fCellsData.fSaturation[i]<<std::endl;
    }
    
    fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(m_sim_data->mTPetroPhysics.mKrModel);
    AssembleResidual();
    REAL res_norm = Norm(Rhs());
    std::cout << "Initial guess residue norm : " <<  res_norm << std::endl;
    
}

bool TMRSTransportAnalysis::QuasiNewtonSteps(TPZFMatrix<STATE> &x, int n){
    TPZMultiphysicsCompMesh * cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
    if (!cmesh) {
        DebugStop();
    }
    REAL res_tol = m_sim_data->mTNumerics.m_res_tol_transport;
    std::cout << "Quasi-Newton process : " <<  std::endl;
    for(m_k_iteration = 1; m_k_iteration <= n; m_k_iteration++){
        
        NewtonIteration();
        
        x += Solution();
        
#ifdef PZDEBUG
        {
            REAL norm = Norm(x);
            if(std::isnan(norm))
            {
                DebugStop();
            }
        }
#endif
        LoadSolution(x);
        cmesh->LoadSolutionFromMultiPhysics();
        
        
        if (m_sim_data->mTPetroPhysics.mKrModel == 0) {
            fAlgebraicTransport.fCellsData.UpdateSaturations(x);
            fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(0);
        }else{
            fAlgebraicTransport.fCellsData.UpdateSaturations(x);
            fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambdaQuasiNewton();
        }
        
        AssembleResidual();
        REAL res_norm = Norm(Rhs());
        std::cout << " Residue norm : " <<  res_norm << std::endl;
        
        res_norm = Norm(Rhs());
        
#ifdef PZDEBUG
        if(std::isnan(res_norm))
        {
            DebugStop();
        }
#endif
        
        bool stop_criterion_Q = (res_norm < res_tol);
        if (stop_criterion_Q) {
            std::cout << "Transport operator: Converged" << std::endl;
            std::cout << "Quasi-Newton iterations = " << m_k_iteration << std::endl;
            std::cout << "Residue norm = " << res_norm << std::endl;
            return true;
        }
    }
    
    return false;
}

void TMRSTransportAnalysis::NewtonIteration(){
    
    if (m_parallel_execution_Q) {
        NewtonIteration_Eigen();
    }else{
        NewtonIteration_serial();
    }
}

void TMRSTransportAnalysis::NewtonIteration_serial(){

    fTransportSpMatrix->Assemble();
    fTransportSpMatrix->Solve();

    auto ds = fTransportSpMatrix->Solution();
    TPZFMatrix<STATE> &sol = Solution();

    for (int i = 0; i < ds.rows(); i++)
    {
        sol(i, 0) = ds(i, 0);
    }
#ifdef PZDEBUG
    STATE norm = Norm(Solution());
    if (std::isnan(norm)) DebugStop();
#endif
}

void TMRSTransportAnalysis::AnalyzePattern(){
//    Assemble_mass_parallel();
//    Assemble_parallel();
//    m_transmissibility += m_mass;
//    m_analysis.analyzePattern(m_transmissibility);
    fTransportSpMatrix = new TPZAnalysisAuxEigen(&fAlgebraicTransport);
    fTransportSpMatrix->SetAlgebraicTransport(&fAlgebraicTransport);
    fTransportSpMatrix->AnalyzePattern();
    // Because in some parts this objects are needed.
    Solution().Resize(fTransportSpMatrix->NRhsRows(), 1);
    Rhs().Resize(fTransportSpMatrix->NRhsRows(), 1);
    

}

void TMRSTransportAnalysis::NewtonIteration_Eigen(){
    
//    Assemble_parallel();
    fTransportSpMatrix->Assemble();
    
#ifdef PZDEBUG
//    {
//        auto norm = fTransportSpMatrix->RhsNorm();
//        if(std::isnan(norm))
//        {
//            DebugStop();
//        }
//    }
#endif
    std::cout<<"\n ---------------------- Solve Transport----------------------"<<std::endl;
    fTransportSpMatrix->Solve();
    
    Eigen::Matrix<REAL, Eigen::Dynamic, 1> ds = fTransportSpMatrix->Solution();
    assert(Solution().Rows() == ds.rows());
    TPZFMatrix<STATE> &sol = Solution();
#ifdef USING_TBB
    tbb::parallel_for(size_t(0), size_t(ds.rows()), size_t(1),
                      [&ds, &sol] (size_t & i){
        sol(i,0) = ds(i,0);
        }
    );
#else
    for (int i = 0; i < ds.rows(); i++) {
        sol(i,0) = ds(i,0);
    }
#endif
#ifdef PZDEBUG
    STATE norm = Norm(Solution());
    if(std::isnan(norm)) DebugStop();
#endif
}

void TMRSTransportAnalysis::AssembleResidual_serial()
{

    if (0)
    {
        int ncells = fAlgebraicTransport.fCellsData.fVolume.size();
        if (!this->fSolver)
        {
            DebugStop();
        }
        Rhs().Resize(ncells, 1);
        Rhs().Zero();
        TPZFMatrix<STATE> &rhs = Rhs();

        // Volumetric Elements
        for (int ivol = 0; ivol < ncells; ivol++)
        {
            int eqindex = fAlgebraicTransport.fCellsData.fEqNumber[ivol];
            TPZFMatrix<double> ef;
            ef.Resize(1, 1);
            fAlgebraicTransport.ContributeResidual(ivol, ef);
            rhs.AddSub(eqindex, 0, ef);
        }

        // Interface Elements
        int interID = m_sim_data->mTGeometry.mInterface_material_id;
        int ninterfaces = fAlgebraicTransport.fInterfaceData[interID].fFluxSign.size();
        if (ninterfaces < 1)
        {
            DebugStop();
        }
        for (int interf = 0; interf < ninterfaces; interf++)
        {
            std::pair<int64_t, int64_t> lrindex = fAlgebraicTransport.fInterfaceData[interID].fLeftRightVolIndex[interf];
            int left = lrindex.first;
            int right = lrindex.second;
            int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
            int righteq = fAlgebraicTransport.fCellsData.fEqNumber[right];
            TPZVec<int64_t> destinationindex(2);
            destinationindex[0] = lefteq;
            destinationindex[1] = righteq;
            TPZFMatrix<double> ef;
            ef.Resize(2, 1);
            fAlgebraicTransport.ContributeInterfaceResidual(interf, ef);
            rhs.AddFel(ef, destinationindex);
        }

        int inlet_mat_id = -2;
        // INLET
        ninterfaces = fAlgebraicTransport.fInterfaceData[inlet_mat_id].fFluxSign.size();
        for (int interf = 0; interf < ninterfaces; interf++)
        {
            std::pair<int64_t, int64_t> lrindex = fAlgebraicTransport.fInterfaceData[inlet_mat_id].fLeftRightVolIndex[interf];
            int left = lrindex.first;
            int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
            TPZVec<int64_t> destinationindex(1);
            destinationindex[0] = lefteq;
            TPZFMatrix<double> ef;
            ef.Resize(1, 1);
            fAlgebraicTransport.ContributeBCInletInterface(interf, ef);
            rhs.AddFel(ef, destinationindex);
        }

        int outlet_mat_id = -4;
        // outlet
        ninterfaces = fAlgebraicTransport.fInterfaceData[outlet_mat_id].fFluxSign.size();
        for (int interf = 0; interf < ninterfaces; interf++)
        {
            std::pair<int64_t, int64_t> lrindex = fAlgebraicTransport.fInterfaceData[outlet_mat_id].fLeftRightVolIndex[interf];
            int left = lrindex.first;
            int lefteq = fAlgebraicTransport.fCellsData.fEqNumber[left];
            TPZVec<int64_t> destinationindex(1);
            destinationindex[0] = lefteq;
            TPZFMatrix<REAL> ef;
            ef.Resize(1, 1);
            fAlgebraicTransport.ContributeBCOutletInterfaceResidual(interf, ef);
            rhs.AddFel(ef, destinationindex);
        }
    }

    fTransportSpMatrix->AssembleResidual();
    TPZFMatrix<STATE> &rhs = fRhs;
    auto r = fTransportSpMatrix->Rhs().toDense();
    for (int i = 0; i < r.rows(); i++)
    {
        rhs(i, 0) = r(i, 0);
    }
}

void TMRSTransportAnalysis::AssembleResidual_Eigen(){
 
    fTransportSpMatrix->AssembleResidual();
    TPZFMatrix<STATE> &rhs = fRhs;
    Eigen::Matrix<REAL, Eigen::Dynamic, 1> r = fTransportSpMatrix->Rhs().toDense();
    assert(Rhs().Rows() == r.rows());
    #ifdef USING_TBB2
        tbb::parallel_for(size_t(0), size_t(r.rows()), size_t(1),
            [this,&r] (size_t & i){
             rhs(i,0) = r(i,0);
            }
        );
    #else
        for (int i = 0; i < r.rows(); i++) {
            rhs(i,0) = r(i,0);
        }
//    std::cout<<"RHS= "<<Rhs()<<std::endl;
    #endif
    
}



void TMRSTransportAnalysis::PostProcessTimeStep(){
    
    TPZStack<std::string,10> scalnames, vecnames;
    scalnames = m_sim_data->mTPostProcess.m_scalnamesTransport;
    constexpr int vtkRes{0}; //resolucao do vtk
    int dim = Mesh()->Reference()->Dimension();
    std::string file = m_sim_data->mTPostProcess.m_file_name_transport;
//    std::ofstream file2("transport"+std::to_string(fpostprocessindex)+".vtk");
    std::ofstream file2("transport"+std::to_string(fpostprocessindex)+".vtk");

    std::set<int> matidsToPost;
    matidsToPost.insert(1);
    
    if(m_sim_data->mTNumerics.m_UseGradRec){
        DefineGraphMesh(2,matidsToPost,scalnames,vecnames,file);
        PostProcess(vtkRes,2);
    }
    else{
            int nels = fGeoMesh->NElements();
            fCompMesh->LoadReferences();
            TPZVec<REAL> elData(nels,0);
            int ncells = fAlgebraicTransport.fCellsData.fVolume.size();
            for (int icell = 0; icell< ncells; icell++) {
                int indexGeo = fAlgebraicTransport.fCellsData.fGeoIndex[icell];
                REAL sat= fAlgebraicTransport.fCellsData.fSaturation[icell];
                auto center = fAlgebraicTransport.fCellsData.fCenterCoordinate[icell];
                TPZGeoEl *gel = fCompMesh->Reference()->Element(indexGeo);
                if(!gel) DebugStop();
        
        #ifdef PZDEBUG
                TPZManVector<REAL,3> xcenter(3,0.), xi(gel->Dimension(),0.);
                gel->CenterPoint(gel->NSides()-1, xi);
                gel->X(xi, xcenter);
                REAL diff = 0;
                std::vector<REAL> celcenter = fAlgebraicTransport.fCellsData.fCenterCoordinate[icell];
                for(int i=0; i<3; i++)
                {
                    diff += fabs(xcenter[i]-celcenter[i]);
                }
                if(diff > 1.e-16)
                {
                    DebugStop();
                }
        # endif
                if(fabs(sat) < 1e-20) sat = 0.;
                elData[indexGeo]=sat;
            }
        
            std::string fileAdjusted = file.substr(0,file.find(".vtk")) + std::to_string(fpostprocessindex) + ".vtk";
            std::ofstream file_ofstream(fileAdjusted);
            TPZVTKGeoMesh::PrintGMeshVTK(fCompMesh->Reference(), file_ofstream, elData);
           
    }
    fpostprocessindex++;
    return;
    std::map<int, TMRSDataTransfer::TFracProperties::FracProp>::iterator it;
    if (m_sim_data->mTGeometry.isThereFracture()) {
        for (it = m_sim_data->mTFracProperties.m_fracprops.begin(); it != m_sim_data->mTFracProperties.m_fracprops.end(); it++)
        {
            int matfracid = it->first;
            matidsToPost.insert(matfracid);
        }
        std::string file_frac("fracture_s.vtk");
        
#ifdef USENEWVTK2
        const std::string plotfile = file_frac.substr(0, file.find("."));
        for (auto nm : vecnames) {
            scalnames.Push(nm);
        }
        auto vtk = TPZVTKGenerator(fCompMesh, matidsToPost, scalnames, plotfile, vtkRes);
        vtk.SetNThreads(8);
        vtk.Do();
#else
        DefineGraphMesh(2,matidsToPost,scalnames,vecnames,file_frac);
        PostProcess(vtkRes,2);
#endif
    }

    
//    
//    std::set<int> matidsToPost;
//    std::map<int, TMRSDataTransfer::TFracProperties::FracProp>::iterator it;
//    if (m_sim_data->mTGeometry.isThereFracture()) {
//        for (it = m_sim_data->mTFracProperties.m_fracprops.begin(); it != m_sim_data->mTFracProperties.m_fracprops.end(); it++)
//        {
//            int matfracid = it->first;
//            matidsToPost.insert(matfracid);
//        }
//        std::string file_frac("fracture_s.vtk");
//        
//#ifdef USENEWVTK
//        const std::string plotfile = file_frac.substr(0, file.find("."));
//        for (auto nm : vecnames) {
//            scalnames.Push(nm);
//        }
//        auto vtk = TPZVTKGenerator(fCompMesh, matidsToPost, scalnames, plotfile, vtkRes);
//        vtk.SetNThreads(8);
//        vtk.Do();
//#else
//        DefineGraphMesh(2,matidsToPost,scalnames,vecnames,file_frac);
//        PostProcess(vtkRes,2);
//#endif
//    }
//    
//    if(m_sim_data->mTFracIntersectProperties.isThereFractureIntersection()){
//        matidsToPost.clear();
//        matidsToPost.insert(m_sim_data->mTGeometry.m_pressureMatId);
//        std::string file_frac2("fracture_s1d.vtk");
//        
//#ifdef USENEWVTK
//        const std::string plotfile = file_frac2.substr(0, file.find("."));
//        for (auto nm : vecnames) {
//            scalnames.Push(nm);
//        }
//        auto vtk = TPZVTKGenerator(fCompMesh, matidsToPost, scalnames, plotfile, vtkRes);
//        vtk.SetNThreads(8);
//        vtk.Do();
//#else
//        DefineGraphMesh(1,matidsToPost,scalnames,vecnames,file_frac2);
//        PostProcess(vtkRes,1);
//#endif
//    }
//    
//#ifdef USENEWVTK
//    const std::string plotfile = file.substr(0, file.find(".")); //sem o .vtk no final
//    for (auto nm : vecnames) {
//        scalnames.Push(nm);
//    }
//    
//    auto vtk = TPZVTKGenerator(fCompMesh, scalnames, plotfile, vtkRes, dim);
//    vtk.SetNThreads(8);
//    vtk.Do();
//#else
//    DefineGraphMesh(dim,scalnames,vecnames,file);
//    PostProcess(vtkRes,dim);
//#endif
    
//    int nels = fGeoMesh->NElements();
//    fCompMesh->LoadReferences();
//    TPZVec<REAL> elData(nels,0);
//    int ncells = fAlgebraicTransport.fCellsData.fVolume.size();
//    for (int icell = 0; icell< ncells; icell++) {
//        int indexGeo = fAlgebraicTransport.fCellsData.fGeoIndex[icell];
//        REAL sat= fAlgebraicTransport.fCellsData.fSaturation[icell];
//        auto center = fAlgebraicTransport.fCellsData.fCenterCoordinate[icell];
//        TPZGeoEl *gel = fCompMesh->Reference()->Element(indexGeo);
//        if(!gel) DebugStop();
//
//#ifdef PZDEBUG
//        TPZManVector<REAL,3> xcenter(3,0.), xi(gel->Dimension(),0.);
//        gel->CenterPoint(gel->NSides()-1, xi);
//        gel->X(xi, xcenter);
//        REAL diff = 0;
//        std::vector<REAL> celcenter = fAlgebraicTransport.fCellsData.fCenterCoordinate[icell];
//        for(int i=0; i<3; i++)
//        {
//            diff += fabs(xcenter[i]-celcenter[i]);
//        }
//        if(diff > 1.e-16)
//        {
//            DebugStop();
//        }
//# endif
//        if(fabs(sat) < 1e-20) sat = 0.;
//        elData[indexGeo]=sat;
//    }
//
//    std::string fileAdjusted = file.substr(0,file.find("vtk")) + std::to_string(fpostprocessindex) + ".vtk";
//    std::ofstream file_ofstream(fileAdjusted);
//    TPZVTKGeoMesh::PrintGMeshVTK(fCompMesh->Reference(), file_ofstream, elData);
//    fpostprocessindex++;
}

void TMRSTransportAnalysis::UpdateInitialSolutionFromCellsData(){
    TPZMultiphysicsCompMesh * cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
    TPZCompMesh * cmesh2 = dynamic_cast<TPZCompMesh *>(Mesh());
    if (!cmesh && !cmesh2) {
        DebugStop();
    }
    fAlgebraicTransport.fCellsData.UpdateSaturationsTo(Solution());
    fAlgebraicTransport.fCellsData.UpdateFractionalFlowsAndLambda(m_sim_data->mTPetroPhysics.mKrModel);
    LoadSolution();
    //REVISAR JOSE EN EL CASO DE NO MULTIFISICO
//    cmesh->LoadSolutionFromMultiPhysics();
}
//-------------------- FUNCTIONS  --------------------
void TMRSTransportAnalysis:: dxdydsFinder( int problemdim, int indexcel, TPZVec<int> GeoToCellIndex, TPZStack<TPZGeoElSide> ElNeighs,TPZFMatrix< REAL> &dxdyel, TPZFMatrix< REAL> &dsel){
    if (problemdim==1) {
        int nneighsdim = ElNeighs.size();
        
        double satEl = fAlgebraicTransport.fCellsData.fSaturation[indexcel];
        std::vector<double> CenterEl = fAlgebraicTransport.fCellsData.fCenterCoordinate[indexcel];
        int nels = fAlgebraicTransport.fCellsData.fSaturation.size();
         dxdyel.Resize(nneighsdim,1);
         dsel.Resize(nneighsdim,1);
        dsel.Zero();
        for (int ineigh=0; ineigh<nneighsdim; ineigh++) {
            auto igel = ElNeighs[ineigh].Element();
            int nsides = igel->NSides();
            int Eldim = igel->Dimension();
            TPZVec<REAL> masscent(Eldim,0);
            TPZVec<REAL> CenterNeigh(3,0);
            igel->CenterPoint(nsides-1, masscent);
            igel->X(masscent, CenterNeigh);
            dxdyel (ineigh,0) = CenterNeigh[0]-CenterEl[0];
            
            int geoindexNeigh = igel->Index();
            
            int icellNeighindex = GeoToCellIndex[geoindexNeigh];
            
            double satNeigh =fAlgebraicTransport.fCellsData.fSaturation[icellNeighindex];

            dsel(ineigh,0) = satNeigh - satEl;
            
            int ok=0;
            
        }
    }
    if (problemdim==2) {
//------------------------------------ PROBLEM 2D --------------------------------------------------
        int nneighsdim = ElNeighs.size();
        
        double satEl = fAlgebraicTransport.fCellsData.fSaturation[indexcel];
        std::vector<double> CenterEl = fAlgebraicTransport.fCellsData.fCenterCoordinate[indexcel];
        int nels = fAlgebraicTransport.fCellsData.fSaturation.size();
        dxdyel.Resize(nneighsdim,2);
        dsel.Resize(nneighsdim,1);
        dsel.Zero();
        for (int ineigh=0; ineigh<nneighsdim; ineigh++) {
            auto igel = ElNeighs[ineigh].Element();
            int nsides = igel->NSides();
            int Eldim = igel->Dimension();
            TPZVec<REAL> masscent(Eldim,0);
            TPZVec<REAL> CenterNeigh(3,0);
            igel->CenterPoint(nsides-1, masscent);
            igel->X(masscent, CenterNeigh);
            dxdyel (ineigh,0) = CenterNeigh[0]-CenterEl[0];
            dxdyel (ineigh,1) = CenterNeigh[1]-CenterEl[1];
            
            int geoindexNeigh = igel->Index();//get geoIndex
            
            int icellNeighindex = GeoToCellIndex[geoindexNeigh];
            
            double satNeigh =fAlgebraicTransport.fCellsData.fSaturation[icellNeighindex];

            dsel(ineigh,0) = satNeigh - satEl;
            
        }
    }
    
   
}

TPZStack<TPZGeoElSide > TMRSTransportAnalysis::FindElNeighs(TPZGeoEl *gel, bool sameDim ){
    int eldim = gel->Dimension();
    int nnodesel = gel->NNodes();
    int nsidesel = gel->NSides();
    int geoindexel = gel->Index();
    

    TPZStack<TPZGeoElSide > ElNeighs;
    TPZVec<int> IdNeighs;
    
    for (int iside = nnodesel; iside < nsidesel-1; iside++) {
        TPZGeoElSide gelside (gel,iside);
        //NNeighbours(dim) => returns the number of elements with the same dimension (eldim) on the side "iside"
        bool condition = sameDim;
        int neighs_same_dim = gelside.NNeighbours(eldim);
        if (neighs_same_dim == 1 && sameDim) {
            TPZStack<TPZGeoElSide> Allneigh;
            
            gelside.AllNeighbours(Allneigh);
            
            int nneighs = Allneigh.size();
            int stop = 0;
        
            for (int ineigh=0; ineigh<nneighs; ineigh++) {
                TPZGeoElSide neighdim = Allneigh[ineigh];
                TPZGeoElSide elsideNeig = Allneigh[ineigh];
                TPZGeoEl *gelNeigh = elsideNeig.Element();
                int dimNeigh = gelNeigh->Dimension();
                if(sameDim){
                    if (dimNeigh== eldim) {
                        ElNeighs.push_back(elsideNeig);
                        int NeighsIndex = gelNeigh->Index();
                        IdNeighs.push_back(NeighsIndex);
                        break;
                    }
                }
                else{
                    ElNeighs.push_back(elsideNeig);
                    int NeighsIndex = gelNeigh->Index();
                    IdNeighs.push_back(NeighsIndex);
                   
                }
               
            }
        }
    }

    return ElNeighs;
}

TPZStack<TPZGeoElSide > TMRSTransportAnalysis::FindAllElNeighs(TPZGeoEl *gel ){
    int eldim = gel->Dimension();
    int nnodesel = gel->NNodes();
    int nsidesel = gel->NSides();
    int geoindexel = gel->Index();
    

    TPZStack<TPZGeoElSide > ElNeighs;
    TPZVec<int> IdNeighs;
    
    for (int iside = nnodesel; iside < nsidesel-1; iside++) {
        TPZGeoElSide gelside (gel,iside);
        //NNeighbours(dim) => returns the number of elements with the same dimension (eldim) on the side "iside"
       
      
            TPZStack<TPZGeoElSide> Allneigh;
            
            gelside.AllNeighbours(Allneigh);
            
            int nneighs = Allneigh.size();
            int stop = 0;
        
            for (int ineigh=0; ineigh<nneighs; ineigh++) {
                TPZGeoElSide neighdim = Allneigh[ineigh];
                TPZGeoElSide elsideNeig = Allneigh[ineigh];
                TPZGeoEl *gelNeigh = elsideNeig.Element();
                int dimNeigh = gelNeigh->Dimension();
                
                    ElNeighs.push_back(elsideNeig);
                    int NeighsIndex = gelNeigh->Index();
                    IdNeighs.push_back(NeighsIndex);
                   
                
               
            }
        
    }

    return ElNeighs;
}

void TMRSTransportAnalysis::FindElNeigsUp(TPZGeoEl *gel, TPZStack<TPZGeoElSide> &ElNeighs, TPZVec<int> &NeighsId){
    //NOTE: Returns the saturation of the gel element and the saturations of its neighbors with the following order [iel,iel+1,iel-1]/[cell,cellRight,cellLeft]
    int eldim = gel->Dimension();
    int nnodesel = gel->NNodes();
    int nsidesel = gel->NSides();
    int geoindexel = gel->Index();
    
    NeighsId.push_back(geoindexel);
    
    
    for (int iside = nnodesel; iside < nsidesel-1; iside++) {
        TPZGeoElSide gelside (gel,iside);
        //NNeighbours(dim) => returns the number of elements with the same dimension (eldim) on the side "iside"
        int neighs_same_dim = gelside.NNeighbours(eldim);
        if (neighs_same_dim == 1) {
            TPZStack<TPZGeoElSide> Allneigh;
            
            gelside.AllNeighbours(Allneigh);
            
            int nneighs = Allneigh.size();
            int stop = 0;
            
//            std::cout<<"Lenght of AllNeigs= "<<nneighs<<std::endl;
            for (int ineigh=0; ineigh<nneighs; ineigh++) {
                TPZGeoElSide elsideNeig = Allneigh[ineigh];
                TPZGeoEl *gelNeigh = elsideNeig.Element();
                int dimNeigh = gelNeigh->Dimension();
                if (dimNeigh== eldim) {
                    ElNeighs.push_back(elsideNeig);
                    int NeighsIndex = gelNeigh->Index();
                    NeighsId.push_back(NeighsIndex);
                    break;
                }
            }
           
        }
        else{
            // inlet outlet cond
        }
    }
}

void TMRSTransportAnalysis::MaxMinMatrix(TPZCompMesh *transportMesh, TPZFMatrix<double> &result){
    int nels = fAlgebraicTransport.fCellsData.fSaturation.size();
    result.Resize(nels, 2);
    
    TPZGeoMesh *gmsh = transportMesh->Reference();
    int geonels = gmsh->NElements();
    TPZVec<int> GeoToIndex(geonels,-1);
    
    for (int jel = 0; jel<nels; jel++) {
        auto idGeoIndex = fAlgebraicTransport.fCellsData.fGeoIndex[jel];
       
        GeoToIndex[idGeoIndex] = jel;
    }
//----------------- ADD SATURATIONS (REFERENCES)------------------------------------------

    for (int iel=0; iel<nels; iel++) {
//        fAlgebraicTransport.fCellsData.fSaturation[0]=0.8;
//        fAlgebraicTransport.fCellsData.fSaturation[1]=0.2;
//        fAlgebraicTransport.fCellsData.fSaturation[2]=0.6;

        
//----------------- 2D SATURATIONS (REFERENCES)------------------------------------------
//        fAlgebraicTransport.fCellsData.fSaturation[0]=0.85;//el 1
//        fAlgebraicTransport.fCellsData.fSaturation[1]= 0.35;//el 3
//        fAlgebraicTransport.fCellsData.fSaturation[2]= 0.25;//el 2
//        fAlgebraicTransport.fCellsData.fSaturation[3]= 0.14;//el 4
//---------------------------------------------------------------------------------------

        
        auto solEl = fAlgebraicTransport.fCellsData.fSaturation[iel];
        auto volEl = fAlgebraicTransport.fCellsData.fVolume[iel];
        auto MatId= fAlgebraicTransport.fCellsData.fMatId[iel];
        auto idGeoIndex = fAlgebraicTransport.fCellsData.fGeoIndex[iel];
        
        TPZStack<TPZGeoElSide> ElNeighs;
        TPZVec<int> NeighsId;
        
        TPZGeoEl *gel = gmsh->Element(idGeoIndex);
        
        FindElNeigsUp(gel, ElNeighs, NeighsId);
        
////------------------------ XLS ---------------------------------------------
//        int nnodes = gel->NNodes();
//        int xlsdim = 2; // dimxls is a function of the problem, in this 1D case is 2
//        TPZVec<REAL> xls (xlsdim,0);
//        for (int i=0; i<nnodes-2; i++) {
//            auto elnode = gel->NodePtr(i);
//            REAL xcoord = elnode->Coord(0);
//            xls[i] = xcoord;
//        }
//------------------------ find neigh sol ---------------------------------------------
        int dimNeighs = NeighsId.size();
        TPZVec<double> AllSat(dimNeighs,0.0);
        for (int i=0; i<dimNeighs; i++) {
            int idGeoElNeig = NeighsId[i];
            int idCompEl = GeoToIndex[idGeoElNeig];
            auto satEl = fAlgebraicTransport.fCellsData.fSaturation[idCompEl];
            AllSat[i] = satEl;
        }
        //get max and min
        auto maxSat = AllSat[0];
        auto minSat = AllSat[0];
        for (int i=0; i<dimNeighs; i++) {
            auto val = AllSat[i];
            if (val>maxSat) {
                maxSat=val;
            }
            if (val<minSat) {
                minSat=val;
            }
        }
        result(iel,0) = maxSat;
        result(iel,1) = minSat;
    }
}
void TMRSTransportAnalysis::GradientLimiter1D(TPZCompMesh *transportMesh, int xlsdim, TPZFMatrix<double>&result){
//    transportMesh->LoadReferences();
// dimxls is a function of the problem, in this 1D case is 2
    TPZGeoMesh *gmesh = transportMesh->Reference();
    int pdim = 1;//dimension of problem
    
    int nels = fAlgebraicTransport.fCellsData.fSaturation.size();
    int ngelmesh = gmesh->NElements();
    result.Resize(nels, 3);
    
    TPZVec<int> GeoToIndex(ngelmesh, -1);
    TPZFMatrix<REAL> GeoGradB(ngelmesh, 3, -1);
    
    for (int jel = 0; jel<nels; jel++) {
        auto idGeoIndex = fAlgebraicTransport.fCellsData.fGeoIndex[jel];
        GeoToIndex[idGeoIndex] = jel;
    }
   
    for(int iel=0; iel<nels;iel++){
        
        auto solEl = fAlgebraicTransport.fCellsData.fSaturation[iel];
        auto volEl = fAlgebraicTransport.fCellsData.fVolume[iel];
        auto idGeoEl = fAlgebraicTransport.fCellsData.fGeoIndex[iel];
        std::vector<std::vector<REAL>> centerEl = fAlgebraicTransport.fCellsData.fCenterCoordinate;
        
        TPZGeoEl *gel = gmesh->Element(idGeoEl);
        int nnodes = gel->NNodes();
        int nsides = gel->NSides();
        TPZCompEl *cel = gel->Reference();
        
        TPZFMatrix<REAL> dxdyelement,test1;
        TPZFMatrix<REAL> dselement;
        auto ElNeighsTest = FindElNeighs(gel);
        TPZStack<TPZGeoElSide> ElNeighs;
        TPZVec<int> NeighsId;
        FindElNeigsUp(gel, ElNeighs, NeighsId);
        dxdydsFinder(pdim,iel, GeoToIndex, ElNeighs, dxdyelement, dselement);
    
        TPZFNMatrix<4, REAL>  AtA, B, AtAInv, reslt;
        auto dxdyelementTemp =dxdyelement;
        dxdyelement.Transpose();
        dxdyelement.Multiply(dxdyelementTemp, AtA);
        dxdyelement.Multiply(dselement, B);
        AtA.Inverse(AtAInv, ELU);
        AtAInv.Multiply(B, reslt);
        
        const std::string varname = "int";
        std::set<int> matids;
        matids.insert(1);
        TPZVec<STATE> value;
        value= cel->IntegrateSolution(2); //Integral de x en el volumen;
        double intx =value[0];
        double m =reslt(0,0);
        double b = solEl - (1/volEl)*(intx*m);//m is the gradient
        
//------------------------ GET ALL XLS ---------------------------------------------
//        int xlsdim = 2; // dimxls is a function of the problem, in this 1D case is 2
        TPZVec<REAL> xls (xlsdim,0);
        for (int i=0; i<nnodes-2; i++) {
            auto elnode = gel->NodePtr(i);
            REAL xcoord = elnode->Coord(0);
            xls[i] = xcoord;
            int stop = 1;
        }
        
        auto MaxMin = FindMaxMin(GeoToIndex, NeighsId);
        double Skmax = MaxMin[0];
        double Skmin = MaxMin[1];
//------------------------ Evaluate xls in SL ---------------------------------------------
        auto Sk = solEl;
        double yval=0.0;
        int stop = 1;
        TPZVec<double> Allyvar (xlsdim);
        for (int j=0 ; j<xlsdim ;j++) {
            double xl = xls[j];
            double Slxl = (m*xl)+b;
            if (Slxl>Sk) {
                yval = (Skmax-Sk)/(Slxl-Sk);//compute yval
            }
            if (Slxl<Sk) {
                yval = (Skmin-Sk)/(Slxl-Sk);
            }
            if (Slxl==Sk) {
                yval = 1.0;
            }
            double yvar = ((yval*yval)+2*yval)/((yval*yval)+yval+2);
            Allyvar[j] = yvar;
        }
        double AlphaK = Allyvar[0];//Find AlphaK
        double AlphaKVal;
        for (int j=0; j<xlsdim; j++) {
            AlphaKVal = Allyvar[j];
            if (AlphaKVal<AlphaK) {
                AlphaK=AlphaKVal;
            }
        }
        double bvar;//compute bvar
        bvar = solEl - (AlphaK/volEl)*(intx*m);
        result(iel,0)= AlphaK;//Store solutions
        result(iel,1)= bvar;
        result(iel,2)= m;
        
//------------------------ STORE SATURATIONS EVALUATED AT THE EDGES  ------------------------
        for (int iside= nnodes; iside<nsides-1; iside++) {
            TPZGeoElSide gelside (gel,iside);
            TPZVec<REAL> Xcenter(3);
            gelside.CenterX(Xcenter);
            REAL xCoordNeigh = Xcenter[0];
            REAL yCoordNeigh = Xcenter[1];
            REAL xCoordCell = centerEl[0][0];
            REAL yCoordCell = centerEl[0][1];
            REAL satEdge = solEl+(AlphaK*m*(xCoordNeigh - xCoordCell));
            fAlgebraicTransport.fCellsData.fSaturationEdges[iel][iside-nnodes] = satEdge;
        }

    }
 
    
}
//------------- TEST GRADIENT LIMITER 2D ---------
void TMRSTransportAnalysis::GradientLimiter2D(TPZCompMesh *transportMesh, TPZFMatrix<double>&result){
//    transportMesh->LoadReferences();
// dimxls is a function of the problem, in this 1D case is 2
    TPZGeoMesh *gmesh = transportMesh->Reference();
    int pdim = 2;//dimension of problem
    
    int nels = fAlgebraicTransport.fCellsData.fSaturation.size();
    int ngelmesh = gmesh->NElements();
    result.Resize(nels, 3);
    
    TPZVec<int> GeoToIndex(ngelmesh, -1);
    TPZFMatrix<REAL> GeoGradB(ngelmesh, 3, -1);
    
    for (int jel = 0; jel<nels; jel++) {
        auto idGeoIndex = fAlgebraicTransport.fCellsData.fGeoIndex[jel];
        GeoToIndex[idGeoIndex] = jel;
    }
 
    for(int iel=0; iel<nels;iel++){
        fAlgebraicTransport.fCellsData.fSaturation[0]=0.85;//el 1
        fAlgebraicTransport.fCellsData.fSaturation[1]= 0.35;//el 3
        fAlgebraicTransport.fCellsData.fSaturation[2]= 0.25;//el 2
        fAlgebraicTransport.fCellsData.fSaturation[3]= 0.14;//el 4
        
        auto solEl = fAlgebraicTransport.fCellsData.fSaturation[iel];
        auto volEl = fAlgebraicTransport.fCellsData.fVolume[iel];
        auto idGeoEl = fAlgebraicTransport.fCellsData.fGeoIndex[iel];
        
        TPZGeoEl *gel = gmesh->Element(idGeoEl);
        TPZCompEl *cel = gel->Reference();
        
        TPZFMatrix<REAL> dxdyelement,test1;
        TPZFMatrix<REAL> dselement;
        auto ElNeighsTest = FindElNeighs(gel);
        TPZStack<TPZGeoElSide> ElNeighs;
        TPZVec<int> NeighsId;
        FindElNeigsUp(gel, ElNeighs, NeighsId);
        dxdydsFinder(pdim,iel, GeoToIndex, ElNeighs, dxdyelement, dselement);
        
    
        TPZFNMatrix<4, REAL>  AtA, B, AtAInv, reslt;
        auto dxdyelementTemp =dxdyelement;
        dxdyelement.Transpose();
        dxdyelement.Multiply(dxdyelementTemp, AtA);
        dxdyelement.Multiply(dselement, B);
        AtA.Inverse(AtAInv, ELU);
        AtAInv.Multiply(B, reslt); //OK
        double gradx =reslt(0,0);
        double grady = reslt(1,0);
        double m =reslt(0,0);
//------------------------ GET ALL XLS ---------------------------------------------
        TPZFMatrix<REAL> AllXlsEl = GetXlsEl(gel);//OK
        int xlsdim = AllXlsEl.Rows();
        
        TPZVec<double> MaxMinSat = FindMaxMin(GeoToIndex,NeighsId);//OK
        double Skmax = MaxMinSat[0];
        double Skmin = MaxMinSat[1];
        
//------------------------ Evaluate xls in SL ---------------------------------------------
        auto Sk = solEl;
        double yval=0.0;
        
        int nsides=gel->NSides();
        TPZGeoElSide gelside (gel,nsides);
        TPZVec<REAL> Xcenter(3);
        gelside.CenterX(Xcenter);//find center of the quadrilateral element
        REAL xcenter = Xcenter[0];
        REAL ycenter = Xcenter[1];
        
//------------------------ Evaluate xls in SL ---------------------------------------------
        int stop = 1;
        TPZVec<double> Allyvar (xlsdim);
        for (int j=0 ; j<xlsdim ;j++) {
            double xlx = AllXlsEl(j,0);
            double xly = AllXlsEl(j,1);
            double xterm = gradx*(xlx-xcenter);
            double yterm = grady*(xly-ycenter);
            double Slxl = Sk+xterm+yterm;//OK
            
            if (Slxl>Sk) {
                yval = (Skmax-Sk)/(Slxl-Sk);//compute yval
            }
            if (Slxl<Sk) {
                yval = (Skmin-Sk)/(Slxl-Sk);
            }
            if (Slxl==Sk) {
                yval = 1.0;
            }
            double yvar = ((yval*yval)+2*yval)/((yval*yval)+yval+2);
            Allyvar[j] = yvar;
        }
       
        double AlphaK = Allyvar[0];//Find AlphaK
        double AlphaKVal;
        for (int j=0; j<xlsdim; j++) {
            AlphaKVal = Allyvar[j];
            if (AlphaKVal<AlphaK) {
                AlphaK=AlphaKVal;
            }
        }
        if (AlphaK>1.0) {
            std::cout<<"AlphaK greater than 1"<<std::endl;
            DebugStop();
        }
      
        result(iel,0)= AlphaK;//Store solutions
        result(iel,1)= gradx;
        result(iel,2)= grady;
    }
}
//----------------------- END  -------------------

TPZFMatrix<REAL>TMRSTransportAnalysis::GetXlsEl(TPZGeoEl *gel){
    int nnodes = gel->NNodes();
    int nsides = gel->NSides();
    TPZFMatrix<REAL> AllCentroids(2*nnodes,2,0.0);
    
    for(int iside=nnodes; iside<nsides-1; iside++){
        TPZGeoElSide gelside (gel,iside);
        TPZVec<REAL> Xcenter(3);
        gelside.CenterX(Xcenter);
        AllCentroids(iside-nnodes,0)= Xcenter[0];
        AllCentroids(iside-nnodes,1)= Xcenter[1];
    }
    for (int i=0; i<nnodes; i++) {
        auto inode = gel->NodePtr(i);
        REAL xcoord = inode->Coord(0);
        REAL ycoord = inode->Coord(1);
        AllCentroids(i+nnodes,0) = xcoord;
        AllCentroids(i+nnodes,1) = ycoord;
    }
    return AllCentroids;//Returns the centroids of the edges and nodes of the element
}
TPZVec<double>TMRSTransportAnalysis::FindMaxMin(TPZVec<int> GeoToIndex, TPZVec<int> NeighsId){
    int dimNeighs = NeighsId.size();
    TPZVec<double> MaxMinSat(2,0.0);
    TPZVec<double> AllSat(dimNeighs,0.0);
    for (int i=0; i<dimNeighs; i++) {
        int idGeoElNeig = NeighsId[i];
        int idCompEl = GeoToIndex[idGeoElNeig];//Warning
        auto satEl = fAlgebraicTransport.fCellsData.fSaturation[idCompEl];
        AllSat[i] = satEl;
    }
    //get max and min
    auto maxSat = AllSat[0];
    auto minSat = AllSat[0];
    for (int i=0; i<dimNeighs; i++) {
        auto val = AllSat[i];
        if (val>maxSat) {
            maxSat=val;
        }
        if (val<minSat) {
            minSat=val;
        }
    }
    MaxMinSat[0]= maxSat;
    MaxMinSat[1]= minSat;
    return MaxMinSat;
   
}
    
//-------------------- END  --------------------

void TMRSTransportAnalysis::GradientReconstruction1D(TPZCompMesh * transportMesh){
    int nels = fAlgebraicTransport.fCellsData.fVolume.size();
    TPZGeoMesh *gmesh = transportMesh->Reference();
    int dim = gmesh->Dimension();
    
    int matVol = 1;
    TPZMaterial * material = transportMesh->FindMaterial(matVol);
    
    TPZTracerFlow *matCool= dynamic_cast<TPZTracerFlow *>(material);
    {
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
    
    int ngelmesh = gmesh->NElements();
    TPZVec<int> GeoToIndex(ngelmesh, -1);
    TPZFMatrix<REAL> GeoGradB(ngelmesh, 4, -1);
    for (int jel = 0; jel<nels; jel++) {
        auto idGeoIndex = fAlgebraicTransport.fCellsData.fGeoIndex[jel];
        GeoToIndex[idGeoIndex] = jel;
        }
    TPZFMatrix<double> MxMn;
    MaxMinMatrix(transportMesh, MxMn);
    
    transportMesh->LoadReferences();
    TPZVec<double> AllSaturation(nels,1);
    for (int jel=0; jel<nels; jel++) {
        double satEl = fAlgebraicTransport.fCellsData.fSaturation[jel];
        AllSaturation[jel] = satEl;
    }
    
    
    int stop = 1;
    TPZFMatrix<double> MatResult;
    GetCellData(transportMesh, MatResult);
    
    
    
    for (int iel =0; iel<nels; iel++) {

        
        auto CenterEl = fAlgebraicTransport.fCellsData.fCenterCoordinate[iel];
        auto solEl = fAlgebraicTransport.fCellsData.fSaturation[iel];
        auto volEl = fAlgebraicTransport.fCellsData.fVolume[iel];
        auto MatId= fAlgebraicTransport.fCellsData.fMatId[iel];
        auto idGeoIndex = fAlgebraicTransport.fCellsData.fGeoIndex[iel];
        
        auto LRIndex=fAlgebraicTransport.fInterfaceData[4].fLeftRightGelIndex[0];
        int Lindex = LRIndex.first;
        int Rindex = LRIndex.first;
        auto flux=fAlgebraicTransport.fInterfaceData[100].fIntegralFlux[0];
        TPZGeoEl *gel = gmesh->Element(idGeoIndex);
        auto XlsEl = GetXlsEl(gel);
        
    
        TPZVec<double> AllSaturation(ngelmesh,0);
        int nvols = fAlgebraicTransport.fCellsData.fSaturationEdges.size();
        

//----------------------- GRADIENT LIMITER FUNCTION  --------------------
        int xlsdim = 2;
        
        TPZFMatrix<double> AlphaKB;
        GradientLimiter1D(transportMesh, xlsdim, AlphaKB);
        double AlphaK_GL= AlphaKB(iel,0);
        double bvar_GL= AlphaKB(iel,1);
        double m = AlphaKB(iel,2);
        int side = gel->NSides()-1;
        TPZVec<REAL> masscent(gel->Dimension(),0.0);
        gel->CenterPoint(side, masscent);
        TPZVec<REAL> result(3,0.0);
        gel->X(masscent, result);
//-------------------------------------------------------------------------
        GeoGradB(idGeoIndex,0)=AlphaK_GL*m ;
        GeoGradB(idGeoIndex,1)=bvar_GL;
        GeoGradB(idGeoIndex,2)=solEl;
        GeoGradB(idGeoIndex,3)=result[0];
     }
    matCool->SetGradB(GeoGradB);
//    GeoGradB.Print(std::cout);
    int test= 1;
    
}

void TMRSTransportAnalysis::GradientReconstruction2D(TPZCompMesh *transportMesh){
    TPZGeoMesh *gmesh = transportMesh->Reference();
    int ngeoels = gmesh->NElements();
    
    {
        std::ofstream out("cmesh.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(transportMesh,out);
        std::ofstream outG("gmsh_2D.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,outG);
    }
    
//------------------------ WARNING TEST ----------------------------
    int matVol = 1;
    TPZMaterial * material = transportMesh->FindMaterial(matVol);
    
    TPZTracerFlow *matCool= dynamic_cast<TPZTracerFlow *>(material);
    TPZFMatrix<REAL> GeoGradB(ngeoels, 5, -1);
//------------------------ END TEST --------------------------------
    
    TPZFMatrix<double> result;
    GradientLimiter2D(transportMesh,result);
    TPZFMatrix<double> M_Result;
    GetCellData(transportMesh, M_Result);
    
    int nels = fAlgebraicTransport.fCellsData.fSaturation.size();
    
    TPZVec<int> GeoToIndex(ngeoels,-1);
    for (int icel=0; icel<nels; icel++) {
        int gelId = fAlgebraicTransport.fCellsData.fGeoIndex[icel];
        GeoToIndex[gelId] = icel;
    }
    
    fAlgebraicTransport.fCellsData.fSaturation[0]=0.85;//el 1
    fAlgebraicTransport.fCellsData.fSaturation[1]= 0.35;//el 3
    fAlgebraicTransport.fCellsData.fSaturation[2]= 0.35;//el 2
    fAlgebraicTransport.fCellsData.fSaturation[3]= 0.14;//el 4
    
    for (int iel = 0; iel<nels; iel++) {
        auto volEl = fAlgebraicTransport.fCellsData.fVolume[iel];
        auto satEl = fAlgebraicTransport.fCellsData.fSaturation[iel];
        int geoId = fAlgebraicTransport.fCellsData.fGeoIndex[iel];
        std::vector<REAL> centerCell = fAlgebraicTransport.fCellsData.fCenterCoordinate[iel];
        
        TPZGeoEl *gel = gmesh->Element(geoId);
        int eldim = gel->Dimension();
        TPZCompEl *cel = gel->Reference();
        
//************************** TEST 19_05 **************************
        TPZVec<double> NeighSols;
        TPZFMatrix<REAL> NeighsCenter;
        GetNeigsSol(gel, GeoToIndex,NeighSols);//OK
        GetNeigsCentr(gel, NeighsCenter);//OK
        
        TPZVec<double>GradCell= GradientCell(centerCell, satEl, NeighsCenter, NeighSols);
        auto AlphaKCell = GradientLimiter(gel, satEl, GradCell, NeighSols);//OK
        if (iel==0) {
            int ncase = IdentifyCase(gel);
            int stop = 1;
        }
//************************** TEST 19_05 **************************
        
        int nsides=gel->NSides();
        TPZGeoElSide gelside (gel,nsides);
        TPZVec<REAL> Xcenter(3);
        gelside.CenterX(Xcenter);//find center of the quadrilateral element
        REAL xcenter = Xcenter[0];
        REAL ycenter = Xcenter[1];
        
        double AlphaK = result(iel,0);
        double gradx = result(iel,1);
        double grady =result(iel,2);
        double AlphaK_times_gradx=AlphaK*gradx;
        double AlphaK_times_grady=AlphaK*grady;
        
//-------------------------- Store solutions --------------------------------
        GeoGradB(geoId,0)=satEl;
        GeoGradB(geoId,1)=AlphaK_times_gradx;
        GeoGradB(geoId,2)=AlphaK_times_grady;
        GeoGradB(geoId,3)=xcenter;
        GeoGradB(geoId,4)=ycenter;
    }
    matCool->SetGradB(GeoGradB);
}

//----------------------- UPDATE SOLUTION  -----------------------
void TMRSTransportAnalysis::UPDATERGK1(TPZFMatrix<STATE> &sn, TPZFMatrix<STATE> &snp1, TPZCompMesh *TransportMesh){
    
    
    
    //----------------------- DATA  -----------------------
    TPZGeoMesh *gmesh = TransportMesh->Reference();
    {
    std::ofstream out("geomeshiMRS.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
    
    int geoNels = gmesh->NElements();
    int nels = fAlgebraicTransport.fCellsData.fSaturation.size();
    TPZVec<double>grad(geoNels,1000);
    TPZFMatrix<double>allSlxs(geoNels,2,1000);
    double sInlet = 1.0;
    int pdim= 1;//Problem Dimension
    int matVol = 1;
    TPZMaterial * material = TransportMesh->FindMaterial(matVol);
    TPZTracerFlow *matCool= dynamic_cast<TPZTracerFlow *>(material);
    TPZFMatrix<REAL> GeoGradB(geoNels, 4, -1);
    TPZFMatrix<double>  Result;
    double dt = fAlgebraicTransport.fdt;
    double dx = 1.0/nels;
    //----------------------- DATA  -----------------------
    
//    fAlgebraicTransport.fCellsData.fSaturation[0]=0.392558;//el 0
//    fAlgebraicTransport.fCellsData.fSaturation[1]= 0.0;//el 4
//    fAlgebraicTransport.fCellsData.fSaturation[2]= 0.007441;//el 1
//    fAlgebraicTransport.fCellsData.fSaturation[3]= 0.0;//el 3
    
    std::cout<<"**** SN  *****"<<std::endl;
//    sn.Print(std::cout);
    
    for (int ivol =0; ivol <nels; ivol++) {
        double Vol = fAlgebraicTransport.fCellsData.fVolume[ivol];
        int volindex = fAlgebraicTransport.fCellsData.fEqNumber[ivol];
        double satlast =fAlgebraicTransport.fCellsData.fSaturationLastState[ivol];
        snp1(volindex) = satlast;
    }
    
    std::cout<<"**** SN before  *****"<<std::endl;
//    snp1.Print(std::cout);
    
    fAlgebraicTransport.fCellsData.UpdateSaturations(snp1);
    
    GetCellData(TransportMesh,Result);
//    Result.Print(std::cout);
    
    for (int iel =0; iel<nels; iel++) {
        int geoindex = fAlgebraicTransport.fCellsData.fGeoIndex[iel];
        int eqNumber = fAlgebraicTransport.fCellsData.fEqNumber[iel];
        snp1(eqNumber,0) -= (dt/dx)*Result(geoindex, 3);
        
    }
    
    std::cout<<"**** SN after  *****"<<std::endl;
//    snp1.Print(std::cout);

    int stop =1;
//    fAlgebraicTransport.fCellsData.UpdateSaturations(snp1);
   // fAlgebraicTransport.fCellsData.UpdateSaturationsLastState(sn);
    
    TPZVec<int> geoToIndex(geoNels,-1);
    for (int jcell = 0; jcell<nels; jcell++) {
        auto geoId = fAlgebraicTransport.fCellsData.fGeoIndex[jcell];
        geoToIndex[geoId]=jcell;
    }
    
    //----- Reconstruction of the  Updated Solution Using Gradient Reconstruction -----
    auto GradPostProcess = m_sim_data->mTNumerics.m_UseGradRec;
    TPZFMatrix<double> Grad_AlphaK_Result;
//    if (GradPostProcess) {
//        UpdatedGradientReconstruction1D(snp1, TransportMesh, Grad_AlphaK_Result);
//    }
    
// ********************  TEST 19_05  ********************
//    if (GradPostProcess) {
//        GradientReconstruction2D(TransportMesh);
//    }
// ********************  TEST 19_05  ********************
 
    //------------------------ PostProcess ------------------------
    for (int icell = 0; icell<nels; icell++) {
        auto geoIndexCell = fAlgebraicTransport.fCellsData.fGeoIndex[icell];
        auto satCell = fAlgebraicTransport.fCellsData.fSaturation[icell];
        
        auto gel = gmesh->Element(geoIndexCell);
        TPZFMatrix<REAL> MatResult;
        GetNeigsCentr(gel, MatResult);
        
        int side = gel->NSides()-1;
        TPZVec<REAL> masscent(gel->Dimension(),0.0);
        gel->CenterPoint(side, masscent);
        TPZVec<REAL> result(3,0.0);
        gel->X(masscent, result);
        //    Result(geoIndexCell,0) = gradCell;
        //    Result(geoIndexCell,1) = AlphaKResult;
        //    Result(geoIndexCell,2) = SR;
        double grady = 0.0;
        
        double AlphaK_times_gradx = Result(geoIndexCell,1)*Result(geoIndexCell,0);
                   //AlphaK_times_gradx = 0.0;
        GeoGradB(geoIndexCell,0)=AlphaK_times_gradx;
        GeoGradB(geoIndexCell,1)=grady;
        GeoGradB(geoIndexCell,2)=satCell;
        GeoGradB(geoIndexCell,3)=result[0];

//        if (GradPostProcess) {
//            double gradCell = Grad_AlphaK_Result(geoIndexCell,0);
//            double alphaKCell = Grad_AlphaK_Result(geoIndexCell,1);
//            double AlphaK_times_gradx = gradCell*alphaKCell;
//
//            GeoGradB(geoIndexCell,0)=AlphaK_times_gradx;
//            GeoGradB(geoIndexCell,1)=grady;
//            GeoGradB(geoIndexCell,2)=satCell;
//            GeoGradB(geoIndexCell,3)=result[0];
//        }
//        else{
//            //------ without Gradient Reconstruction-----
//            double AlphaK_times_gradx = Result(geoIndexCell,1)*Result(geoIndexCell,0);
//            //AlphaK_times_gradx = 0.0;
//
//            GeoGradB(geoIndexCell,0)=AlphaK_times_gradx;
//            GeoGradB(geoIndexCell,1)=grady;
//            GeoGradB(geoIndexCell,2)=satCell;
//            GeoGradB(geoIndexCell,3)=result[0];
//        }
    }
    matCool->SetGradB(GeoGradB);
    
}
//----------------------- UPDATE SOLUTION TEST -----------------------

void TMRSTransportAnalysis::FindAlphaK(double SlxsR,double SlxsL, double Skmax, double Skmin, double cellSat, double &AlphaKResult){
//    transportMesh->LoadReferences();
// dimxls is a function of the problem, in this 1D case is 2

    int pdim = 1;//dimension of problem
    auto Sk = cellSat;
    double yval=0.0;
    int stop = 1;
    
    TPZVec<double> Allyvar(2);
    TPZVec<double>xls(2);
    xls[0]=SlxsR;
    xls[1]=SlxsL;
    
    for (int j=0 ; j<xls.size() ;j++) {
        double Slxl = xls[j];
        if (Slxl>Sk) {
            yval = (Skmax-Sk)/(Slxl-Sk);//compute yval
        }
        if (Slxl<Sk) {
            yval = (Skmin-Sk)/(Slxl-Sk);
        }
        if (Slxl==Sk) {
            yval = 1.0;
        }
        double yvar = ((yval*yval)+2*yval)/((yval*yval)+yval+2);
            Allyvar[j] = yvar;
        }
        double AlphaK = Allyvar[0];//Find AlphaK
        double AlphaKVal;
        for (int j=0; j<xls.size(); j++) {
            AlphaKVal = Allyvar[j];
            if (AlphaKVal<AlphaK) {
                AlphaK=AlphaKVal;
            }
        }
    AlphaKResult=AlphaK;
}

void TMRSTransportAnalysis::GetCellData(TPZCompMesh *TransportMesh, TPZFMatrix<double> &Result){
    //----------------------- DATA  -----------------------
    TPZGeoMesh *gmesh = TransportMesh->Reference();
    int geoNels = gmesh->NElements();
    int nels = fAlgebraicTransport.fCellsData.fSaturation.size();
    TPZVec<double>grad(geoNels,1000);
    TPZFMatrix<double>allSlxs(geoNels,2,1000);
    double sInlet = 1.0;
    int pdim= 1;//Problem Dimension
    TPZFMatrix<REAL> GeoGradB(geoNels, 4, -1);
    Result.Resize(geoNels, 4);
    TPZFMatrix<double> LeftRightCon(nels,4);
    for (int i=0; i<nels; i++) {
        LeftRightCon(i,0)=-1;
        LeftRightCon(i,1)=-1;
        LeftRightCon(i,2)=-1;
        LeftRightCon(i,3)=-1;
    }
    
    //----------------------- DATA  -----------------------
    TPZVec<int> geoToIndex(geoNels,-1);
    
    for (int jcell = 0; jcell<nels; jcell++) {
        auto geoId = fAlgebraicTransport.fCellsData.fGeoIndex[jcell];
        geoToIndex[geoId]=jcell;
    }
    
//    for (int jcell = 0; jcell<geoNels; jcell++) {
//        std::cout<<"geocell= "<<jcell<<" geoid= "<<geoToIndex[jcell]<<std::endl;
//
//    }
    
    for (int icell = 0; icell<nels; icell++) {
        
//        std::cout<<"icell=  "<<icell<<std::endl;
        
        auto geoIndexCell = fAlgebraicTransport.fCellsData.fGeoIndex[icell];
        auto satCell = fAlgebraicTransport.fCellsData.fSaturation[icell];
        
    
        TPZStack<TPZGeoElSide> cellNeighs;//Find Neighbors
        TPZFMatrix<REAL> dxdyelement;
        TPZFMatrix<REAL> dselement;
        TPZVec<int> neighsId;
        TPZGeoEl *gel = gmesh->Element(geoIndexCell);
        FindElNeigsUp(gel, cellNeighs, neighsId);
        dxdydsFinder(pdim,icell, geoToIndex, cellNeighs, dxdyelement, dselement);
        auto dx = std::abs(dxdyelement(0));
    
        int dimNeighs = neighsId.size();//Find Saturations of the Neighbors
        TPZVec<double> allSats(dimNeighs,0.0);
       
        for (int ineigh = 0; ineigh<dimNeighs; ineigh++) {
            int idGeoNeigh = neighsId[ineigh];
            if(idGeoNeigh==10){
                int ok=0;
            }
            int idCell = geoToIndex[idGeoNeigh];
            auto satNeigh = fAlgebraicTransport.fCellsData.fSaturation[idCell];
            
//            std::cout<<"idNeigh= "<<idCell<<" satNeigh= "<<satNeigh<<std::endl;
            
            allSats[ineigh] = satNeigh;//TEST THIS SOLUTION
            LeftRightCon(icell,ineigh) = idGeoNeigh;
            
        }
        
        //-------- Compute the gradient ---------------
    
        int nnodesel = gel->NNodes();
        int nsidesel = gel->NSides();
    
        TPZVec<int> allMatIdNeighs;
        int type = 0;
        for (int iside = nnodesel; iside<nsidesel-1; iside++) {
            TPZGeoElSide gelside (gel,iside);
            int gelSideDim = gelside.Dimension();
            int neighs_same_dim = gelside.NNeighbours(gelSideDim);
            if (neighs_same_dim != 0) {
                TPZStack<TPZGeoElSide> allNeigh;
                gelside.AllNeighbours(allNeigh);
    
                int nNeighs = allNeigh.size();
                double stop = 1.0;
                for (int iNeig=0; iNeig<nNeighs; iNeig++) {
                    TPZGeoElSide elSideNeigh = allNeigh[iNeig];
                    auto elNeigh = elSideNeigh.Element();
                    auto matIdNeigh = elNeigh->MaterialId();
                    allMatIdNeighs.push_back(matIdNeigh);
                }
            }
        }
        
        for (int jneigh = 0; jneigh<allMatIdNeighs.size(); jneigh++) {
            int value = allMatIdNeighs[jneigh];
            if (value==2) {
                type=value;
                break;
            }
            if (value==3) {
                type=value;
                break;
            }
        }
//        std::cout<<" cell type= "<<type<<std::endl;
        LeftRightCon(icell,3)=type;
        if (type==2) {//InletElement
            auto neighSat=allSats[1];//CAUTION WITH THE INDEX
            grad[geoIndexCell]=(neighSat-sInlet)/(2*dx);
            auto gradientval =(neighSat-sInlet)/(2*dx);
            int stop =1;
        }
        if (type==3) {//OutletElement
            auto cellSat = allSats[0];//CAUTION WITH THE INDEX
            auto neighSat = allSats[1];
            grad[geoIndexCell]=(cellSat-neighSat)/(dx);
            auto gradientval =(cellSat-neighSat)/(dx);
            int stop =1;
        }
        if (type==0) {
            auto neighSatR = allSats[1];//CAUTION WITH THE INDEX
            auto neighSatL = allSats[2];
            grad[geoIndexCell] = (neighSatR-neighSatL)/(2*dx);
            auto gradientval =(neighSatR-neighSatL)/(2*dx);
            int stop =1;
        }
//        std::cout<<"Type Element= "<<type<<" IdIndexEl= "<<geoIndexCell<<std::endl;
//        std::cout<<"GRADIENT= "<<grad[geoIndexCell]<<std::endl;
        
        //-------- End Compute the gradient ---------------
    
        TPZVec<double> MaxMinSat = FindMaxMin(geoToIndex,neighsId);//Find Max&Min
        double Skmax = MaxMinSat[0];
        double Skmin = MaxMinSat[1];
        if (type==2) {
            Skmax=sInlet;
        }
    
        auto gradCell = grad[geoIndexCell];//Compute Edges Values
        allSlxs(icell,0)=satCell+(gradCell*dx*0.5);
        allSlxs(icell,1)=satCell-(gradCell*dx*0.5);
        auto SlxsR = allSlxs(icell,0);
        auto SlxsL = allSlxs(icell,1);
    
        double AlphaKResult;
        FindAlphaK(SlxsR, SlxsL, Skmax,  Skmin, satCell, AlphaKResult);//Get AlphaK
        auto AlphaK_times_gradx= AlphaKResult*gradCell;
    
    
        auto SR = satCell + (AlphaKResult*gradCell*dx*0.5);//SR is the updated solution of the cell
        
        Result(geoIndexCell,0) = gradCell;
        Result(geoIndexCell,1) = AlphaKResult;
        Result(geoIndexCell,2) = SR;
       
    }
    
    for (int icell=0; icell<nels; icell++) {
        int myGeoIndex = LeftRightCon(icell,0);
        int RightGeoIndex = LeftRightCon(icell,1);
        int LeftGeoIndex = LeftRightCon(icell,2);
        
        
        double MySat = Result(myGeoIndex, 2);
        double SRcell = Result(myGeoIndex,2);
        double RightSat = Result(RightGeoIndex, 2);
        
        //Note: SRGlob = SR[icell]-SR[icell-1]
        
        if(LeftRightCon(icell,3)==2){//Type 2 => Inlet Element
            double SRLeftCell = sInlet;
            Result(myGeoIndex, 3) = SRcell - SRLeftCell; //MySat - 1.0;
            auto value = SRcell - SRLeftCell;
            int stop =1;
        }
        if(LeftRightCon(icell,3)==3){//Type 3 => OutLet Element
            double SRLeftCell = Result(RightGeoIndex,2);//Use RightGeoIndex (second index) because type 3 only have one neighbor
            Result(myGeoIndex, 3) = SRcell - SRLeftCell; //MySat - RightSat;
            double value = SRcell - SRLeftCell;
            int stop = 1;
        }
        if(LeftRightCon(icell,3)==0){//Type 0 =>  (!= InletElement && != OutletElemet)
            double SRLeftCell = Result(LeftGeoIndex, 2);
            Result(myGeoIndex, 3) = SRcell - SRLeftCell;
            double value = SRcell-SRLeftCell;
        }

    }
    
}

void TMRSTransportAnalysis::UpdatedGradientReconstruction1D(TPZFMatrix<double> snp1,TPZCompMesh *TransportMesh,TPZFMatrix<double> &Result){
    //----------- DATA--------
    int nCells = snp1.Rows();
    TPZVec<double> grad(nCells,1000);
    TPZGeoMesh *gmesh = TransportMesh->Reference();
    int geoNels = gmesh->NElements();
    Result.Resize(geoNels, 2);
    double sInlet = 1.0;
    int pdim = 1;
    //----------- DATA--------
    
    TPZVec<int> geoToIndex(geoNels,-1);
    
    for (int jcell = 0; jcell<nCells; jcell++) {
        auto geoId = fAlgebraicTransport.fCellsData.fGeoIndex[jcell];
        geoToIndex[geoId]=jcell;
    }
    
    for (int icell=0; icell<nCells; icell++) {
        
        int geoIndexCell = fAlgebraicTransport.fCellsData.fGeoIndex[icell];
        auto satCell = fAlgebraicTransport.fCellsData.fSaturation[icell];
        double gradCell = -100;
        TPZGeoEl *gel = gmesh->Element(geoIndexCell);
        TPZStack<TPZGeoElSide> cellNeighs;//Find Neighbors
        TPZFMatrix<REAL> dxdyelement;
        TPZFMatrix<REAL> dselement;
        TPZVec<int> neighsId;
        FindElNeigsUp(gel, cellNeighs, neighsId);
        dxdydsFinder(pdim,icell, geoToIndex, cellNeighs, dxdyelement, dselement);
        auto dx = std::abs(dxdyelement(0));
        
        int dimNeighs = neighsId.size();//Find Saturations of the Neighbors
        TPZVec<double> allSats(dimNeighs,0.0);
           
        for (int ineigh = 0; ineigh<dimNeighs; ineigh++) {
            int idGeoNeigh = neighsId[ineigh];
            if(idGeoNeigh==10){
                int ok=0;
            }
            int idCell = geoToIndex[idGeoNeigh];
            auto satNeigh = fAlgebraicTransport.fCellsData.fSaturation[idCell];
            allSats[ineigh] = satNeigh;
        }
        
        //-------- Compute the gradient ---------------
    
        int nnodesel = gel->NNodes();
        int nsidesel = gel->NSides();
    
        TPZVec<int> allMatIdNeighs;
        int type = 0;
        for (int iside = nnodesel; iside<nsidesel-1; iside++) {
            TPZGeoElSide gelside (gel,iside);
            int gelSideDim = gelside.Dimension();
            int neighs_same_dim = gelside.NNeighbours(gelSideDim);
            if (neighs_same_dim != 0) {
                TPZStack<TPZGeoElSide> allNeigh;
                gelside.AllNeighbours(allNeigh);
    
                int nNeighs = allNeigh.size();
                double stop = 1.0;
                for (int iNeig=0; iNeig<nNeighs; iNeig++) {
                    TPZGeoElSide elSideNeigh = allNeigh[iNeig];
                    auto elNeigh = elSideNeigh.Element();
                    auto matIdNeigh = elNeigh->MaterialId();
                    allMatIdNeighs.push_back(matIdNeigh);
                }
            }
        }
        
        
        for (int jneigh = 0; jneigh<allMatIdNeighs.size(); jneigh++) {
            int value = allMatIdNeighs[jneigh];
            if (value==2) {
                type=value;
                break;
            }
            if (value==3) {
                type=value;
                break;
            }
        }
        if (type==2) {//InletElement
            auto neighSat=allSats[1];//CAUTION WITH THE INDEX
            auto gradCell=(neighSat-sInlet)/(2*dx);
        }
        if (type==3) {//OutletElement
            auto cellSat = allSats[0];//CAUTION WITH THE INDEX
            auto neighSat = allSats[1];
            auto gradCell=(cellSat-neighSat)/(dx);
        }
        if (type==0) {
            auto neighSatR = allSats[1];//CAUTION WITH THE INDEX
            auto neighSatL = allSats[2];
            auto gradCell = (neighSatR-neighSatL)/(2*dx);
        }
        
        //-------- End Compute the gradient ---------------
        
        TPZVec<double> MaxMinSat = FindMaxMin(geoToIndex,neighsId);//Find Max&Min
        double Skmax = MaxMinSat[0];
        double Skmin = MaxMinSat[1];
        if (type==2) {
            Skmax=sInlet;
        }
    
        auto SlxsR = satCell+(gradCell*dx*0.5);//Compute Edges Values
        auto SlxsL = satCell-(gradCell*dx*0.5);
    
        double AlphaKResult;
        FindAlphaK(SlxsR, SlxsL, Skmax,  Skmin, satCell, AlphaKResult);//Get AlphaK
        
        Result(geoIndexCell,0) = gradCell;
        Result(geoIndexCell,1) = AlphaKResult;
        
    }
    
}

void TMRSTransportAnalysis::GetNeigsCentr(TPZGeoEl *gel , TPZFMatrix<REAL> &MatResult){
    int eldim = gel->Dimension();
    int nnodes = gel->NNodes();
    int nsides = gel->NSides();
    int geoindexel = gel->Index();
    TPZVec<int> NeighIds;
    TPZStack<TPZGeoElSide> ElNeighs;
    
    for (int iside = nnodes; iside< nsides-1 ; iside++) {
        TPZGeoElSide gelside (gel, iside);
        int neighs_same_dim = gelside.NNeighbours(eldim);
        if (neighs_same_dim==1) {
            TPZStack<TPZGeoElSide> allneigs;
            gelside.AllNeighbours(allneigs);
            int nneighs = allneigs.size();
            
            for (int ineig =0; ineig<nneighs; ineig++) {
                TPZGeoElSide elsideNeigh = allneigs[ineig];
                TPZGeoEl *gelNeigh = elsideNeigh.Element();
                int neighDim = gelNeigh->Dimension();
                if (neighDim == eldim) {
                    ElNeighs.push_back(elsideNeigh);
                    int neighIndex = gelNeigh->Index();
                    NeighIds.push_back(neighIndex);
                }
            }
        }
    }
    int nneighs = NeighIds.size();
    MatResult.Resize(nneighs, 2);//Resize Matrix
    for (int ineigh = 0; ineigh<nneighs; ineigh++) {
        auto igel = ElNeighs[ineigh].Element();
        int nsides = igel->NSides();
        int Neighdim = igel->Dimension();
        TPZVec<REAL> masscent(Neighdim,0);
        TPZVec<REAL> centerNeigh (3,0);
        igel->CenterPoint(nsides-1, masscent);
        igel->X(masscent, centerNeigh);
        MatResult(ineigh,0) = centerNeigh[0];
        MatResult(ineigh,1) = centerNeigh[1];
    }
}

void TMRSTransportAnalysis::GetNeigsSol(TPZGeoEl *gel, TPZVec<int> GeoToCellIndex,TPZVec<double> &Result){
    int eldim = gel->Dimension();
    int nnodes = gel->NNodes();
    int nsides = gel->NSides();
    int geoindexel = gel->Index();
    TPZVec<int> NeighIds;
    TPZStack<TPZGeoElSide> ElNeighs;
    
    for (int iside = nnodes; iside< nsides-1 ; iside++) {
        TPZGeoElSide gelside (gel, iside);
        int neighs_same_dim = gelside.NNeighbours(eldim);
        if (neighs_same_dim==1) {
            TPZStack<TPZGeoElSide> allneigs;
            gelside.AllNeighbours(allneigs);
            int nneighs = allneigs.size();
            
            for (int ineig =0; ineig<nneighs; ineig++) {
                TPZGeoElSide elsideNeigh = allneigs[ineig];
                TPZGeoEl *gelNeigh = elsideNeigh.Element();
                int neighDim = gelNeigh->Dimension();
                if (neighDim == eldim) {
                    ElNeighs.push_back(elsideNeigh);
                    int neighIndex = gelNeigh->Index();
                    NeighIds.push_back(neighIndex);
                    break;
                }
            }
        }
    }
    int nneighs = NeighIds.size();
    Result.resize(nneighs);
    for (int ineigh = 0; ineigh<nneighs; ineigh++) {
        auto igel = ElNeighs[ineigh].Element();
        int geoIndexNeigh = igel->Index();
        int cellNeighIndex = GeoToCellIndex[geoIndexNeigh];
        double satNeigh = fAlgebraicTransport.fCellsData.fSaturation[cellNeighIndex];
        Result[ineigh] = satNeigh;
    }
}

TPZVec<double> TMRSTransportAnalysis::GradientCell(std::vector<REAL> CenterCell, double CellSol, TPZFMatrix<REAL> CenterNeighs, TPZVec<double> NeighsSols){
    
    int nneighs = NeighsSols.size();
    int dimcmesh;
    TPZFMatrix<REAL> dxdy(nneighs,2);
    TPZFMatrix<REAL> ds(nneighs);
    for (int ineigh=0; ineigh< nneighs; ineigh++) {
        auto xcoordNeigh = CenterNeighs(ineigh,0);
        auto ycoordNeigh = CenterNeighs(ineigh,1);
        dxdy(ineigh,0) = xcoordNeigh - CenterCell[0];
        dxdy(ineigh,1) = ycoordNeigh - CenterCell[1];
        
        ds(ineigh) = NeighsSols[ineigh] - CellSol;
    }
    TPZFNMatrix<4, REAL>  AtA, B, AtAInv, GradSol;
    auto dxdyTemp =dxdy;
    dxdy.Transpose();
    dxdy.Multiply(dxdyTemp, AtA);
    dxdy.Multiply(ds, B);
    AtA.Inverse(AtAInv, ELU);
    AtAInv.Multiply(B, GradSol);
    GradSol.Print(std::cout);
    double gradx = GradSol(0,0);
    double grady = GradSol(1,0);
    TPZVec<double>GradientSol(2);
    GradientSol[0]= gradx;
    GradientSol[1]= grady;
    return GradientSol;
}

double TMRSTransportAnalysis::GradientLimiter(TPZGeoEl *gel, double cellSol, TPZVec<double> gradCell, TPZVec<double> neighsSols){
    
    int nnodes = gel->NNodes();
    int nsides = gel->NSides();
    double tolerance = 0.0000001;
//    std::cout<<"tolerance= "<<tolerance<<std::endl;
    
    TPZFMatrix<REAL> allElNodes(nnodes,2);//consider mesh 2D
    TPZVec<double> allSlxl(nnodes);
    
    TPZGeoElSide gelside (gel,nsides);//find the centroid of the element
    TPZVec<REAL> Xcenter(3);
    gelside.CenterX(Xcenter);
    REAL xcenter = Xcenter[0];
    REAL ycenter = Xcenter[1];
    
    for (int i=0; i<nnodes; i++) {//compute the solution at the vertices
        TPZGeoNode *inode = gel->NodePtr(i);
        REAL xcoord = inode->Coord(0);
        REAL ycoord = inode->Coord(1);
        allSlxl[i] = cellSol+((gradCell[0]*(xcoord-xcenter))+(gradCell[1]*(ycoord-ycenter)));
    }
    
    //get max and min
    neighsSols.push_back(cellSol);//Add cellSol to allneighsSol
    double maxSat = neighsSols[0];
    double minSat = neighsSols[0];
    int nneighs = neighsSols.size();
    
    for (int j=0; j<nneighs; j++) {
        double satNeigh = neighsSols[j];
        if (satNeigh>maxSat) {
            maxSat=satNeigh;
        }
        if (satNeigh<minSat) {
            minSat=satNeigh;
        }
    }
    
    TPZVec<double> allYs;
    double yval;
    
    for (int i=0; i<nnodes; i++) {
        auto slxlVal = allSlxl[i];
        double checkval = std::abs(slxlVal-cellSol);
        
        if (checkval<tolerance) {
            yval = 1.0;
        }
        else{
            if (slxlVal>cellSol) {
                yval = (maxSat-cellSol)/(slxlVal-cellSol);
            }
            if (slxlVal<cellSol) {
                yval = (minSat-cellSol)/(slxlVal-cellSol);
            }
        }
        double ybar = ((yval*yval)+2*yval)/((yval*yval)+yval+2);
        allYs.push_back(ybar);
    }
    allYs.push_back(3.0/4);
    
    int nyvals = allYs.size();
    
    //find limiter val AlphaK
    double alphaK = allYs[0];
    for (int j = 0; j< nyvals; j++) {
        double alphaKVal = allYs[j];
        if (alphaKVal<alphaK) {
            alphaK=alphaKVal;
        }
    }
    return alphaK;
}

int TMRSTransportAnalysis::IdentifyCase(TPZGeoEl *gel){
    int nnodes = gel->NNodes();
    int nsides = gel->NSides();
    int elDimen = gel->Dimension();
    int counter = 0;
    
    for (int iside = nnodes; iside<nsides-1; iside++) {
        TPZGeoElSide gelSide (gel, iside);
        int neighs_same_dim = gelSide.NNeighbours(elDimen);
        if (neighs_same_dim==1) {
            counter++;
        }
    }
    
}



void TMRSTransportAnalysis::GetNeigsCenterAndSolTest(TPZGeoEl *gel, TPZStack<std::pair<TPZManVector<REAL,3>, REAL>> &CenterAndSol, std::function<void(const TPZVec<REAL> &, TPZVec<STATE> &, TPZFMatrix<STATE> &)> SatFunction, int istep, TPZVec<double> AllCellSols) {
  int dim = gel->Mesh()->Dimension();
    
    int eldim = gel->Dimension();
    int nnodes = gel->NNodes();
    int nsides = gel->NSides();
    int geoindexel = gel->Index();
    TPZVec<int> NeighIds;
    TPZStack<TPZGeoElSide> ElNeighs;
    int mistep = m_sim_data->mTNumerics.m_istep;
    for (int iside = nnodes; iside < nsides - 1; iside++)
    {
      TPZGeoElSide gelside(gel, iside);
      int neighs_same_dim = gelside.NNeighbours(eldim);
      if (neighs_same_dim == 1)
      {
        TPZStack<TPZGeoElSide> allneigs;
        gelside.AllNeighbours(allneigs);
        int nneighs = allneigs.size();

        for (int ineig = 0; ineig < nneighs; ineig++)
        {
          TPZGeoElSide elsideNeigh = allneigs[ineig];
          TPZGeoEl *gelNeigh = elsideNeigh.Element();
          int neighDim = gelNeigh->Dimension();
          if (neighDim == eldim)
          {
            ElNeighs.push_back(elsideNeigh);
            int neighIndex = gelNeigh->Index();
            NeighIds.push_back(neighIndex);
          }
        }
      }
    }
    int nneighs = NeighIds.size();

//    CenterAndSol.resize(nneighs);

    for (int ineigh = 0; ineigh < nneighs; ineigh++)
    {
      auto igel = ElNeighs[ineigh].Element();
      int geoIndexNeigh = NeighIds[ineigh];

      int nsides = igel->NSides();
      int Neighdim = igel->Dimension();
      TPZVec<REAL> masscent(Neighdim, 0);
      TPZVec<REAL> centerNeigh(3, 0);
      igel->CenterPoint(nsides - 1, masscent);
      igel->X(masscent, centerNeigh);

      if (eldim == 2)
      {
        centerNeigh.resize(eldim);
      }

      REAL neighSol = AllCellSols[geoIndexNeigh];
        
      CenterAndSol.push_back(std::make_pair(centerNeigh, neighSol));
        int AllNeighs=CenterAndSol.size();
        auto pair = CenterAndSol[0];
        auto centertest = pair.first;
//        std::cout<<centertest[0]<<centertest[1]<<std::endl;
        int stop=1;

    }
      int firstside = gel->FirstSide(dim-1);
      int lastside = gel->NSides()-1;
      for (int iside = firstside; iside<lastside; iside++) {
          TPZGeoElSide gelside (gel,iside);
          int neigs_same_dim = gelside.NNeighbours(gelside.Dimension());
          if (neigs_same_dim!=0) {
              TPZStack<TPZGeoElSide> allneigh;
              gelside.AllNeighbours(allneigh);
              for (int ineig = 0; ineig<allneigh.size(); ineig++) {
                  TPZGeoElSide neighElSide = allneigh[ineig];
                  auto elNeigh = neighElSide.Element();
                  auto matIdNeigh = elNeigh->MaterialId();
                  if (matIdNeigh==2) {//Inlet Element
//                      TPZGeoElSide GelSideTest(gel,lastside);
//                      TPZVec<REAL> XcenterTest(3);
//                      GelSideTest.CenterX(XcenterTest);
//                      std::cout<<XcenterTest[0]<<" "<<XcenterTest[1]<<std::endl;
                      
                      TPZVec<REAL> Xcenter(3);
                      neighElSide.CenterX(Xcenter);
                      double dt = fAlgebraicTransport.fdt;
                      Xcenter[2] = dt*istep;
                      auto timevalue = dt*istep;
                      TPZVec<STATE> rhsVal(1);
                      TPZFMatrix<STATE> matVal;
                      SatFunction(Xcenter,rhsVal,matVal);
                      auto x= Xcenter[0];
                      auto y= Xcenter[1];
                      Xcenter.resize(2);
                      auto boundaryval = rhsVal[0];
                      CenterAndSol.push_back(std::make_pair(Xcenter, rhsVal[0]));
                      
                      break;
                  }
              }
          }
      }
  }
void TMRSTransportAnalysis:: GradientCell(TPZVec<REAL> &CenterCell, REAL CellSol, TPZStack<std::pair<TPZManVector<REAL,3>, REAL>> &CenterAndSol, TPZVec<REAL> &GradientSol, int dim)
{

  int nneighs = CenterAndSol.size();

  TPZFMatrix<REAL> dxdy(nneighs, dim);
  TPZFMatrix<REAL> ds(nneighs);
  for (int ineigh = 0; ineigh < nneighs; ineigh++)
  {
    auto pair = CenterAndSol[ineigh];
    TPZVec<REAL> CenterNeigh = pair.first;
    REAL NeighSol = pair.second;

    if (dim == 2)
    {
      auto xcoordNeigh = CenterNeigh[0];
      auto ycoordNeigh = CenterNeigh[1];
      dxdy(ineigh, 0) = xcoordNeigh - CenterCell[0];
      dxdy(ineigh, 1) = ycoordNeigh - CenterCell[1];
    }
    if (dim == 3)
    {
      auto xcoordNeigh = CenterNeigh[0];
      auto ycoordNeigh = CenterNeigh[1];
      auto zcoordNeigh = CenterNeigh[2];
      dxdy(ineigh, 0) = xcoordNeigh - CenterCell[0];
      dxdy(ineigh, 1) = ycoordNeigh - CenterCell[1];
      dxdy(ineigh, 2) = zcoordNeigh - CenterCell[2];
    }

    ds(ineigh) = NeighSol - CellSol;
  }
  TPZFNMatrix<4, REAL> AtA, B, AtAInv, GradSol;
  auto dxdyTemp = dxdy;
  dxdy.Transpose();
  dxdy.Multiply(dxdyTemp, AtA);
  dxdy.Multiply(ds, B);
  AtA.Inverse(AtAInv, ELU);
  AtAInv.Multiply(B, GradSol);
  //    GradSol.Print(std::cout);
  if (dim == 2)
  {
    GradientSol.resize(dim);
    GradientSol[0] = GradSol(0, 0);
    GradientSol[1] = GradSol(1, 0);
  }
  if (dim == 3)
  {
    GradientSol.resize(dim);
    GradientSol[0] = GradSol(0, 0);
    GradientSol[1] = GradSol(1, 0);
    GradientSol[2] = GradSol(2, 0);
  }
}

