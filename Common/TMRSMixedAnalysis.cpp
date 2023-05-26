//
//  TMRSMixedAnalysis.cpp
//
//  Created by Omar Durán on 10/15/19.

#include "TMRSMixedAnalysis.h"
#include "TPZSpStructMatrix_Eigen.h"
#include "TPZSpMatrixEigen.h"
#include <pzshapequad.h>
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapetetra.h"
#include "pzshapecube.h"
#include "pzshapeprism.h"
#include "pzshapepiram.h"
#include "pzshapepoint.h"
#include "TPZCompElHDivCollapsed.h"

// Uses the new vtk function developed by Fran
#define USENEWVTK

#ifdef USENEWVTK
#include "TPZVTKGenerator.h"
#endif

using namespace std;

TMRSMixedAnalysis::TMRSMixedAnalysis(){
    fmixed_report_data = new std::ofstream("Report_Mixed.txt");
    
}

TMRSMixedAnalysis::~TMRSMixedAnalysis(){
    
}

TMRSMixedAnalysis::TMRSMixedAnalysis(TPZMultiphysicsCompMesh * cmesh_mult, bool must_opt_band_width_Q) : TPZLinearAnalysis(cmesh_mult, must_opt_band_width_Q){
    fsoltransfer.BuildTransferData(cmesh_mult);
    fmixed_report_data = new std::ofstream("Report_Mixed.txt");
}

void TMRSMixedAnalysis::SetDataTransfer(TMRSDataTransfer * sim_data){
    m_sim_data = sim_data;
    fmixed_report_data = new std::ofstream("Report_Mixed.txt");
    
}

TMRSDataTransfer * TMRSMixedAnalysis::GetDataTransfer(){
    return m_sim_data;
}

int TMRSMixedAnalysis::GetNumberOfIterations(){
    return m_k_iteration;
}

void TMRSMixedAnalysis::Configure(int n_threads, bool UsePardiso_Q,bool UsePZ){
    
    if (UsePardiso_Q) {
        if(UsePZ){
            TPZSSpStructMatrix<STATE> matrix(Mesh());
            matrix.SetNumThreads(n_threads);
            SetStructuralMatrix(matrix);
            TPZStepSolver<STATE> step;
            step.SetDirect(ELDLt);
            SetSolver(step);
        }
        else{
            TPZSymetricSpStructMatrixEigen matrix(Mesh());
            matrix.SetNumThreads(n_threads);
            SetStructuralMatrix(matrix);
            TPZStepSolver<STATE> step;
            step.SetDirect(ELDLt);
            SetSolver(step);
        }


    }else{
        TPZSkylineStructMatrix<STATE> matrix(Mesh());
        matrix.SetNumThreads(n_threads);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        SetSolver(step);
        SetStructuralMatrix(matrix);
    }
    
   
    //    Assemble();
}

void TMRSMixedAnalysis::RunTimeStep(){
    
    TPZMultiphysicsCompMesh * cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
    if (!cmesh)
        DebugStop();
    
    int n = m_sim_data->mTNumerics.m_max_iter_mixed;
    bool stop_criterion_Q = false;
    bool stop_criterion_corr_Q = false;
    REAL res_norm = 1.0;
    REAL corr_norm = 1.0;
    REAL res_tol = m_sim_data->mTNumerics.m_res_tol_mixed;
    REAL corr_tol = m_sim_data->mTNumerics.m_corr_tol_mixed;
    TPZFMatrix<STATE> dx,x(Solution());
    for(m_k_iteration = 1; m_k_iteration <= n; m_k_iteration++){
        
        dx = Solution();
        corr_norm = Norm(dx);
        res_norm = Norm(Rhs());
        
        NewtonIteration();
        REAL fistAssemTime = flastAssembleTime;
        dx = Solution();
        corr_norm = Norm(dx);
        res_norm = Norm(Rhs());
        x +=dx;
//        LoadSolution(x);
        cmesh->UpdatePreviousState(1.0);
        fsoltransfer.TransferFromMultiphysics();
        
        
        Assemble();
        
//        TPZMatrixSolver<STATE> *matsol = dynamic_cast<TPZMatrixSolver<STATE> *>(fSolver);
//        matsol->Matrix()->SetIsDecomposed(0);
        
//        this->Solver()->ResetMatrix();
    
//        this->Solver()->
          
//        Solver()->Matrix()->SetIsDecomposed(0);
//        Solve();
        res_norm = Norm(Rhs());
        REAL normsol = Norm(Solution());
        if(m_k_iteration == 1){
            (*fmixed_report_data) <<"         "<< m_k_iteration<<"     "<<res_norm<<"     "<<normsol<< "         "<<fistAssemTime<<"         "<< flastSolveTime<<"         " << flastAssembleTime<<std::endl;
        }
        else{
            (*fmixed_report_data)<<"   M            "<< m_k_iteration<<"     "<<res_norm<<"     "<<normsol<< "         "<<fistAssemTime<<"         "<< flastSolveTime<<"         " << flastAssembleTime<<std::endl;
        }
        
#ifdef PZDEBUG
        {
            if(std::isnan(corr_norm) || std::isnan(res_norm))
                DebugStop();
        }
#endif

        stop_criterion_Q = res_norm < res_tol;
        stop_criterion_corr_Q = corr_norm < corr_tol;
        if (stop_criterion_Q || stop_criterion_corr_Q) {
            std::cout << "\n\n\t================================================" << std::endl;
            std::cout << "Mixed operator: " << std::endl;
            std::cout << "Iterative method converged with res_norm = " << res_norm << std::endl;
            std::cout << "Number of iterations = " << m_k_iteration << std::endl;
            break;
        }
        else if (m_k_iteration == n) {
            std::cout << "Mixed operator: " << std::endl;
            std::cout<<"Mixed not converge =( "<<std::endl;
            std::cout << "Iterative method converged with res_norm = " << res_norm << std::endl;
            std::cout << "Number of iterations = " << m_k_iteration << std::endl;
//            DebugStop();
        }
    }
}


void TMRSMixedAnalysis::NewtonIteration(){
    
    if(mIsFirstAssembleQ == true)
    {
        fStructMatrix->SetNumThreads(0);
        int64_t nel = fCompMesh->NElements();
        for(int64_t el = 0; el<nel; el++)
        {
            TPZCompEl *cel = fCompMesh->Element(el);
            TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
            if(sub)
            {
            int numthreads = 0;
//                sub->SetAnalysisSparse(0); sub->Analysis()->StructMatrix()->SetNumThreads(m_sim_data->mTNumerics.m_nThreadsMixedProblem);
                TPZSSpStructMatrix<STATE> matrix(sub);
                matrix.SetNumThreads(numthreads);
                sub->Analysis()->SetStructuralMatrix(matrix);
//                TPZStepSolver<STATE> step;
//                step.SetDirect(ELDLt);
//                sub->Analysis()->SetSolver(step);
            }
        }
        mIsFirstAssembleQ=false;
    }
    else{
        fStructMatrix->SetNumThreads(m_sim_data->mTNumerics.m_nThreadsMixedProblem);
        int64_t nel = fCompMesh->NElements();
        for(int64_t el = 0; el<nel; el++)
        {
            TPZCompEl *cel = fCompMesh->Element(el);
            TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
            if(sub)
            {
            int numthreads = m_sim_data->mTNumerics.m_nThreadsMixedProblem;
                TPZSSpStructMatrix<STATE> matrix(sub);
                matrix.SetNumThreads(numthreads);
                sub->Analysis()->SetStructuralMatrix(matrix);
            }
        }
    }
//    std::string outputFolder("NOSE");
//    FilterZeroNeumann( outputFolder, m_sim_data, this->StructMatrix(), Mesh());
    Assemble();
    

    
    Solve();
    //   this->MatrixSolver<REAL>().ResetMatrix();

    
//    this->PostProcessTimeStep();
}

void TMRSMixedAnalysis::PostProcessTimeStep(int dimToPost){
    
    const int dim = Mesh()->Dimension();
    auto start_time_pp = std::chrono::steady_clock::now();
    cout << "\n--------------------- Post process dim = " << dimToPost << " ---------------------\n" << endl;
        
    TPZStack<std::string,10> scalnames, vecnames;
    
    scalnames = m_sim_data->mTPostProcess.m_scalnamesDarcy;
    vecnames = m_sim_data->mTPostProcess.m_vecnamesDarcy;
    
    int div = 0;
    if (dimToPost < 0){
        dimToPost = Mesh()->Reference()->Dimension();
    }
//    dim = 2;
//    std::set<int> mat_id_2D;
//    mat_id_2D.insert(10);
//    std::string file_frac("fractureFlux_s.vtk");
//    DefineGraphMesh(2,mat_id_2D,scalnames,vecnames,file_frac);
//    PostProcess(div,2);
    std::string file = m_sim_data->mTPostProcess.m_file_name_mixed;

    constexpr int vtkRes{0}; //resolucao do vtk
    if (dimToPost == dim-1){
        std::set<int> matids;
        map<int, TMRSDataTransfer::TFracProperties::FracProp>::iterator it;
        for (it = m_sim_data->mTFracProperties.m_fracprops.begin(); it != m_sim_data->mTFracProperties.m_fracprops.end(); it++)
        {
            int matfracid = it->first;
            matids.insert(matfracid);
        }

#ifdef USENEWVTK
        const std::string plotfile = file.substr(0, file.find(".")) + "_frac";
        for (auto nm : vecnames) {
            scalnames.Push(nm);
        }
        auto vtk = TPZVTKGenerator(fCompMesh, matids, scalnames, plotfile, vtkRes);
        vtk.SetNThreads(8);
        vtk.Do();
#else
        file = file.substr(0, file.find(".")) + "_frac.vtk";
        DefineGraphMesh(dimToPost, matids, scalnames, vecnames, file);
        PostProcess(div,dimToPost);
#endif
    }
    else{
#ifdef USENEWVTK
        std::string plotfile = file.substr(0, file.find(".")); //sem o .vtk no final
        for (auto nm : vecnames) {
            scalnames.Push(nm);
        }

        auto vtk = TPZVTKGenerator(fCompMesh, scalnames, plotfile, vtkRes, dimToPost);
        vtk.SetStep(fpostprocessindex);
        vtk.SetNThreads(8);
        vtk.Do();
        fpostprocessindex++;
#else
        DefineGraphMesh(dimToPost,scalnames,vecnames,file);
        PostProcess(div,dimToPost);
#endif
    }
    
    
    
    auto total_time_pp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_pp).count()/1000.;
    cout << "Total time post process = " << total_time_pp << " seconds" << endl;

}

void TMRSMixedAnalysis::Assemble(){

    auto start_time_ass = std::chrono::steady_clock::now();

    cout << "\n---------------------- Assemble Flux Problem ----------------------" << endl;
    cout << "Number of equations: " << fCompMesh->NEquations() << endl;
    cout << "Number of elements: " << fCompMesh->NElements() << endl;
    TPZLinearAnalysis::Assemble();


    auto total_time_ass = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_ass).count()/1000.;
    cout << "\nTotal time assemble = " << total_time_ass << " seconds" << endl;
    flastAssembleTime=total_time_ass;
}

void TMRSMixedAnalysis::Solve(){

    auto start_time_solve = std::chrono::steady_clock::now();
    
    cout << "\n---------------------- Solve Flux Problem ----------------------" << endl;
    TPZLinearAnalysis::Solve();
    
    auto total_time_solve = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_solve).count()/1000.;
    cout << "Total time solve = " << total_time_solve << " seconds" << endl;
    flastSolveTime = total_time_solve;
}

void TMRSMixedAnalysis::VerifyElementFluxes(){
    const REAL tol = 1.e-7;
	TPZMultiphysicsCompMesh *mixedmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh()) ;
	TPZCompMesh *cmesh = mixedmesh->MeshVector()[0];
//	std::ofstream file("fuxmesh.txt");
	const TPZFMatrix<STATE> &meshSol = cmesh->Solution();
//	mixedmesh->MeshVector()[0]->Print(file);
	int nels =mixedmesh->MeshVector()[0]->NElements();
	for (int iel =0; iel<nels-1; iel++) {
		TPZCompEl *cel = mixedmesh->MeshVector()[0]->Element(iel);
		if (!cel) continue;
//            if(cel->Dimension() != cmesh->Dimension()) continue;
		TPZCompElHDiv<pzshape::TPZShapeCube> *hdivel = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeCube> *>(cel);
		TPZCompElHDiv<pzshape::TPZShapeTetra> *hdiveltet = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeTetra> *>(cel);
		TPZCompElHDiv<pzshape::TPZShapeQuad> *hdivelq = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeQuad> *>(cel);
		TPZCompElHDiv<pzshape::TPZShapeTriang> *hdivelqTri = dynamic_cast<TPZCompElHDiv<pzshape::TPZShapeTriang> *>(cel);
		TPZCompElHDivCollapsed<pzshape::TPZShapeQuad> *hdivCollaps = dynamic_cast<TPZCompElHDivCollapsed<pzshape::TPZShapeQuad> *>(hdivelq);
		TPZCompElHDivCollapsed<pzshape::TPZShapeTriang> *hdivCollapsTri = dynamic_cast<TPZCompElHDivCollapsed<pzshape::TPZShapeTriang> *>(hdivelqTri);
		if (!hdivel && !hdiveltet && !hdivCollaps && !hdivCollapsTri) {
			continue;
		}
		TPZInterpolatedElement* intel = dynamic_cast<TPZInterpolatedElement*>(cel);
		
		int ncon = cel->NConnects();
		int nCorners = cel->Reference()->NCornerNodes();
		int nsides1 = cel->Reference()->NSides(1);
		int nsides = cel->Reference()->NSides();
		REAL sumel=0.0;
		for (int icon=0; icon<ncon-1; icon++) {
			TPZConnect &con = cel->Connect(icon);
			int sideOrient =0;
			const bool isColaps = hdivCollaps || hdivCollapsTri;
			const int iconNum = hdivCollaps ? 5 : 4;
			const int sideForEl = hdivCollaps ? 8 : 6;
			if(!isColaps){
				 sideOrient = intel->GetSideOrient(nCorners+nsides1+icon);
			}
			else{
				if (isColaps && icon < iconNum) {
					sideOrient = intel->GetSideOrient(nCorners+icon);
				}
				 
			}
			if (isColaps && icon == iconNum-1) {
				ncon++;
				continue;
			}
			if (isColaps && icon == iconNum) {
				sideOrient = intel->GetSideOrient(sideForEl);
			}
			if (isColaps && icon == iconNum+1) {
				sideOrient = intel->GetSideOrient(sideForEl+1);
			}
			
		   
			int sequence =con.fSequenceNumber;
			int64_t pos = cmesh->Block().Position(sequence);
			if(sequence==-1) continue;
			int blocksize = cmesh->Block().Size(sequence);
			for(int ieq=0; ieq< blocksize; ieq++)
			{
                sumel += sideOrient*meshSol.GetVal(cmesh->Block().Index(sequence,ieq),0);
			}
		}
		if(std::abs(sumel)> tol ){
            std::cout << "\n\nERROR! Conservation of element index " << cel->Reference()->Index() << " is " << sumel << std::endl;
//			DebugStop();
		}
	}
    std::cout << "\n\n===> Nice! All flux elements satisfy conservation up to tolerance " << tol << std::endl;
}
void TMRSMixedAnalysis::AllZero(TPZCompMesh *cmesh){
    
//    std::ofstream filetoprint("malha_antes.txt");
//    cmesh->Print(filetoprint);
    
    cmesh->Solution().Zero();
    cmesh->SolutionN().Zero();
    int count=0;
    TPZMultiphysicsCompMesh *multcmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(cmesh);
    if(multcmesh){
        multcmesh->Solution().Zero();

        multcmesh->ComputeNodElCon();
        multcmesh->CleanUpUnconnectedNodes();
        fsoltransfer.TransferFromMultiphysics();
        for (auto atomicmesh: multcmesh->MeshVector()) {
//            std::cout<<count<<std::endl;
            count++;
            atomicmesh->Solution().Zero();
        }
    }
    int nels = cmesh->NElements();
    for (int iel=0; iel<nels; iel++) {
        TPZCompEl *cel = cmesh->Element(iel);
        if(!cel) continue;
        TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *> (cel);
        if(!sub) continue;
        AllZero(sub);
    }
}
void TMRSMixedAnalysis::FilterZeroNeumann(std::string& outputFolder, TMRSDataTransfer* sim_data, TPZAutoPointer<TPZStructMatrix> strmat, TPZCompMesh* cmesh) {

    std::cout << "\n---------------------- Filtering zero neumann equations ----------------------" << std::endl;
//    TPZSimpleTimer timer_filter("Timer Filter Equations");
    
   
    std::set<int64_t> matidset;
    
    // First find all the zero neumann in in the 3d domain
    std::cout << "Domain BC matids: ";
    for(auto &chunk : sim_data->mTBoundaryConditions.mBCFlowMatIdToTypeValue) {
        const int bc_id = chunk.first;
        const std::pair<int,REAL>& typeAndVal = chunk.second;
        const int bc_type = typeAndVal.first;
        const REAL val = typeAndVal.second;
        if(bc_type == 1 && fabs(val) < ZeroTolerance()){
            std::cout << bc_id << " ";
            matidset.insert(bc_id);
        }
    }
    
    // Then all the zero neumann in the fractures
    std::cout << "\nFracture BC matids: ";
    for (auto& chunk : sim_data->mTBoundaryConditions.mBCFlowFracMatIdToTypeValue) {
        const int bc_id   = chunk.first;
        const std::pair<int,REAL>& typeAndVal = chunk.second;
        const int bc_type = typeAndVal.first;
        const REAL val = typeAndVal.second;
        if(bc_type == 1 && fabs(val) < ZeroTolerance()){
            std::cout << bc_id << " ";
            matidset.insert(bc_id);
        }
    }
    
    // Set equations to be filted
    std::set<int64_t> eqset;
    cmesh->GetEquationSetByMat(matidset, eqset);
    if (eqset.size()) {
        strmat->EquationFilter().ExcludeEquations(eqset);
    }
    
    int count = 0;
    for(auto cel : cmesh->ElementVec()){
        TPZSubCompMesh* subcmesh = dynamic_cast<TPZSubCompMesh*>(cel);
        if (subcmesh) {
//            std::cout << "\n\t------- Submesh " << count << " -------" << std::endl;
//            if(count == 1 || count == 3){
//                count++;
//                std::cout << "===> Skipping filter" << std::endl;
//                continue;
//            }
            count++;
#ifdef PZDEBUG
            {
                std::string filename = outputFolder + "submesh_" + to_string(count) + ".vtk";
                std::ofstream out(filename);
//                TPZVTKGeoMesh::PrintCMeshVTK(subcmesh, out);
            }
#endif
            std::set<int64_t> eqsetsub;
            subcmesh->GetEquationSetByMat(matidset, eqsetsub);
            const int64_t ninteq = subcmesh->NumInternalEquations();
            const int64_t neq = subcmesh->TPZCompMesh::NEquations();
            if (eqsetsub.size()) {
                subcmesh->Analysis()->StructMatrix()->EquationFilter().SetNumEq(neq); // Setting again bcz els could have been condensed
                subcmesh->Analysis()->StructMatrix()->EquationFilter().ExcludeEquations(eqsetsub); // will set as active all eqs (internal and external) that are not zero neumann
                auto& activeeq = subcmesh->Analysis()->StructMatrix()->EquationFilter().GetActiveEquations();
                
//                std::cout << "size eqsetsub = " << eqsetsub.size() << " | eqsetsub = " << eqsetsub;
//                std::cout << "size activeeq = " << activeeq.size() << " | active eq = " << activeeq << std::endl;
                auto& excludedeq = subcmesh->Analysis()->StructMatrix()->EquationFilter().GetExcludedEquations();
//                std::cout << "size excludedeq = " << excludedeq.size() << " | excluded eq = " << excludedeq << std::endl;
                const int64_t nexteq = neq - ninteq;
#ifdef PZDEBUG
                // All the external equations must be present in the active equations
                int64_t count = 0;
                for(auto& eq : activeeq){
                    if(eq > ninteq-1) count++;
                }
                if(count < nexteq) DebugStop();
#endif
                // Now we want to set as active only the ones that are internal and not zero Neumann
                // Since activeeq is ordered, we simply pick equations until we reach the first non internal equation
                TPZVec<int64_t> activeinternal(ninteq);
                int64_t i = 0;
                for(auto& eq : activeeq){
                    if(eq > ninteq-1) break;
                    activeinternal[i++] = eq;
                }
                activeinternal.Resize(i);
                subcmesh->Analysis()->StructMatrix()->EquationFilter().Reset();
                subcmesh->Analysis()->StructMatrix()->EquationFilter().SetActiveEquations(activeinternal); // sets new filter
                // Note that, in subcmesh, after we create the data structure of the matred matrices, we need to the
                // set the external equations as active again
            }
        }
    }
    
//    std::cout << "\n==> Total Filter time: " << timer_filter.ReturnTimeDouble()/1000. << " seconds" << std::endl;
}
