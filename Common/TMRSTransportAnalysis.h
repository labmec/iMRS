//
//  TMRSTransportAnalysis.h
//
//  Created by Omar Durán on 10/15/19.
//

#ifndef TMRSTransportAnalysis_h
#define TMRSTransportAnalysis_h

#include <stdio.h>
#include "TPZLinearAnalysis.h"

#include "TPZMultiphysicsCompMesh.h"
#include "TPZSpStructMatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TMRSDataTransfer.h"
#include "TPZMFSolutionTransfer.h"
#include "TPZAlgebraicTransport.h"

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
//#include <Eigen/PardisoSupport>
#include "TPZAnalysisAuxEigen.h"

template<typename StorageIndex = typename Eigen::SparseMatrix<REAL>::StorageIndex >
class Triplet
{
public:
  
 Triplet() : m_row(0), m_col(0), m_value(0) {}

  Triplet(const StorageIndex& i, const StorageIndex& j, const REAL& v = REAL(0))
    : m_row(i), m_col(j), m_value(v)
  {}

    const StorageIndex& row() const { return m_row; }
    const StorageIndex& col() const { return m_col; }
    const REAL & value() const { return m_value; }
    
protected:
    
  StorageIndex m_row, m_col;
    
  REAL m_value;
};

class TMRSTransportAnalysis : public TPZLinearAnalysis {
    
private:
    
    TPZAnalysisAuxEigen *fTransportSpMatrix;
    /// Data transfer object
    TMRSDataTransfer * m_sim_data;
    
    int fpostprocessindex=1;
    
    /// Number of iterations
    int m_k_iteration = 0;
    
    REAL m_current_time;
    
    TPZFMatrix<REAL> F_inlet ;
    
    TPZFMatrix<STATE>  M_diag;
    
    bool m_parallel_execution_Q = false;
    
    Eigen::SparseMatrix<REAL> m_mass;
    Eigen::SparseMatrix<REAL> m_transmissibility;
    Eigen::SparseMatrix<REAL> m_rhs;
    
    std::vector< Triplet<REAL> >           m_trans_triplets;
    std::vector< Triplet<REAL> >           m_rhs_triplets;
    std::vector< Triplet<REAL> >           m_mass_triplets;
//    Eigen::PardisoLU<Eigen::SparseMatrix<REAL>>  m_analysis;
    Eigen::SparseLU<Eigen::SparseMatrix<REAL>>  m_analysis;
public:
    
    TPZAlgebraicTransport fAlgebraicTransport;
    
    TPZMFSolutionTransfer m_soltransportTransfer;

    /// Default constructor
    TMRSTransportAnalysis();
    
    /// Default destructor
    ~TMRSTransportAnalysis();
    
    /// Constructor based on a cmesh and optimization band directive
    TMRSTransportAnalysis(TPZCompMesh * cmesh_mult,
                          const RenumType& renumtype = RenumType::EDefault);
    
    /// Configurates iternal members
    void Configure(int n_threads, bool UsePardiso_Q);
    void Assemble();
    void AssembleResidual();
    
    /// Set data transfer object
    void SetDataTransfer(TMRSDataTransfer * sim_data);
    
    //set Current time
    void SetCurrentTime(REAL time){m_current_time = time;};
    
    //Get Current time
    REAL GetCurrentTime(){return m_current_time;};
    
    /// Get data transfer object
    TMRSDataTransfer * GetDataTransfer();
    
    /// Get the number of iterations
    int GetNumberOfIterations();
    
    /// Run a time step
    void RunTimeStep();
    
    /// Render a vtk file with requested variables for a time step
    void PostProcessTimeStep();
    
    /// Perform a Newton iteration
    void NewtonIteration();
    
    /// Perform a Newton iteration pz based.
    void NewtonIteration_serial();
    
    void Assemble_serial();
    
    void AssembleResidual_serial();
    
    /// Perform a Newton iteration eigen based.
    void NewtonIteration_Eigen();
    
    void Assemble_mass_parallel();
    void Assemble_parallel();
    void AssembleResidual_Eigen();
    void AnalyzePattern();
    
    void ComputeInitialGuess(TPZFMatrix<STATE> &x);
    
    bool QuasiNewtonSteps(TPZFMatrix<STATE> &x, int n);
    
    void UpdateInitialSolutionFromCellsData();
    
    //-------------------- FUNCTIONS  --------------------

    void RGK1(TPZFMatrix<STATE> &sn, TPZFMatrix<STATE> &snp1);
    void dxdydsFinder( int problemdim, int indexcel, TPZVec<int> GeoToCellIndex, TPZStack<TPZGeoElSide> ElNeighs,TPZFMatrix< REAL> &dxdyel, TPZFMatrix< REAL> &dsel);
    TPZStack<TPZGeoElSide> FindElNeighs(TPZGeoEl *gel, bool sameDim = true);
    
    TPZStack<TPZGeoElSide > FindAllElNeighs(TPZGeoEl *gel );
    void FindElNeigsUp(TPZGeoEl *gel, TPZStack<TPZGeoElSide> &ElNeighs, TPZVec<int> &NeighsId);
    void MaxMinMatrix(TPZCompMesh *transportMesh, TPZFMatrix<double> &result);
    void GradientLimiter1D(TPZCompMesh *transportMesh, int xlsdim, TPZFMatrix<double>&result);
    void GradientLimiter2D(TPZCompMesh *transportMesh, TPZFMatrix<double>&result);
    TPZFMatrix<REAL>GetXlsEl(TPZGeoEl *gel);
    TPZVec<double>FindMaxMin(TPZVec<int> GeoToIndex, TPZVec<int> NeighsId);
    void GradientReconstruction1D(TPZCompMesh * transportMesh);
    void GradientReconstruction2D(TPZCompMesh *transportMesh);
    void UPDATERGK1(TPZFMatrix<STATE> &sn, TPZFMatrix<STATE> &snp1, TPZCompMesh *TransportMesh);
    void FindAlphaK(double SlxsR,double SlxsL, double Skmax, double Skmin, double cellSat, double &AlphaKResult);
    void GetCellData(TPZCompMesh *TransportMesh, TPZFMatrix<double> &Result);
    void UpdatedGradientReconstruction1D(TPZFMatrix<double> snp1,TPZCompMesh *TransportMesh,TPZFMatrix<double> &Result);
    void GetNeigsCentr(TPZGeoEl *gel , TPZFMatrix<REAL> &MatResult);
    void GetNeigsSol(TPZGeoEl *gel, TPZVec<int> GeoToCellIndex,TPZVec<double> &Result);
    TPZVec<double>GradientCell(std::vector<REAL> CenterCell, double CellSol, TPZFMatrix<REAL> CenterNeighs, TPZVec<double> NeighsSols);
    double GradientLimiter(TPZGeoEl *gel, double cellSol, TPZVec<double> gradCell, TPZVec<double> neighsSols);
    int IdentifyCase(TPZGeoEl *gel);
    
    void GetNeigsCenterAndSol(TPZGeoEl *gel, TPZStack<std::pair<TPZManVector<REAL,3>, REAL>> &CenterAndSol, std::function<void(const TPZVec<REAL> &, TPZVec<STATE> &, TPZFMatrix<STATE> &)> SatFunction);
    void GetNeigsCenterAndSolTest(TPZGeoEl *gel, TPZStack<std::pair<TPZManVector<REAL,3>, REAL>> &CenterAndSol, std::function<void(const TPZVec<REAL> &, TPZVec<STATE> &, TPZFMatrix<STATE> &)> SatFunction, int istep, TPZVec<REAL> AllCellSols);
    void GradientCell(TPZVec<REAL> &CenterCell, REAL CellSol, TPZStack<std::pair<TPZManVector<REAL,3>, REAL>> &CenterAndSol, TPZVec<REAL> &GradientSol, int dim);

}; 

#endif /* TMRSTransportAnalysis_h */


