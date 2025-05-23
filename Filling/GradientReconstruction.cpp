// C++ includes
// nothing

// PZ includes
#include <filesystem>

#include "TMRSApproxSpaceGenerator.h"
#include "TPZGenGrid3D.h"
#include "TPZSimpleTimer.h"
#include "TPZVTKGeoMesh.h"
#include "TPZVTKGenerator.h"
#include "imrs_config.h"
#include "json.hpp"
#include "pzintel.h"
#include "pzlog.h"
#include "pzsmanal.h"
#include "TMRSPropertiesFunctions.h"
#include "Material/Projection/TPZL2Projection.h"

// ----- Namespaces -----
using namespace std;
namespace fs = std::filesystem;
typedef std::function<void(const TPZVec<REAL> &, TPZVec<double> &, TPZFMatrix<double> &)>
    ExactSolFunc;

// ----- End of namespaces -----

// ----- Global vars -----
const int glob_n_threads = 0;

// ----- Functions -----
TPZGeoMesh *ReadMeshFromGmsh(TMRSDataTransfer &sim_data);

void GetNeigsCenterAndSol(TPZGeoEl *gel, TPZVec<REAL> AllCellSols, TPZStack<std::pair<TPZManVector<REAL,3>, REAL>> &CenterAndSol);
void GetNeigsCenterAndSol(TPZGeoEl *gel, TPZStack<std::pair<TPZManVector<REAL,3>, REAL>> &CenterAndSol);

void RunProblem(TPZGeoMesh *gmesh, const ExactSolFunc &exsol, TPZVec<REAL> &AllCelSols);

void FillDataTransfer(std::string filename, TMRSDataTransfer &simdata);
void GradientCell(TPZVec<REAL> &CenterCell, REAL CellSol, TPZVec<std::pair<TPZVec<REAL>, REAL>> &CenterAndSol, TPZVec<REAL> &GradientSol, int dim);
REAL GradientLimiter(TPZGeoEl *gel, REAL cellSol, TPZVec<REAL> &gradCell, TPZStack<std::pair<TPZManVector<REAL,3>, REAL>> &CenterAndSol);
TPZCompMesh *CreateCMesh(TPZGeoMesh *gmesh, const ExactSolFunc &ExactSolution);

void ComputeGradients(TPZGeoMesh *gmesh, TPZVec<REAL> &AllCellSols, TPZFMatrix<REAL> &Gradients);

void ComputeGradients(TPZCompMesh *cmesh);
void EvaluateError(TPZGeoMesh *gmesh,
                   std::function<void(const TPZVec<REAL> &, TPZVec<STATE> &, TPZFMatrix<STATE> &)> solution,TPZVec<REAL>& error);
// Definition of the left and right boundary conditions
auto solution = [](const TPZVec<REAL> &coord, TPZVec<STATE> &rhsVal, TPZFMatrix<STATE> &matVal) -> void
{
  REAL x = coord[0];
  REAL y = coord[1];
  REAL z = coord[2];

  rhsVal[0] = cos(M_PI * 0.5 * x) * cos(M_PI * 0.5 * y);
  matVal.Resize(3, 1);
  matVal(0, 0) = -0.5 * M_PI * sin(M_PI * 0.5 * x) * cos(M_PI * 0.5 * y);
  matVal(1, 0) = -0.5 * M_PI * cos(M_PI * 0.5 * x) * sin(M_PI * 0.5 * y);
  matVal(2, 0) = 0.0;
    
    
  int stop = 2;
};
// ----- Logger -----
#ifdef PZ_LOG
static TPZLogger mainlogger("imrs");
static TPZLogger fracIntersectLogger("imrs_fracIntersect");
#endif

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  string basemeshpath(FRACMESHES);
#ifdef PZ_LOG
  string logpath = basemeshpath + "/../Gradient/log4cxx.cfg";
  TPZLogger::InitializePZLOG(logpath);
  if (mainlogger.isDebugEnabled())
  {
    std::stringstream sout;
    sout << "\nLogger for Gradient problem target\n"
         << endl;
    ;
    LOGPZ_DEBUG(mainlogger, sout.str())
  }
#endif
  // =========> Read the json file and fill the data transfer object
  TMRSDataTransfer sim_data;
  FillDataTransfer(basemeshpath + "/../Filling/Gradient_Test", sim_data);

  // =========> Create GeoMesh
  TPZGeoMesh *gmesh = ReadMeshFromGmsh(sim_data);
  // TPZCheckGeom geom
  //  Print gmesh to vtk
  std::ofstream out("gmesh_up.vtk");
  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);

//    int ngels = gmesh->NElements();
//    int dim = gmesh->Dimension();
//    for (int i = 0; i<ngels; i++) {
//        TPZGeoEl *gel = gmesh->Element(i);
//
//        if (gel->Dimension()==dim) {
//            TPZVec<TPZGeoEl *> child;
//            gel->Divide(child);
//            break;
//        }
//    }
//    std::ofstream outE("gmesh_edit.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outE);

    
  auto cmesh = CreateCMesh(gmesh, solution);
  TPZManVector<REAL, 3> errorsum(3, 0.);
  cmesh->ElementSolution().Redim(cmesh->NElements(),3);
  cmesh->EvaluateError(true, errorsum);
  std::cout << "Errorsum constant approx: " << errorsum[0] << " " << errorsum[1] << " " << errorsum[2] << std::endl;
  ComputeGradients(cmesh);
  errorsum.Fill(0.);
  cmesh->EvaluateError(true, errorsum);
  std::cout << "Errorsum reconstructed: " << errorsum[0] << " " << errorsum[1] << " " << errorsum[2] << std::endl;
    
    TPZVec<REAL> error_m_1;
    ofstream out_m("gmesh_refined_0.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out_m);
    EvaluateError(gmesh, solution, error_m_1);
    
//  Refine
    TPZCheckGeom geom(gmesh);
    geom.UniformRefine(1);
    ofstream out_m1("gmesh_refined_1.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out_m1);
    TPZVec<REAL> error_m_2;
    EvaluateError(gmesh, solution, error_m_2);
    
    TPZCheckGeom geom1(gmesh);
    geom1.UniformRefine(1);
    ofstream out_m2("gmesh_refined_2.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out_m2);
    TPZVec<REAL> error_m_3;
    EvaluateError(gmesh, solution, error_m_3);
    
    TPZCheckGeom geom2(gmesh);
    geom2.UniformRefine(1);
    ofstream out_m3("gmesh_refined_3.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out_m3);
    TPZVec<REAL> error_m_4;
    EvaluateError(gmesh, solution, error_m_4);
    
    TPZCheckGeom geom3(gmesh);
    geom3.UniformRefine(1);
    ofstream out_m4("gmesh_refined_4.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out_m4);
    TPZVec<REAL> error_m_5;
    EvaluateError(gmesh, solution, error_m_5);
    
    TPZCheckGeom geom4(gmesh);
    geom4.UniformRefine(1);
    ofstream out_m5("gmesh_refined_5.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out_m5);
    TPZVec<REAL> error_m_6;
    EvaluateError(gmesh, solution, error_m_6);
    
    TPZCheckGeom geom5(gmesh);
    geom5.UniformRefine(1);
    ofstream out_m6("gmesh_refined_5.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out_m6);
    TPZVec<REAL> error_m_7;
    EvaluateError(gmesh, solution, error_m_7);
    
    ofstream errortxt("error.txt");
    errortxt<<error_m_1[0]<< " " <<error_m_2[0]<<" "<<error_m_3[0]<<" "<<error_m_4[0]<<" "<<error_m_5[0]<<" "<<" "<<error_m_4[0]<<" "<<error_m_5[0]<<" ";
    errortxt.close();
    ofstream GRerrortxt("GRerror.txt");
    GRerrortxt<<error_m_1[1]<< " " <<error_m_2[1]<<" "<<error_m_3[1]<<" "<<error_m_4[1]<<" "<<error_m_5[1]<<" "<<" "<<error_m_4[1]<<" "<<error_m_5[1]<<" ";
    GRerrortxt.close();
   
    
    
    
  return 0;
    
  TPZVec<REAL> AllCellSols;
  TPZFMatrix<REAL> AllCellSolGradients;
  RunProblem(gmesh, solution, AllCellSols);
  TPZCompMesh* cmsh = CreateCMesh(gmesh, solution);
  ComputeGradients(gmesh, AllCellSols, AllCellSolGradients);

delete gmesh;
cout << "-------------------- Simulation Finished --------------------" << endl;
return 0;
}
// ----------------- End of Main -----------------
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZGeoMesh *ReadMeshFromGmsh(TMRSDataTransfer &sim_data)
{
  // read mesh from gmsh
  TPZGeoMesh *gmesh;
  gmesh = new TPZGeoMesh();
  std::string file_name = std::string(FRACMESHES) + "/../Filling/" + sim_data.mTGeometry.mGmeshFileName;
  {
    TPZGmshReader reader;
    TPZManVector<std::map<std::string, int>, 4> stringtoint(5);
    stringtoint[2]["dom"] = sim_data.mTGeometry.mDomainNameAndMatId["dom"];

    stringtoint[1]["bcl"] = sim_data.mTBoundaryConditions.mDomainNameAndMatId["bcl"];
    stringtoint[1]["bcr"] = sim_data.mTBoundaryConditions.mDomainNameAndMatId["bcr"];
    stringtoint[1]["bct"] = sim_data.mTBoundaryConditions.mDomainNameAndMatId["bct"];
    stringtoint[1]["bcb"] = sim_data.mTBoundaryConditions.mDomainNameAndMatId["bcb"];
    reader.SetDimNamePhysical(stringtoint);
    reader.GeometricGmshMesh(file_name, gmesh);
  }

  return gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void FillDataTransfer(std::string filenameBase, TMRSDataTransfer &sim_data)
{
  using json = nlohmann::json;
  std::string filenamejson = filenameBase + ".json";

  std::ifstream filejson(filenamejson);
  json input = json::parse(filejson, nullptr, true, true); // to ignore comments in json file

  // ------------------------ Getting number of domains and fractures ------------------------
  if (input.find("Domains") == input.end())
    DebugStop();
  if (input.find("Mesh") == input.end())
    DebugStop();
  const int ndom = input["Domains"].size();
std:
  string mesh = input["Mesh"];
  sim_data.mTGeometry.mGmeshFileName = mesh;
  sim_data.mTReservoirProperties.mPorosityAndVolumeScale.resize(ndom + 1);
  int countPhi = 0;

  // ------------------------ Reading 3D Domain matids ------------------------
  std::map<int, REAL> idPerm;
  for (auto &domain : input["Domains"])
  {
    if (domain.find("matid") == domain.end())
      DebugStop();
    if (domain.find("name") == domain.end())
      DebugStop();
    if (domain.find("K") == domain.end())
      DebugStop();
    if (domain.find("phi") == domain.end())
      DebugStop();
    const int matid = domain["matid"];
    const string name = domain["name"];
    const REAL permeability = domain["K"];
    const REAL phi = domain["phi"];
    sim_data.mTGeometry.mDomainNameAndMatId[name] = matid;
    idPerm[matid] = permeability;
    sim_data.mTReservoirProperties.mPorosityAndVolumeScale[countPhi++] = std::make_tuple(matid, phi, 1.0);
  }

  sim_data.mTReservoirProperties.m_permeabilitiesbyId = idPerm;

  // ------------------------ Reading 3D Domain BC matids ------------------------
  //  if (input.find("Boundary") == input.end()) DebugStop();
  //  std::map<int, std::pair<int, REAL>>& BCFlowMatIdToTypeValue = sim_data.mTBoundaryConditions.mBCFlowMatIdToTypeValue;
  //  std::map<int, std::pair<int, ForcingFunctionBCType<REAL>>>& BCFlowMatIdToFunctionId = sim_data.mTBoundaryConditions.mBCFlowMatIdToFunctionId;
  //  std::map<int, std::pair<int, REAL>>& BCTransportMatIdToTypeValue = sim_data.mTBoundaryConditions.mBCTransportMatIdToTypeValue;
  //  std::map<int, std::pair<int, ForcingFunctionBCType<REAL>>>& BCTransportMatIdToFunctionId = sim_data.mTBoundaryConditions.mBCTransportMatIdToFunctionId;
  //  for (auto& bc : input["Boundary"]) {
  //    if (bc.find("matid") == bc.end()) DebugStop();
  //    if (bc.find("type") == bc.end()) DebugStop();
  //    if (bc.find("value") == bc.end()) DebugStop();
  //    if (bc.find("name") == bc.end()) DebugStop();
  //    const int matid = bc["matid"];
  //    const int type = bc["type"];
  //    const REAL value = bc["value"];
  //    const std::string name = bc["name"];
  //    int functionID = 0;
  //    if (bc.find("functionID") != bc.end()) {
  //      functionID = bc["functionID"];
  //    }
  //    REAL external_saturation = 0.0;
  //    int saturation_functionID = 0;
  //    if (input["Numerics"]["RunWithTransport"])
  //    {
  //      if (bc.find("ExternalSaturation") == bc.end()) DebugStop();
  //      external_saturation = bc["ExternalSaturation"];
  //      if (bc.find("SaturationFunctionID") != bc.end()) {
  //        saturation_functionID = bc["SaturationFunctionID"];
  //      }
  //
  //
  //
  //
  //    }
  //
  //
  //
  //    sim_data.mTBoundaryConditions.mDomainNameAndMatId[name] = matid;
  //    if (BCFlowMatIdToTypeValue.find(matid) != BCFlowMatIdToTypeValue.end()) DebugStop();
  //    BCFlowMatIdToTypeValue[matid] = std::make_pair(type, value);
  //    BCFlowMatIdToFunctionId[matid] = std::make_pair(functionID, forcingfunctionBC[functionID]);
  //    BCTransportMatIdToTypeValue[matid] = std::make_pair(type, external_saturation);
  //    BCTransportMatIdToFunctionId[matid] = std::make_pair(saturation_functionID, forcingfunctionBC[saturation_functionID]);
  //  }

  // ------------------------ Numerics Parameters ------------------------
  if (input.find("Numerics") != input.end())
  {
    auto numerics = input["Numerics"];
    sim_data.mTNumerics.m_run_with_transport = numerics["RunWithTransport"];

    if (sim_data.mTNumerics.m_run_with_transport)
    {
      if (numerics.find("UseGradientRec") != numerics.end())
        sim_data.mTNumerics.m_UseGradRec = numerics["UseGradientRec"];
      if (numerics.find("SolverTransportMethod") != numerics.end())
        sim_data.mTNumerics.m_TransportSolMethod = numerics["SolverTransportMethod"];
      if (numerics.find("DeltaT") == numerics.end())
        DebugStop();
      sim_data.mTNumerics.m_dt = numerics["DeltaT"];
      if (numerics.find("NSteps") == numerics.end())
        DebugStop();
      sim_data.mTNumerics.m_n_steps = numerics["NSteps"];
    }
    if (numerics.find("Gravity") == numerics.end())
      DebugStop();
    std::vector<REAL> grav(3, 0.0);
    for (int i = 0; i < 3; i++)
    {
      grav[i] = numerics["Gravity"][i];
    }
    sim_data.mTNumerics.m_gravity = grav;

    if (numerics.find("IsAxisymmetric") != numerics.end())
    {
      sim_data.mTNumerics.m_is_axisymmetric = numerics["IsAxisymmetric"];
    }
    if (numerics.find("IsLinearTrace") != numerics.end())
    {
      sim_data.mTNumerics.m_is_linearTrace = numerics["IsLinearTrace"];
    }
  }

  // ------------------------ Fluids Properties ------------------------
  if (input.find("FluidProperties") != input.end())
  {
    auto properties = input["FluidProperties"];
    if (properties.find("WaterDensity") == properties.end())
      DebugStop();
    sim_data.mTFluidProperties.mWaterDensityRef = properties["WaterDensity"];
    if (properties.find("WaterViscosity") == properties.end())
      DebugStop();
    sim_data.mTFluidProperties.mWaterViscosity = properties["WaterViscosity"];
    if (properties.find("WaterCompressibility") == properties.end())
      DebugStop();
    sim_data.mTFluidProperties.mWaterCompressibility = properties["WaterCompressibility"];
    if (properties.find("OilDensity") == properties.end())
      DebugStop();
    sim_data.mTFluidProperties.mOilDensityRef = properties["OilDensity"];
    if (properties.find("OilViscosity") == properties.end())
      DebugStop();
    sim_data.mTFluidProperties.mOilViscosity = properties["OilViscosity"];
    if (properties.find("OilCompressibility") == properties.end())
      DebugStop();
    sim_data.mTFluidProperties.mOilCompressibility = properties["OilCompressibility"];
    if (properties.find("DensityModel") == properties.end())
      DebugStop();
    if (properties["DensityModel"] == 0)
    {
      sim_data.mTFluidProperties.CreateLinearDensityFunction();
    }
    else
    {
      sim_data.mTFluidProperties.CreateExponentialDensityFunction();
    }
  }

  // ------------------------ Petro Physics ------------------------
  if (input.find("PetroPhysics") != input.end())
  {
    auto petro = input["PetroPhysics"];
    if (petro.find("KrModel") == petro.end())
      DebugStop();
    sim_data.mTPetroPhysics.mKrModel = petro["KrModel"];
    if (petro["KrModel"] == 2)
    {
      if (petro.find("Swr") == petro.end())
        DebugStop();
      if (petro.find("Sor") == petro.end())
        DebugStop();
      sim_data.mTPetroPhysics.mSwr = petro["Swr"];
      sim_data.mTPetroPhysics.mSor = petro["Sor"];
      sim_data.mTPetroPhysics.CreateQuadraticResidualKrModel(); // It is necessary to call this method after the residual saturations are set
    }
    sim_data.mTPetroPhysics.mWaterViscosity = sim_data.mTFluidProperties.mWaterViscosity;
    sim_data.mTPetroPhysics.mOilViscosity = sim_data.mTFluidProperties.mOilViscosity;
  }

  // ------------------------ Reservoir Properties ------------------------
  if (input.find("ReservoirProperties") != input.end())
  {
    auto reservoir = input["ReservoirProperties"];
    if (reservoir.find("s0") != reservoir.end())
    {
      auto s0 = reservoir["s0"];
      TMRSPropertiesFunctions::EFunctionType s0_functionType = s0["functionType"];
      REAL constant_val = s0["value"];
      TMRSPropertiesFunctions reservoir_properties;
      reservoir_properties.set_function_type_s0(s0_functionType, constant_val);
      auto s0_function = reservoir_properties.Create_s0();
      sim_data.mTReservoirProperties.s0 = s0_function;
    }
  }

  // ------------------------ Setting extra stuff that is still not in JSON ------------------------
  const int D_Type = 0, N_Type = 1, Mixed_Type = 2;
  // sim_data.mTGeometry.mInterface_material_id = 100;
  // sim_data.mTGeometry.mInterface_material_idFracInf = 102;
  // sim_data.mTGeometry.mInterface_material_idFracSup = 101;
  // sim_data.mTGeometry.mInterface_material_idFracFrac = 103;
  // sim_data.mTGeometry.mInterface_material_idFracBound = 104;

  // sim_data.mTGeometry.mSkeletonDiv = 0;
  sim_data.mTNumerics.m_sfi_tol = 0.00000001;
  sim_data.mTNumerics.m_res_tol_transport = 0.00000001;
  sim_data.mTNumerics.m_corr_tol_transport = 0.00000001;
  sim_data.mTNumerics.m_four_approx_spaces_Q = true;
  sim_data.mTNumerics.m_nThreadsMixedProblem = glob_n_threads;
  sim_data.mTNumerics.m_max_iter_sfi = 20;
  sim_data.mTNumerics.m_max_iter_mixed = 20;
  sim_data.mTNumerics.m_max_iter_transport = 20;

  sim_data.mTPostProcess.m_file_name_mixed = "postdarcy.vtk";
  sim_data.mTPostProcess.m_file_name_transport = "posttransport.vtk";
  TPZStack<std::string, 10> scalnames, vecnames, scalnamesTransport;
  vecnames.Push("Flux");
  scalnames.Push("Pressure");
  scalnames.Push("div_q");
  if (sim_data.mTNumerics.m_four_approx_spaces_Q)
  {
    scalnames.Push("g_average");
    scalnames.Push("p_average");
  }
  scalnamesTransport.Push("Sw");
  scalnamesTransport.Push("So");

  if (sim_data.mTNumerics.m_UseGradRec)
  {
    scalnamesTransport.Push("GradPos");
  }

  sim_data.mTPostProcess.m_vecnamesDarcy = vecnames;
  sim_data.mTPostProcess.m_scalnamesDarcy = scalnames;
  sim_data.mTPostProcess.m_scalnamesTransport = scalnamesTransport;

  int n_steps = sim_data.mTNumerics.m_n_steps;
  sim_data.mTPostProcess.m_file_time_step = sim_data.mTNumerics.m_dt;
  REAL dt = sim_data.mTNumerics.m_dt;
  TPZStack<REAL, 100> reporting_times;
  REAL time = sim_data.mTPostProcess.m_file_time_step;
  // int n_reporting_times = (n_steps) / (time * 100 / dt) + 1;
  int shift_time = 1;
  int n_reporting_times = (n_steps) / (time * shift_time / dt) + 1;
  REAL r_time = 0.0;
  int j = 1;
  for (int i = 1; i <= n_reporting_times; i++)
  {

    r_time = j * dt * (time / dt);
    reporting_times.push_back(r_time);
    j += shift_time;
  }
  sim_data.mTPostProcess.m_vec_reporting_times = reporting_times;
}

void GetNeigsCenterAndSol(TPZGeoEl *gel, TPZStack<std::pair<TPZManVector<REAL,3>, REAL>> &CenterAndSol) {
  int dim = gel->Mesh()->Dimension();
  TPZCompEl *cel = gel->Reference();
  if(!cel) DebugStop();
  TPZCompMesh *cmesh = cel->Mesh();
  TPZBlock &block = cmesh->Block();
  TPZFMatrix<STATE> &sol = cmesh->Solution();
  int firstside = gel->FirstSide(dim-1);
  int lastside = gel->NSides() - 1;
  for (int side = firstside; side < lastside; side++)
  {
    TPZGeoElSide gelside(gel, side);
    TPZStack<TPZCompElSide> celstack;
    gelside.EqualorHigherCompElementList3(celstack, 0, 0);
    int ncel = celstack.size();
    if(ncel == 0) {
      TPZCompElSide celside = gelside.LowerLevelCompElementList2(0);
      if(celside) celstack.Push(celside);
    }
    std::set<int64_t> elindices;
    for(auto it : celstack) {
      TPZCompEl *cel = it.Element();
      if(!cel) DebugStop();
      TPZGeoEl *gel = cel->Reference();
      if(!gel) DebugStop();
      int64_t index = gel->Index();
      if(elindices.find(index) != elindices.end()) continue;
      elindices.insert(index);
      TPZVec<REAL> center(3,0.);
      TPZVec<REAL> centerNeigh(3,0.);
      TPZGeoElSide gelside = it.Reference();
      TPZVec<REAL> masscent(gelside.Dimension(),0.);
      gelside.CenterPoint(masscent);
      gelside.X(masscent, centerNeigh);
      TPZConnect &connect = cel->Connect(0);
      int64_t seqnum = connect.SequenceNumber();
      int64_t pos = block.Position(seqnum);
      int nshape = connect.NShape();
      REAL cellSol = sol(pos+nshape-1,0);
      CenterAndSol.push_back(std::make_pair(centerNeigh, cellSol));
    }
  }
}


void GetNeigsCenterAndSol(TPZGeoEl *gel, TPZVec<REAL> AllCellSols, TPZStack<std::pair<TPZManVector<REAL,3>, REAL>> &CenterAndSol)
{
  int eldim = gel->Dimension();
  int nnodes = gel->NNodes();
  int nsides = gel->NSides();
  int geoindexel = gel->Index();
  TPZVec<int> NeighIds;
  TPZStack<TPZGeoElSide> ElNeighs;

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

  CenterAndSol.resize(nneighs);

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

    CenterAndSol[ineigh] = std::make_pair(centerNeigh, neighSol); // store solution
  }
}

void GradientCell(TPZVec<REAL> &CenterCell, REAL CellSol, TPZStack<std::pair<TPZManVector<REAL,3>, REAL>> &CenterAndSol, TPZVec<REAL> &GradientSol, int dim)
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

REAL GradientLimiter(TPZGeoEl *gel, REAL cellSol, TPZVec<REAL> &gradCell, TPZStack<std::pair<TPZManVector<REAL,3>, REAL>> &CenterAndSol)
{

  int nnodes = gel->NNodes();
  int nsides = gel->NSides();
  int geldim = gel->Dimension();
  double tolerance = 0.0000001;
  int nneighs = CenterAndSol.size();
  TPZVec<REAL> neighsSols(nneighs);

  for (int ineigh = 0; ineigh < nneighs; ineigh++)
  {
    auto pair = CenterAndSol[ineigh];
    REAL NeighSol = pair.second;
    neighsSols[ineigh] = NeighSol;
  }

  TPZFMatrix<REAL> allElNodes(nnodes, 2); // consider mesh 2D
  TPZVec<double> allSlxl(nnodes);

  TPZGeoElSide gelside(gel, nsides); // find the centroid of the element
  TPZVec<REAL> Xcenter(3);
  gelside.CenterX(Xcenter);
  REAL xcenter = Xcenter[0];
  REAL ycenter = Xcenter[1];
  REAL zcenter = Xcenter[2];

  for (int i = 0; i < nnodes; i++)
  { // compute the solution at the vertices
    TPZGeoNode *inode = gel->NodePtr(i);
    REAL xcoord = inode->Coord(0);
    REAL ycoord = inode->Coord(1);
    REAL zcoord = inode->Coord(2);
    if (geldim == 2)
    {
      allSlxl[i] = cellSol + ((gradCell[0] * (xcoord - xcenter)) + (gradCell[1] * (ycoord - ycenter)));
    }
    if (geldim == 3)
    {
      allSlxl[i] = cellSol + ((gradCell[0] * (xcoord - xcenter)) + (gradCell[1] * (ycoord - ycenter)) + (gradCell[2] * (zcoord - zcenter))); // Verify
    }
  }

  // get max and min
  neighsSols.push_back(cellSol); // Add cellSol to allneighsSol
  REAL maxSat = neighsSols[0];
  REAL minSat = neighsSols[0];

  for (int j = 0; j < neighsSols.size(); j++)
  {
    REAL satNeigh = neighsSols[j];
    if (satNeigh > maxSat)
    {
      maxSat = satNeigh;
    }
    if (satNeigh < minSat)
    {
      minSat = satNeigh;
    }
  }

  TPZVec<REAL> allYs;
  REAL yval;

  for (int i = 0; i < nnodes; i++)
  {
    auto slxlVal = allSlxl[i];
    REAL checkval = std::abs(slxlVal - cellSol);

    if (checkval < tolerance)
    {
      yval = 1.0;
    }
    else
    {
      if (slxlVal > cellSol)
      {
        yval = (maxSat - cellSol) / (slxlVal - cellSol);
      }
      if (slxlVal < cellSol)
      {
        yval = (minSat - cellSol) / (slxlVal - cellSol);
      }
    }
    REAL ybar = ((yval * yval) + 2 * yval) / ((yval * yval) + yval + 2);
    allYs.push_back(ybar);
  }

  int nyvals = allYs.size();

  // find limiter val AlphaK
  REAL alphaK = allYs[0];
  for (int j = 0; j < nyvals; j++)
  {
    REAL alphaKVal = allYs[j];
    if (alphaKVal < alphaK)
    {
      alphaK = alphaKVal;
    }
  }
  return alphaK;
}

void RunProblem(TPZGeoMesh *gmesh, const ExactSolFunc &sol, TPZVec<REAL> &AllCelAverage)
{
  int ngels = gmesh->NElements();
  int meshdim = gmesh->Dimension();
  AllCelAverage.resize(ngels);

  for (int igel = 0; igel < ngels; igel++)
  {
    TPZGeoEl *gel = gmesh->Element(igel);
    int eldim = gel->Dimension();
    if (eldim == meshdim)
    {
      int nsides = gel->NSides();
      TPZGeoElSide gelside(gel, nsides);
      TPZVec<REAL> Xcenter(3);
      gelside.CenterX(Xcenter);
      TPZVec<STATE> rhsVal(1);
      TPZFMatrix<STATE> matVal;
      sol(Xcenter, rhsVal, matVal);
      auto cellAverage = rhsVal[0];
      AllCelAverage[igel] = cellAverage;
    }
  }
}

TPZCompMesh *CreateCMesh(TPZGeoMesh *gmesh, const ExactSolFunc &ExactSolution)
{
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  const int dim = gmesh->Dimension();
  cmesh->SetDimModel(dim);
  cmesh->SetDefaultOrder(0);
  cmesh->SetAllCreateFunctionsDiscontinuous();
  int matid = 1;

  TPZL2Projection<STATE> *mat = new TPZL2Projection(matid, dim, 1);

  mat->SetExactSol(ExactSolution, 2);

  cmesh->InsertMaterialObject(mat);

  // Constructs mesh
  cmesh->AutoBuild();

  {
    int64_t ncel = cmesh->NElements();
    for (int64_t i = 0; i < ncel; i++)
    {
      TPZCompEl *cel = cmesh->Element(i);
      if (cel)
      {
        TPZCompElDisc *celDisc = dynamic_cast<TPZCompElDisc *>(cel);
        if (celDisc)
        {
          celDisc->SetFalseUseQsiEta();
          celDisc->SetTotalOrderShape();
          celDisc->SetConstC(1.);
          celDisc->PRefine(1);          
        }
        else
        {
          DebugStop();
        }
      }
    }
  }

  cmesh->ExpandSolution();
  TPZVec<REAL> AllCelAverage;
  RunProblem(gmesh, ExactSolution, AllCelAverage);
  {
    int64_t ncel = cmesh->NElements();
    TPZBlock &block = cmesh->Block();
    TPZFMatrix<REAL> &sol = cmesh->Solution();

    for (int64_t i = 0; i < ncel; i++)
    {
      TPZCompEl *cel = cmesh->Element(i);
      if (cel)
      {
        TPZGeoEl *gel = cel->Reference();
        int gelindex = gel->Index();
        REAL cellSol = AllCelAverage[gelindex];

        int ncon = cel->NConnects();
        if(ncon != 1)
        {
          DebugStop();
        }
        TPZConnect &con = cel->Connect(0);
        int64_t seqnum = con.SequenceNumber();
        int64_t pos = block.Position(seqnum);
        int blsize = block.Size(seqnum);
        if (blsize != 3)
        {
          DebugStop();
        }
        sol(pos+blsize-1,0) = cellSol;
      }
    }
  }
  TPZStack<std::string> fields;
  fields.Push("Solution");
  fields.Push("Derivative");
  TPZVTKGenerator vtk(cmesh, fields, "cmesh.vtk" , 0, gmesh->Dimension());
  vtk.Do();
  return cmesh;
}

STATE CellAverage(TPZCompElDisc *disc) {
  TPZCompEl *cel = disc;
  TPZGeoEl *gel = cel->Reference();
  int gelindex = gel->Index();
  int64_t seqnum = cel->Connect(0).SequenceNumber();
  TPZBlock &block = cel->Mesh()->Block();
  TPZFMatrix<REAL> &sol = cel->Mesh()->Solution();
  int pos = block.Position(seqnum);
  int blsize = block.Size(seqnum);
  REAL cellSol = sol(pos + blsize - 1, 0);
  return cellSol;
}

void ComputeGradients(TPZCompMesh *cmesh) {
  int64_t ncels = cmesh->NElements();
  int mdim = cmesh->Dimension();
  TPZBlock &block = cmesh->Block();
  TPZFMatrix<REAL> &sol = cmesh->Solution();
  for (int64_t i = 0; i < ncels; i++)
  {
    TPZCompEl *cel = cmesh->Element(i);
    if (!cel) {
      continue;
    }
    TPZGeoEl *gel = cel->Reference();
    int geldim = gel->Dimension();
    if (geldim != mdim) {
      continue;
    }
    TPZCompElDisc *celDisc = dynamic_cast<TPZCompElDisc *>(cel);
    if (!celDisc) {
      DebugStop();
    }
    int nsides = gel->NSides();
    REAL CellSol = CellAverage(celDisc);
    TPZGeoElSide gelside(gel, nsides);
    TPZVec<REAL> centerCell(3);
    gelside.CenterX(centerCell);
    TPZStack<std::pair<TPZManVector<REAL,3>, REAL>> CenterAndSol;
    GetNeigsCenterAndSol(gel, CenterAndSol);
    TPZVec<REAL> GradientSol;
    GradientCell(centerCell, CellSol, CenterAndSol, GradientSol, mdim);
    REAL AlphaK = GradientLimiter(gel, CellSol, GradientSol, CenterAndSol);
    TPZConnect &con = cel->Connect(0);
    int64_t seqnum = con.SequenceNumber();
    int blsize = block.Size(seqnum);
    if (blsize != 3) {
      DebugStop();
    }
        int64_t pos = block.Position(seqnum);
    sol(pos,0) = GradientSol[0];
    sol(pos+1,0) = GradientSol[1];

//
//    std::cout << "gelId = " << i << "  AlphaK= " << AlphaK << "   Gradient= " << GradientSol[0] << "  " << GradientSol[1] << std::endl;
  }
  TPZStack<std::string> fields;
  fields.Push("Solution");
  fields.Push("Derivative");
  TPZVTKGenerator vtk(cmesh, fields, "cmeshafter.vtk" , 0, cmesh->Dimension());
  vtk.Do();
}

void ComputeGradients(TPZGeoMesh *gmesh, TPZVec<REAL> &AllCellSols, TPZFMatrix<REAL> &Gradients)
{
  int ngels = gmesh->NElements();
  int mdim = gmesh->Dimension();

  for (int i = 0; i < ngels; i++)
  {
    TPZGeoEl *gel = gmesh->Element(i);
    int geldim = gel->Dimension();

    if (geldim == mdim)
    {
      int nsides = gel->NSides();
      REAL CellSol = AllCellSols[i];
      TPZGeoElSide gelside(gel, nsides);
      TPZVec<REAL> centerCell(3);
      gelside.CenterX(centerCell);
      TPZStack<std::pair<TPZManVector<REAL,3>, REAL>> CenterAndSol;
      GetNeigsCenterAndSol(gel, AllCellSols, CenterAndSol);
      TPZVec<REAL> GradientSol;
      GradientCell(centerCell, CellSol, CenterAndSol, GradientSol, mdim);
      REAL AlphaK = GradientLimiter(gel, CellSol, GradientSol, CenterAndSol);

      std::cout << "gelId = " << i << "  AlphaK= " << AlphaK << "   Gradient= " << GradientSol[0] << "  " << GradientSol[1] << std::endl;
      int stop = 1;
    }
  }
}


void EvaluateError(TPZGeoMesh *gmesh,
                           std::function<void(const TPZVec<REAL> &, TPZVec<STATE> &, TPZFMatrix<STATE> &)> solution,TPZVec<REAL>& error)
{
    TPZCompMesh *cmesh =CreateCMesh(gmesh, solution);
   
    TPZManVector<REAL, 3> errorsum(3, 0.);
    cmesh->ElementSolution().Redim(cmesh->NElements(),3);
    cmesh->EvaluateError(true, errorsum);
   
    REAL errorL2 = errorsum[1];
    
    ComputeGradients(cmesh);
    errorsum.Fill(0.);
    cmesh->EvaluateError(true, errorsum);
    REAL errorReconstL2 =errorsum[1];
  
    error.resize(2);
    error[0]=errorL2;
    error[1]=errorReconstL2;
    cout<<"error L2= "<<errorL2<<"  error reconstru L2= "<<errorReconstL2<<endl;
 
    
}

  // ---------------------------------------------------------------------
  // ---------------------------------------------------------------------
