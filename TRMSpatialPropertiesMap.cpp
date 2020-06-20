//
//  TRMSpatialPropertiesMap.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMSpatialPropertiesMap.h"
#include <boost/math/distributions/normal.hpp>

/** @brief default constructor */
TRMSpatialPropertiesMap::TRMSpatialPropertiesMap(){
 
    fSPE10Cmesh = NULL;
    
    /** @brief SPE10 fields file */
    fPermPorFields.first  = "";
    fPermPorFields.second = "";
    
    /** @brief number of blocks i, j and k  */
    fNBlocks.resize(0);
    
    /** @brief size of blocks dx, dy and dz  */
    fBlocks_sizes.resize(0);
}

/** @brief default destructor */
TRMSpatialPropertiesMap::~TRMSpatialPropertiesMap(){
    
}

/** @brief Absolute Permeability m2  $\kappa$ */
void TRMSpatialPropertiesMap::Kappa(TPZManVector<STATE,3> &x, TPZFMatrix<STATE> &kappa){

    REAL phi;
//    this->ComputePropertieSPE10Map(x, kappa, phi);
}

/** @brief Porosity fraction  $\phi$ */
void TRMSpatialPropertiesMap::phi(TPZManVector<STATE,3> &x, TPZManVector<STATE,10> &phi){
    
    TPZFMatrix<STATE> kappa, inv_kappa;
    long index = 0;
    phi.Resize(10, 0.0);
    REAL val;
//    this->ComputePropertieSPE10Map(index, x, kappa, inv_kappa, val);
    phi[0] = val;
}

void TRMSpatialPropertiesMap::LoadSPE10Map(bool PrintMapQ)
{

    // Cartesian mesh for SPE10 spatial properties
    // http://www.spe.org/web/csp/datasets/set02.htm
    
    int nel_x = fNBlocks[0];
    int nel_y = fNBlocks[1];
    int nel_z = fNBlocks[2];
    TPZManVector<REAL,2> dx(2,nel_x), dy(2,nel_y), dz(2,nel_z);
    dx[0] = fBlocks_sizes[0];
    dy[0] = fBlocks_sizes[1];
    dz[0] = fBlocks_sizes[2];
    TPZGeoMesh * gmesh =this->CreateGeometricBoxMesh(dx, dy, dz);
    
    REAL angle = 90.0;
    int axis = 3; // z -axis;
    this->RotateGeomesh(gmesh, angle, axis);
    
    REAL s = 1.0;
    this->ExpandGeomesh(gmesh, s, s, s);
    
    TPZVec<REAL> t_vec(3,0.0);
    t_vec[0] = 0.0;
    t_vec[1] = 0.0;
    t_vec[2] = 0.0;
    this->TraslateGeomesh(gmesh, t_vec);
    
    if(!gmesh)
    {
        std::cout<< "SpatialPropertiesMap:: Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    
    int dim = gmesh->Dimension();
    int order = 0;
    int rock_id = 1;
    
    // Malha computacional
    fSPE10Cmesh = new TPZCompMesh(gmesh);

    TRMSpatialMap * mat = new TRMSpatialMap(rock_id,dim);
    fSPE10Cmesh->InsertMaterialObject(mat);
    fSPE10Cmesh->SetDimModel(dim);
    fSPE10Cmesh->SetDefaultOrder(order);
    fSPE10Cmesh->SetAllCreateFunctionsDiscontinuous();
    fSPE10Cmesh->AutoBuild();
    
    
//#ifdef PZDEBUG
//    std::ofstream out("Cmesh_SPE10.txt");
//    fSPE10Cmesh->Print(out);
//
//    //  Print Geometrical Base Mesh
//    std::ofstream planefile("SPE10Geometry.txt");
//    gmesh->Print(planefile);
//    std::ofstream file("SPE10Geometry.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh,file, true);
//#endif
    

    int n_data = nel_x * nel_y * nel_z;
    this->Insert_Inside_Map(n_data);
    
    if (PrintMapQ) {
        TPZAnalysis * an_map = new TPZAnalysis(fSPE10Cmesh,false);
        an_map->LoadSolution();
        
        int div = 0;
        TPZStack<std::string> scalnames, vecnames;
        std::string plotfile =  "Spatial_map.vtk";
        scalnames.Push("phi");
        vecnames.Push("kappa");
        
        an_map->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
        an_map->PostProcess(div);
    }
    
}

/** @brief Insert spatial properties from SPE10 on pz mesh of order zero */
bool TRMSpatialPropertiesMap::Insert_Inside_Map(int n_data){

    bool IsLoadedQ = false;
    
#ifdef PZDEBUG
    
    if(!fSPE10Cmesh || n_data > 1122000 )
    {
        std::cout<< "SpatialPropertiesMap:: Computational interpolation mesh doesn't exist" << std::endl;
        DebugStop();
    }
    
#endif
    
    bool IsConstantQ = false;
    bool UseKozenyCarmanQ = false;
    REAL mdTom2 = 9.869233e-16;
    
    TPZFMatrix<REAL> properties(n_data,4);
    properties.Zero();
    std::string dirname = PZSOURCEDIR;
    std::string perm_file, phi_file;
    perm_file = dirname + "/Projects/iRMS/" + fPermPorFields.first;
    phi_file  = dirname + "/Projects/iRMS/" + fPermPorFields.second;
    
    //////////////////////////////////////////////// perm data /////////////////////////////////////
    
    REAL kab = 1.0e-13;
    REAL phi = 0.25;
    
    if (IsConstantQ) {
        for (int i = 0; i < n_data; i++) {
            properties(i,0) = kab;
            properties(i,1) = kab;
            properties(i,2) = kab;
            properties(i,3) = phi;
        }
    }
    else{
        
        // reading a cartesian porosity file
        std::ifstream read_por (phi_file.c_str());
        REAL phi_val;
        int count = 0;
        for (int k = 0; k < fNBlocks[2]; k++) {
            for (int j = 0; j < fNBlocks[1]; j++) {
                for (int i = 0; i < fNBlocks[0]; i++) {
                    read_por >> phi_val;
                    properties(count,3) = phi_val + 0.001;
                    if (UseKozenyCarmanQ) {
                        double k_val, c = 0.0001;
                        k_val = phi_val*phi_val*phi_val/(c*(1.0-phi_val)*(1.0-phi_val));
                        properties(count,0) = k_val*mdTom2;
                        properties(count,1) = k_val*mdTom2;
                        properties(count,2) = k_val*mdTom2;
                    }
                    count++;
                }
            }
        }
        
        read_por.close();
        
        if (!UseKozenyCarmanQ) {
            // reading a cartesian permeability file
            std::ifstream read_perm(perm_file.c_str());
            
            REAL k_val;
            int count = 0;
            for (int k = 0; k < fNBlocks[2]; k++) {
                for (int j = 0; j < fNBlocks[1]; j++) {
                    for (int i = 0; i < fNBlocks[0]; i++) {
                        read_perm >> k_val;
                        properties(count,0) = k_val*mdTom2 + kab;
                        count++;
                    }
                }
            }
            
            count = 0;
            for (int k = 0; k < fNBlocks[2]; k++) {
                for (int j = 0; j < fNBlocks[1]; j++) {
                    for (int i = 0; i < fNBlocks[0]; i++) {
                        read_perm >> k_val;
                        properties(count,1) = k_val*mdTom2 + kab;
                        count++;
                    }
                }
            }
            
            count = 0;
            for (int k = 0; k < fNBlocks[2]; k++) {
                for (int j = 0; j < fNBlocks[1]; j++) {
                    for (int i = 0; i < fNBlocks[0]; i++) {
                        read_perm >> k_val;
                        properties(count,2) = k_val*mdTom2 + kab;
                        count++;
                    }
                }
            }
            
            read_perm.close();
        }

    }
    

    
    
#ifdef PZDEBUG
    if(fSPE10Cmesh->NElements() != n_data) {
        DebugStop();
    }
#endif
    
    int n_elements = fSPE10Cmesh->NElements();
    TPZVec<REAL> x_c(3,0.0);
    TPZVec<REAL> qsi(3,0.0);
    TPZManVector<long,4> dof_indexes;
    for (int icel = 0; icel < n_elements; icel++) {
        
        TPZCompEl * cel = fSPE10Cmesh->Element(icel);
        
#ifdef PZDEBUG
        if(!cel) {
            DebugStop();
        }
#endif
        
        TPZInterpolationSpace * int_vol = dynamic_cast<TPZInterpolationSpace *>(cel);
        this->ElementDofIndexes(int_vol, dof_indexes);
        
        for (int k = 0; k < dof_indexes.size(); k++) {
            fSPE10Cmesh->Solution()(dof_indexes[k],0) = properties(icel,k);
        }

    }
    
    IsLoadedQ = true;
    return IsLoadedQ;
}

void TRMSpatialPropertiesMap::ElementDofIndexes(TPZInterpolationSpace * &intel, TPZVec<long> &dof_indexes){
    
#ifdef PZDEBUG
    if (!intel) {
        DebugStop();
    }
#endif
    
    int n_state = 4;
    
    TPZStack<long> index(0,0);
    int nconnect = intel->NConnects();
    for (int icon = 0; icon < nconnect; icon++) {
        TPZConnect  & con = intel->Connect(icon);
        long seqnumber = con.SequenceNumber();
        long position = intel->Mesh()->Block().Position(seqnumber);
        int nshape = con.NShape() * n_state;
        for (int ish=0; ish < nshape; ish++) {
            index.Push(position+ ish);
        }
    }
    
    dof_indexes = index;
    return;
}

bool TRMSpatialPropertiesMap::ComputePropertieSPE10Map(TPZVec<STATE> &x, TPZFMatrix<STATE> &kappa, REAL & phi){
    
    int index = 0;
#ifdef PZDEBUG
    
    if(!fSPE10Cmesh)
    {
        std::cout<< "SpatialPropertiesMap:: Computational interpolation mesh doesn't exist" << std::endl;
        DebugStop();
    }
    
#endif
    fSPE10Cmesh->LoadReferences();
    TPZGeoMesh * geometry = fSPE10Cmesh->Reference();
    
#ifdef PZDEBUG
    if(!geometry)
    {
        DebugStop();
    }
#endif
    TPZVec<REAL> qsi;
    int target_dim = 3;
//    TPZGeoEl * gel = geometry->FindApproxElement(x, qsi, index, target_dim);
    TPZGeoEl * gel = geometry->FindElement(x, qsi, index, target_dim);
    
    if(!gel)
    {
        std::cout<< "SpatialPropertiesMap:: Warning Filling map with constant values " << std::endl;
        TPZManVector<STATE,3> x_c(x);
        TPZManVector<STATE, 10> state_vars(4,0.0);
        this->Kappa_c(x_c, kappa, inv_kappa, state_vars);
        TPZManVector<STATE, 10> phi_v;
        this->phi_c(x_c, phi_v, state_vars);
        phi = phi_v[0];
        return false;
    }
    
    kappa.Resize(3,3);
    kappa.Zero();
    inv_kappa.Resize(3,3);
    inv_kappa.Zero();
    
    TPZCompEl * cel = gel->Reference();
    
#ifdef PZDEBUG
    if(!cel)
    {
        DebugStop();
    }
#endif
    
    int k_index = 0;
    TPZManVector<STATE,3> sol;
    cel->Solution(qsi, k_index, sol);
    
    REAL s = 1.0e+6;
    kappa(0,0) = sol[0]*s;
    kappa(1,1) = sol[1]*s;
    kappa(2,2) = sol[2]*s;
    
    int phi_index = 1;
    cel->Solution(qsi, phi_index, sol);
    phi = sol[0];
    
    return true;
    
}

/** @brief Create a reservoir-box geometry i-j-k element ordering  */
TPZGeoMesh *  TRMSpatialPropertiesMap::CreateGeometricBoxMesh(TPZManVector<REAL,2> dx, TPZManVector<REAL,2> dy, TPZManVector<REAL,2> dz){
    
    REAL t=0.0;
    REAL dt;
    int n;
    bool IsTetrahedronMeshQ = false;
    
    int rock =  1;
    int bcs =  -1;
    
    // Creating a 0D element to be extruded
    TPZGeoMesh * GeoMesh0D = new TPZGeoMesh;
    GeoMesh0D->NodeVec().Resize(1);
    TPZGeoNode Node;
    TPZVec<REAL> coors(3,0.0);
    Node.SetCoord(coors);
    Node.SetNodeId(0);
    GeoMesh0D->NodeVec()[0]=Node;
    
    TPZVec<long> Topology(1,0);
    int elid=0;
    
    new TPZGeoElRefPattern < pzgeom::TPZGeoPoint >(elid,Topology,rock,*GeoMesh0D);
    GeoMesh0D->BuildConnectivity();
    GeoMesh0D->SetDimension(0);
    
    TPZHierarquicalGrid CreateGridFrom0D(GeoMesh0D);
    TPZAutoPointer<TPZFunction<STATE> > ParFuncX = new TPZDummyFunction<STATE>(ParametricfunctionZ);
    CreateGridFrom0D.SetParametricFunction(ParFuncX);
    CreateGridFrom0D.SetFrontBackMatId(bcs,bcs);
    
    dt = dz[0];
    t  = -dz[1]*(dz[0]/2.0);
    n = int(dz[1]);
    // Computing Mesh extruded along the parametric curve Parametricfunction
    TPZGeoMesh * GeoMesh1D = CreateGridFrom0D.ComputeExtrusion(t, dt, n);
    
    TPZHierarquicalGrid CreateGridFrom1D(GeoMesh1D);
    TPZAutoPointer<TPZFunction<STATE> > ParFuncY = new TPZDummyFunction<STATE>(ParametricfunctionY);
    CreateGridFrom1D.SetParametricFunction(ParFuncY);
    CreateGridFrom1D.SetFrontBackMatId(bcs,bcs);
    if(IsTetrahedronMeshQ){
        CreateGridFrom1D.SetTriangleExtrusion();
    }
    
    
    dt = dy[0];
    t  = -dy[1]*(dy[0]/2.0);
    n = int(dy[1]);
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh2D = CreateGridFrom1D.ComputeExtrusion(t, dt, n);
    
    TPZHierarquicalGrid CreateGridFrom2D(GeoMesh2D);
    TPZAutoPointer<TPZFunction<STATE> > ParFuncZ = new TPZDummyFunction<STATE>(ParametricfunctionX);
    CreateGridFrom2D.SetParametricFunction(ParFuncZ);
    CreateGridFrom2D.SetFrontBackMatId(bcs,bcs);
    if(IsTetrahedronMeshQ){
        CreateGridFrom2D.SetTriangleExtrusion();
        CreateGridFrom2D.SetTetrahedonExtrusion();
    }
    
    
    dt = dx[0];
    t  = -dx[1]*(dx[0]/2.0);
    n = int(dx[1]);
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh3D = CreateGridFrom2D.ComputeExtrusion(t, dt, n);
    
    long last_node = GeoMesh3D->NNodes() - 1;
    long last_element = GeoMesh3D->NElements() - 1;
    long node_id = GeoMesh3D->NodeVec()[last_node].Id();
    long element_id = GeoMesh3D->Element(last_element)->Id();
    const std::string name("SPE10 reservoir box ");
    GeoMesh3D->SetName(name);
    GeoMesh3D->SetMaxNodeId(node_id);
    GeoMesh3D->SetMaxElementId(element_id);
    return GeoMesh3D;
}


void TRMSpatialPropertiesMap::ParametricfunctionX(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = par[0];
    X[1] = 0.0*sin(0.01*par[0]);
    X[2] = 0.0;
}

void TRMSpatialPropertiesMap::ParametricfunctionY(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = 0.0;
    X[1] = par[0];
    X[2] = 0.0;
}

void TRMSpatialPropertiesMap::ParametricfunctionZ(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = 0.0;
    X[1] = 0.0;
    X[2] = par[0];
}

void TRMSpatialPropertiesMap::ExpandGeomesh(TPZGeoMesh *gmesh, REAL sx, REAL  sy, REAL  sz)
{
    
    TPZVec<STATE> iCoords(3,0.0);
    TPZVec<STATE> iCoordsExpanded(3,0.0);
    
    int NumberofGeoNodes = gmesh->NNodes();
    for (int inode = 0; inode < NumberofGeoNodes; inode++)
    {
        TPZGeoNode GeoNode = gmesh->NodeVec()[inode];
        GeoNode.GetCoordinates(iCoords);
        // Apply rotation
        iCoordsExpanded[0] = iCoords[0]*sx;
        iCoordsExpanded[1] = iCoords[1]*sy;
        iCoordsExpanded[2] = iCoords[2]*sz;
        GeoNode.SetCoord(iCoordsExpanded);
        gmesh->NodeVec()[inode] = GeoNode;
    }
    
}

void TRMSpatialPropertiesMap::TraslateGeomesh(TPZGeoMesh *gmesh, TPZVec<REAL> t_vec)
{
    
    TPZVec<STATE> iCoords(3,0.0);
    TPZVec<STATE> iCoordsTraslated(3,0.0);
    
    int NumberofGeoNodes = gmesh->NNodes();
    for (int inode = 0; inode < NumberofGeoNodes; inode++)
    {
        TPZGeoNode GeoNode = gmesh->NodeVec()[inode];
        GeoNode.GetCoordinates(iCoords);
        // Apply rotation
        iCoordsTraslated[0] = iCoords[0] + t_vec[0];
        iCoordsTraslated[1] = iCoords[1] + t_vec[1];
        iCoordsTraslated[2] = iCoords[2] + t_vec[2];
        GeoNode.SetCoord(iCoordsTraslated);
        gmesh->NodeVec()[inode] = GeoNode;
    }
    
}

void TRMSpatialPropertiesMap::RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis)
{
    REAL theta =  (M_PI/180.0)*CounterClockwiseAngle;
    // It represents a 3D rotation around the axis -> i.
    TPZFMatrix<STATE> RotationMatrix(3,3,0.0);
    
    switch (Axis) {
        case 1:
        {
            RotationMatrix(0,0) = 1.0;
            RotationMatrix(1,1) =   +cos(theta);
            RotationMatrix(1,2) =   -sin(theta);
            RotationMatrix(2,1) =   +sin(theta);
            RotationMatrix(2,2) =   +cos(theta);
        }
            break;
        case 2:
        {
            RotationMatrix(0,0) =   +cos(theta);
            RotationMatrix(0,2) =   +sin(theta);
            RotationMatrix(1,1) = 1.0;
            RotationMatrix(2,0) =   -sin(theta);
            RotationMatrix(2,2) =   +cos(theta);
        }
            break;
        case 3:
        {
            RotationMatrix(0,0) =   +cos(theta);
            RotationMatrix(0,1) =   -sin(theta);
            RotationMatrix(1,0) =   +sin(theta);
            RotationMatrix(1,1) =   +cos(theta);
            RotationMatrix(2,2) = 1.0;
        }
            break;
        default:
        {
            RotationMatrix(0,0) =   +cos(theta);
            RotationMatrix(0,1) =   -sin(theta);
            RotationMatrix(1,0) =   +sin(theta);
            RotationMatrix(1,1) =   +cos(theta);
            RotationMatrix(2,2) = 1.0;
        }
            break;
    }
    
    TPZVec<STATE> iCoords(3,0.0);
    TPZVec<STATE> iCoordsRotated(3,0.0);
    
    //RotationMatrix.Print("Rotation = ");
    
    int NumberofGeoNodes = gmesh->NNodes();
    for (int inode = 0; inode < NumberofGeoNodes; inode++)
    {
        TPZGeoNode GeoNode = gmesh->NodeVec()[inode];
        GeoNode.GetCoordinates(iCoords);
        // Apply rotation
        iCoordsRotated[0] = RotationMatrix(0,0)*iCoords[0]+RotationMatrix(0,1)*iCoords[1]+RotationMatrix(0,2)*iCoords[2];
        iCoordsRotated[1] = RotationMatrix(1,0)*iCoords[0]+RotationMatrix(1,1)*iCoords[1]+RotationMatrix(1,2)*iCoords[2];
        iCoordsRotated[2] = RotationMatrix(2,0)*iCoords[0]+RotationMatrix(2,1)*iCoords[1]+RotationMatrix(2,2)*iCoords[2];
        GeoNode.SetCoord(iCoordsRotated);
        gmesh->NodeVec()[inode] = GeoNode;
    }
    
}
