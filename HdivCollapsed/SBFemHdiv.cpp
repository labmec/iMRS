#ifdef HAVE_CONFIG_H
#include "pz_config.h"
#endif

#include "SBFemHdiv.h"
#include "TPZLagrangeMultiplier.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include <stack> 

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.sbfem"));
#endif


#ifdef _AUTODIFF
TLaplaceExample1 LaplaceExact;
#endif

bool gHybrid, gSbfemhdiv, gHdivcollapsed;

TPZGeoMesh * SetupSquareMesh(int nelx)
{
    TPZManVector<REAL, 4> x0(3, -1.), x1(3, 1.);
    x0[0] = -1, x0[1] = -5, x0[2] = 0.;
    x1[0] = 2, x1[1] = 1, x1[2] = 0.;

    TPZManVector<int, 4> nx(2, nelx);
    TPZGenGrid2D gengrid(nx, x0, x1);
    gengrid.SetElementType(MMeshType::EQuadrilateral);
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gengrid.Read(gmesh, Emat1);

    // Setting boundary conditions
    gengrid.SetBC(gmesh, 4, Ebc1);
    gengrid.SetBC(gmesh, 5, Ebc2);
    gengrid.SetBC(gmesh, 6, Ebc3);
    gengrid.SetBC(gmesh, 7, Ebc4);

    gmesh->BuildConnectivity();

    if (gHdivcollapsed)
    {
        TPZGeoMesh *gmeshcollapsed = new TPZGeoMesh();
        gmeshcollapsed->SetDimension(1);
        auto nnodes = gmesh->NNodes();
        for (auto node : gmesh->NodeVec())
        {
            auto newnode0 = gmeshcollapsed->NodeVec().AllocateNewElement();
            TPZManVector<REAL,3> xco(3,0.);
            node.GetCoordinates(xco);
            gmeshcollapsed->NodeVec()[newnode0].Initialize(xco, *gmeshcollapsed);
        }
        
        for (auto gel : gmesh->ElementVec())
        {
            if (!gel) continue;
            if (gel->MaterialId() != Emat1) continue;

            for (auto isides = 0; isides < gel->NSides(); isides++)
            {
                TPZGeoElSide gelside = gel->Neighbour(isides);
                if (!gelside || gelside.Dimension() != 1) continue;
                if (gelside.Element()->Index() < gel->Index()) continue; // if the neighbour is a previous element, continue (the interface element was already created)
                
                TPZManVector<int64_t> nodes(2,0);
                nodes[0] = gelside.SideNodeIndex(0);
                nodes[1] = gelside.SideNodeIndex(1);
                int64_t index;
                gmeshcollapsed->CreateGeoElement(EOned, nodes, ESkeleton, index);
            }
        }
        
        gmeshcollapsed->BuildConnectivity();
        
        for (auto gel : gmeshcollapsed->ElementVec())
        {
            if(!gel) continue;
            if (gel->MaterialId() != ESkeleton) continue;
            
            TPZManVector<REAL,3> xco1(3), xco2(3);
            gel->Node(0).GetCoordinates(xco1);
            gel->Node(1).GetCoordinates(xco2);
            if (IsZero(xco1[1] + 1) && IsZero(xco2[1] + 1))
            {
                TPZGeoElBC(gel,2,Ebc1);
            }
            else if(IsZero(xco1[0] - 1) && IsZero(xco2[0] - 1))
            {
                TPZGeoElBC(gel,2,Ebc2);
            }
            else if(IsZero(xco1[1] - 1) && IsZero(xco2[1] - 1))
            {
                TPZGeoElBC(gel,2,Ebc3);
            }
            else if(IsZero(xco1[0] + 1) && IsZero(xco2[0] + 1))
            {
                TPZGeoElBC(gel,2,Ebc4);
            }
        }
        std::ofstream out("gmesh.txt");
        gmeshcollapsed->Print(out);
        return gmeshcollapsed;
    }
    
    return gmesh;
}

TPZGeoMesh * SetupCollapsedMesh(int nelx)
{
    TPZManVector<REAL, 4> x0(3, -1.), x1(3, 1.);
    x0[0] = -1, x0[1] = -1, x0[2] = 0.;
    x1[0] = 1, x1[1] = 1, x1[2] = 0.;

    TPZManVector<int, 4> nx(2, nelx);
    TPZGenGrid2D gengrid(nx, x0, x1);
    gengrid.SetElementType(MMeshType::EQuadrilateral);
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    gengrid.Read(gmesh, Emat1);

    // Setting boundary conditions
    gengrid.SetBC(gmesh, 4, Ebc1);
    gengrid.SetBC(gmesh, 5, Ebc2);
    gengrid.SetBC(gmesh, 6, Ebc3);
    gengrid.SetBC(gmesh, 7, Ebc4);

    gmesh->BuildConnectivity();

    TPZGeoMesh *gmeshcollapsed = new TPZGeoMesh();
    gmeshcollapsed->SetDimension(2);
    auto nnodes = gmesh->NNodes();
    for (auto node : gmesh->NodeVec())
    {
        auto newnode0 = gmeshcollapsed->NodeVec().AllocateNewElement();
        TPZManVector<REAL,3> xco(3,0.);
        node.GetCoordinates(xco);
        gmeshcollapsed->NodeVec()[newnode0].Initialize(xco, *gmeshcollapsed);
    }
    
    for (auto gel : gmesh->ElementVec())
    {
        if (!gel) continue;
        if (gel->MaterialId() != Emat1) continue;

        // creating node in the element's center
        TPZManVector<REAL,3> xco(3,0.);
        auto nnodes = gel->NNodes();
        TPZManVector<int64_t,4> nodeindices(nnodes,0);
        gel->GetNodeIndices(nodeindices);

        for (auto inode : nodeindices)
        {
            auto node = gel->Node(inode);
            TPZManVector<REAL,3> xnode(3,0.);
            node.GetCoordinates(xnode);
            xco[0] += xnode[0]; xco[1] += xnode[1]; xco[2] += xnode[2];
        }
        xco[0] = xco[0]/nnodes; xco[1] = xco[1]/nnodes; xco[2] = xco[2]/nnodes;

        auto newnode0 = gmeshcollapsed->NodeVec().AllocateNewElement();
        gmeshcollapsed->NodeVec()[newnode0].Initialize(xco, *gmeshcollapsed);

        for (auto isides = 0; isides < gel->NSides(); isides++)
        {
            TPZGeoElSide gelside = gel->Neighbour(isides);
            if (!gelside || gelside.Dimension() != 1) continue;
            if (gelside.Element()->Index() < gel->Index()) continue; // if the neighbour is a previous element, continue (the interface element was already created)
            
            TPZManVector<int64_t> nodes(2,0);
            nodes[0] = gelside.SideNodeIndex(0);
            nodes[1] = gelside.SideNodeIndex(1);
            int64_t index;
            gmeshcollapsed->CreateGeoElement(EOned, nodes, ESkeleton, index);
            nodes.Resize(4);
            nodes[2] = newnode0; nodes[3] = newnode0;
            gmeshcollapsed->CreateGeoElement(EQuadrilateral, nodes, Emat0, index);
        }
    }
    
    gmeshcollapsed->BuildConnectivity();
    
    // Boundary conditions
    if (!gHybrid)
    {
        for (auto gel : gmeshcollapsed->ElementVec())
        {
            if(!gel) continue;
            if (gel->MaterialId() != ESkeleton) continue;
            TPZGeoElBC(gel,2,Ebc1);
            
            TPZManVector<REAL,3> xco1(3), xco2(3);
            gel->Node(0).GetCoordinates(xco1);
            gel->Node(1).GetCoordinates(xco2);
            if (IsZero(xco1[1] + 1) && IsZero(xco2[1] + 1))
            {
                TPZGeoElBC(gel,2,Ebc1);
            }
            else if(IsZero(xco1[0] - 1) && IsZero(xco2[0] - 1))
            {
                TPZGeoElBC(gel,2,Ebc2);
            }
            else if(IsZero(xco1[1] - 1) && IsZero(xco2[1] - 1))
            {
                TPZGeoElBC(gel,2,Ebc3);
            }
            else if(IsZero(xco1[0] + 1) && IsZero(xco2[0] + 1))
            {
                TPZGeoElBC(gel,2,Ebc4);
            }
        }
    }
    
    gmeshcollapsed->BuildConnectivity();

    if (gHybrid)
    {
        for (auto gel : gmeshcollapsed->ElementVec())
        {
            if (!gel || gel->MaterialId() != ESkeleton)
            {
                continue;
            }
            TPZGeoElBC(gel,2,Eleftpressure);
            TPZGeoElBC(gel,2,Erightpressure);
            TPZGeoElBC(gel,2,Eleftflux);
            TPZGeoElBC(gel,2,Erightflux);
        }
    }
    
    std::ofstream out("gmesh.txt");
    gmeshcollapsed->Print(out);
    return gmeshcollapsed;
}

TPZCompMesh * flux(TPZAutoPointer<TPZGeoMesh> gmesh, int POrder)
{
    auto dim = gmesh->Dimension();
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(POrder);
    
    TPZVecL2 *mat = new TPZVecL2(Emat1);
    mat->SetDimension(dim);
    mat->SetForcingFunction(LaplaceExact.ForcingFunction());
    mat->SetForcingFunctionExact(LaplaceExact.Exact());
    
    cmesh->InsertMaterialObject(mat); //Insere material na malha
    cmesh->SetDimModel(dim);

    int nstate = 1;
    TPZFMatrix<STATE> val1(2, 2, 0.), val2(2, 1, 0.);
    {
        auto bcond = mat->CreateBC(mat, Ebc1, 0, val1, val2);
        cmesh->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, Ebc2, 0, val1, val2);
        cmesh->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, Ebc3, 0, val1, val2);
        cmesh->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, Ebc4, 0, val1, val2);
        cmesh->InsertMaterialObject(bcond);
    }
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
}

TPZCompMesh * fluxhdivcollapsed(TPZAutoPointer<TPZGeoMesh> gmesh, int POrder)
{
    auto dim = 2; auto nstate = 1;
    TPZCompMesh * cmeshcollapsed = new TPZCompMesh(gmesh);
    cmeshcollapsed->SetDefaultOrder(POrder);
    cmeshcollapsed->CleanUp();
    cmeshcollapsed->SetDimModel(dim);

    auto mat = new TPZNullMaterial(ESkeleton, dim, nstate);
    cmeshcollapsed->InsertMaterialObject(mat); //Insere material na malha
    for (auto gel : gmesh->ElementVec())
    {
        if (!gel) continue;
        if (gel->Dimension() != 1) continue;
        if (gel->MaterialId() != ESkeleton) continue;
        int64_t index;
        auto celhdivc = new TPZCompElHDivCollapsed<pzshape::TPZShapeLinear>(*cmeshcollapsed, gel, index);
    }

    cmeshcollapsed->SetDimModel(2);
    cmeshcollapsed->SetAllCreateFunctionsHDiv();
    
    TPZFMatrix<STATE> val1(2, 2, 0.), val2(2, 1, 0.);
    {
        auto bcond = mat->CreateBC(mat, Ebc1, 0, val1, val2);
        cmeshcollapsed->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, Ebc2, 0, val1, val2);
        cmeshcollapsed->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, Ebc3, 0, val1, val2);
        cmeshcollapsed->InsertMaterialObject(bcond);
    }
    {
        auto bcond = mat->CreateBC(mat, Ebc4, 0, val1, val2);
        cmeshcollapsed->InsertMaterialObject(bcond);
    }
    std::set<int> matids = {Ebc1, Ebc2, Ebc3, Ebc4};
    cmeshcollapsed->Reference()->ResetReference();
    cmeshcollapsed->AutoBuild(matids);
    cmeshcollapsed->AdjustBoundaryElements();
    cmeshcollapsed->CleanUpUnconnectedNodes();
    
    cmeshcollapsed->ExpandSolution();
    return cmeshcollapsed;
}

using namespace std;
TPZCompMesh * fluxhdivsbfem(TPZAutoPointer<TPZGeoMesh> gmesh, int POrder)
{
    auto dim = 1; auto nstate = 1;
    auto cmeshcollapsed = new TPZCompMesh(gmesh);
    cmeshcollapsed->SetDefaultOrder(POrder);
    cmeshcollapsed->SetDimModel(dim);
    cmeshcollapsed->CleanUp();

    auto mat = new TPZNullMaterial(ESkeleton, dim, nstate);
    if (gHybrid)
    {
        auto matleft = new TPZNullMaterial(Eleftflux, dim, nstate);
        cmeshcollapsed->InsertMaterialObject(matleft);

        auto matright = new TPZNullMaterial(Erightflux, dim, nstate);
        cmeshcollapsed->InsertMaterialObject(matright);
    }
    
    cmeshcollapsed->InsertMaterialObject(mat); //Insere material na malha

    map<int64_t,TPZCompEl *> geltocel;
    for (auto gel : gmesh->ElementVec())
    {
        if (!gel) continue;
        if (gel->MaterialId() != Emat0) continue;
        int64_t index;
        TPZGeoElSide gelside(gel,4);
        auto gel1dside = gelside.HasNeighbour(ESkeleton);
        if (!gelside)
        {
            DebugStop();
        }
        auto gel1d = gel1dside.Element();
        auto celhdivc = new TPZCompElHDivSBFem<pzshape::TPZShapeLinear>(*cmeshcollapsed, gel1d, gelside, index);
        geltocel[gel1d->Index()] = celhdivc;
    }
    cmeshcollapsed->SetDimModel(dim);
    cmeshcollapsed->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmeshcollapsed->ExpandSolution();
    gmesh->ResetReference();

    if (gHybrid)
    {
        for (auto gel : gmesh->ElementVec())
        {
            if (!gel) continue;
            if (gel->MaterialId() != Eleftflux) continue;

            TPZGeoElSide gelside(gel,2);
            auto intfluxside = gelside.HasNeighbour(ESkeleton);
            auto cel  = geltocel[intfluxside.Element()->Index()];
            TPZCompElHDivSBFem<pzshape::TPZShapeLinear> * celhdivc = dynamic_cast<TPZCompElHDivSBFem<pzshape::TPZShapeLinear> * >(cel);
            if (!celhdivc)
            {
                DebugStop();
            }

            int64_t index;
            auto hdivboundleft = new TPZCompElHDivBound2<pzshape::TPZShapeLinear>(*cmeshcollapsed,gel,index);
            hdivboundleft->SetConnectIndex(0,celhdivc->ConnectIndex(3));
            gel->ResetReference();

            auto rightfluxside = gelside.HasNeighbour(Erightflux);
            auto gel1dright = rightfluxside.Element();
            auto hdivboundright = new TPZCompElHDivBound2<pzshape::TPZShapeLinear>(*cmeshcollapsed,gel1dright,index);
            hdivboundright->SetConnectIndex(0,celhdivc->ConnectIndex(4));
            gel1dright->ResetReference();
        }
    }

    
    if (!gHybrid)
    {
        int nstate = 1;
        TPZFMatrix<STATE> val1(2, 2, 0.), val2(2, 1, 0.);
        {
            auto bcond = mat->CreateBC(mat, Ebc1, 0, val1, val2);
            cmeshcollapsed->InsertMaterialObject(bcond);
        }
        {
            auto bcond = mat->CreateBC(mat, Ebc2, 0, val1, val2);
            cmeshcollapsed->InsertMaterialObject(bcond);
        }
        {
            auto bcond = mat->CreateBC(mat, Ebc3, 0, val1, val2);
            cmeshcollapsed->InsertMaterialObject(bcond);
        }
        {
            auto bcond = mat->CreateBC(mat, Ebc4, 0, val1, val2);
            cmeshcollapsed->InsertMaterialObject(bcond);
        }

        std::set<int> matids = {Ebc1, Ebc2, Ebc3, Ebc4};
        cmeshcollapsed->AutoBuild(matids);
        cmeshcollapsed->Reference()->ResetReference();
        cmeshcollapsed->AdjustBoundaryElements();
        cmeshcollapsed->CleanUpUnconnectedNodes();
    }
    
    cmeshcollapsed->LoadReferences();
    cmeshcollapsed->Reference()->ResetReference();
    cmeshcollapsed->CleanUpUnconnectedNodes();
    cmeshcollapsed->ExpandSolution();
    
    return cmeshcollapsed;
}


TPZCompMesh * pressure(TPZAutoPointer<TPZGeoMesh> gmesh, int POrder)
{        
    auto fmat = Emat1;
    auto dim = 2; auto nstate = 1;
    if (gHdivcollapsed || gSbfemhdiv)
    {
        fmat = ESkeleton;
        dim = 1;
    }

    auto cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(POrder);
    cmesh->SetDimModel(dim);
    
    if (gHdivcollapsed || gSbfemhdiv)
    {
        auto mat = new TPZNullMaterial(ESkeleton, dim, nstate);
        cmesh->InsertMaterialObject(mat); //Insere material na malha
        if (gHybrid)
        {
            auto matleft = new TPZNullMaterial(Eleftpressure, dim, nstate);
            cmesh->InsertMaterialObject(matleft);

            auto matright = new TPZNullMaterial(Erightpressure, dim, nstate);
            cmesh->InsertMaterialObject(matright);
        }
        
    }
    else
    {
        auto mat = new TPZMixedPoisson(fmat, dim);
        cmesh->InsertMaterialObject(mat); //Insere material na malha
    }

    std::set<int> matid = {ESkeleton, Eleftpressure, Erightpressure};
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->SetDimModel(dim);
    cmesh->AutoBuild(matid);

    for(auto newnod : cmesh->ConnectVec())
    {
        newnod.SetLagrangeMultiplier(1);
    }

    return cmesh;
}

TPZMultiphysicsCompMesh *  multiphysics(TPZAutoPointer<TPZGeoMesh> gmesh, TPZCompMesh * cmeshp, TPZCompMesh * cmeshf, int POrder)
{
    gmesh->ResetReference();
    auto dim = gmesh->Dimension();
    auto fmat = Emat1;

    auto cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDefaultOrder(POrder);
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    
    auto mat = new TPZMixedPoisson(fmat, dim);
    
    mat->SetForcingFunction(LaplaceExact.ForcingFunction());
    mat->SetForcingFunctionExact(LaplaceExact.Exact());
    cmesh->InsertMaterialObject(mat);
    
    TPZFMatrix<STATE> val1(2,2,0.), val2(3,1,0.);
    {
        TPZMaterial * bcond = mat->CreateBC(mat, Ebc1, 0, val1, val2); 
        bcond->SetForcingFunction(LaplaceExact.TensorFunction());
        cmesh->InsertMaterialObject(bcond);
    }
    {
        TPZMaterial * bcond = mat->CreateBC(mat, Ebc2, 0, val1, val2); 
        bcond->SetForcingFunction(LaplaceExact.TensorFunction());
        cmesh->InsertMaterialObject(bcond);
    }
    {
        TPZMaterial * bcond = mat->CreateBC(mat, Ebc3, 0, val1, val2); 
        bcond->SetForcingFunction(LaplaceExact.TensorFunction());
        cmesh->InsertMaterialObject(bcond);
    }
    {
        TPZMaterial * bcond = mat->CreateBC(mat, Ebc4, 0, val1, val2); 
        bcond->SetForcingFunction(LaplaceExact.TensorFunction());
        cmesh->InsertMaterialObject(bcond);
    }
    std::cout << "Creating multiphysics mesh\n";
    TPZManVector< TPZCompMesh *, 2> meshvector(2);
    meshvector[0] = cmeshf;
    meshvector[1] = cmeshp;

    TPZManVector<int> active(2,1);
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->LoadReferences();
    
    return cmesh;
}

TPZMultiphysicsCompMesh *  multiphysicscollapsed(TPZAutoPointer<TPZGeoMesh> gmesh, TPZCompMesh * cmeshp, TPZCompMesh * cmeshf, int POrder)
{
    // gmesh->ResetReference();
    auto fmat = ESkeleton;
    int dim = 2;
    int nstate = 1;

    auto cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDefaultOrder(POrder);
    cmesh->SetDimModel(dim);
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();

    {
        auto mat = new TPZHybridPoissonCollapsed(ESkeleton,dim);
        cmesh->InsertMaterialObject(mat);
    }
    {
        auto mat = new TPZNullMaterial(Eleftpressure, dim, nstate);
        cmesh->InsertMaterialObject(mat);
    }
    {
        auto mat = new TPZNullMaterial(Erightpressure, dim, nstate);
        cmesh->InsertMaterialObject(mat);
    }
    {
        auto mat = new TPZNullMaterial(Eleftflux, dim, nstate);
        cmesh->InsertMaterialObject(mat);
    }
    {
        auto mat = new TPZNullMaterial(Erightflux, dim, nstate);
        cmesh->InsertMaterialObject(mat);
    }
    {
        auto mat = new TPZLagrangeMultiplier(Eint, dim, nstate);
        cmesh->InsertMaterialObject(mat);
    }

    std::cout << "Creating multiphysics mesh\n";
    TPZManVector< TPZCompMesh *, 2> meshvector(2);
    meshvector[0] = cmeshf;
    meshvector[1] = cmeshp;

    TPZManVector<int> active(2,1);
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->LoadReferences();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
}

void AddInterfaceElements(TPZCompMesh * cmesh)
{
    auto gmesh = cmesh->Reference();
    for (auto gel : gmesh->ElementVec())
    {
        if (!gel || gel->MaterialId() != Eleftflux)
        {
            continue;
        }
        auto nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        auto gelsidepr = gelside.HasNeighbour(Eleftpressure);
        if (!gelsidepr)
        {
            DebugStop();
        }

        TPZCompElSide celside = gelside.Reference();
        TPZCompElSide celneigh = gelsidepr.Reference();
        if (!celside || !celneigh) {
            DebugStop();
        }
        TPZGeoElBC gelbc(gelside, Eint);
        int64_t index;
        TPZMultiphysicsInterfaceElement *intf = new TPZMultiphysicsInterfaceElement(*cmesh,gelbc.CreatedElement(),index,celneigh,celside);
    }

    for (auto gel : gmesh->ElementVec())
    {
        if (!gel || gel->MaterialId() != Erightflux)
        {
            continue;
        }
        auto nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        auto gelsidepr = gelside.HasNeighbour(Erightpressure);
        if (!gelsidepr)
        {
            DebugStop();
        }

        TPZCompElSide celside = gelside.Reference();
        TPZCompElSide celneigh = gelsidepr.Reference();
        if (!celside || !celneigh) {
            DebugStop();
        }
        TPZGeoElBC gelbc(gelside, Eint);
        int64_t index;
        TPZMultiphysicsInterfaceElement *intf = new TPZMultiphysicsInterfaceElement(*cmesh,gelbc.CreatedElement(),index,celneigh,celside);
    }
}

void GroupandCondense(TPZCompMesh * cmesh)
{
    auto nel = cmesh->NElements();

    int64_t index;
    TPZElementGroup * elgr = new TPZElementGroup(*cmesh, index);

    TPZCompEl *celgr = cmesh->Element(index);
    if(!celgr) DebugStop();
    
    for (auto cel : cmesh->ElementVec())
    {
        if (!cel)
        {
            continue;
        }
        std::cout << "Element index " << cel->Index() << " ";
        if (cel->Reference())
        {
            std::cout << "matid " << cel->Reference()->MaterialId();
            if (cel->Reference()->MaterialId() == Eleftpressure || cel->Reference()->MaterialId() == Erightpressure)
            {
                std::cout << " not added\n";
                continue;
            }
        }
        if(cel == elgr)
        {
            std::cout << " group not added\n";
            continue;
        }
        elgr->AddElement(cel);

    }
    cmesh->ComputeNodElCon();
    
    {
        std::ofstream out("cmesh.txt");
        cmesh->Print(out);
    }
    bool keepmatrix = false;
    auto cond = new TPZCondensedCompEl(elgr, keepmatrix);
    cmesh->ExpandSolution();
    
}

void AdjustExtPressureConnectivity(TPZCompMesh * cmeshm, TPZCompMesh * cmeshf, TPZManVector<int64_t> &perm)
{ 
    auto newid = -1;
    auto nsides = 4;
    auto dim = cmeshm->Reference()->Dimension();

    TPZStack<int64_t> internalprcon;

    auto ncon = cmeshm->NConnects() - nsides * dim;
    perm.resize(ncon);

    for (auto con : cmeshm->ConnectVec())
    {
        if (con.SequenceNumber() == -1)
        {
            continue;
        }
        perm[con.SequenceNumber()] = con.SequenceNumber();
    }

    int64_t nf = cmeshf->NConnects() - 2*nsides;
    auto id = nf+3*nsides;

    for (int is = 0; is < nsides; is++)
    {
        for (int ic = 0; ic < 3; ic++)
        {
            auto pos = nf + 3*nsides + is*6 + ic;
            perm[pos] = id;
            id++;
        }
    }
    for (int is = 0; is < nsides; is++)
    {
        for (int ic = 0; ic < 3; ic++)
        {
            auto pos = nf + 3*nsides + is*6 + ic + 3;
            perm[pos] = id;
            id++;
        }
    }

}

void ComputeMatrices(TPZFMatrix<STATE> &E0fMat, TPZFMatrix<STATE> &E1fMat, TPZFMatrix<STATE> &E2fMat, TPZMatrix<REAL> &KCond)
{
    auto n = KCond.Rows()/2;
    E0fMat.Resize(n,n);
    E1fMat.Resize(n,n);
    E2fMat.Resize(n,n);

    for (auto i = 0; i < n; i++)
    {
        for (auto j = 0; j < n; j++)
        {
            E0fMat(i,j) = KCond(i,j)*1/4;
            E1fMat(i,j) = KCond(i+n,j)*1/2;
            E2fMat(i,j) = KCond(i+n,j+n);
        }
    }
    
}

void ComputeEigenvalues(TPZFMatrix<STATE> &E0fMat, TPZFMatrix<STATE> &E1fMat, TPZFMatrix<STATE> &E2fMat)
{
    int n = E0fMat.Rows();
    auto dim = 2;
    // int dim = Mesh()->Dimension();
    
    TPZFMatrix<STATE> E0Inv(E0fMat);

    TPZVec<int> pivot(E0Inv.Rows(),0);
    int nwork = 4*n*n + 2*n;
    TPZVec<STATE> work(2*nwork,0.);
    int info=0;
#ifdef STATEdouble
    dgetrf_(&n, &n, &E0Inv(0,0), &n, &pivot[0], &info);
#endif
#ifdef STATEfloat
    sgetrf_(&n, &n, &E0Inv(0,0), &n, &pivot[0], &info);
#endif
    if (info != 0) {
        DebugStop();
    }
#ifdef STATEdouble
    dgetri_(&n, &E0Inv(0,0), &n, &pivot[0], &work[0], &nwork, &info);
#endif
#ifdef STATEfloat
    sgetri_(&n, &E0Inv(0,0), &n, &pivot[0], &work[0], &nwork, &info);
#endif
    if (info != 0) {
        DebugStop();
    }
    
    TPZFMatrix<STATE> globmat(2*n,2*n,0.);
    
#ifdef STATEdouble
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1., &E0Inv(0,0), n, &E1fMat(0,0), n, 0., &globmat(0,0), 2*n);
#elif defined STATEfloat
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1., &E0Inv(0,0), n, &E1fMat(0,0), n, 0., &globmat(0,0), 2*n);
#else
    std::cout << "SBFem does not execute for this configuration\n";
    DebugStop();
#endif
    
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            globmat(i,j+n) = -E0Inv(i,j);
        }
    }
#ifdef STATEdouble
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1., &E1fMat(0,0), n, &globmat(0,0), 2*n, 0., &globmat(n,0), 2*n);
#elif defined STATEfloat
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1., &E1fMat(0,0), n, &globmat(0,0), 2*n, 0., &globmat(n,0), 2*n);
#else
    std::cout << "SBFem does not execute for this configuration\n";
    DebugStop();
#endif
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            globmat(i+n,j) -= E2fMat(i,j);
        }
    }

#ifdef STATEdouble
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, -1., &E1fMat(0,0), n, &E0Inv(0,0), n, 0., &globmat(n,n), 2*n);
#elif defined STATEfloat
    cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, -1., &E1fMat(0,0), n, &E0Inv(0,0), n, 0., &globmat(n,n), 2*n);
#else
    std::cout << "SBFem does not execute for this configuration\n";
    DebugStop();
#endif

    for (int i=0; i<n; i++) {
        globmat(i,i) -= (dim-2)*0.5;
        globmat(i+n,i+n) += (dim-2)*0.5;
    }
    
    TPZFMatrix<STATE> globmatkeep(globmat);
    TPZFMatrix< std::complex<double> > eigenVectors;
    TPZManVector<std::complex<double> > eigenvalues;
    globmatkeep.SolveEigenProblem(eigenvalues, eigenVectors);

    std::cout << eigenvalues << std::endl;

}

// Criar um novo ElementGroup?
// O SBFemElementGroup vai conter a implementação do main? 
// Como que eu incorporo meu projeto externo no RefactorGeom?