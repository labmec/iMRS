//
//  ReservoirTools.hpp
//  Buckley_Levertt
//
//  Created by Jose on 5/26/20.
//
#include "TPZCompMeshTools.h"
#include "pzelchdiv.h"
#include "pzshapepiram.h"
#include "TPZOneShapeRestraint.h"
#include "pzshapetriang.h"
#include "pzshapetetra.h"
#include "pzelementgroup.h"
#include "pzsubcmesh.h"
#include "pzcondensedcompel.h"
#include "pzmultiphysicselement.h"
#include "TPZMeshSolution.h"
#include "TPZFastCondensedElement.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZCompMeshTools.h"
#ifndef ReservoirTools_hpp
#define ReservoirTools_hpp

#include <stdio.h>

class TPZReservoirTools
{
public:
    static void CreatedCondensedElements(TPZCompMesh *cmesh, bool KeepOneLagrangian, bool keepmatrix);
    static void CondenseElements(TPZCompMesh *cmesh, char LagrangeLevelNotCondensed, bool keepmatrix);
    static void CondenseElements(TPZCompMesh *cmesh, char LagrangeLevelNotCondensed, bool keepmatrix, std::set<int> matids);
    /// push the connect with LagrangeLevel with highest sequence number to LagrangeDest and apply saddle permute
    static void PushConnectBackward(TPZCompMesh *cmesh, char LagrangeStart, char LagrangeDest);
    
    static void AddDependency(std::vector<std::pair<TPZGeoEl*, std::vector<TPZGeoEl*>>> &fatherAndSons);
    static void TakeFatherSonsCorrespondence(TPZCompMesh *fluxCmesh, TPZVec<int64_t> &subdomain, std::vector<std::pair<TPZGeoEl*, std::vector<TPZGeoEl*>>> &fatherAndSons);
    // duplicate the flux elements and put each of them in a separate subdomain
    static void PutFluxElementsinSubdomain(TPZCompMesh *fluxCmesh, TPZVec<int64_t> &subdomain, std::vector<std::pair<TPZGeoEl*, std::vector<TPZGeoEl*>>> &fatherAndSons);
    static void TakeElementsbyID(TPZGeoMesh *mGeometry, std::map<int, std::vector<TPZGeoEl* >> &interfaces, std::vector<int> &matIds);
    static void FindCondensed(TPZCompEl *cel, TPZStack<TPZFastCondensedElement *> &condensedelements);
    static void GroupNeighbourElements(TPZCompMesh *cmesh, const std::set<int64_t> &seed_elements, std::set<int64_t> &groupindexes);
    static void TakeSeedElements(TPZCompMesh *cmesh,  std::set<int64_t> &seed_elements);
    
};





#endif /* ReservoirTools_hpp */
