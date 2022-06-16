#ifndef TRULYSEAMLESS3D_H
#define TRULYSEAMLESS3D_H

#include <iostream>
#include <math.h>
#include <string>

#define HEXEX_TESTING
#include "HexEx/HexExtractor.hh"
#undef HEXEX_TESTING

#include <algorithm>
#include <iomanip>
#include <queue>
#include <tuple>

#include "HexEx/ExactPredicates.hh"
#include "HexEx/Utils.hh"

#include "HexEx/DerivedExactPredicates.hh"
#include "HexEx/Direction.hh"
#include "HexEx/FileAccessor.hh"

#include <Eigen/Sparse>

namespace TS3D
{
using namespace Eigen;
using namespace HexEx;

class TrulySeamless3D : public HexExtractor
{
  private:
    enum SheetType
    {
        SHEET_NONE = -1,
        SHEET_ALIGNX = 0,
        SHEET_ALIGNY = 1,
        SHEET_ALIGNZ = 2
    };

  private:
    bool m_failFlag = false;

    std::string m_fileName;
    bool m_withBoundaryAlignment;
    bool m_refillIfSheetEmpty;

    // Tracking
    CellProperty<bool> m_cellVisited;
    VertexProperty<bool> m_vertexUpdated;
    HalfFaceProperty<bool> m_orientationType;

    // Sheet Info
    FaceProperty<SheetType> m_alignmentType; // u, v or w : 0,1,2
    FaceProperty<int> m_sheet;               //-2 -> usual face, >=0 -> correspoing sheet id
    EdgeProperty<int> m_branches;
    VertexProperty<bool> m_node;

    // Variables
    PerCellVertexProperty<int> m_nodeSector;        // cell to vertex to sector id
    PerCellVertexProperty<bool> m_cellVisitedCheck; // ensure all sectors are filled
    PerCellVertexProperty<Parameter> m_paramOld;    // For Debugging

    std::vector<std::vector<CellHandle>> m_sectorCells; // id to cells
    std::vector<VertexHandle> m_sectorVh;               // sector id to node vertex handle

    // sheet Information
    std::vector<std::vector<int>> m_sheetNodeSectors; // sheet -> Node sectors (in pairs)
    std::vector<HalfFaceHandle> m_sheetType;          // sheet: one face to determine cut or alignment

    int m_sheetCountCut;
    int m_sheetCountAlign;

    int m_nodeSectorCount;
    int m_branchID;

    std::vector<Vec3d> m_nodeSectorValues;

    std::map<int, Transition> m_transitions;
    std::map<int, Transition> m_transitionsInv;

    int m_totalSheetID;

  public:
    TrulySeamless3D();

    template <typename MeshT>
    TrulySeamless3D(const MeshT& tetmesh) : TrulySeamless3D()
    {
        convertToHexExTetrahedralMesh(tetmesh, inputMesh);
        for (auto ch : inputMesh.cells())
            cellVertices[ch] = inputMesh.get_cell_vertices(ch);
    }

    template <typename MeshT, typename ParameterT>
    TrulySeamless3D(const MeshT& tetmesh, PerCellVertexProperty<ParameterT>& parameters) : TrulySeamless3D(tetmesh)
    {
        for (auto ch : tetmesh.cells())
            for (auto cv_it = tetmesh.cv_iter(ch); cv_it.valid(); ++cv_it)
                vertexParameters[ch][*cv_it] = toVec3d(parameters[ch][*cv_it]);
        init();
    }

    TrulySeamless3D(std::string& fileName) : TrulySeamless3D()
    {
        HexEx::readFromFile(fileName, inputMesh, vertexParameters);
        for (auto ch : inputMesh.cells())
            cellVertices[ch] = inputMesh.get_cell_vertices(ch);
        init();
    }

    bool sanitize(double perturb = 0.0, bool keepOriginalTransitions = true);

    Parameter& parameter(CellHandle ch, VertexHandle vh)
    {
        return vertexParameters[ch][vh];
    }

    bool init();

    void writeToFile(const std::string& fileName)
    {
        HexEx::writeToFile(fileName, inputMesh, vertexParameters);
    }

  private:
    // Compute Transition Functions
    void extractTransitionFunction(FaceHandle fh);
    void extractTransitionFunctions();
    const Transition& getTransitionFunction(HalfFaceHandle hfh);

    void calculateEdgeSingularity(EdgeHandle eh);
    void calculateEdgeSingularities();

    void doTransition(HalfFaceHandle hfh, CellHandle& ch);
    void doTransition(HalfFaceHandle hfh, std::vector<Parameter>& params);
    void doTransition(HalfFaceHandle hfh, Parameter& parameter);
    void doTransition(HalfFaceHandle hfh, Direction& dir);
    void doTransition(HalfFaceHandle hfh, Transition& tranFun);

    template <typename T, typename... Rest>
    void doTransition(HalfFaceHandle hfh, T& target, Rest&... rest)
    {
        doTransition(hfh, target);
        doTransition(hfh, rest...);
    }

    template <class T>
    void resizeVec(std::vector<T>& A);
    std::vector<HalfFaceHandle> halffacesAroundHalfedge(HalfEdgeHandle he);
    bool faceContainsEdge(HalfFaceHandle f, VertexHandle& v1, VertexHandle& v2);
    HalfFaceHandle otherEdgeFace(HalfFaceHandle& hf1, HalfEdgeHandle he);
    bool sameRotation(Transition& tranFun1, Transition& tranFun2);
    bool sameRotation(HalfFaceHandle& hf1, HalfFaceHandle& hf2);

    std::pair<int, int> identityTransitionCount();
    int orientationCount();

    void pushTransitionsOntoCutgraph();

    void perturbParametrization(double perturb);

    void mark();
    void markSheets();
    void markBranches();
    void markNodes();
    void markNode(VertexHandle vnode);

    void trace();
    void traceBranches();
    void traceSheets();

    void addNodeToSheet(VertexHandle& v, HalfFaceHandle& hf1, std::set<std::pair<int, int>>& sheet_nodes);

    void transformID(Vec3i& v, Vec3i& inv, Transition t);
    void computeSeamlessnessVariables();

    double maxUVW();
    void fixPrecision(double& uv_max);
    std::vector<HalfFaceHandle> getSheetLoop(HalfEdgeHandle he);
    bool containsBranchVertex(HalfFaceHandle hf, VertexHandle node, std::set<VertexHandle>& bv);
    bool fillBranchAxis(VertexHandle node, HalfFaceHandle finp, std::set<VertexHandle>& bv);
    void fillSeamlessParameterization(VectorXd& x, double& uv_max);

    std::list<HalfFaceHandle> fillSector(HalfFaceHandle hfx, VertexHandle v, Parameter p, bool skip_filled = false);
    void fillSector(HalfFaceHandle hfx, VertexHandle v, bool skip_filled);
    void fillSector(std::list<HalfFaceHandle> hlist, VertexHandle v, bool skip_filled);

    bool checkSeamlessness();
};

} // namespace TS3D

#endif // TRULYSEAMLESS3D_H
