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
        SHEET_ALIGNZ = 2,
    };

    enum BranchType
    {
        BRANCH_NONE = -1,
        BRANCH_ALIGNX = 0,
        BRANCH_ALIGNY = 1,
        BRANCH_ALIGNZ = 2
    };

  private:
    bool m_failFlag = false;
    bool m_algorithmFinished = false;

    std::string m_fileName;
    bool m_withBoundaryAlignment;
    bool m_refillIfSheetEmpty;

    // Tracking
    CellProperty<bool> m_cellVisited;
    VertexProperty<bool> m_vertexUpdated;
    HalfFaceProperty<bool> m_orientationType;

    // Sheet Info
    FaceProperty<SheetType> m_alignmentType; // u, v or w : 0,1,2
    EdgeProperty<BranchType> m_branchType; // u, v or w : 0,1,2
    EdgeProperty<CellHandle> m_branchCell; // u, v or w : 0,1,2
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

    EdgeProperty<bool> m_edgeFeature;
    VertexProperty<bool> m_vertexFeature;
    FaceProperty<bool> m_faceFeature;

    // Mapping to original
    std::vector<VertexHandle> m_vertexLocal;
    std::vector<CellHandle> m_cellLocal;
    std::vector<EdgeHandle> m_edgeLocal;
    std::vector<FaceHandle> m_faceLocal;

  public:
    TrulySeamless3D();

    template <typename MeshT>
    TrulySeamless3D(const MeshT& tetmesh) : TrulySeamless3D()
    {
        using inPoint = typename MeshT::PointT;

        auto toVec3d = [&](const inPoint& point)
        {
            return Vec3d(point[0], point[1], point[2]);
        };

        m_vertexLocal = std::vector<VertexHandle>(tetmesh.n_vertices());
        m_cellLocal = std::vector<CellHandle>(tetmesh.n_cells());
        m_edgeLocal = std::vector<EdgeHandle>(tetmesh.n_edges());
        m_faceLocal = std::vector<FaceHandle>(tetmesh.n_faces());

        // add vertices
        inputMesh.clear(false);

        for (auto v_it = tetmesh.vertices_begin(); v_it != tetmesh.vertices_end(); ++v_it)
        {
            auto v = inputMesh.add_vertex(toVec3d(tetmesh.vertex(*v_it)));
            m_vertexLocal[v_it->idx()] = v;
        }

        for (auto e: tetmesh.edges())
        {
            auto vs = tetmesh.edge_vertices(e);
            auto eNew = inputMesh.add_edge(m_vertexLocal[vs[0].idx()], m_vertexLocal[vs[1].idx()]);
            m_edgeLocal[e.idx()] = eNew;
        }

        for (auto f: tetmesh.faces())
        {
            std::vector<VertexHandle> vs;
            for (auto v: tetmesh.get_halfface_vertices(tetmesh.halfface_handle(f, 0)))
                vs.push_back(m_vertexLocal[v.idx()]);
            auto fNew = inputMesh.add_face(vs);
            m_faceLocal[f.idx()] = fNew;
        }

        // add tets
        for (auto c_it = tetmesh.cells_begin(); c_it != tetmesh.cells_end(); ++c_it)
        {
            auto vertices = tetmesh.get_cell_vertices(*c_it);
            auto tet = inputMesh.add_cell(vertices);
            m_cellLocal[c_it->idx()] = tet;
        }

        for (auto ch : inputMesh.cells())
            cellVertices[ch] = inputMesh.get_cell_vertices(ch);
    }

    template <typename MeshT, typename ParameterT>
    TrulySeamless3D(const MeshT& tetmesh, PerCellVertexProperty<ParameterT>& parameters) : TrulySeamless3D(tetmesh)
    {
        for (auto ch : tetmesh.cells())
            for (auto cv_it = tetmesh.cv_iter(ch); cv_it.valid(); ++cv_it)
                setParam(ch, *cv_it, toVec3d(parameters[ch][*cv_it]));
        init();
    }

    TrulySeamless3D(std::string& fileName) : TrulySeamless3D()
    {
        HexEx::readFromFile(fileName, inputMesh, vertexParameters);
        for (auto ch : inputMesh.cells())
            cellVertices[ch] = inputMesh.get_cell_vertices(ch);
        init();
    }

    void setFeature(FaceHandle orig_face)
    {
        m_faceFeature[m_faceLocal[orig_face.idx()]] = true;
    }

    void setFeature(EdgeHandle orig_edge)
    {
        m_edgeFeature[m_edgeLocal[orig_edge.idx()]] = true;
    }

    void setFeature(VertexHandle orig_vertex)
    {
        m_vertexFeature[m_vertexLocal[orig_vertex.idx()]] = true;
    }

    void setParam(CellHandle orig_cell, VertexHandle orig_vertex, const Parameter& param)
    {
        vertexParameters[m_cellLocal[orig_cell.idx()]][m_vertexLocal[orig_vertex.idx()]] = param;
    }

    Parameter getParam(CellHandle orig_cell, VertexHandle orig_vertex)
    {
        return vertexParameters[m_cellLocal[orig_cell.idx()]][m_vertexLocal[orig_vertex.idx()]];
    }

    bool sanitize(double perturb = 0.0, bool keepOriginalTransitions = true);

    bool init();

    void writeToFile(const std::string& fileName)
    {
        HexEx::writeToFile(fileName, inputMesh, vertexParameters);
    }

  private:

    Parameter& parameter(CellHandle ch, VertexHandle vh)
    {
        return vertexParameters[ch][vh];
    }

    // Compute Transition Functions
    void extractTransitionFunction(FaceHandle fh);
    void extractTransitionFunctions();
    const Transition& getTransitionFunction(HalfFaceHandle hfh);

    void calculateEdgeSingularity(EdgeHandle eh);
    void calculateEdgeSingularities();
    bool isSingularEdge(EdgeHandle eh);

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

    void debugExportFace(const FaceHandle& f) const;
    void debugExportBranch(int branchID) const;

    void exportOVMFile(const std::set<CellHandle>& cells) const;
    void exportOVMFile(const std::set<FaceHandle>& faces) const;
    void exportOVMFile(const std::set<EdgeHandle>& edges) const;
    void exportOVMFile(const std::set<VertexHandle>& vertices) const;
};

} // namespace TS3D

#endif // TRULYSEAMLESS3D_H
