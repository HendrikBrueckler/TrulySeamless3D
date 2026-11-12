#ifndef TRULYSEAMLESS3D_H
#define TRULYSEAMLESS3D_H

#include <math.h>
#include <string>

#define HEXEX_TESTING
#include "HexEx/HexExtractor.hh"
#undef HEXEX_TESTING

#include <stdexcept>

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
  public:
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

  public:
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
    HalfFaceProperty<SheetType> m_alignmentType; // u, v or w : 0,1,2
    FaceProperty<bool> m_isCut;                  // u, v or w : 0,1,2
    EdgeProperty<BranchType> m_branchType;       // u, v or w : 0,1,2
    EdgeProperty<CellHandle> m_branchCell;       // u, v or w : 0,1,2
    FaceProperty<int> m_sheet;                   //-2 -> usual face, >=0 -> correspoing sheet id
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
    std::vector<CellHandle> m_cellLocal;                   // tet cell -> submesh cell (a cell has at most one clone)
    std::vector<std::vector<VertexHandle>> m_vertexLocal;  // tet vertex -> clones in the submesh
    std::vector<std::vector<EdgeHandle>> m_edgeLocal;      // tet edge -> clones in the submesh
    std::vector<std::vector<FaceHandle>> m_faceLocal;      // tet face -> clones in the submesh
    PerCellVertexProperty<VertexHandle> m_cellVertexLocal; // submesh cell -> (submesh vertex -> tet vertex)

  public:
    TrulySeamless3D();

    template <typename MeshT>
    TrulySeamless3D(const MeshT& tetmesh)
        : TrulySeamless3D()
    {
        using inPoint = typename MeshT::PointT;

        auto toVec3d = [&](const inPoint& point) { return Vec3d(point[0], point[1], point[2]); };

        m_vertexLocal = std::vector<std::vector<VertexHandle>>(tetmesh.n_vertices());
        m_cellLocal = std::vector<CellHandle>(tetmesh.n_cells());
        m_edgeLocal = std::vector<std::vector<EdgeHandle>>(tetmesh.n_edges());
        m_faceLocal = std::vector<std::vector<FaceHandle>>(tetmesh.n_faces());

        // add vertices
        inputMesh.clear(false);

        for (auto v_it = tetmesh.vertices_begin(); v_it != tetmesh.vertices_end(); ++v_it)
        {
            auto v = inputMesh.add_vertex(toVec3d(tetmesh.vertex(*v_it)));
            m_vertexLocal[v_it->idx()].push_back(v);
        }

        for (auto e : tetmesh.edges())
        {
            auto vs = tetmesh.edge_vertices(e);
            auto eNew = inputMesh.add_edge(m_vertexLocal[vs[0].idx()][0], m_vertexLocal[vs[1].idx()][0]);
            m_edgeLocal[e.idx()].push_back(eNew);
        }

        for (auto f : tetmesh.faces())
        {
            std::vector<VertexHandle> vs;
            for (auto v : tetmesh.get_halfface_vertices(tetmesh.halfface_handle(f, 0)))
                vs.push_back(m_vertexLocal[v.idx()][0]);
            auto fNew = inputMesh.add_face(vs);
            m_faceLocal[f.idx()].push_back(fNew);
        }

        // add tets
        for (auto c_it = tetmesh.cells_begin(); c_it != tetmesh.cells_end(); ++c_it)
        {
            auto vertices = tetmesh.get_cell_vertices(*c_it);
            auto tet = inputMesh.add_cell(vertices);
            m_cellLocal[c_it->idx()] = tet;
        }

        std::vector<VertexHandle> vertexToTet(inputMesh.n_vertices());
        for (size_t i = 0; i < m_vertexLocal.size(); ++i)
            vertexToTet[m_vertexLocal[i].front().idx()] = VertexHandle((int)i);

        for (auto ch : inputMesh.cells())
        {
            cellVertices[ch] = inputMesh.get_cell_vertices(ch);
            for (auto vh : cellVertices[ch])
                m_cellVertexLocal[ch][vh] = vertexToTet[vh.idx()];
        }
    }

    // cut faces become literal boundary faces; only explicitly marked feature faces should
    // receive boundary alignment
    template <typename MeshT>
    TrulySeamless3D(const MeshT& tetmesh,
                    CellProperty<bool>& includeCell,
                    FaceProperty<bool>& nonCutFace,
                    EdgeProperty<int>& edgeParametricValence)
        : TrulySeamless3D()
    {
        using inPoint = typename MeshT::PointT;

        auto toVec3d = [&](const inPoint& point) { return Vec3d(point[0], point[1], point[2]); };

        m_withBoundaryAlignment = false;

        m_vertexLocal = std::vector<std::vector<VertexHandle>>(tetmesh.n_vertices());
        m_edgeLocal = std::vector<std::vector<EdgeHandle>>(tetmesh.n_edges());
        m_faceLocal = std::vector<std::vector<FaceHandle>>(tetmesh.n_faces());
        m_cellLocal = std::vector<CellHandle>(tetmesh.n_cells());

        inputMesh.clear(false);

        // corner (cell, vertex) -> new vertex of the sector containing the corner
        std::vector<std::array<VertexHandle, 4>> cornerVtx(tetmesh.n_cells());

        auto cornerOf = [&](CellHandle c, VertexHandle v)
        {
            auto cvs = tetmesh.get_cell_vertices(c);
            for (int i = 0; i < 4; ++i)
                if (cvs[i] == v)
                    return i;
            return -1;
        };

        auto sectorBFS = [&](CellHandle cSeed, const std::function<bool(HalfFaceHandle)>& containsX)
        {
            std::set<CellHandle> sector;
            std::list<CellHandle> q;
            q.push_back(cSeed);
            sector.insert(cSeed);
            while (!q.empty())
            {
                auto cCurr = q.front();
                q.pop_front();
                for (auto chf : tetmesh.cell_halffaces(cCurr))
                {
                    if (!containsX(chf) || !nonCutFace[tetmesh.face_handle(chf)])
                        continue;
                    auto hfOpp = tetmesh.opposite_halfface_handle(chf);
                    auto cNext = tetmesh.incident_cell(hfOpp);
                    if (!cNext.is_valid() || !includeCell[cNext] || sector.contains(cNext))
                        continue;
                    sector.insert(cNext);
                    q.push_back(cNext);
                }
            }
            return sector;
        };

        for (auto v : tetmesh.vertices())
        {
            std::set<CellHandle> visited;
            for (auto tetSeed : tetmesh.vertex_cells(v))
            {
                if (visited.count(tetSeed) || !includeCell[tetSeed])
                    continue;
                auto sector = sectorBFS(tetSeed,
                                        [&](HalfFaceHandle chf)
                                        {
                                            for (auto hfv : tetmesh.get_halfface_vertices(chf))
                                                if (hfv == v)
                                                    return true;
                                            return false;
                                        });
                visited.insert(sector.begin(), sector.end());
                auto newV = inputMesh.add_vertex(toVec3d(tetmesh.vertex(v)));
                m_vertexLocal[v.idx()].push_back(newV);
                for (auto c : sector)
                    cornerVtx[c.idx()][cornerOf(c, v)] = newV;
            }
        }

        for (auto e : tetmesh.edges())
        {
            bool sing = edgeParametricValence[e] != 0;
            auto vs = tetmesh.edge_vertices(e);
            auto v1 = vs[0], v2 = vs[1];
            std::set<CellHandle> visited;
            for (auto tetSeed : tetmesh.edge_cells(e))
            {
                if (visited.count(tetSeed) || !includeCell[tetSeed])
                    continue;
                auto sector = sectorBFS(tetSeed,
                                        [&](HalfFaceHandle chf)
                                        {
                                            for (auto hfe : tetmesh.halfface_edges(chf))
                                                if (hfe == e)
                                                    return true;
                                            return false;
                                        });
                visited.insert(sector.begin(), sector.end());
                auto eNew = inputMesh.add_edge(cornerVtx[tetSeed.idx()][cornerOf(tetSeed, v1)],
                                               cornerVtx[tetSeed.idx()][cornerOf(tetSeed, v2)]);
                m_edgeLocal[e.idx()].push_back(eNew);
                edgeSingularity[eNew] = sing;
            }
        }

        auto addFaceOf = [&](CellHandle c, FaceHandle f)
        {
            for (auto chf : tetmesh.cell(c).halffaces())
            {
                if (tetmesh.face_handle(chf) != f)
                    continue;
                std::vector<VertexHandle> vs;
                for (auto hv : tetmesh.get_halfface_vertices(chf))
                    vs.push_back(cornerVtx[c.idx()][cornerOf(c, hv)]);
                return inputMesh.add_face(vs);
            }
            return FaceHandle();
        };
        for (auto f : tetmesh.faces())
        {
            auto cs = tetmesh.face_cells(f);
            auto c0 = cs[0], c1 = cs[1];
            bool incl0 = c0.is_valid() && includeCell[c0];
            bool incl1 = c1.is_valid() && includeCell[c1];
            if (!incl0 && !incl1)
                continue;
            if (incl0 && incl1 && !nonCutFace[f])
            {
                m_faceLocal[f.idx()].push_back(addFaceOf(c0, f));
                m_faceLocal[f.idx()].push_back(addFaceOf(c1, f));
            }
            else
                m_faceLocal[f.idx()].push_back(addFaceOf(incl0 ? c0 : c1, f));
        }

        for (auto c : tetmesh.cells())
        {
            if (!includeCell[c])
                continue;
            auto cvs = tetmesh.get_cell_vertices(c);
            auto cNew = inputMesh.add_cell(
                cornerVtx[c.idx()][0], cornerVtx[c.idx()][1], cornerVtx[c.idx()][2], cornerVtx[c.idx()][3]);
            m_cellLocal[c.idx()] = cNew;
            for (int i = 0; i < 4; ++i)
                m_cellVertexLocal[cNew][cornerVtx[c.idx()][i]] = cvs[i];
        }

        for (auto ch : inputMesh.cells())
            cellVertices[ch] = inputMesh.get_cell_vertices(ch);

        edgeSingularitiesCalculated = true;
    }

    template <typename MeshT, typename ParameterT>
    TrulySeamless3D(const MeshT& tetmesh, PerCellVertexProperty<ParameterT>& parameters)
        : TrulySeamless3D(tetmesh)
    {
        for (auto ch : tetmesh.cells())
            for (auto cv_it = tetmesh.cv_iter(ch); cv_it.valid(); ++cv_it)
                setParam(ch, *cv_it, toVec3d(parameters[ch][*cv_it]));
        init();
    }

    TrulySeamless3D(std::string& fileName)
        : TrulySeamless3D()
    {
        HexEx::readFromFile(fileName, inputMesh, vertexParameters);
        for (auto ch : inputMesh.cells())
            cellVertices[ch] = inputMesh.get_cell_vertices(ch);
        init();
    }

    void setFeature(FaceHandle orig_face)
    {
        for (auto f : m_faceLocal[orig_face.idx()])
            m_faceFeature[f] = true;
    }

    void setFeature(EdgeHandle orig_edge)
    {
        for (auto e : m_edgeLocal[orig_edge.idx()])
            m_edgeFeature[e] = true;
    }

    void setFeature(VertexHandle orig_vertex)
    {
        for (auto v : m_vertexLocal[orig_vertex.idx()])
            m_vertexFeature[v] = true;
    }

    void setParam(CellHandle orig_cell, VertexHandle orig_vertex, const Parameter& param)
    {
        // orig_cell, orig_vertex are tet handles; the submesh vertex of this tet's cell that
        // corresponds to the tet vertex is unique (a corner belongs to exactly one vertex class)
        auto ch = m_cellLocal[orig_cell.idx()];
        for (auto& e : m_cellVertexLocal[ch])
            if (e.second == orig_vertex)
            {
                vertexParameters[ch][e.first] = param;
                return;
            }
        throw std::logic_error("setParam: vertex not found in cell");
    }

    Parameter getParam(CellHandle orig_cell, VertexHandle orig_vertex)
    {
        auto ch = m_cellLocal[orig_cell.idx()];
        for (auto& e : m_cellVertexLocal[ch])
            if (e.second == orig_vertex)
                return vertexParameters[ch][e.first];
        throw std::logic_error("getParam: vertex not found in cell");
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
        return vertexParameters[ch].at(vh);
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
