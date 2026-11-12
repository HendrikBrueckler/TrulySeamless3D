#ifndef TRULYSEAMLESS3D_H
#define TRULYSEAMLESS3D_H

#include <iomanip>
#include <math.h>
#include <stdexcept>
#include <string>

#define HEXEX_TESTING
#include "HexEx/HexExtractor.hh"
#undef HEXEX_TESTING

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

    // -- temporary debug write tracing --
    std::string dbgPhase;
    std::set<VertexHandle> dbgVertices;

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
    std::vector<CellHandle> m_cellLocal;                  // tet cell -> submesh cell (a cell has at most one clone)
    std::vector<std::vector<VertexHandle>> m_vertexLocal; // tet vertex -> clones in the submesh
    std::vector<std::vector<EdgeHandle>> m_edgeLocal;     // tet edge -> clones in the submesh
    std::vector<std::vector<FaceHandle>> m_faceLocal;     // tet face -> clones in the submesh
    std::vector<std::map<VertexHandle, VertexHandle>> m_tetVertexLocal; // tet + tet vertex -> local vtx

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
            auto eNew = inputMesh.add_edge(m_vertexLocal[vs[0].idx()].front(), m_vertexLocal[vs[1].idx()].front());
            m_edgeLocal[e.idx()].push_back(eNew);
        }

        for (auto f : tetmesh.faces())
        {
            std::vector<VertexHandle> vs;
            for (auto v : tetmesh.get_halfface_vertices(tetmesh.halfface_handle(f, 0)))
                vs.push_back(m_vertexLocal[v.idx()].front());
            auto fNew = inputMesh.add_face(vs);
            m_faceLocal[f.idx()].push_back(fNew);
        }

        // add tets
        for (auto c_it = tetmesh.cells_begin(); c_it != tetmesh.cells_end(); ++c_it)
        {
            std::vector<VertexHandle> vs;
            for (auto v : tetmesh.tet_vertices(*c_it))
            {
                m_tetVertexLocal[c_it->idx()][v] = m_vertexLocal[v.idx()].front();
                vs.push_back(m_vertexLocal[v.idx()].front());
            }
            auto tet = inputMesh.add_cell(vs);
            m_cellLocal[c_it->idx()] = tet;
        }

        for (auto ch : inputMesh.cells())
            cellVertices[ch] = inputMesh.get_cell_vertices(ch);
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
        m_tetVertexLocal = std::vector<std::map<VertexHandle, VertexHandle>>(tetmesh.n_cells());
        std::vector<std::map<EdgeHandle, EdgeHandle>> tetE2newE(tetmesh.n_cells());
        std::vector<std::map<FaceHandle, FaceHandle>> tetF2newF(tetmesh.n_cells());

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
                    m_tetVertexLocal[c.idx()][v] = newV;
            }
        }

        size_t nE = 0;
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
                auto eNew = inputMesh.add_edge(
                    m_tetVertexLocal[tetSeed.idx()].at(v1), m_tetVertexLocal[tetSeed.idx()].at(v2), true);
                if (!eNew.is_valid() || ++nE != inputMesh.n_logical_edges())
                    throw std::logic_error("Impossible to add edge; topological badness");
                m_edgeLocal[e.idx()].push_back(eNew);
                edgeSingularity[eNew] = sing;
                for (auto c : sector)
                    tetE2newE[c.idx()][e] = eNew;
            }
        }

        auto addFaceOf = [&](CellHandle c, FaceHandle f)
        {
            for (auto chf : tetmesh.cell_halffaces(c))
            {
                if (tetmesh.face_handle(chf) != f)
                    continue;
                std::vector<HalfEdgeHandle> hes;
                for (auto hfhe : tetmesh.halfface_halfedges(chf))
                {
                    EdgeHandle eNew = tetE2newE[c.idx()][hfhe.edge_handle()];
                    HalfEdgeHandle heNew = eNew.halfedge_handle(0);
                    if (m_tetVertexLocal[c.idx()].at(tetmesh.to_vertex_handle(hfhe))
                        != inputMesh.to_vertex_handle(heNew))
                        heNew = heNew.opposite_handle();
                    hes.push_back(heNew);
                }
                return inputMesh.add_face(hes);
            }
            return FaceHandle();
        };
        size_t nF = 0;
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
                auto fNew = addFaceOf(c0, f);
                if (!fNew.is_valid() || ++nF != inputMesh.n_logical_faces())
                    throw std::logic_error("Impossible to add face; topological badness");
                m_faceLocal[f.idx()].push_back(fNew);
                tetF2newF[c0.idx()][f] = fNew;

                fNew = addFaceOf(c1, f);
                if (!fNew.is_valid() || ++nF != inputMesh.n_logical_faces())
                    throw std::logic_error("Impossible to add face; topological badness");
                m_faceLocal[f.idx()].push_back(fNew);
                tetF2newF[c1.idx()][f] = fNew;
            }
            else
            {
                auto fNew = addFaceOf(incl0 ? c0 : c1, f);
                m_faceLocal[f.idx()].push_back(fNew);
                if (!fNew.is_valid() || ++nF != inputMesh.n_logical_faces())
                    throw std::logic_error("Impossible to add face; topological badness");
                if (incl0)
                    tetF2newF[c0.idx()][f] = fNew;
                if (incl1)
                    tetF2newF[c1.idx()][f] = fNew;
            }
        }

        size_t nT = 0;
        for (auto c : tetmesh.cells())
        {
            if (!includeCell[c])
                continue;
            std::vector<HalfFaceHandle> hfs;
            for (auto hf : tetmesh.cell_halffaces(c))
            {
                FaceHandle fNew = tetF2newF[c.idx()][hf.face_handle()];
                HalfFaceHandle hfNew = fNew.halfface_handle(0);
                auto vsRef = tetmesh.halfedge_vertices(*tetmesh.hfhe_iter(hf));
                bool hasMatchingHe = false;
                for (auto heNew : inputMesh.halfface_halfedges(hfNew))
                    if (inputMesh.from_vertex_handle(heNew) == m_tetVertexLocal[c.idx()].at(vsRef[0])
                        && inputMesh.to_vertex_handle(heNew) == m_tetVertexLocal[c.idx()].at(vsRef[1]))
                    {
                        hasMatchingHe = true;
                        break;
                    }
                if (!hasMatchingHe)
                    hfNew = hfNew.opposite_handle();
                hfs.push_back(hfNew);
            }
            auto cNew = inputMesh.add_cell(hfs);
            if (!cNew.is_valid() || ++nT != inputMesh.n_logical_cells())
                throw std::logic_error("Impossible to add cell; topological badness");
            m_cellLocal[c.idx()] = cNew;
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
        vertexParameters[ch][m_tetVertexLocal[orig_cell.idx()].at(orig_vertex)] = param;
    }

    Parameter getParam(CellHandle orig_cell, VertexHandle orig_vertex)
    {
        auto ch = m_cellLocal[orig_cell.idx()];
        return vertexParameters[ch][m_tetVertexLocal[orig_cell.idx()].at(orig_vertex)];
    }

    bool sanitize(double perturb = 0.0, bool keepOriginalTransitions = true);

    bool init();

    void writeToFile(const std::string& fileName)
    {
        HexEx::writeToFile(fileName, inputMesh, vertexParameters);
    }

  private:
    const Parameter& parameter(CellHandle ch, VertexHandle vh) const
    {
        return vertexParameters[ch].at(vh);
    }

    void set_parameter(CellHandle ch, VertexHandle vh, const Parameter& pNew)
    {
        auto& p = vertexParameters[ch].at(vh);
        if (!dbgPhase.empty() && dbgVertices.count(vh))
        {
            std::ostringstream oss;
            oss << "WRIT [" << dbgPhase << "] cell=" << ch << " v=" << vh << " old=(" << std::setprecision(17) << p[0]
                << "," << p[1] << "," << p[2] << ")"
                << " new=(" << std::setprecision(17) << pNew[0] << "," << pNew[1] << "," << pNew[2] << ")";
            std::cout << oss.str() << std::endl;
        }
        p = pNew;
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
