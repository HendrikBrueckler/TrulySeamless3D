#include <fstream>
#include <iomanip>
#include <iostream>

#include <list>

#include "TS3D/helpers.h"
#include "TS3D/timer.h"
#include "TS3D/trulyseamless.h"
#include <Eigen/Dense>

#include <OpenVolumeMesh/FileManager/FileManager.hh>

const double M = 1000000000;

namespace TS3D
{
using namespace OpenVolumeMesh;
using Vec3d = HexEx::Vec3d;
using Vec3i = HexEx::Vec3i;

int nonzero(SparseMatrixi& m)
{
    int count = 0;
    for (int i = 0; i < m.outerSize(); i++)
    {
        for (SparseMatrixi::InnerIterator it(m, i); it; ++it)
        {
            if (it.value() != 0)
                count++;
        }
    }
    return count;
}

Vec3d absolute(Vec3d x)
{
    for (int i = 0; i < 3; i++)
    {
        if (x[i] < 0.0)
            x[i] *= -1.0;
    }
    return x;
}

inline double randomDouble(double a, double b)
{
    double range = b - a;
    return a + (double(std::rand()) / double(RAND_MAX)) * range;
}

inline double randomDouble(double a)
{
    return randomDouble(-a, a);
}

double volume(Vec3d u, Vec3d v, Vec3d w)
{
    Eigen::Matrix3d m;
    m << u[0], v[0], w[0], u[1], v[1], w[1], u[2], v[2], w[2];
    return m.determinant();
}

Vec3d randomVec(double a)
{
    return Vec3d(randomDouble(a), randomDouble(a), randomDouble(a));
}

TrulySeamless3D::TrulySeamless3D()
    : HexExtractor(), m_cellVisited(inputMesh.request_cell_property<bool>()),
      m_vertexUpdated(inputMesh.request_vertex_property<bool>()),
      m_orientationType(inputMesh.request_halfface_property<bool>()),
      m_alignmentType(inputMesh.request_face_property<SheetType>()),
      m_branchType(inputMesh.request_edge_property<BranchType>()),
      m_branchCell(inputMesh.request_edge_property<CellHandle>()), m_sheet(inputMesh.request_face_property<int>()),
      m_branches(inputMesh.request_edge_property<int>()), m_node(inputMesh.request_vertex_property<bool>()),
      m_nodeSector(inputMesh.request_cell_property<VertexMapProp<int>>()),
      m_cellVisitedCheck(inputMesh.request_cell_property<VertexMapProp<bool>>()),
      m_paramOld(inputMesh.request_cell_property<VertexMapProp<Parameter>>()),
      m_edgeFeature(inputMesh.request_edge_property<bool>()),
      m_vertexFeature(inputMesh.request_vertex_property<bool>()), m_faceFeature(inputMesh.request_face_property<bool>())
{
    m_withBoundaryAlignment = true;
    m_refillIfSheetEmpty = true;
}

bool TrulySeamless3D::init()
{
    m_withBoundaryAlignment = true;
    m_refillIfSheetEmpty = true;
    for (auto ch : inputMesh.cells())
    {
        cellVertices[ch] = inputMesh.get_cell_vertices(ch);
        for (VertexHandle v : cellVertices[ch])
            m_nodeSector[ch][v] = -1;
    }

    extractTransitionFunctions();

    m_failFlag = false;
    for (auto ch : inputMesh.cells())
    {
        if (isCellDegenerate(ch))
        {
            m_failFlag = true;
        }
    }

#ifndef TRULYSEAMLESS_SILENT
    if (m_failFlag)
        printf("Degenerate triangles!\n");
#endif
    return !m_failFlag;
}

void TrulySeamless3D::extractTransitionFunction(FaceHandle fh)
{
    auto hf0 = inputMesh.halfface_handle(fh, 0);
    auto hf1 = inputMesh.halfface_handle(fh, 1);

    if (inputMesh.is_boundary(fh)) // no transition for boundary faces
    {
        setTransitionFunction(hf0, identity);
        setTransitionFunction(hf1, identity);
    }
    else
    {
        auto c0 = inputMesh.incident_cell(hf0);
        auto c1 = inputMesh.incident_cell(hf1);

        auto vertices0 = inputMesh.get_cell_vertices(hf0);

        auto u0 = parameter(c0, vertices0[0]);
        auto v0 = parameter(c0, vertices0[1]);
        auto w0 = parameter(c0, vertices0[2]);

        auto u1 = parameter(c1, vertices0[0]);
        auto v1 = parameter(c1, vertices0[1]);
        auto w1 = parameter(c1, vertices0[2]);

        // first check if transition function is identity
        if ((u0 == u1) && (v0 == v1) && (w0 == w1))
        {
            setTransitionFunction(hf0, identity);
            setTransitionFunction(hf1, identity);
        }
        else
        {
            auto v0_ = v0 - u0;
            auto w0_ = w0 - u0;

            auto v1_ = v1 - u1;
            auto w1_ = w1 - u1;

            auto min_dist = std::numeric_limits<double>::max();
            auto min_transition = identity;
            for (auto i = 0u; i < all24Transitions.size(); ++i)
            {
                auto v0_transformed = all24Transitions[i].transform_point(v0_);
                auto dist1 = (v1_ - v0_transformed) | (v1_ - v0_transformed);
                auto w0_transformed = all24Transitions[i].transform_point(w0_);
                auto dist2 = (w1_ - w0_transformed) | (w1_ - w0_transformed);

                if (dist1 + dist2 < min_dist)
                {
                    min_dist = dist1 + dist2;
                    min_transition = all24Transitions[i];
                }
            }

            auto u0_t = min_transition.transform_point(u0);
            auto t = u1 - u0_t;

            setTranslation(min_transition, t);

            setTransitionFunction(hf0, min_transition);
            min_transition.invert();
            setTransitionFunction(hf1, min_transition);
        }
    }
}

void TrulySeamless3D::extractTransitionFunctions()
{
    auto n = inputMesh.n_faces();

    //#pragma omp parallel for
    for (auto i = 0u; i < n; ++i)
    {
        auto fh = FaceHandle(i);
        if (inputMesh.is_deleted(fh))
            continue;
        extractTransitionFunction(fh);
    }
    transitionFunctionsComputed = true;
}

const Transition& TrulySeamless3D::getTransitionFunction(HalfFaceHandle hfh)
{
    if (!transitionFunctionsComputed)
        extractTransitionFunctions();
    return transitionFunctions[hfh];
}

void TrulySeamless3D::doTransition(HalfFaceHandle hfh, CellHandle& ch)
{
    ch = inputMesh.incident_cell(inputMesh.opposite_halfface_handle(hfh));
}

void TrulySeamless3D::doTransition(HalfFaceHandle hfh, std::vector<Parameter>& params)
{
    auto ch = inputMesh.incident_cell(inputMesh.opposite_halfface_handle(hfh));
    if (ch.is_valid())
        params = getParameters(ch);
}

void TrulySeamless3D::doTransition(HalfFaceHandle hfh, Parameter& parameter)
{
    auto tranFun = transitionFunctions[hfh];
    setTranslation(tranFun, Vec3d(0, 0, 0));
    parameter = tranFun.transform_point(parameter);
}

void TrulySeamless3D::doTransition(HalfFaceHandle hfh, Direction& dir)
{
    auto tranFun = transitionFunctions[hfh];
    setTranslation(tranFun, Vec3d(0, 0, 0));
    dir.transform(tranFun);
}

void TrulySeamless3D::doTransition(HalfFaceHandle hfh, Transition& tranFun)
{
    auto additionalTranFun = getTransitionFunction(hfh);
    setTranslation(additionalTranFun, Vec3d(0, 0, 0));
    tranFun = additionalTranFun * tranFun;
}

void TrulySeamless3D::calculateEdgeSingularity(EdgeHandle eh)
{
    if (inputMesh.is_boundary(eh))
    {
        // for safety, all boundary edges incident to a degenerate cell are considered singular
        for (auto hec_it = inputMesh.hec_iter(inputMesh.halfedge_handle(eh, 0)); hec_it.valid(); ++hec_it)
        {
            if (isCellDegenerate(*hec_it))
            {
                edgeSingularity[eh] = true;
                return;
            }
        }

        // get first boundary halfface
        auto heh = inputMesh.halfedge_handle(eh, 0);
        auto boundaryHalfFace1 = *inputMesh.hehf_iter(heh);
        while (!inputMesh.is_boundary(inputMesh.opposite_halfface_handle(boundaryHalfFace1)))
        {
            auto currentCell = inputMesh.incident_cell(inputMesh.opposite_halfface_handle(boundaryHalfFace1));
            auto transitionFace = rotateAroundHalfedge(currentCell, heh, false);
            boundaryHalfFace1 = transitionFace;
        }

        // find other boundary face and transition
        auto tranFun = identity;
        auto boundaryHalfFace2 = inputMesh.adjacent_halfface_in_cell(boundaryHalfFace1, heh);

        for (auto i = 0; i < (int)inputMesh.valence(eh) - 2; ++i)
        {
            auto transitionFace = boundaryHalfFace2;
            if (isFaceDegenerate(transitionFace))
            {
                edgeSingularity[eh] = true;
                return;
            }
            doTransition(transitionFace, tranFun);
            boundaryHalfFace2
                = inputMesh.adjacent_halfface_in_cell(inputMesh.opposite_halfface_handle(boundaryHalfFace2), heh);
        }

        auto n0 = getParameterNormal(boundaryHalfFace1);
        auto n1 = getParameterNormal(boundaryHalfFace2);
        n0 = tranFun.transform_vector(n0);

        auto res = (n0 | n1) < 0.5;

        if ((res == true) && (edgeValences[eh] == 2))
        {
            edgeSingularity[eh] = false;
            return;
        }

        edgeSingularity[eh] = res;
        return;
    }
    else
    {
        auto tranFun = identity;
        auto heh = inputMesh.halfedge_handle(eh, 0);
        auto currentCell = *inputMesh.hec_iter(heh);
        for (auto i = 0u; i < inputMesh.valence(eh); ++i)
        {
            auto transitionFace = rotateAroundHalfedge(currentCell, heh);
            if (isFaceDegenerate(transitionFace))
            {
                edgeSingularity[eh] = true;
                return;
            }
            doTransition(transitionFace, currentCell, tranFun);
        }

        if (tranFun == identity)
        {
            edgeSingularity[eh] = edgeValences[eh] > 6;
            if (edgeValences[eh] > 6)
            {
                edgeSingularity[eh] = false;
            }
            return;
        }
        else
        {
            edgeSingularity[eh] = true;
            return;
        }
    }
}

void TrulySeamless3D::calculateEdgeSingularities()
{
    for (auto e_it = inputMesh.edges_begin(); e_it != inputMesh.edges_end(); ++e_it)
        edgeValences[*e_it] = edgeValence(*e_it);

    for (auto eh : inputMesh.edges())
        calculateEdgeSingularity(eh);
    edgeSingularitiesCalculated = true;
}

bool TrulySeamless3D::isSingularEdge(EdgeHandle eh)
{
    if (!edgeSingularitiesCalculated)
        calculateEdgeSingularities();
    return edgeSingularity[eh];
}

template <class T>
void TrulySeamless3D::resizeVec(std::vector<T>& A)
{
    T x;
    A.push_back(x);
}

void TrulySeamless3D::pushTransitionsOntoCutgraph()
{
    extractTransitionFunctions();

    // Mark all cells unvisited
    for (auto ch : inputMesh.cells())
        m_cellVisited[ch] = false;

    std::list<HalfFaceHandle> bfs_list;

    // Initialization
    CellHandle ch0 = *(inputMesh.cells_begin());
    for (auto hf : inputMesh.cell(ch0).halffaces())
        bfs_list.push_back(hf);
    m_cellVisited[ch0] = true;

    while (!bfs_list.empty())
    {
        auto hf0 = bfs_list.front();
        bfs_list.pop_front();
        auto hf1 = inputMesh.opposite_halfface_handle(hf0);

        if (inputMesh.is_boundary(hf1)) // Stop at boundary
            continue;

        ch0 = inputMesh.incident_cell(hf0);
        CellHandle ch1 = inputMesh.incident_cell(hf1);

        if (m_cellVisited[ch1]) // Stop if OPPOSITE Cell is already visited
            continue;
        m_cellVisited[ch1] = true; // Mark visited

        // Push Transition to its neighbour =>
        // Bring opposite tet to this side of the "seam"

        // Copy parameterization for 3 vertices of the common face
        // Apply Transition to 4th vertex
        extractTransitionFunction(inputMesh.face_handle(hf0));
        auto vf = inputMesh.get_halfface_vertices(hf0);
        for (VertexHandle v : cellVertices[ch1])
        {
            auto p = parameter(ch0, v);
            if (v == vf[0] || v == vf[1] || v == vf[2])
                parameter(ch1, v) = p;
            else
            {
                auto T = transitionFunctions[hf1];
                parameter(ch1, v) = T.transform_point(parameter(ch1, v));
            }
        }

        // Push faces of this cell in the list
        for (auto hf : inputMesh.cell(ch1).halffaces())
            bfs_list.push_back(hf);
    }

    extractTransitionFunctions();
}

bool TrulySeamless3D::sameRotation(Transition& tranFun1, Transition& tranFun2)
{
    Vec3d t1(1, 2, 3);
    Vec3d t2(1, 2, 3);
    t1 = tranFun1.transform_vector(t1);
    t2 = tranFun2.transform_vector(t2);
    return (t1 == t2);
}

bool TrulySeamless3D::sameRotation(HalfFaceHandle& hf1, HalfFaceHandle& hf2)
{
    auto tranFun1 = transitionFunctions[hf1];
    auto tranFun2 = transitionFunctions[hf2];
    return sameRotation(tranFun1, tranFun2);
}

bool TrulySeamless3D::faceContainsEdge(HalfFaceHandle hf, VertexHandle& v1, VertexHandle& v2)
{
    auto v = inputMesh.get_halfface_vertices(hf);
    bool is_v1 = false, is_v2 = false;
    for (int i = 0; i < 3; i++)
    {
        if (v1 == v[i])
            is_v1 = true;
        else if (v2 == v[i])
            is_v2 = true;
    }

    if (is_v1 && is_v2) // Other face contains same edge in this cell
        return true;
    return false;
}

std::vector<HalfFaceHandle> TrulySeamless3D::halffacesAroundHalfedge(HalfEdgeHandle he)
{
    std::vector<HalfFaceHandle> f_vec;

    for (auto hehfIt = inputMesh.hehf_iter(he); hehfIt.valid(); ++hehfIt)
        f_vec.push_back(*hehfIt);
    return f_vec;

    auto e = inputMesh.halfedge(he);
    auto v1 = e.from_vertex();
    auto v2 = e.to_vertex();

    HalfFaceHandle f_start;
    for (auto hehfIt = inputMesh.hehf_iter(he); hehfIt.valid(); ++hehfIt)
    {
        f_start = *hehfIt;
        auto fh = inputMesh.face_handle(*hehfIt);
        if (inputMesh.is_boundary(fh))
            break;
    }

    if (inputMesh.is_boundary(f_start))
        f_start = inputMesh.opposite_halfface_handle(f_start);

    HalfFaceHandle f_end = f_start;
    do
    {
        f_vec.push_back(f_start);
        if (inputMesh.is_boundary(f_start))
            break;
        bool found = false;
        auto ch = inputMesh.incident_cell(f_start);
        for (auto f : inputMesh.cell(ch).halffaces())
        {
            if (f == f_start)
                continue;
            if (faceContainsEdge(f, v1, v2))
            {
                found = true;
                f_start = inputMesh.opposite_halfface_handle(f);
                break;
            }
        }
        if (!found)
        {
#ifndef TRULYSEAMLESS_SILENT
            std::cout << "ERROR: halffacesAroundHalfedge() - No neighbouring face?" << std::endl;
#endif
            m_failFlag = true;
            return {};
        }
    } while (f_start != f_end);
    return f_vec;
}

HalfFaceHandle TrulySeamless3D::otherEdgeFace(HalfFaceHandle& hf, HalfEdgeHandle he)
{
    if (inputMesh.is_boundary(hf))
    {
#ifndef TRULYSEAMLESS_SILENT
        std::cout << "ERROR: should not be boundary" << std::endl;
#endif
        m_failFlag = true;
        return HalfFaceHandle(-1);
    }

    auto e = inputMesh.halfedge(he);
    auto v1 = e.from_vertex();
    auto v2 = e.to_vertex();
    HalfFaceHandle hf_new = hf;
    while (true)
    {
        bool found = false;
        auto ch = inputMesh.incident_cell(hf_new);
        for (auto f : inputMesh.cell(ch).halffaces())
        {
            if (f == hf_new)
                continue;

            if (faceContainsEdge(f, v1, v2))
            {
                found = true;
                if (m_sheet[inputMesh.face_handle(f)] > -2)
                    return f;
                hf_new = inputMesh.opposite_halfface_handle(f);
                break;
            }
        }
        if (!found)
        {
#ifndef TRULYSEAMLESS_SILENT
            std::cout << "ERROR: No neighbouring face?" << std::endl;
#endif
            m_failFlag = true;
            return HalfFaceHandle(-1);
        }
    }
    return hf;
}

void TrulySeamless3D::markSheets()
{
    for (auto f_it = inputMesh.faces_begin(); f_it != inputMesh.faces_end(); ++f_it)
    {
        auto hf1 = inputMesh.halfface_handle(*f_it, 0);
        auto hf2 = inputMesh.halfface_handle(*f_it, 1);

        // Initialize
        m_orientationType[hf1] = false;
        m_orientationType[hf2] = false;
        m_alignmentType[*f_it] = SHEET_NONE;
        m_sheet[*f_it] = -2;

        if (inputMesh.is_boundary(*f_it)) // alignment Sheets
        {
            if (!m_withBoundaryAlignment && !m_faceFeature[*f_it])
                continue;

            m_sheet[*f_it] = -1;
            if (inputMesh.is_boundary(hf1)) // Non border half-face
                hf1 = hf2;
            auto ch = inputMesh.incident_cell(hf1);
            auto vertices = inputMesh.get_halfface_vertices(hf1);
            auto a = parameter(ch, vertices[0]);
            auto b = parameter(ch, vertices[1]);
            auto c = parameter(ch, vertices[2]);
            auto error = absolute(b - a) + absolute(c - b) + absolute(c - a); // direction for minimum alignment error
            for (int i = 0; i < 3; i++)
            {
                if (error[i] <= error[(i + 1) % 3] && error[i] <= error[(i + 2) % 3])
                {
                    m_alignmentType[*f_it] = (SheetType)i;
                    break;
                }
            }
        }
        else if (m_faceFeature[*f_it])
        {
            m_sheet[*f_it] = -1;
            auto ch = inputMesh.incident_cell(hf1);
            auto vertices = inputMesh.get_halfface_vertices(hf1);
            auto a = parameter(ch, vertices[0]);
            auto b = parameter(ch, vertices[1]);
            auto c = parameter(ch, vertices[2]);
            auto error = absolute(b - a) + absolute(c - b) + absolute(c - a); // direction for minimum alignment error
            for (int i = 0; i < 3; i++)
            {
                if (error[i] <= error[(i + 1) % 3] && error[i] <= error[(i + 2) % 3])
                {
                    m_alignmentType[*f_it] = (SheetType)i;
                    break;
                }
            }
        }
        else
        {
            auto ch1 = inputMesh.incident_cell(hf1);
            auto ch2 = inputMesh.incident_cell(hf2);
            auto vertices = inputMesh.get_halfface_vertices(hf1);
            if (parameter(ch1, vertices[0]) != parameter(ch2, vertices[0])
                || parameter(ch1, vertices[1]) != parameter(ch2, vertices[1])
                || parameter(ch1, vertices[2]) != parameter(ch2, vertices[2]))
                m_sheet[*f_it] = -1;
        }
    }
}

void TrulySeamless3D::markBranches()
{
    // Mark Branches: Singularity OR ...
    // Boundary: #non-identity-faces > 2 OR different alignment faces
    // Otherwise: #non-identity-faces != 0, 2
    for (auto e_it = inputMesh.edges_begin(); e_it != inputMesh.edges_end(); ++e_it)
    {
        m_branches[*e_it] = -2;
        if (isSingularEdge(*e_it))
        {
            m_branches[*e_it] = -1;
            continue;
        }

        if (m_edgeFeature[*e_it])
        {
            m_branches[*e_it] = -1;
            continue;
        }

        int sheet_count = 0;
        std::set<int> alignments;
        auto he = inputMesh.halfedge_handle(*e_it, 0);
        std::vector<HalfFaceHandle> fvec = halffacesAroundHalfedge(he);
        for (auto hfx : fvec)
        {
            auto fh = inputMesh.face_handle(hfx);
            if (-1 == m_sheet[fh])
            {
                sheet_count++;
                alignments.insert(m_alignmentType[fh]);
            }
        }

        if (1 == sheet_count || sheet_count > 2 || alignments.size() > 1)
            m_branches[*e_it] = -1;
    }

    for (auto e : inputMesh.edges())
    {
        if (-1 == m_branches[e])
        {
            m_branchType[e] = BRANCH_NONE;
            m_branchCell[e] = *inputMesh.ec_iter(e);

            auto vertices = inputMesh.edge_vertices(e);
            auto a = parameter(m_branchCell[e], vertices[0]);
            auto b = parameter(m_branchCell[e], vertices[1]);
            auto error = absolute(b - a); // direction for minimum alignment error
            for (int i = 0; i < 3; i++)
            {
                if (error[i] >= error[(i + 1) % 3] && error[i] >= error[(i + 2) % 3])
                {
                    m_branchType[e] = (BranchType)i;
                    break;
                }
            }
        }
    }
}

void TrulySeamless3D::markNodes()
{
    // Mark nodes: 1 or more than 2 branch edges are incident
    // If node => mark its variables in its cells
    m_nodeSectorCount = 0;
    int node_count = 0;
    for (auto v_it = inputMesh.vertices_begin(); v_it != inputMesh.vertices_end(); ++v_it)
    {
        m_node[*v_it] = false;
        int n_branches = 0;
        bool singularity_branch = false;
        bool non_singularity_branch = false;
        int internal_branches = 0;
        std::set<int> alignments;

        std::set<CH> tets;
        for (auto tet : inputMesh.vertex_cells(*v_it))
            tets.insert(tet);
        std::map<CH, Transition> tet2rotation;
        std::list<std::pair<CH, Transition>> tetQ;
        auto tetSeed = *tets.begin();
        tet2rotation.insert({tetSeed, Transition()});
        tetQ.push_back({tetSeed, Transition()});

        while (!tetQ.empty())
        {
            auto tetRotation = tetQ.front();
            tetQ.pop_front();

            for (auto hf : inputMesh.cell_halffaces(tetRotation.first))
            {
                auto tetNext = inputMesh.incident_cell(inputMesh.opposite_halfface_handle(hf));
                if (!tetNext.is_valid() || tets.count(tetNext) == 0 || tet2rotation.count(tetNext) != 0)
                    continue;

                Transition rotationNext
                    = Transition(transitionFunctions[hf].toMatrix() * tetRotation.second.toMatrix());
                rotationNext.setTranslation(Vec3d(0, 0, 0));
                tetQ.push_back({tetNext, rotationNext});
                tet2rotation.insert({tetNext, rotationNext});
            }
        }

        bool hasFeatureBranch = false;
        bool hasAlignedSheet = false;
        bool hasCutSheet = false;
        set<FH> sheetFaces;
        set<int> nSheetFaces;
        for (auto voh_it = inputMesh.voh_iter(*v_it); voh_it.valid(); ++voh_it)
        {
            auto e = inputMesh.edge_handle(*voh_it);
            int n = 0;
            for (auto f : inputMesh.edge_faces(e))
                if (m_sheet[f] > -2)
                {
                    if (m_alignmentType[f] > -1)
                        hasAlignedSheet = true;
                    else
                        hasCutSheet = true;
                    sheetFaces.insert(f);
                    n++;
                }
            nSheetFaces.insert(n);

            if (-1 == m_branches[e])
            {
                n_branches++;
                if (isSingularEdge(e))
                    singularity_branch = true;
                else
                    non_singularity_branch = true;
                if (m_edgeFeature[e])
                    hasFeatureBranch = true;
                if (inputMesh.is_boundary(*v_it) && !inputMesh.is_boundary(*voh_it))
                    internal_branches++;
                if (m_edgeFeature[e] && !isSingularEdge(e))
                {
                    int alignment = m_branchType[e];
                    Vec3d in(0, 0, 0);
                    in[alignment] = 1;
                    Vec3d out = tet2rotation.at(m_branchCell[e]).inverted().transform_vector(in);
                    for (int coord = 0; coord < 3; coord++)
                        if (out[coord] != 0)
                            alignment = coord;
                    alignments.insert(alignment);
                }
            }
        }

        // To debug: mark all feature vs as nodes
        if (1 == n_branches || n_branches > 2 || (singularity_branch && non_singularity_branch) || internal_branches > 0
            || alignments.size() > 1)
        {
            node_count++;
            markNode(*v_it);
        }
        // TODO this may still be overly strict, test
        else if (!inputMesh.is_boundary(*v_it) && hasAlignedSheet && hasCutSheet && n_branches == 0)
        {
            node_count++;
            markNode(*v_it);
        }
        // TODO this may be overly strict (or possibly not strict enough), test
        else if (!inputMesh.is_boundary(*v_it) && hasFeatureBranch && nSheetFaces.size() > 1)
        {
            node_count++;
            markNode(*v_it);
        }
    }

#ifndef TRULYSEAMLESS_SILENT
    printf("Nodes = %d\n", node_count);
    printf("Sectors = %d\n", m_nodeSectorCount);
#endif
}

void TrulySeamless3D::markNode(VertexHandle vnode)
{
    if (m_node[vnode])
        return;
    m_node[vnode] = true;

    // Mark Variables within each sector
    std::set<HalfFaceHandle> f_set; // Collect all local faces for easy vertex check
    for (auto f_it = inputMesh.vf_iter(vnode); f_it.valid(); ++f_it)
    {
        auto hf1 = inputMesh.halfface_handle(*(f_it), 0);
        auto hf2 = inputMesh.halfface_handle(*(f_it), 1);
        if (!inputMesh.is_boundary(hf1))
            f_set.insert(hf1);
        if (!inputMesh.is_boundary(hf2))
            f_set.insert(hf2);
    }

    std::list<HalfFaceHandle> circ_list;
    auto hf0 = inputMesh.halfface_handle(*(inputMesh.vf_iter(vnode)), 0);
    if (!inputMesh.is_boundary(hf0))
        circ_list.push_back(hf0);
    else
        circ_list.push_back(inputMesh.opposite_halfface_handle(hf0));

    while (!circ_list.empty())
    {
        auto circ_hf = circ_list.front();
        circ_list.pop_front();
        if (inputMesh.is_boundary(circ_hf)) // ignore boundary: no sectors there
            continue;

        CellHandle ch = inputMesh.incident_cell(circ_hf);
        if (-1 < m_nodeSector[ch][vnode]) // already visited cell
            continue;

        std::list<HalfFaceHandle> sector_list;
        sector_list.push_back(circ_hf);
        resizeVec(m_sectorCells);
        m_nodeSectorValues.push_back(parameter(ch, vnode)); // Save Original value

        while (!sector_list.empty())
        {
            auto s_hf = sector_list.front();
            auto s_ch = inputMesh.incident_cell(s_hf);
            sector_list.pop_front();

            if (-1 < m_nodeSector[s_ch][vnode])
                continue;

            m_nodeSector[s_ch][vnode] = m_nodeSectorCount;
            m_sectorCells[m_nodeSectorCount].push_back(s_ch);

            // IF Sheet => push opposite face in circ_list ELSE in sector_list
            for (auto hf : inputMesh.cell(s_ch).halffaces())
            {
                if (f_set.end() == f_set.find(hf))
                    continue;
                auto hf_opp = inputMesh.opposite_halfface_handle(hf);
                if (inputMesh.is_boundary(hf_opp) || -1 < m_nodeSector[inputMesh.incident_cell(hf_opp)][vnode])
                    continue;

                if (m_sheet[inputMesh.face_handle(hf_opp)] > -2)
                    circ_list.push_back(hf_opp);
                else
                    sector_list.push_back(hf_opp);
            }
        }

        // store id node relation
        m_sectorVh.push_back(vnode);
        m_nodeSectorCount++;
    }
}

void TrulySeamless3D::traceBranches()
{
    // Trace and mark branches ids
    m_branchID = 0;
    for (auto e_it = inputMesh.edges_begin(); e_it != inputMesh.edges_end(); ++e_it)
    {
        if (-1 != m_branches[*e_it])
            continue;

        // Trace this branch
        std::list<VertexHandle> bfs_list;
        m_branches[*e_it] = m_branchID;

        auto e = inputMesh.edge(*e_it);
        if (!m_node[e.from_vertex()])
            bfs_list.push_back(e.from_vertex());
        if (!m_node[e.to_vertex()])
            bfs_list.push_back(e.to_vertex());

        while (!bfs_list.empty())
        {
            auto v = bfs_list.front();
            bfs_list.pop_front();

            for (auto voh_it = inputMesh.voh_iter(v); voh_it.valid(); ++voh_it)
            {
                if (-1 != m_branches[inputMesh.edge_handle(*voh_it)])
                    continue;
                m_branches[inputMesh.edge_handle(*voh_it)] = m_branchID;

                auto e2 = inputMesh.halfedge(*voh_it);
                if (v != e2.to_vertex() && !m_node[e2.to_vertex()])
                    bfs_list.push_back(e2.to_vertex());
                if (v != e2.from_vertex() && !m_node[e2.from_vertex()])
                    bfs_list.push_back(e2.from_vertex());
            }
        }
        m_branchID++;
    }
#ifndef TRULYSEAMLESS_SILENT
    printf("Total Branches = %d\n", m_branchID);
#endif
}

void TrulySeamless3D::addNodeToSheet(VertexHandle& v, HalfFaceHandle& hf1, std::set<std::pair<int, int>>& sheet_nodes)
{
    if (!m_node[v])
        return;
    int sheet_id = m_sheet[inputMesh.face_handle(hf1)];
    bool boundary = inputMesh.is_boundary(inputMesh.face_handle(hf1));

    int id1 = m_nodeSector[inputMesh.incident_cell(hf1)][v];
    int id2 = -1;
    if (!boundary)
    {
        auto hf2 = inputMesh.opposite_halfface_handle(hf1);
        id2 = m_nodeSector[inputMesh.incident_cell(hf2)][v];
    }
    if (sheet_nodes.insert(std::make_pair(id1, id2)).second)
    {
        m_sheetNodeSectors[sheet_id].push_back(id1);
        if (!boundary) // Save other sector only for CUT-sheets
        {
            m_sheetNodeSectors[sheet_id].push_back(id2);
            sheet_nodes.insert(std::make_pair(id2, id1));
        }
    }
}

void TrulySeamless3D::traceSheets()
{
    // Trace and mark sheet ids
    m_totalSheetID = 0;
    m_sheetCountCut = 0;
    m_sheetCountAlign = 0;
    int invalid_sheets = 0;
    std::set<VertexHandle> invalid_sheet_branch_vertices;

    for (auto f_it = inputMesh.faces_begin(); f_it != inputMesh.faces_end(); ++f_it)
    {
        if (-1 != m_sheet[*f_it])
            continue;

        std::set<VertexHandle> sheet_branch_vertices;
        int f_count = 1;

        // Trace this sheet
        std::list<HalfFaceHandle> bfs_list;
        auto hfSeed = inputMesh.halfface_handle(*f_it, 0);
        if (inputMesh.is_boundary(hfSeed)) // Save only non-boundary halfface
            hfSeed = inputMesh.halfface_handle(*f_it, 1);
        bfs_list.push_back(hfSeed);

        std::vector<HalfFaceHandle> sheet_faces;
        sheet_faces.push_back(hfSeed);

        std::set<std::pair<int, int>> sheet_nodes; // Avoid adding multiple node equations

        m_orientationType[hfSeed] = true;
        resizeVec(m_sheetNodeSectors);
        m_sheet[*f_it] = m_totalSheetID;
        m_sheetType.push_back(hfSeed);
        bool boundary = inputMesh.is_boundary(*f_it);
        int align = m_alignmentType[*f_it];
        if (align > SHEET_NONE)
            m_sheetCountAlign++;
        else
            m_sheetCountCut++;

        while (!bfs_list.empty())
        {
            auto hf1 = bfs_list.front();
            bfs_list.pop_front();
            m_orientationType[hf1] = true;

            for (HalfFaceHalfEdgeIter hfhe_iter = inputMesh.hfhe_iter(hf1); hfhe_iter.valid(); ++hfhe_iter)
            {
                auto eh = inputMesh.edge_handle(*hfhe_iter);
#ifndef TRULYSEAMLESS_SILENT
                if (m_branches[eh] == -1)
                    std::cout << "ERROR: cutBranch id not set?" << std::endl;
#endif

                // Insert if its a new node of sheet
                auto v1 = inputMesh.halfedge(*hfhe_iter).to_vertex();
                auto v2 = inputMesh.halfedge(*hfhe_iter).from_vertex();

                addNodeToSheet(v1, hf1, sheet_nodes);
                addNodeToSheet(v2, hf1, sheet_nodes);
                if (m_branches[eh] > -1)
                {
                    if (!m_node[v1])
                        sheet_branch_vertices.insert(v1);
                    if (!m_node[v2])
                        sheet_branch_vertices.insert(v2);
                }
                else // Edge is not a branch
                {
                    auto f_opp = otherEdgeFace(hf1, *hfhe_iter);
                    auto f_opp_handle = inputMesh.face_handle(f_opp);

                    assert(align == m_alignmentType[f_opp_handle]);

                    if (-1 == m_sheet[f_opp_handle] && align == m_alignmentType[f_opp_handle])
                    {
                        f_count++;
#ifndef TRULYSEAMLESS_SILENT
                        if (!boundary)
                        {
                            if (!sameRotation(hf1, f_opp))
                                printf("ERROR: Different Rotations for neighbor faces\n");
                        }
#endif
                        m_sheet[f_opp_handle] = m_totalSheetID;
                        m_orientationType[f_opp] = true;
                        bfs_list.push_back(f_opp);
                        sheet_faces.push_back(f_opp);
                    }
                }
            }
        }
        if ((boundary && m_sheetNodeSectors[m_totalSheetID].size() < 2)
            || (!boundary && m_sheetNodeSectors[m_totalSheetID].size() < 4))
        {
            invalid_sheets++;

            int node_count = m_sheetNodeSectors[m_totalSheetID].size();
            if (!boundary)
                node_count /= 2;
            for (auto v : sheet_branch_vertices)
            {
                if (m_node[v])
                    continue;

                invalid_sheet_branch_vertices.insert(v);
                node_count++;
                if (2 == node_count)
                    break;
            }
        }
        m_totalSheetID++;
    }
#ifndef TRULYSEAMLESS_SILENT
    printf("Total Sheets = %d (align = %d and cut = %d)\n\n", m_totalSheetID, m_sheetCountAlign, m_sheetCountCut);
#endif

    if (invalid_sheets > 0)
    {
#ifndef TRULYSEAMLESS_SILENT
        printf("Invalid sheets (not enough nodes on them) = %d\n", invalid_sheets);
#endif
        if (m_refillIfSheetEmpty)
        {
            m_refillIfSheetEmpty = false;
#ifndef TRULYSEAMLESS_SILENT
            printf("Adding extra nodes and Refilling\n");
#endif
            // mark new nodes
            for (auto v : invalid_sheet_branch_vertices)
                markNode(v);

            // cleanup tracing
            for (auto e_it = inputMesh.edges_begin(); e_it != inputMesh.edges_end(); ++e_it)
            {
                if (m_branches[*e_it] > -1)
                    m_branches[*e_it] = -1;
            }
            for (auto f_it = inputMesh.faces_begin(); f_it != inputMesh.faces_end(); ++f_it)
            {
                if (m_sheet[*f_it] > -1)
                    m_sheet[*f_it] = -1;
            }
            m_sheetType.clear();
            m_sheetNodeSectors.clear();
            return trace();
        }
    }
}

void TrulySeamless3D::mark()
{
    calculateEdgeSingularities();
    markSheets();
    markBranches();
    markNodes();
}

void TrulySeamless3D::trace()
{
    traceBranches();
    traceSheets();
}

void TrulySeamless3D::transformID(Vec3i& v, Vec3i& inv, Transition t)
{
    Vec3d d(1, 2, 3);
    d = t.transform_vector(d);

    v = Vec3i((int)d[0], (int)d[1], (int)d[2]);
    for (int i = 0; i < 3; i++)
    {
        inv[i] = 1;
        if (v[i] < 0)
        {
            inv[i] = -1;
            v[i] *= -1;
        }
        v[i] -= 1;
    }
}

void TrulySeamless3D::computeSeamlessnessVariables()
{
    typedef Eigen::Triplet<int, int> Triplets_int;

#ifndef TRULYSEAMLESS_SILENT
    std::cout << "Setting EQUATION SYSTEM" << std::endl;
#endif
    // Fill equations
    std::vector<Triplets_int> triplets;
    int eq = 0;
    int sheets_ignored = 0;
    for (unsigned int i = 0; i < m_sheetNodeSectors.size(); i++)
    {
        std::vector<int>& nodeSectors = m_sheetNodeSectors[i];
        if (nodeSectors.size() < 2)
        {
            sheets_ignored++;
            continue;
        }
        auto hf = m_sheetType[i];
        auto f = inputMesh.face_handle(hf);

        bool boundary = inputMesh.is_boundary(f);
        bool align = m_alignmentType[f] > SHEET_NONE;
        if (align)
        {
            int alignmentCoord = m_alignmentType[f];
            if (boundary)
            {
                for (unsigned int j = 1; j < nodeSectors.size(); j++)
                {
                    triplets.push_back(Triplets_int(eq, 3 * nodeSectors[0] + alignmentCoord, 1));
                    triplets.push_back(Triplets_int(eq, 3 * nodeSectors[j] + alignmentCoord, -1));
                    eq++;
                }
            }
            else
            {
                for (unsigned int j = 2; j < nodeSectors.size(); j += 2)
                {
                    triplets.push_back(Triplets_int(eq, 3 * nodeSectors[0] + alignmentCoord, 1));
                    triplets.push_back(Triplets_int(eq, 3 * nodeSectors[j] + alignmentCoord, -1));
                    eq++;
                }
            }
        }
        if (!boundary) // cut sheet
        {
            // EQUATION: - pi*u0_p + u0_n + pi*uj_p - uj_n = 0

            Vec3i inv, ids(1, 2, 3);
            transformID(ids, inv, transitionFunctions[hf]);

            for (unsigned int j = 2; j < nodeSectors.size(); j += 2)
            {
                for (int k = 0; k < 3; k++)
                {
                    triplets.push_back(Triplets_int(eq, 3 * nodeSectors[0] + ids[k], -1 * inv[k]));
                    triplets.push_back(Triplets_int(eq, 3 * nodeSectors[1] + k, 1));
                    triplets.push_back(Triplets_int(eq, 3 * nodeSectors[j] + ids[k], 1 * inv[k]));
                    triplets.push_back(Triplets_int(eq, 3 * nodeSectors[j + 1] + k, -1));
                    eq++;
                }
            }
        }
    }

    std::vector<std::vector<VertexHandle>> featureBranchEndpoints(m_branchID);
    std::vector<std::vector<EdgeHandle>> featureBranchEndEdges(m_branchID);
    for (auto v : inputMesh.vertices())
        if (m_node[v])
            for (auto e : inputMesh.vertex_edges(v))
                if (m_edgeFeature[e] && !isSingularEdge(e))
                {
                    assert(m_branches[e] > -1);
                    featureBranchEndpoints[m_branches[e]].emplace_back(v);
                    featureBranchEndEdges[m_branches[e]].emplace_back(e);
                }

    for (int i = 0; i < m_branchID; i++)
    {
        assert(featureBranchEndpoints[i].size() == 0 || featureBranchEndpoints[i].size() == 2);
        assert(featureBranchEndEdges[i].size() == 0 || featureBranchEndEdges[i].size() == 2);
        if (featureBranchEndEdges[i].size() == 0)
            continue;

        auto v1 = featureBranchEndpoints[i].front();
        auto v2 = featureBranchEndpoints[i].back();
        auto e1 = featureBranchEndEdges[i].front();
        auto e2 = featureBranchEndEdges[i].back();
        int alignment1 = m_branchType[e1];
        auto cell1 = m_branchCell[e1];
        auto cell2 = m_branchCell[e2];

        bool onSheet = false;
        for (auto f : inputMesh.edge_faces(e1))
            if (m_sheet[f] > -1)
            {
                onSheet = true;
                break;
            }
#ifndef NDEBUG
        bool onSheet2 = false;
        for (auto f : inputMesh.edge_faces(e2))
            if (m_sheet[f] > -1)
            {
                onSheet2 = true;
                break;
            }
        assert(onSheet2 == onSheet);

        set<EH> eVisited({e1});
        list<EH> es({e1});
        bool foundNext = true;
        while (foundNext)
        {
            foundNext = false;
            for (auto v : inputMesh.edge_vertices(es.back()))
            {
                if (foundNext)
                    break;
                for (auto e : inputMesh.vertex_edges(v))
                    if (eVisited.count(e) == 0 && m_branches[e] == i)
                    {
                        eVisited.insert(e);
                        es.push_back(e);
                        foundNext = true;
                        break;
                    }
            }
        }
        set<int> sheetsCheck;
        for (auto f : inputMesh.edge_faces(e1))
            if (m_sheet[f] > -1)
                sheetsCheck.insert(m_sheet[f]);
        for (auto e : es)
        {
            set<int> sheets;
            for (auto f : inputMesh.edge_faces(e))
                if (m_sheet[f] > -1)
                    sheets.insert(m_sheet[f]);
            if (sheets != sheetsCheck)
                debugExportBranch(i);
            assert(sheets == sheetsCheck);
        }
#endif
        if (onSheet)
        {
            set<CellHandle> tetVisited({cell1});
            list<CellHandle> tetQ({cell1});

            while (!tetQ.empty())
            {
                auto tet = tetQ.front();
                tetQ.pop_front();

                for (auto hf : inputMesh.cell_halffaces(tet))
                {
                    if (m_sheet[inputMesh.face_handle(hf)] > -1)
                        continue;
                    auto tetNext = inputMesh.incident_cell(inputMesh.opposite_halfface_handle(hf));
                    if (!tetNext.is_valid() || tetVisited.count(tetNext) != 0)
                        continue;
                    bool onBranch = false;
                    for (auto v : inputMesh.cell_vertices(tetNext))
                    {
                        for (auto e : inputMesh.vertex_edges(v))
                        {
                            if (m_branches[e] == i)
                            {
                                onBranch = true;
                                break;
                            }
                        }
                        if (onBranch)
                            break;
                    }
                    if (!onBranch)
                        continue;
                    tetQ.push_back(tetNext);
                    tetVisited.insert(tetNext);
                }
            }
            bool found = false;
            for (auto tet : inputMesh.edge_cells(e2))
                if (tetVisited.count(tet) != 0)
                {
                    found = true;
                    cell2 = tet;
                    break;
                }
            if (!found)
                throw std::logic_error("Branch start end sectors not connected");
        }

        int sector1 = m_nodeSector[cell1][v1];
        int sector2 = m_nodeSector[cell2][v2];

        assert(sector1 < m_nodeSectorCount);
        assert(sector2 < m_nodeSectorCount);
        assert(sector1 >= 0);
        assert(sector2 >= 0);
        assert(sector1 != sector2);

        int alignCoord1 = (alignment1 + 1) % 3;
        int alignCoord2 = (alignment1 + 2) % 3;

        triplets.push_back(Triplets_int(eq, 3 * sector1 + alignCoord1, 1));
        triplets.push_back(Triplets_int(eq, 3 * sector2 + alignCoord1, -1));
        eq++;
        triplets.push_back(Triplets_int(eq, 3 * sector1 + alignCoord2, 1));
        triplets.push_back(Triplets_int(eq, 3 * sector2 + alignCoord2, -1));
        eq++;
    }

    // rows: number of equations
    // cols: number of variables
    int n = 3 * m_nodeSectorCount;

    // Resize matrices
    SparseMatrixi C;
    SparseMatrixi C_orig;
    VectorXd X_bar, X, b;
    VectorXi indexC, indexR;

    C.resize(eq, n);
    C_orig.resize(eq, n);

    X_bar.resize(n);
    X_bar.setZero();

    X.resize(n);
    X.setZero();

    b.resize(eq);
    b.setZero();

#ifndef TRULYSEAMLESS_SILENT
    printf("Overall equations: %ld x %ld\n", C.rows(), C.cols());
#endif

    double m_uv_max = maxUVW();
    for (unsigned int i = 0; i < m_nodeSectorValues.size(); i++)
    {
        for (int k = 0; k < 3; k++)
            X_bar(3 * i + k) = m_nodeSectorValues[i][k];
    }

    C.setFromTriplets(triplets.begin(), triplets.end());
    C_orig = C;

#ifndef TRULYSEAMLESS_SILENT
    cout << "entries(non-zero elements):" << (double)100 * nonzero(C) / (C.outerSize() * C.innerSize())
         << "% : " << C.nonZeros() << endl
         << endl;
#endif

    double total_time = 0.0;

#ifndef TRULYSEAMLESS_SILENT
    cout << "IREF:" << endl;
#endif
    IREF(C, b, indexC, indexR, true, total_time);

#ifndef TRULYSEAMLESS_SILENT
    cout << "IRREF:" << endl;
#endif
    IRREF(C, b, indexC, true, total_time);

    statistic(C, indexC);

#ifndef TRULYSEAMLESS_SILENT
    cout << "evaluation:" << endl;
#endif
    X = evaluate(C, X_bar, b, M, m_uv_max, indexC, indexR);

    VectorXd b_orig = C_orig.cast<double>() * X_bar;
    double max_b = 0;
    for (int i = 0; i < b_orig.size(); i++)
    {
        if (max_b < fabs(b_orig(i)))
        {
            max_b = fabs(b_orig(i));
        }
    }
    VectorXd b_final = C_orig.cast<double>() * X;
    int non_zero = 0;
    for (int i = 0; i < b_final.size(); i++)
    {
        if (0.0 != b_final(i))
            non_zero++;
    }

#ifndef TRULYSEAMLESS_SILENT
    printf("C_ij max input error = %0.17f\n", max_b);
    printf("non zero output entries: %d\n", non_zero);
#endif

    fillSeamlessParameterization(X, m_uv_max);
}

std::list<HalfFaceHandle> TrulySeamless3D::fillSector(HalfFaceHandle hfx, VertexHandle v, Parameter p, bool skip_filled)
{
    std::list<HalfFaceHandle> sector_incident_sheethf;
    std::set<int> sector_incident_sheets;
    sector_incident_sheets.insert(-1);
    sector_incident_sheets.insert(-2);

    if (m_node[v])
        return sector_incident_sheethf;

    std::set<HalfFaceHandle> f_set; // All faces for easy check
    for (auto f_it = inputMesh.vf_iter(v); f_it.valid(); ++f_it)
    {
        auto hf1 = inputMesh.halfface_handle(*(f_it), 0);
        auto hf2 = inputMesh.halfface_handle(*(f_it), 1);
        if (!inputMesh.is_boundary(hf1))
            f_set.insert(hf1);
        if (!inputMesh.is_boundary(hf2))
            f_set.insert(hf2);
    }

    std::set<CellHandle> c_set; // Visited cells
    std::list<HalfFaceHandle> sector_list;
    sector_list.push_back(hfx);
    while (!sector_list.empty())
    {
        auto s_hf = sector_list.front();

        auto s_ch = inputMesh.incident_cell(s_hf);
        sector_list.pop_front();

        if (c_set.end() != c_set.find(s_ch))
            continue;
        c_set.insert(s_ch);

        if (skip_filled && m_cellVisitedCheck[s_ch][v])
            continue;

        parameter(s_ch, v) = p; // Update value in this cell
        m_cellVisitedCheck[s_ch][v] = true;

        // IF !sheet => push opposite face in sector_list
        for (auto hf : inputMesh.cell(s_ch).halffaces())
        {
            if (f_set.end() == f_set.find(hf))
                continue;
            auto hf_opp = inputMesh.opposite_halfface_handle(hf);
            if (inputMesh.is_boundary(hf_opp))
                continue;

            if (-2 == m_sheet[inputMesh.face_handle(hf_opp)])
                sector_list.push_back(hf_opp);
            else if (!m_cellVisitedCheck[inputMesh.incident_cell(hf_opp)][v] || !skip_filled)
                sector_incident_sheethf.push_back(hf_opp);
        }
    }
    return sector_incident_sheethf;
}

void TrulySeamless3D::fillSector(HalfFaceHandle hfx, VertexHandle v, bool skip_filled)
{
    if (m_node[v])
        return;
    if (inputMesh.is_boundary(hfx))
        return;
    // Assume we want to fully fill the cells around vertex v
    for (auto c = inputMesh.vc_iter(v); c.valid(); ++c)
        m_cellVisitedCheck[*c][v] = false;
    auto p_old = parameter(inputMesh.incident_cell(hfx), v);
    std::list<HalfFaceHandle> hlist = fillSector(hfx, v, p_old, skip_filled);

    fillSector(hlist, v, skip_filled);
}

void TrulySeamless3D::fillSector(std::list<HalfFaceHandle> hlist, VertexHandle v, bool skip_filled)
{
    if (m_node[v])
        return;
    while (!hlist.empty())
    {
        auto s_hf = hlist.front();
        hlist.pop_front();

        if (inputMesh.is_boundary(inputMesh.face_handle(s_hf)))
            continue;
        auto s_ch = inputMesh.incident_cell(s_hf);
        if (m_cellVisitedCheck[s_ch][v])
            continue;

        // Get value from opposite side and fill sector
        auto hf_opp = inputMesh.opposite_halfface_handle(s_hf);
        auto p_opp = parameter(inputMesh.incident_cell(hf_opp), v);

        int sheet_id = m_sheet[inputMesh.face_handle(s_hf)];

        auto T = m_transitions[sheet_id];
        if (!m_orientationType[hf_opp])
            T = m_transitionsInv[sheet_id];
        auto p = T.transform_point(p_opp);
        std::list<HalfFaceHandle> hlist2 = fillSector(s_hf, v, p, skip_filled);
        hlist.insert(hlist.end(), hlist2.begin(), hlist2.end());
    }
}

double TrulySeamless3D::maxUVW()
{
    double m = -1.0;
    for (auto c_it = inputMesh.cells_begin(); c_it != inputMesh.cells_end(); ++c_it)
    {
        auto vertices = cellVertices[*c_it];
        for (int v = 0; v < 4; v++)
        {
            auto p = parameter(*c_it, vertices[v]);
            for (int k = 0; k < 3; k++)
                m = std::max(m, fabs(p[k]));
        }
    }
    return 2.0 * m;
}

void TrulySeamless3D::fixPrecision(double& uv_max)
{
    double m = pow(2.0, floor(log2(uv_max)));

    double diff = 0.0;
    for (auto ch : inputMesh.cells())
    {
        for (VertexHandle v : cellVertices[ch])
        {
            Vec3d p = parameter(ch, v);
            Vec3d sign(p[0] < 0 ? -m : m, p[1] < 0 ? -m : m, p[2] < 0 ? -m : m);
            Vec3d q = (p + sign) - sign;
            parameter(ch, v) = q;

            p = p - q;
            diff = std::max(diff, fabs(p[0]));
            diff = std::max(diff, fabs(p[1]));
            diff = std::max(diff, fabs(p[2]));
        }
    }
#ifndef TRULYSEAMLESS_SILENT
    if (0 == diff)
        printf("\nNO CHANGE via precision fixing\n");
    else
        printf("\nmax change via precision fixing = %.17f\n", diff);
#endif
}

std::vector<HalfFaceHandle> TrulySeamless3D::getSheetLoop(HalfEdgeHandle he)
{
    auto e = inputMesh.halfedge(he);
    auto v1 = e.from_vertex();
    auto v2 = e.to_vertex();

    std::vector<HalfFaceHandle> sheet_loop;

    // Get starting face
    HalfFaceHandle hf_start;
    for (auto hehfIt = inputMesh.hehf_iter(he); hehfIt.valid(); ++hehfIt)
    {
        auto f = inputMesh.face_handle(*hehfIt);
        if (m_sheet[f] < 0 && !inputMesh.is_boundary(f))
            continue;

        hf_start = *hehfIt;
        if (inputMesh.is_boundary(hf_start))
            break;
        hf_start = inputMesh.opposite_halfface_handle(hf_start);
        if (inputMesh.is_boundary(hf_start))
            break;
    }
    if (!hf_start.is_valid())
        return sheet_loop;
    hf_start = inputMesh.opposite_halfface_handle(hf_start); // not pointing to boundary
    sheet_loop.push_back(hf_start);

    HalfFaceHandle hf_new = hf_start;
    do
    {
        auto ch = inputMesh.incident_cell(hf_new);
        for (auto hf : inputMesh.cell(ch).halffaces())
        {
            if (hf == hf_new)
                continue;
            if (faceContainsEdge(hf, v1, v2)) // Other face with same edge in this cell
            {
                hf_new = hf;
                break;
            }
        }
        hf_new = inputMesh.opposite_halfface_handle(hf_new); // Jump to next cell around edge
        auto f = inputMesh.face_handle(hf_new);

        if (hf_new != hf_start && m_sheet[f] >= 0)
            sheet_loop.push_back(hf_new);

        if (inputMesh.is_boundary(f)) // Reached another boundary face
            break;
    } while (hf_start != hf_new);

    return sheet_loop;
}

bool TrulySeamless3D::containsBranchVertex(HalfFaceHandle hf, VertexHandle node, std::set<VertexHandle>& bv)
{
    auto fv = inputMesh.get_halfface_vertices(hf);
    for (VertexHandle v : fv)
    {
        if (v != node && bv.end() != bv.find(v))
            return true;
    }
    return false;
}

bool TrulySeamless3D::fillBranchAxis(VertexHandle node, HalfFaceHandle finp, std::set<VertexHandle>& bv)
{
    auto p0 = parameter(inputMesh.incident_cell(finp), node);
    for (auto v : bv)
    {
        for (auto c = inputMesh.vc_iter(v); c.valid(); ++c)
            m_cellVisited[*c] = false;
    }

    std::vector<std::pair<CellHandle, VertexHandle>> s;
    std::map<VertexHandle, HalfFaceHandle> vertex_sectors;
    Vec3d pn = p0;
    bool other_node_found = false;
    VertexHandle node2;

    std::list<HalfFaceHandle> sector_list;
    sector_list.push_back(finp);
    while (!sector_list.empty())
    {
        auto s_hf = sector_list.front();
        auto s_ch = inputMesh.incident_cell(s_hf);
        sector_list.pop_front();

        if (m_cellVisited[s_ch])
            continue;
        m_cellVisited[s_ch] = true;

        // Fill vertices of this cell
        auto scv = inputMesh.get_cell_vertices(s_ch);
        for (VertexHandle v : scv)
        {
            if (bv.end() == bv.find(v)) // Is part of branch?
                continue;
            if (!other_node_found && m_node[v] && pn != parameter(s_ch, v)) // v != node )
            {
                pn = parameter(s_ch, v);
                other_node_found = true;
                node2 = v;
            }
            else if (!m_node[v])
                vertex_sectors[v] = s_hf;         // Fill in the sector axis later
            s.push_back(std::make_pair(s_ch, v)); // Debug
        }

        // Move to neighboring cells
        for (auto chf : inputMesh.cell(s_ch).halffaces())
        {
            auto f = inputMesh.face_handle(chf);
            if (inputMesh.is_boundary(f) || m_sheet[f] >= 0)
                continue;
            auto hf_opp = inputMesh.opposite_halfface_handle(chf);
            if (inputMesh.is_boundary(chf) || inputMesh.is_boundary(hf_opp))
                continue;

            auto ch_opp = inputMesh.incident_cell(hf_opp);
            if (m_cellVisited[ch_opp])
                continue;

            if (containsBranchVertex(hf_opp, node, bv)
                && (!other_node_found || containsBranchVertex(hf_opp, node2, bv)))
                sector_list.push_back(hf_opp);
        }
    }

    if (vertex_sectors.size() <= 0)
        return true;
    if (!other_node_found)
    {
        std::vector<HalfFaceHandle> vec_bad_faces;
        std::set<int> vec_sheets;

        // Fill data
        for (auto f_itxx = inputMesh.vf_iter(node); f_itxx.valid(); ++f_itxx)
        {
            if (m_sheet[*f_itxx] > -1)
            {
                vec_sheets.insert(m_sheet[*f_itxx]);
                auto hf1 = inputMesh.halfface_handle(*f_itxx, 0);
                if (inputMesh.is_boundary(hf1))
                    hf1 = inputMesh.halfface_handle(*f_itxx, 1);
                vec_bad_faces.clear();
                vec_bad_faces.push_back(hf1);
            }
        }

#ifndef TRULYSEAMLESS_SILENT
        printf("Could not find other node along the branch. Is it an issue?\n");
#endif
        return false;
    }

    // Fill the values
    int axis = 0;
    auto er0 = fabs(p0[0] - pn[0]);
    auto er1 = fabs(p0[1] - pn[1]);
    auto er2 = fabs(p0[2] - pn[2]);
    if (er0 < er1 && er2 < er1)
        axis = 1;
    else if (er0 < er2 && er1 < er2)
        axis = 2;

    for (auto vf : vertex_sectors)
    {
        if (m_node[vf.first])
            continue;
        for (auto c = inputMesh.vc_iter(vf.first); c.valid(); ++c)
            m_cellVisitedCheck[*c][vf.first] = false;
    }

    std::map<VertexHandle, std::list<HalfFaceHandle>> hlist;
    for (auto vf : vertex_sectors)
    {
        if (m_node[vf.first])
            continue;
        auto p = parameter(inputMesh.incident_cell(vf.second), vf.first);
        auto pnew = p0;
        pnew[axis] = p[axis];
        std::list<HalfFaceHandle> hl2 = fillSector(vf.second, vf.first, pnew, true);

        if (hlist.end() != hlist.find(vf.first))
        {
            auto hl1 = hlist[vf.first];
            hl1.insert(hl1.end(), hl2.begin(), hl2.end());
            hlist[vf.first] = hl1;
        }
        else
            hlist[vf.first] = hl2;
    }

    // Fill rest of sectors
    for (auto vh : hlist)
        fillSector(vh.second, vh.first, true);

    return true;
}

void TrulySeamless3D::fillSeamlessParameterization(VectorXd& X, double& uv_max)
{
    // Fix max precision for whole data
    fixPrecision(uv_max);

    for (auto ch : inputMesh.cells())
    {
        auto cV = inputMesh.get_cell_vertices(ch);
        for (VertexHandle v : cV)
            m_cellVisitedCheck[ch][v] = false;
    }

    for (auto v : inputMesh.vertices())
        m_vertexUpdated[v] = false;

        // Compute Unique translation per cut-sheet
#ifndef TRULYSEAMLESS_SILENT
    printf("Compute Final Transitions\n");
#endif
    m_transitions.clear();
    m_transitionsInv.clear();

    std::map<int, Vec3d> alignments;
    for (unsigned int i = 0; i < m_sheetNodeSectors.size(); i++)
    {
        auto hf = m_sheetType[i];
        assert(m_orientationType[hf]);
        auto f = inputMesh.face_handle(hf);
        auto vf = inputMesh.get_halfface_vertices(hf);
        auto hf_opp = inputMesh.opposite_halfface_handle(hf);

        Vec3d a, b;
        int u_p = 3 * m_sheetNodeSectors[i][0];
        a = Vec3d(X(u_p), X(u_p + 1), X(u_p + 2));

        if (m_alignmentType[f] != SHEET_NONE)
            alignments[i] = a;
        if (!inputMesh.is_boundary(f))
        {
            int u_n = 3 * m_sheetNodeSectors[i][1];
            b = Vec3d(X(u_n), X(u_n + 1), X(u_n + 2));

            Transition T = transitionFunctions[hf];
            Vec3d a1 = T.transform_vector(a);
            Transition T2 = transitionFunctions[hf_opp];
            Vec3d b1 = T2.transform_vector(b);

            setTranslation(T, Vec3d(0, 0, 0));
            setTranslation(T2, Vec3d(0, 0, 0));
            setTranslation(T, b - a1);
            setTranslation(T2, a - b1);
            m_transitions[i] = T;
            m_transitionsInv[i] = T2;
        }
    }

    // Fill values in Nodes
#ifndef TRULYSEAMLESS_SILENT
    printf("Fill values in Nodes\n");
#endif
    for (unsigned int i = 0; i < m_sectorVh.size(); i++)
    {
        VertexHandle v = m_sectorVh[i];
        m_vertexUpdated[v] = true;
        Vec3d p(X(3 * i), X(3 * i + 1), X(3 * i + 2));
        for (CellHandle ch : m_sectorCells[i])
        {
            parameter(ch, v) = p;
            m_cellVisitedCheck[ch][v] = true;
        }
    }

    // Collect All branch vertices
    std::set<VertexHandle> all_branch_set;
    // Get Branch Vertices
    std::vector<std::set<VertexHandle>> branch_vertices;
    std::vector<int> branch_edge_count;
    std::vector<bool> branch_processed; // Tag to check which branch already processed
    for (int i = 0; i < m_branchID; i++)
    {
        std::set<VertexHandle> s;
        branch_vertices.push_back(s);

        branch_edge_count.push_back(0);
        branch_processed.push_back(false);
    }
    for (auto e_it = inputMesh.edges_begin(); e_it != inputMesh.edges_end(); ++e_it)
    {
        int branch_id = m_branches[*e_it];
        if (branch_id < 0)
            continue;
        auto e = inputMesh.edge(*e_it);
        auto v1 = e.from_vertex();
        auto v2 = e.to_vertex();
        all_branch_set.insert(v1);
        all_branch_set.insert(v2);

        branch_vertices[branch_id].insert(v1);
        branch_vertices[branch_id].insert(v2);

        branch_edge_count[branch_id]++;
    }

    for (int i = 0; i < m_branchID; i++)
    {
        int node_count = 0;
        for (auto v : branch_vertices[i])
            if (m_node[v])
                node_count++;
        if (2 != node_count)
        {
#ifndef TRULYSEAMLESS_SILENT
            printf("WARNING: Number of nodes in a branch = %d\n", node_count);
#endif
        }
    }

#ifndef TRULYSEAMLESS_SILENT
    printf("Fill values in boundary Branches\n");
#endif
    // Fill alignment constraints along the branches
    for (auto e_it = inputMesh.edges_begin(); e_it != inputMesh.edges_end(); ++e_it)
    {
        if (-2 == m_branches[*e_it] || !inputMesh.is_boundary(*e_it)
            || (!isSingularEdge(*e_it) && m_edgeFeature[*e_it]))
            continue;

        // If Both vertices already filled then continue
        auto e = inputMesh.edge(*e_it);
        auto v1 = e.from_vertex();
        auto v2 = e.to_vertex();
        if (m_vertexUpdated[v1] && m_vertexUpdated[v2])
            continue;

        std::vector<HalfFaceHandle> sheet_loop = getSheetLoop(inputMesh.halfedge_handle(*e_it, 0));
        auto f = inputMesh.face_handle(sheet_loop[sheet_loop.size() - 1]);
        if (!inputMesh.is_boundary(f)) // Ignore non boundary branches
            continue;

        int align = m_alignmentType[f];
        int sheet_id = m_sheet[f];
        auto xx = alignments[sheet_id][align];

        // second last: non-boundary halfface
        auto ch_opp = inputMesh.incident_cell(inputMesh.opposite_halfface_handle(sheet_loop[sheet_loop.size() - 1]));
        auto p1 = parameter(ch_opp, v1);
        auto p2 = parameter(ch_opp, v2);
        p1[align] = xx;
        p2[align] = xx;

        // Loop one side
        for (unsigned int i = sheet_loop.size() - 2; i > 0; i--)
        {
            sheet_id = m_sheet[inputMesh.face_handle(sheet_loop[i])];
            auto hf_opp = inputMesh.opposite_halfface_handle(sheet_loop[i]);

            Transition T = m_transitionsInv[sheet_id];
            if (!m_orientationType[hf_opp])
                T = m_transitions[sheet_id];
            p1 = T.transform_point(p1);
            p2 = T.transform_point(p2);
        }

        // Fill final alignment
        // Loop back and fill sectors
        f = inputMesh.face_handle(sheet_loop[0]);
        sheet_id = m_sheet[f];
        align = m_alignmentType[f];
        xx = alignments[sheet_id][align];
        p1[align] = xx;
        p2[align] = xx;

        if (!m_vertexUpdated[v1])
        {
            parameter(inputMesh.incident_cell(sheet_loop[0]), v1) = p1;
            fillSector(sheet_loop[0], v1, true);
        }
        if (!m_vertexUpdated[v2])
        {
            parameter(inputMesh.incident_cell(sheet_loop[0]), v2) = p2;
            fillSector(sheet_loop[0], v2, true);
        }

        m_vertexUpdated[v1] = true;
        m_vertexUpdated[v2] = true;
    }

    // Fill Boundary Sheets first
#ifndef TRULYSEAMLESS_SILENT
    printf("Fill values in non-branch Boundary-Sheet vertices\n");
#endif
    for (auto f_it = inputMesh.faces_begin(); f_it != inputMesh.faces_end(); ++f_it)
    {
        if (m_sheet[*f_it] < 0)
            continue;
        if (!inputMesh.is_boundary(*f_it))
            continue;

        auto hf = inputMesh.halfface_handle(*f_it, 0);
        if (inputMesh.is_boundary(hf))
            hf = inputMesh.halfface_handle(*f_it, 1);

        int sheet_id = m_sheet[*f_it];
        int align = m_alignmentType[*f_it];

        CellHandle c1 = inputMesh.incident_cell(hf);
        auto vx = inputMesh.get_halfface_vertices(hf);
        for (auto v : vx)
        {
            if (m_vertexUpdated[v])
                continue;
            if (all_branch_set.end() != all_branch_set.find(v))
                continue;

            auto p = parameter(c1, v);
            if (align >= 0) // Fill at branches too!
            {
                p[align] = alignments[sheet_id][align];
                parameter(c1, v) = p;
            }

            fillSector(hf, v, true);
            m_vertexUpdated[v] = true;
        }
    }

    // Fill sheets
#ifndef TRULYSEAMLESS_SILENT
    printf("Fill values in non-branch Sheet vertices\n");
#endif
    for (auto f_it = inputMesh.faces_begin(); f_it != inputMesh.faces_end(); ++f_it)
    {
        if (m_sheet[*f_it] < 0)
            continue;
        if (inputMesh.is_boundary(*f_it))
            continue;

        auto hf = inputMesh.halfface_handle(*f_it, 0);
        assert(!inputMesh.is_boundary(hf));
        bool swapped = !m_orientationType[hf];
        if (swapped)
            hf = inputMesh.halfface_handle(*f_it, 1);

        int sheet_id = m_sheet[*f_it];
        int align = m_alignmentType[*f_it];
        if (align > SHEET_NONE && swapped)
        {
            Vec3d idx(0, 0, 0);
            idx[align] = 1.0;
            idx = m_transitionsInv[sheet_id].transform_vector(idx);
            for (int i = 0; i < 3; i++)
                if (idx[i] != 0)
                    align = i;
        }

        assert(m_orientationType[hf]);

        CellHandle c1 = inputMesh.incident_cell(hf);
        auto vx = inputMesh.get_halfface_vertices(hf);
        for (auto v : vx)
        {
            if (m_vertexUpdated[v] || inputMesh.is_boundary(v))
                continue;

            auto p = parameter(c1, v);
            if (align >= 0) // Fill at branches too!
            {
                p[align] = alignments[sheet_id][align];

                parameter(c1, v) = p;
            }

            fillSector(hf, v, true);
            if (all_branch_set.end()
                == all_branch_set.find(v)) // But marked only if not a branch to fill internal cuts, if any
                m_vertexUpdated[v] = true;
        }
    }

    // Fill values in Branches
#ifndef TRULYSEAMLESS_SILENT
    printf("Fill values in singular Branches\n");
#endif
    for (auto v_it = inputMesh.vertices_begin(); v_it != inputMesh.vertices_end(); ++v_it)
    {
        if (!m_node[*v_it])
            continue;

        // Get branches
        for (auto voh_it = inputMesh.voh_iter(*v_it); voh_it.valid(); ++voh_it)
        {
            auto e = inputMesh.edge_handle(*voh_it);
            int branch_id = m_branches[e];
            if (branch_id < 0 || branch_processed[branch_id])
                continue;
            branch_processed[branch_id] = true;
            if (branch_edge_count[branch_id] < 2)
                continue;

            if (inputMesh.is_boundary(e) && (isSingularEdge(e) || !m_edgeFeature[e]))
                continue;

            std::vector<HalfFaceHandle> sheet_loop = getSheetLoop(inputMesh.halfedge_handle(e, 0));
            auto tranFun = identity;
            for (auto sheet_f : sheet_loop)
            {
                auto sheet_id = m_sheet[inputMesh.face_handle(sheet_f)];
                auto TT = m_transitions[sheet_id];
                if (!m_orientationType[sheet_f])
                    TT = m_transitionsInv[sheet_id];
                tranFun = tranFun * TT;
            }

            if (tranFun == identity && !m_edgeFeature[e])
                continue;

            // Trace this branch sector by sector
            bool fixed_axis = false;
            if (sheet_loop.empty())
            {
                for (auto hf : inputMesh.halfedge_halffaces(inputMesh.halfedge_handle(e, 0)))
                    if (fillBranchAxis(*v_it, hf, branch_vertices[branch_id]))
                    {
                        fixed_axis = true;
                        break;
                    }
            }
            else
            {
                for (auto sheet_f : sheet_loop)
                {
                    if (fillBranchAxis(*v_it, sheet_f, branch_vertices[branch_id]))
                    {
                        fixed_axis = true;
                        break;
                    }
                }
            }

            if (!fixed_axis)
            {
#ifndef TRULYSEAMLESS_SILENT
                printf("ERROR: In fixing the axis of a singularity\n");
#endif
                m_failFlag = true;
                return;
            }

            for (auto v : branch_vertices[branch_id])
                m_vertexUpdated[v] = true;
        }
    }

#ifndef TRULYSEAMLESS_SILENT
    printf("Fill values in rest of cut Branches\n");
#endif
    // Fill using Transitions
    for (auto e_it = inputMesh.edges_begin(); e_it != inputMesh.edges_end(); ++e_it)
    {
        if (-2 == m_branches[*e_it])
            continue;

        // If Both vertices already filled then continue
        auto e = inputMesh.edge(*e_it);
        auto v1 = e.from_vertex();
        auto v2 = e.to_vertex();
        if (m_vertexUpdated[v1] && m_vertexUpdated[v2])
            continue;

        std::vector<HalfFaceHandle> sheet_loop = getSheetLoop(inputMesh.halfedge_handle(*e_it, 0));
        for (unsigned int i = 0; i < sheet_loop.size(); i++)
        {
            auto sheet_f = sheet_loop[i];
            if (inputMesh.is_boundary(sheet_f))
                sheet_f = inputMesh.opposite_halfface_handle(sheet_f);
            if (!m_vertexUpdated[v1])
                fillSector(sheet_f, v1, true);
            if (!m_vertexUpdated[v2])
                fillSector(sheet_f, v2, true);
            break;
        }
        m_vertexUpdated[v1] = true;
        m_vertexUpdated[v2] = true;
    }

#ifndef TRULYSEAMLESS_SILENT
    printf("Done : Filled new values\n");
#endif
}

int TrulySeamless3D::orientationCount()
{
    int inverted = 0;
    int degen = 0;
    double min_volume = 1.0;
    double max_volume = 0.0;
    computeCellTypes();
    for (auto c_it = inputMesh.cells_begin(); c_it != inputMesh.cells_end(); ++c_it)
    {
        auto vertices = cellVertices[*c_it];

        auto u = parameter(*c_it, vertices[1]) - parameter(*c_it, vertices[0]);
        auto v = parameter(*c_it, vertices[2]) - parameter(*c_it, vertices[0]);
        auto w = parameter(*c_it, vertices[3]) - parameter(*c_it, vertices[0]);

        double d = volume(u, v, w);
        min_volume = std::min(min_volume, d);
        max_volume = std::max(max_volume, d);

        if (isCellFlipped(*c_it))
            inverted++;
        if (isCellDegenerate(*c_it))
            degen++;
    }
#ifndef TRULYSEAMLESS_SILENT
    std::cout << "Mesh Info: " << inverted << " inverted and " << degen
              << " degenerate cells with volume min = " << min_volume << " max = " << max_volume << std::endl;
#endif
    return degen + inverted;
}

std::pair<int, int> TrulySeamless3D::identityTransitionCount()
{
    extractTransitionFunctions();
    int count1 = 0;
    int count2 = 0;
    for (auto f_it = inputMesh.faces_begin(); f_it != inputMesh.faces_end(); ++f_it)
    {
        auto hf0 = inputMesh.halfface_handle(*f_it, 0);
        if (transitionFunctions[hf0] == identity)
            count1++;
        else
            count2++;
    }
    return std::make_pair(count1, count2);
}

bool TrulySeamless3D::checkSeamlessness()
{
    extractTransitionFunctions();

    int f_align = 0;
    int f_cut = 0;
    int f_identity = 0;
    int v_count = 0;
    int bad_alignment = 0;
    int bad_cut = 0;
    int bad_identity = 0;

    std::set<VertexHandle> error_vset;
    for (auto f_it = inputMesh.faces_begin(); f_it != inputMesh.faces_end(); ++f_it)
    {
        auto hf1 = inputMesh.halfface_handle(*f_it, 0);
        auto hf2 = inputMesh.opposite_halfface_handle(hf1);

        if (inputMesh.is_boundary(*f_it) || m_faceFeature[*f_it])
        {
            f_align++;
            if (inputMesh.is_boundary(hf1)) // Get Non border half-face
                hf1 = hf2;

            auto ch1 = inputMesh.incident_cell(hf1);
            auto vx = inputMesh.get_halfface_vertices(hf1);
            bool face_aligned = false;
            Vec3d a = parameter(ch1, vx[0]);
            Vec3d b = parameter(ch1, vx[1]);
            Vec3d c = parameter(ch1, vx[2]);
            for (int i = 0; i < 3; i++)
            {
                if (a[i] == b[i] && a[i] == c[i])
                {
                    face_aligned = true;
                    break;
                }
            }
            if (!face_aligned)
            {
                bad_alignment++;
                v_count += 3;
            }
        }
        else
        {
            if (m_sheet[*f_it] > -1)
            {
                if (!m_orientationType[hf1])
                    std::swap(hf1, hf2);
            }
            auto T = transitionFunctions[hf1];

            auto ch1 = inputMesh.incident_cell(hf1);
            auto ch2 = inputMesh.incident_cell(hf2);

            bool er = false;
            auto vx = inputMesh.get_halfface_vertices(hf1);

            std::vector<Vec3d> error;
            error.push_back(parameter(ch2, vx[0]) - T.transform_vector(parameter(ch1, vx[0])));
            error.push_back(parameter(ch2, vx[1]) - T.transform_vector(parameter(ch1, vx[1])));
            error.push_back(parameter(ch2, vx[2]) - T.transform_vector(parameter(ch1, vx[2])));
            if (error[0] != error[1] || error[1] != error[2] || error[2] != error[0])
            {

                std::vector<bool> br = {false, false, false};
                for (int i = 0; i < 3; i++)
                {
                    for (auto voh_it = inputMesh.voh_iter(vx[i]); voh_it.valid(); ++voh_it)
                    {
                        auto e = inputMesh.edge_handle(*voh_it);
                        if (m_branches[e] >= 0)
                        {
                            br[i] = true;
                            break;
                        }
                    }
                }

                er = true;
                v_count++;
            }

            if (T == identity)
            {
                f_identity++;
                if (er)
                    bad_identity++;
            }
            else
            {
                f_cut++;
                if (er)
                    bad_cut++;
            }
        }
    }

    double diff = 0.0;
    for (auto ch : inputMesh.cells())
    {
        auto cV = inputMesh.get_cell_vertices(ch);
        for (VertexHandle v : cV)
        {
            auto p = m_paramOld[ch][v] - parameter(ch, v);
            diff = std::max(diff, fabs(p[0]));
            diff = std::max(diff, fabs(p[1]));
            diff = std::max(diff, fabs(p[2]));
        }
    }
#ifndef TRULYSEAMLESS_SILENT
    printf("\nMax change in Parametrization = %f\n", diff);
#endif

    int e_singularAlign = 0;
    int e_featureAlign = 0;
    int bad_singularAlign = 0;
    int bad_featureAlign = 0;

    for (auto e : inputMesh.edges())
    {
        if (m_edgeFeature[e] || isSingularEdge(e))
        {
            if (isSingularEdge(e))
                e_singularAlign++;
            else
                e_featureAlign++;

            auto vs = inputMesh.edge_vertices(e);
            for (auto c : inputMesh.edge_cells(e))
            {
                Vec3d a = parameter(c, vs[0]);
                Vec3d b = parameter(c, vs[1]);
                bool edge_aligned = (a[0] == b[0] && (a[1] == b[1] || a[2] == b[2])) || (a[1] == b[1] && a[2] == b[2]);

                if (!edge_aligned)
                {
                    if (isSingularEdge(e))
                        bad_singularAlign++;
                    else
                    {
                        bad_featureAlign++;
                    }
                    v_count += 2;
                    break;
                }
            }
        }
    }

    if (0 != v_count)
    {
#ifndef TRULYSEAMLESS_SILENT
        std::cout << "Total Alignment Faces: " << f_align << " and cut Faces = " << f_cut << " and identity "
                  << f_identity << endl;
        std::cout << "Total singular Edges: " << e_singularAlign << " and feature edges: " << e_featureAlign << endl;
        std::cout << (m_algorithmFinished ? "ERROR" : "INFO" ) << ": Misaligned singular edges: " << bad_singularAlign << endl;
        std::cout << (m_algorithmFinished ? "ERROR" : "INFO" ) << ": Misaligned feature edges: " << bad_featureAlign << endl;
        std::cout << (m_algorithmFinished ? "ERROR" : "INFO" ) << ": Alignment Faces: " << bad_alignment << endl;
        std::cout << (m_algorithmFinished ? "ERROR" : "INFO" ) << ": " << v_count << " vertices not seamless across " << bad_cut << " cut faces and "
                  << bad_identity << " identity faces\n"
                  << endl;
#endif
        return false;
    }
#ifndef TRULYSEAMLESS_SILENT
    std::cout << "Output is Truly Seamless!" << endl;
#endif
    return true;
}

void TrulySeamless3D::perturbParametrization(double perturb)
{
    if (0.0 == perturb)
        return;
#ifndef TRULYSEAMLESS_SILENT
    printf("\nPerturb the parametrization with range %f\n", perturb);
#endif
    int f_perturbed = 0;
    for (auto f_it = inputMesh.faces_begin(); f_it != inputMesh.faces_end(); ++f_it)
    {
        auto hf1 = inputMesh.halfface_handle(*f_it, 0);
        auto hf2 = inputMesh.halfface_handle(*f_it, 1);

        if (inputMesh.is_boundary(*f_it)) // alignment sheets
        {
            if (!m_withBoundaryAlignment)
                continue;
            if (inputMesh.is_boundary(hf1)) // Non border halfface
                hf1 = hf2;
        }
        else if (transitionFunctions[hf1] == identity) // No need to perturb the non-sheets
            continue;

        auto ch = inputMesh.incident_cell(hf1);
        auto vertices = inputMesh.get_halfface_vertices(hf1);

        for (auto v : vertices)
        {
            auto p1 = parameter(ch, v);
            auto p2 = p1 + randomVec(perturb);
            for (auto c = inputMesh.vc_iter(v); c.valid(); ++c)
            {
                if (p1 == parameter(*c, v))
                    parameter(*c, v) = p2;
            }
        }
        f_perturbed++;
    }
#ifndef TRULYSEAMLESS_SILENT
    std::cout << "Perturbed Faces = " << f_perturbed << std::endl;
#endif

#ifndef TRULYSEAMLESS_SILENT
    std::pair<int, int> c2 = identityTransitionCount();
    std::cout << "Perturbation: #Identity Transitions = " << c2.first << " and #Non-Identity = " << c2.second
              << std::endl
              << std::endl;
#endif

    // Save perturbed mesh
    std::string ff = m_fileName + "-perturbed.hexex";
    HexEx::writeToFile(ff, inputMesh, vertexParameters);
}

bool TrulySeamless3D::sanitize(double perturb, bool keepOriginalTransitions)
{
    if (m_failFlag)
    {
#ifndef TRULYSEAMLESS_SILENT
        std::cout << "ERROR: Not properly initialized" << endl;
#endif
        return false;
    }
    if (orientationCount() > 0)
    {
#ifndef TRULYSEAMLESS_SILENT
        std::cout << "ERROR: Algorithm cannot work with flipped or degenerate tets" << endl;
#endif
        return false;
    }

#ifndef TRULYSEAMLESS_SILENT
    auto t_start = std::chrono::high_resolution_clock::now();
#endif

#ifndef TRULYSEAMLESS_SILENT
    std::pair<int, int> c1 = identityTransitionCount();
    std::cout << "INPUT: #Identity Transitions = " << c1.first << " and #Non-Identity = " << c1.second << std::endl;
#endif

    if (!keepOriginalTransitions)
    {
        pushTransitionsOntoCutgraph();
#ifndef TRULYSEAMLESS_SILENT
        std::pair<int, int> c2 = identityTransitionCount();
        std::cout << "CutGraph: #Identity Transitions = " << c2.first << " and #Non-Identity = " << c2.second
                  << std::endl;
#endif
        extractTransitionFunctions();
    }

    // Copy old values
    for (auto ch : inputMesh.cells())
    {
        auto cV = inputMesh.get_cell_vertices(ch);
        for (VertexHandle v : cV)
            m_paramOld[ch][v] = parameter(ch, v);
    }

    // Perturb the mesh?
    if (perturb != 0.0)
        perturbParametrization(perturb);

    if (checkSeamlessness() && orientationCount() == 0)
        return true;

    // CODE FLOW:
    // 1) Mark cut sheets and alignments
    // 2) Mark all branches
    // 3) Mark nodes and their sector variables
    // 4) Trace and mark branch id
    // 5) Trace and mark sheet id
    mark();
    trace();

    computeSeamlessnessVariables();

    if (m_failFlag)
        return false;

    m_algorithmFinished = true;

#ifndef TRULYSEAMLESS_SILENT
    cout << "Overall running time: " << subtractTimes(t_start) << " ms" << endl;
#endif

    return (checkSeamlessness() && orientationCount() == 0);
}

void TrulySeamless3D::debugExportFace(const FaceHandle& f) const
{
    // Export face itself
    exportOVMFile({f});

    // Export nodes
    set<VH> nodes;
    for (auto v : inputMesh.face_vertices(f))
        if (m_node[v])
            nodes.insert(v);
    if (!nodes.empty())
        exportOVMFile(nodes);

    // Gather and export all branches
    set<int> branches;
    for (auto v : inputMesh.face_vertices(f))
        for (auto e : inputMesh.vertex_edges(v))
            if (m_branches[e] > -1)
                branches.insert(m_branches[e]);
    for (auto branch : branches)
    {
        set<EH> otherEdges;
        for (auto e : inputMesh.edges())
            if (m_branches[e] == branch)
                otherEdges.insert(e);
        exportOVMFile(otherEdges);
    }

    // Gather and export all sheets incident on branch
    set<int> sheets;
    for (auto v : inputMesh.face_vertices(f))
        for (auto f2 : inputMesh.vertex_faces(v))
            if (m_sheet[f2] > -1 && sheets.count(m_sheet[f2]) == 0)
                sheets.insert(m_sheet[f2]);
    for (auto sheet : sheets)
    {
        set<FH> faces;
        for (auto e : inputMesh.faces())
            if (m_sheet[e] == sheet)
                faces.insert(e);
        exportOVMFile(faces);
    }
}

void TrulySeamless3D::debugExportBranch(int branchID) const
{
    set<EH> edges;
    for (auto e : inputMesh.edges())
        if (m_branches[e] == branchID)
            edges.insert(e);

    // Export branch itself
    exportOVMFile(edges);

    // Export nodes
    set<VH> nodes;
    for (auto e : edges)
        for (auto v : inputMesh.edge_vertices(e))
            if (m_node[v])
                nodes.insert(v);
    exportOVMFile(nodes);

    // Gather and export all other branches branching off this branch
    set<int> otherBranches;
    for (auto e1 : edges)
        for (auto v : inputMesh.edge_vertices(e1))
            for (auto e2 : inputMesh.vertex_edges(v))
                if (m_branches[e2] > -1 && m_branches[e2] != branchID)
                    otherBranches.insert(m_branches[e2]);
    for (auto branch : otherBranches)
    {
        set<EH> otherEdges;
        for (auto e : inputMesh.edges())
            if (m_branches[e] == branch)
                otherEdges.insert(e);
        exportOVMFile(otherEdges);
    }

    // Gather and export all sheets incident on branch
    set<int> sheets;
    for (auto e : edges)
        for (auto v : inputMesh.edge_vertices(e))
            for (auto f : inputMesh.vertex_faces(v))
                if (m_sheet[f] > -1 && sheets.count(m_sheet[f]) == 0)
                    sheets.insert(m_sheet[f]);
    for (auto sheet : sheets)
    {
        set<FH> faces;
        for (auto e : inputMesh.faces())
            if (m_sheet[e] == sheet)
                faces.insert(e);
        exportOVMFile(faces);
    }
}

void TrulySeamless3D::exportOVMFile(const set<CH>& cells) const
{
    HexEx::TetrahedralMesh exportMesh;
    std::map<VH, VH> v2vExport;
    for (auto tet : cells)
    {
        std::vector<VH> vs;
        for (auto v : inputMesh.tet_vertices(tet))
        {
            auto vIt = v2vExport.find(v);
            if (vIt == v2vExport.end())
            {
                auto vExport = exportMesh.add_vertex(inputMesh.vertex(v));
                vs.emplace_back(vExport);
            }
            else
                vs.emplace_back(vIt->second);
        }
        exportMesh.add_cell(vs);
    }

    static int i = 0;
    OpenVolumeMesh::IO::FileManager f;
    f.writeFile("./cells" + std::to_string(i++) + ".ovm", exportMesh);
}

void TrulySeamless3D::exportOVMFile(const set<FH>& faces) const
{
    HexEx::TetrahedralMesh exportMesh;
    std::map<VH, VH> v2vExport;
    for (auto f : faces)
    {
        std::vector<VH> vs;
        for (auto v : inputMesh.face_vertices(f))
        {
            auto vIt = v2vExport.find(v);
            if (vIt == v2vExport.end())
            {
                auto vExport = exportMesh.add_vertex(inputMesh.vertex(v));
                vs.emplace_back(vExport);
            }
            else
                vs.emplace_back(vIt->second);
        }
        exportMesh.add_face(vs);
    }

    static int i = 0;
    OpenVolumeMesh::IO::FileManager f;
    f.writeFile("./faces" + std::to_string(i++) + ".ovm", exportMesh);
}

void TrulySeamless3D::exportOVMFile(const set<EH>& edges) const
{
    HexEx::TetrahedralMesh exportMesh;
    std::map<VH, VH> v2vExport;
    for (auto e : edges)
    {
        std::vector<VH> vs;
        for (auto v : inputMesh.edge_vertices(e))
        {
            auto vIt = v2vExport.find(v);
            if (vIt == v2vExport.end())
            {
                auto vExport = exportMesh.add_vertex(inputMesh.vertex(v));
                vs.emplace_back(vExport);
            }
            else
                vs.emplace_back(vIt->second);
        }
        exportMesh.add_edge(vs[0], vs[1]);
    }

    static int i = 0;
    OpenVolumeMesh::IO::FileManager f;
    f.writeFile("./edges" + std::to_string(i++) + ".ovm", exportMesh);
}

void TrulySeamless3D::exportOVMFile(const set<VH>& vertices) const
{
    HexEx::TetrahedralMesh exportMesh;
    for (auto v : vertices)
    {
        exportMesh.add_vertex(inputMesh.vertex(v));
    }

    static int i = 0;
    OpenVolumeMesh::IO::FileManager f;
    f.writeFile("./vertices" + std::to_string(i++) + ".ovm", exportMesh);
}

} // namespace TS3D
