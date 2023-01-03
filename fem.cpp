#include "mfem.hpp"
#include "mfem/general/forall.hpp"

#include "calculus.hpp"
#include "fem.hpp"

// We could just setup a GridFunctionCoefficient of the solution in
// H1, and then project that coefficient onto the new L2 space.
// However, this (currently) only runs on the host, and we would like
// to be able to do it on the GPU.
//
// In the projection from H1 to L2, we go from having DoFs in the
// vertices to having them in the middle of the element. Since we have
// a first order elements, the projection from H1 to L2 is simply an
// averaging of the H1 DoFs in element i, to the single L2 dof in the
// center of element i.
//
// First we get the element to dof mappings from the H1 and L2 spaces
// (which we can get read pointers to on the device). Then we can loop
// over the elements (in parallel!), get the 4 DoFs belonging to
// element i in H1, and set the value of the single DoF in L2 for
// element i to be the average of the four.
//
// DISCLAIMER: This is hardcoded for, and will only work for, tetrahedron elements.
void project_h1_to_l2(
        double *l2_vals,          // Write-only pointer to the solution in L2
        const double *h1_vals,    // Read-only pointer to the solution in H1
        const int ne,             // Number of elements
        const int *h1_table_col,  // Read-only pointer to the H1 element-to-DoF table columns
        const int *h1_table_row,  // Read-only pointer to the H1 element-to-DoF table rows
        const int *l2_table_col,  // Read-only pointer to the L2 element-to-DoF table columns
        const int *l2_table_row)  // Read-only pointer to the L2 element-to-DoF table rows
{
    mfem::MFEM_FORALL(i, ne, {
        const int *h1_dofs = &h1_table_row[h1_table_col[i]];
        double avg = h1_vals[h1_dofs[0]]
                   + h1_vals[h1_dofs[1]]
                   + h1_vals[h1_dofs[2]]
                   + h1_vals[h1_dofs[3]];
        avg *= 0.25;
        l2_vals[l2_table_row[l2_table_col[i]]] = avg;
    });
}

// XXX: For some reason amdgcn-link refuses to find and link against the cross
// and vecdot implementations in calculus.{hpp|cpp}. Redefine them here.

// Cross product of two 3D vectors a and b, store in c.
MFEM_HOST_DEVICE
void _vec3_cross(vec3& c, vec3& a, vec3& b)
{
    // c = a x b
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

// Dot product of vectors a and b
MFEM_HOST_DEVICE
double _vec3_dot(vec3& a, vec3& b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void compute_gradient(
        double *gradient,         // Write-only pointer to the gradient in L2
        const double *laplace,    // Read-only pointer to the solution in H1
        const double *vert,       // Read-only pointer to the vertex coordinates
        const int ne,             // Number of elements
        const int nv,             // Number of vertices
        const int *h1_table_col,  // Read-only pointer to the H1 element-to-DoF table columns
        const int *h1_table_row,  // Read-only pointer to the H1 element-to-DoF table rows
        const int *l2_table_col,  // Read-only pointer to the L2 element-to-DoF table columns
        const int *l2_table_row)  // Read-only pointer to the L2 element-to-DoF table rows
{
    mfem::MFEM_FORALL(i, ne, {
        const int *h1_dofs = &h1_table_row[h1_table_col[i]];

        // FIXME: Which vertices are set as i,j,k and h is purely decided by
        // what actually worked. The problem is that we want the vertices in a
        // counter-clockwise orientation: i, j, k, h. If they are not, they
        // won't span a volume at all. If that is the case, it should be as
        // simple as swapping two of the vertices.

        // The vertices are in SoA-format:
        //     x0  x1  x2 ... xn y0 y1 y2 ... yn z0 z1 z2 ... zn
        // where n is the total number of vertices.
        vec3 v_j = { vert[0*nv+h1_dofs[0]], vert[1*nv+h1_dofs[0]], vert[2*nv+h1_dofs[0]] };
        vec3 v_k = { vert[0*nv+h1_dofs[1]], vert[1*nv+h1_dofs[1]], vert[2*nv+h1_dofs[1]] };
        vec3 v_h = { vert[0*nv+h1_dofs[2]], vert[1*nv+h1_dofs[2]], vert[2*nv+h1_dofs[2]] };
        vec3 v_i = { vert[0*nv+h1_dofs[3]], vert[1*nv+h1_dofs[3]], vert[2*nv+h1_dofs[3]] };

        // Get the values of the scalar H1 DoFs at each vertex of the tetrahedron.
        const double f_j = laplace[h1_dofs[0]];
        const double f_k = laplace[h1_dofs[1]];
        const double f_h = laplace[h1_dofs[2]];
        const double f_i = laplace[h1_dofs[3]];

        vec3 v_ik = v_i - v_k;
        vec3 v_hk = v_h - v_k;

        vec3 v_ih = v_i - v_h;
        vec3 v_jh = v_j - v_h;

        vec3 v_ki = v_k - v_i;
        vec3 v_ji = v_j - v_i;

        // The normal vector of the plane adjacent to v_j
        vec3 v_ik_x_v_hk; _vec3_cross(v_ik_x_v_hk, v_ik, v_hk);

        // The normal vector of the plane adjacent to v_k
        vec3 v_ih_x_v_jh; _vec3_cross(v_ih_x_v_jh, v_ih, v_jh);

        // The normal vector of the plane adjacent to v_h
        vec3 v_ki_x_v_ji; _vec3_cross(v_ki_x_v_ji, v_ki, v_ji);

        double six_vol;
        {
            vec3 v_hi = v_h - v_i;
            six_vol = abs(_vec3_dot(v_hi, v_ki_x_v_ji));
        }

        v_ik_x_v_hk *= (f_j - f_i);
        v_ih_x_v_jh *= (f_k - f_i);
        v_ki_x_v_ji *= (f_h - f_i);

        vec3 grad = v_ik_x_v_hk;
        grad += v_ih_x_v_jh;
        grad += v_ki_x_v_ji;
        grad *= (1.0 / six_vol);

        const int l2_dof = l2_table_row[l2_table_col[i]];

        gradient[3*l2_dof+0] = grad[0];
        gradient[3*l2_dof+1] = grad[1];
        gradient[3*l2_dof+2] = grad[2];
    });
}

void interpolate_gradient_to_h1(
        double *h1_vals,          // Write-only pointer to the gradient in H1
        const double *l2_vals,    // Read-only pointer to the gradient in L2
        const int nv,             // Number of vertices
        const int *v2e_table_col, // Read-only pointer to the vertex-to-element table columns
        const int *v2e_table_row, // Read-only pointer to the vertex-to-element table rows
        const int *l2_table_col,  // Read-only pointer to the L2 element-to-DoF table columns
        const int *l2_table_row)  // Read-only pointer to the L2 element-to-DoF table rows
{
    // FIXME (ivhak): So this implementation is not really well defined when
    // using more than one process. Along the mesh cuts, the gradients in the
    // vertex will have different values when coming from the two different
    // sides. E.g. if a vertex v is part of the "cut" in a mesh that is divided
    // into two sub-meshes M_1 and M_2, then the value of the gradient in v
    // when in M_1 will be an average of the elements that have v as a vertex
    // _and_ that are in M_1. The same is true for M_2. So, assuming two
    // processes/sub-meshes, the real value of the gradient in v is
    // 0.5*(grad_v_M_1 + grad_v_M_2).
    // This will probably result on some weird fibers in the cuts.

    // We interpolate the gradients in a DoF i in H1 by taking the average of
    // the gradients in the elements that i is a part of.
    mfem::MFEM_FORALL(i, nv, {
        // Find the elements connected to vertex i
        const int *element_indices = &v2e_table_row[v2e_table_col[i]];
        const int num_elements = v2e_table_col[i+1] - v2e_table_col[i];

        vec3 grad = {0};
        for (int j = 0; j < num_elements; j++) {
            // Find the L2 DoF of this element
            const int l2_dof = l2_table_row[l2_table_col[element_indices[j]]];

            vec3 element_grad = {
                l2_vals[3*l2_dof+0],
                l2_vals[3*l2_dof+1],
                l2_vals[3*l2_dof+2],
            };
            grad += element_grad;
        }
        grad *= (1.0 / (double) num_elements);

        h1_vals[3*i+0] = grad[0];
        h1_vals[3*i+1] = grad[1];
        h1_vals[3*i+2] = grad[2];
    });
}
