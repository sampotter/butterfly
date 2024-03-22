//                                MFEM Example 7
//
// Compile with: make ex7
//
// Sample runs:  ex7 -e 0 -o 2 -r 4
//               ex7 -e 1 -o 2 -r 4 -snap
//               ex7 -e 0 -amr 1
//               ex7 -e 1 -amr 2 -o 2
//
// Description:  This example code demonstrates the use of MFEM to define a
//               triangulation of a unit sphere and a simple isoparametric
//               finite element discretization of the Laplace problem with mass
//               term, -Delta u + u = f.
//
//               The example highlights mesh generation, the use of mesh
//               refinement, high-order meshes and finite elements, as well as
//               surface-based linear and bilinear forms corresponding to the
//               left-hand side and right-hand side of the discrete linear
//               system. Simple local mesh refinement is also demonstrated.
//
//               We recommend viewing Example 1 before viewing this example.

#include <fstream>
#include <iostream>
#include <cstdint>

#include "mfem.hpp"

using namespace std;
using namespace mfem;

// Exact solution and r.h.s., see below for implementation.
double analytic_solution(const Vector &x);
double analytic_rhs(const Vector &x);
void SnapNodes(Mesh &mesh);

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   int elem_type = 1;
   int ref_levels = 2;
   int order = 2;
   bool visualization = 1;

   OptionsParser args(argc, argv);
   args.AddOption(&elem_type, "-e", "--elem",
                  "Type of elements to use: 0 - triangles, 1 - quads.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh uniformly.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   // 2. Generate an initial high-order (surface) mesh on the unit sphere. The
   //    Mesh object represents a 2D mesh in 3 spatial dimensions. We first add
   //    the elements and the vertices of the mesh, and then make it high-order
   //    by specifying a finite element space for its nodes.
   int Nvert = 8, Nelem = 6;
   if (elem_type == 0)
   {
      Nvert = 6;
      Nelem = 8;
   }
   Mesh *mesh = new Mesh(2, Nvert, Nelem, 0, 3);

   if (elem_type == 0) // inscribed octahedron
   {
      const double tri_v[6][3] =
      {
         { 1,  0,  0}, { 0,  1,  0}, {-1,  0,  0},
         { 0, -1,  0}, { 0,  0,  1}, { 0,  0, -1}
      };
      const int tri_e[8][3] =
      {
         {0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4},
         {1, 0, 5}, {2, 1, 5}, {3, 2, 5}, {0, 3, 5}
      };

      for (int j = 0; j < Nvert; j++)
      {
         mesh->AddVertex(tri_v[j]);
      }
      for (int j = 0; j < Nelem; j++)
      {
         int attribute = j + 1;
         mesh->AddTriangle(tri_e[j], attribute);
      }
      mesh->FinalizeTriMesh(1, 1, true);
   }
   else // inscribed cube
   {
      const double quad_v[8][3] =
      {
         {-1, -1, -1}, {+1, -1, -1}, {+1, +1, -1}, {-1, +1, -1},
         {-1, -1, +1}, {+1, -1, +1}, {+1, +1, +1}, {-1, +1, +1}
      };
      const int quad_e[6][4] =
      {
         {3, 2, 1, 0}, {0, 1, 5, 4}, {1, 2, 6, 5},
         {2, 3, 7, 6}, {3, 0, 4, 7}, {4, 5, 6, 7}
      };

      for (int j = 0; j < Nvert; j++)
      {
         mesh->AddVertex(quad_v[j]);
      }
      for (int j = 0; j < Nelem; j++)
      {
         int attribute = j + 1;
         mesh->AddQuad(quad_e[j], attribute);
      }
      mesh->FinalizeQuadMesh(1, 1, true);
   }

   // Set the space for the high-order mesh nodes.
   H1_FECollection fec(order, mesh->Dimension());
   FiniteElementSpace nodal_fes(mesh, &fec, mesh->SpaceDimension());
   mesh->SetNodalFESpace(&nodal_fes);

   // 3. Refine the mesh while snapping nodes to the sphere.
   for (int l = 0; l <= ref_levels; l++)
   {
     mesh->UniformRefinement();
     SnapNodes(*mesh);
   }

   // 4. Define a finite element space on the mesh. Here we use isoparametric
   //    finite elements -- the same as the mesh nodes.
   FiniteElementSpace *fespace = new FiniteElementSpace(mesh, &fec);
   cout << "Number of unknowns: " << fespace->GetTrueVSize() << endl;

   // 5. Set up the linear form b(.) which corresponds to the right-hand side of
   //    the FEM linear system, which in this case is (1,phi_i) where phi_i are
   //    the basis functions in the finite element fespace.
   LinearForm *b = new LinearForm(fespace);
   ConstantCoefficient one(1.0);
   FunctionCoefficient rhs_coef (analytic_rhs);
   FunctionCoefficient sol_coef (analytic_solution);
   b->AddDomainIntegrator(new DomainLFIntegrator(rhs_coef));
   b->Assemble();

   // 7. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
   //    and Mass domain integrators.
   BilinearForm *a = new BilinearForm(fespace);
   a->AddDomainIntegrator(new DiffusionIntegrator(one));
   a->Assemble();
   a->Finalize();

   FILE *fp = NULL;

   BilinearForm *m = new BilinearForm(fespace);
   m->AddDomainIntegrator(new MassIntegrator(one));
   m->Assemble();
   m->Finalize();

   const SparseMatrix& A = a->SpMat();
   const int *A_indptr = A.GetI();
   const int *A_indices = A.GetJ();
   const double *A_data = A.GetData();

   const SparseMatrix& M = m->SpMat();
   const int *M_indptr = M.GetI();
   const int *M_indices = M.GetJ();
   const double *M_data = M.GetData();

   fp = fopen("A_indptr.bin", "w");
   fwrite(A_indptr, sizeof(int), A.NumRows() + 1, fp);
   fclose(fp);

   fp = fopen("A_indices.bin", "w");
   fwrite(A_indices, sizeof(int), A.NumNonZeroElems(), fp);
   fclose(fp);

   fp = fopen("A_data.bin", "w");
   fwrite(A_data, sizeof(double), A.NumNonZeroElems(), fp);
   fclose(fp);

   fp = fopen("M_indptr.bin", "w");
   fwrite(M_indptr, sizeof(int), M.NumRows() + 1, fp);
   fclose(fp);

   fp = fopen("M_indices.bin", "w");
   fwrite(M_indices, sizeof(int), M.NumNonZeroElems(), fp);
   fclose(fp);

   fp = fopen("M_data.bin", "w");
   fwrite(M_data, sizeof(double), M.NumNonZeroElems(), fp);
   fclose(fp);

   const GridFunction *nodes = mesh->GetNodes();
   std::cout << nodes << std::endl;

   fp = fopen("nodes.bin", "w");
   fwrite(nodes->GetData(), sizeof(double), nodes->Size(), fp);
   fclose(fp);

   // 14. Free the used memory.
   delete a;
   delete m;
   delete fespace;
   delete mesh;

   return EXIT_SUCCESS;
}

double analytic_solution(const Vector &x)
{
   double l2 = x(0)*x(0) + x(1)*x(1) + x(2)*x(2);
   return x(0)*x(1)/l2;
}

double analytic_rhs(const Vector &x)
{
   double l2 = x(0)*x(0) + x(1)*x(1) + x(2)*x(2);
   return 7*x(0)*x(1)/l2;
}

void SnapNodes(Mesh &mesh)
{
   GridFunction &nodes = *mesh.GetNodes();
   Vector node(mesh.SpaceDimension());
   for (int i = 0; i < nodes.FESpace()->GetNDofs(); i++)
   {
      for (int d = 0; d < mesh.SpaceDimension(); d++)
      {
         node(d) = nodes(nodes.FESpace()->DofToVDof(i, d));
      }

      node /= node.Norml2();

      for (int d = 0; d < mesh.SpaceDimension(); d++)
      {
         nodes(nodes.FESpace()->DofToVDof(i, d)) = node(d);
      }
   }
   if (mesh.Nonconforming())
   {
      // Snap hanging nodes to the master side.
      Vector tnodes;
      nodes.GetTrueDofs(tnodes);
      nodes.SetFromTrueDofs(tnodes);
   }
}
