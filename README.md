# Finite-Element-Method-ICES-deal.ii

## Description
The Finite Element Method (FEM) is an essential numerical tool for solving boundary value problems of PDE's. FEM is used here to model the properties of materials that exhibit elastic-plastic deformation. The "Mesh" is a 2-D plane which can easily be extrapolated to a 3-D object. 

Petroleum Engineers often run stress tests on core samples where the loading/unloading Young's modulus, Poisson's ratio, and different Lame's parameters can be obtained. Many combinations of axial and confining pressures can be ran on cores which adds to time and labor in the lab. The following method is the beginning of being able to predict rock strain under certain parameters. This method can be potentially be used to predict rock slippage effects and propagation of fractures in hydraulic fracturing modeling. Simulating reservoir conditions on a Limestone core requires the aforemented inputs to simulate a reservoir environment and predict when a rock will yield.

### Assumptions
  * No micro-fractures or planes for the elements to 'slip' on.
  * Young's modulus of 20GPa.
  * Poisson's ratio of 0.25.
  * Overburden stress at 2500m.

## Method and Materials
* **Deal.ii** is a C++ library that supports the implementation of FEM code. 
* When outlining the degrees of freedom on the mesh, the degrees of freedom are only assigned to the mesh's vertices. 
* The figure shows a general description of the iterative process after the inputs have been defined. The mesh grid and stiffness matrix are assembled before the force matrix is generated. Post-processing begins after the final displacement, strain, and stresses are realized.
* Maximum iterations were used as ‘break’ criteria instead of residual calculations because the iterations needed for convergence are relatively small.
* An order of integration of 2 was used for simplicity. Higher orders of integration require more memory and tend to be more appropriate for higher order elements.

![Force Matrix Generator](https://github.com/FracThePermian/Finite-Element-Method-ICES-deal.ii/blob/master/Force_Matrix_Generator.png "Force_Matrix_Generator")

## Results

* The 'stressed' Mesh only reacts to a vertical force on the mesh. The Mesh with Confining Stress enacts a 3-to-2 vertical to horizontal stress ratio.

![Mesh Stress Pic](https://github.com/FracThePermian/Finite-Element-Method-ICES-deal.ii/blob/master/Mesh_Stress2.png "Mesh Stress Visualized")

* The Stress-Strain diagram is qualitatively similar to one seen performed on a core in the lab.

* When the plastic yield conditions are met, the stress-strain enters the plastic zone as indicated by the horizontal line.

* The Stress-Strain responds accordingly when Young’s Modulus is reduced.

![Stress vs. Strain](https://github.com/FracThePermian/Finite-Element-Method-ICES-deal.ii/blob/master/stressVstrain.PNG "Stress versus Strain")

```cpp
void P_mesh()
{
Triangulation<2> tri_locate; //list of vertices & cells
GridGenerator::hyper_cube(tri_locate); //2-space 
tri_locate.refine_global(8); //refined 2, 4, or 8 times
  std::ofstream out("raw_grid.svg");
  GridOut       mesh_out;
  mesh_out.write_svg(tri_locate, out); //write to output.svg
  std::cout << "Write to raw_grid.svg" << std::endl;
}
```

```cpp
//enumerate all dof
template <int dim>
void TopLevel<dim>::assemble_system()
{
  sys_lim    = 0;
  sys_mesh = 0;
  FEValues<dim> fe_values(hv,
                          quad_form,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);
  const unsigned int cell_dof = fe.n_dofs_per_cell();
  const unsigned int NQP    = quadrature_formula.size();
  FullMatrix<double> cell_matrix(dofs_per_cell, cell_dof);
  Vector<double>     cell_rhs(cell_dof);
  std::vector<types::global_dof_index> local_dof_indices(cell_dof);
  BodyForce<dim>              body_force;
  std::vector<Vector<double>> force_value(NQP, Vector<double>(dim));
 for (const auto &cell : dof_handler.active_cell_iterators())
  if (cell->is_locally_owned())
    {
      hv_matrix = 0;
      cell_rhs = 0;
      elem_value.reinit(cell);
 for (unsigned int i = 0; i < cell_dof; ++i)
  for (unsigned int j = 0; j < cell_dof; ++j)
    for (unsigned int q_p = 0; q_p < NQP; ++q_p)
      {
        const SymmetricTensor<2, dim>
          ep_i = fetch_strain(elem_value, i, q_p),
          ep_j = fetch_strain(elem_value, j, q_p);
        hv_matrix(i, j) += (ep_i * stress_strain_tensor * ep_j) * fe_values.JxW(q_p); 
      }
      
      //raw output then translated to MATLAB for visualization*
```


## Conclusions

* The mesh only either is in an elastic state or plastic state.

* Increasing order of integration consumes significantly more memory but leads to a more accurate result.

* The type of element matters when determining accuracy of results. For example, a 9 node square configuration would create more degrees of freedom.

* The boundary conditions could contribute to the non-symmetrical deforming of the Mesh with Confining Stress.

* Future refining will include a 3-D mesh coupled with a confining pressure that does not unrealistically deform the mesh. This could include a Von-Mises criteria for yield failure along with strain hardening, and strain softening. It would be worthwhile to perform stresses on a different geometry such as an annulus which resembles a wellbore. This could also be extended to model fracture propagation. Dynamic control over loading and unloading the object given yielding conditions would also provide a better visualization of permanent deformations.


### References

  * Bangerth, W., & Heister, T. (n.d.). The deal.IIFinite Element Library. Retrieved July 28, 2016, from http://www.dealii.org/
  * Hughes, Thomas J. R.The Finite Element Method: Linear Static and Dynamic Finite Element Analysis. Englewood Cliffs, NJ: Prentice-Hall, 1987. Print.
  * Kolukula, Siva. MatlabMesh Postprocessing. Net.
  * Simo, J. C., and Thomas J. R. Hughes.Computational Inelasticity. New York: Springer, 1998. Print.
  * Zienkiewicz, OlgierdC., and Robert L. Taylor.The Basis. Oxford: Butterworth-Heinemann, 2000. Print.

### Acknowledgements
I.C.E.S (Institute of Computational Engineering) provided this opportunity, stipend, and office for me to conduct my research.

### License
[Deal.ii](https://github.com/dealii/dealii/blob/master/LICENSE.md "Deal.ii Copyright and License")












