# Top2DFVT
Top2DFVT is an open source Matlab implementation for topology optimization of two-dimensional structures based on the finite-volume theory for linear elastic materials.

# Syntax

* **Top2DFVT(*L*, *H*, *nx*, *ny*, *volfrac*, *penal*, *frad*, *ft*)** Performs topology optimization algorithm on the reference domain defined by the horizontal length, *L*, and vertical lenght, *H*, using a discretization of *nx* horizontal subvolumes, and *ny* vertical subvolumes, where *volfrac* is the prescribed volume fraction, *penal* is the penalty factor, *frad* is the filter radius, *ft* is the variable that specifies the adopted filter, sensitivity (*ft = 1*), density filter (*ft = 2*), or no filter (*ft = 0*), and other default parameters.
* **Top2DFVT(*L*, *H*, *nx*, *ny*, *volfrac*, *penal*, *frad*, *ft*, *varargin*)** Performs topology optimization using *L*, *H*, *nx*, *ny*, *volfrac*, *penal*, *neig*, *ft*, and the *fsparse* routine when *varargin = 'fast'*.

 ### Table 1: Inputs parameters' declaration
--------------------------------------------------------------------------------------------------------------------
 |    Name-value     |                      Description          |    Value                          |
 |-------------------|---------------------------------------------------------|-----------------------------------|
 |  P |   applied concentrated load   |           -1                        |
 |  E0     |   Young modulus of solid material                |      1     |
 |  Emin      |    soft (void) material stiffness                         |  1e-9                       |
 |  nu  |  Poisson ratio |          0.3               |
 | model      |  material interpolation method  | 'SIMP'|
 |                   |                                                         | 'RAMP'|
 |  eta     |   damping factor | 1/2|
 |  move             |  move-limit parameter      | 0.2|
 |  F      |  global force vector (resultant forces on the subvolumes' faces)          | |
 |  supp      |  fixed degrees of freedom    | |
 
 ## Examples:

The following examples show how **Top2DFVT** can be used. Some input values are presented in Table 1 for the data initialization parameters.

  *  Example 1:
    Optimize a cantilever beam using default values of optional inputs, as presented in Table 1, and the sensitivity filter taking the adjacent neighborhood

     - **Top2DFVT(100,50,202,101,0.4,1:0.5:4,0.71,1)** or **Top2DFVT(100,50,202,101,0.4,1:0.5:4,0.71,1,'fast')**
     
     where the supporting conditions are set up as
     
     **supp = unique(dof(1:nx:end-nx+1,7:8))**
     
     and the global force vector is
     
     **F = sparse(dof(nx\*(ny+1)/2,4)',1,P,ndof,1)**
     
     - For the density filter approach: **Top2DFVT(100,50,202,101,0.4,1:0.5:4,0.71,2)** or **Top2DFVT(100,50,202,101,0.4,1:0.5:4,0.71,2,'fast')**
     
     - For the no filter approach: **Top2DFVT(100,50,202,101,0.4,1:0.5:4,[],0)** or **Top2DFVT(100,50,202,101,0.4,1:0.5:4,[],0,'fast')**

Obs: For the no filter approach employing the SIMP material interpolation method, the *eta* parameter is adjusted to 1/2.6 to avoid the oscillatory phenomenon.

  *  Example 2:
     Optimize a half-MBB beam using default values of optional inputs and the sensitivity filter, taking the two neighborhoods around the subvolume

     - **Top2DFVT(150,50,240,80,0.5,1:0.5:4,1.42,1)** or **Top2DFVT(150,50,240,80,0.5,1:0.5:4,1.42,1,'fast')**
     
     where the supporting conditions are set up as
     
     **supp = unique([dof(1:nx:end-nx+1,7);dof(nx,2)])**
     
     and the global force vector is
     
     **F = sparse(dof(nx\*ny-nx+1,6)',1,P,ndof,1)**
     
     - For the density filter approach: **Top2DFVT(150,50,240,80,0.5,1:0.5:4,1.42,2)** or **Top2DFVT(150,50,240,80,0.5,1:0.5:4,1.42,2,'fast')**
     
     - For the no filter approach: **Top2DFVT(150,50,240,80,0.5,1:0.5:4,[],0)** or **Top2DFVT(150,50,240,80,0.5,1:0.5:4,[],0,'fast')**

Obs: For the SIMP approaches, the *eta* parameter is adjusted to 1/2.6 to avoid the oscillatory phenomenon.

## Supporting Open-Source Codes
**Top2DFVT** utilizes other open-source codes, such as [fsparse](https://github.com/stefanengblom/stenglib.git) routine by
[Engblom and Lukarski](https://doi.org/10.1016/j.parco.2016.04.001), for fast convergence when medium and large-scale problems are analyzed.

