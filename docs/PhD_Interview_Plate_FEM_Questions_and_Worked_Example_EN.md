---
title: "PhD Interview Preparation: Plate Finite Elements"
subtitle: "28 potential questions, one general-purpose worked example, and detailed mathematical verification"
lang: en
date: 2026-09-10
---

# How to use this guide

This guide follows the mathematical core of **2D-Plate-Project**: small-deformation, linear-elastic, static bending of homogeneous, isotropic plates of constant thickness. The questions are preparation prompts derived from the mathematics, not a record of questions asked by a particular university or supervisor.

Start with a common sequence that also applies to beams, continuum elements, shells, and nonlinear finite element methods:

**Physical assumptions → kinematics → constitutive relations → energy and weak form → discretization and integration → assembly and constraints → solution → internal forces and stresses → verification.**

Use Part I for answers lasting approximately 30–60 seconds, followed by a whiteboard explanation. Part II takes one problem through the entire sequence: first an analytical solution within plate theory, then a finite element approximation for comparison. The analytical solution is exact for the stated plate model and boundary conditions; it is not an exact three-dimensional elasticity solution.

## Notation and sign conventions

| Symbol | Meaning | SI unit |
|---|---|---|
| $x,y,z$ | Midsurface coordinates and thickness coordinate | m |
| $a,b,t$ | Plate length, width, and thickness | m |
| $w$ | Midsurface deflection in the positive $z$ direction | m |
| $\theta_x,\theta_y$ | Slope-type variables describing director rotation | rad; dimensionless |
| $\boldsymbol\kappa$ | Engineering curvature vector; its third component contains two mixed derivatives | m$^{-1}$ |
| $\boldsymbol\gamma$ | Transverse engineering shear strain | Dimensionless |
| $D,\mathbf D_b$ | Bending rigidity coefficient and matrix | N·m |
| $S,\mathbf D_s$ | Shear rigidity coefficient and matrix | N/m |
| $\mathbf M$ | Bending moment resultants work-conjugate to curvature | N·m/m, equivalent to N |
| $\mathbf Q$ | Transverse shear force resultants | N/m |
| $q$ | Transverse surface load, positive in the same direction as $w$ | N/m$^2$ |

The displacement assumption is $\mathbf u=(-z\theta_x,-z\theta_y,w)^T$, consistent with the project. Consequently, $\theta_x$ must not automatically be identified with another program's right-hand rotation about the global $x$ axis. Transform the kinematic convention before comparing results between programs.

# Part I. Potential interview questions

## A. Physical model and kinematics

### 01 | Can you explain the project in one minute?

**Core answer.** The inputs are the midsurface geometry, thickness, elastic properties, transverse loading, and supports. The model solves for deflections and rotations, then recovers curvatures, bending moments, shear forces, and surface stresses. The central task is to discretize continuum energy into $\mathbf K\mathbf u=\mathbf f$ and assess the results through equilibrium, analytical comparisons, and convergence.

**Follow-up.** How do beams, plates, and shells differ? A beam is primarily described along an axis; a plate bends over a two-dimensional midsurface; a general shell also involves midsurface curvature and membrane action. This project does not include the membrane degrees of freedom of a complete shell formulation.

### 02 | Why can a three-dimensional problem be reduced to two dimensions?

**Core answer.** Instead of solving independently for displacement throughout the thickness, we assume its through-thickness distribution and integrate over that direction. Reissner–Mindlin theory assumes that initially straight normals remain straight, but allows them to cease being perpendicular to the deformed midsurface.

**At the whiteboard.** Write $u_x=-z\theta_x$, $u_y=-z\theta_y$, and $u_z=w$, and explain the physical meaning of each term.

**Limitation.** Near supports or concentrated loads, or when detailed through-thickness stresses are needed, a three-dimensional model or local refinement may be necessary. Refining the plate mesh alone cannot remove an error in the physical model assumptions.

### 03 | Where do curvature and shear strain come from?

**Core answer.** Differentiate the displacement field rather than starting from a memorized strain matrix:

$$
\boldsymbol\kappa=
\begin{bmatrix}\theta_{x,x}\\\theta_{y,y}\\\theta_{x,y}+\theta_{y,x}\end{bmatrix},
\qquad
\boldsymbol\gamma=
\begin{bmatrix}w_{,x}-\theta_x\\w_{,y}-\theta_y\end{bmatrix}.
$$

The in-plane strain through the thickness is $-z\boldsymbol\kappa$. Spatial variations of rotation produce bending; the difference between the deflection gradient and the director rotation produces transverse shear.

**Follow-up.** Why is the third curvature component not multiplied by $1/2$? This formulation uses engineering shear strain. If tensorial shear strain is used instead, the constitutive and energy coefficients must be changed consistently.

### 04 | How do Kirchhoff–Love and Reissner–Mindlin theories differ?

**Core answer.** Kirchhoff theory imposes $\boldsymbol\gamma=0$, so $\boldsymbol\theta=\nabla w$, and retains only bending energy. Mindlin theory treats deflection and rotation independently and retains transverse shear energy. The distinction originates in their kinematic assumptions.

**Follow-up.** Can Mindlin theory be used for thin plates? Yes, provided that the discretization handles the thin-plate limit correctly. Otherwise, shear locking can occur.

**Common mistake.** Interpreting negligible physical shear deformation as permission to remove shear stiffness arbitrarily from a finite element.

### 05 | Is a thickness-to-span ratio of $1/20$ a theoretical dividing line?

**Core answer.** It is a default selection rule in the project's high-level interface, not a universal theorem. The choice should also consider boundary conditions, load wavelength, the shear contribution, and the acceptable error in the response of interest. The characteristic length belongs to the whole plate geometry; changing the element size should not change the selected physical theory.

**Connection to the example.** The worked example has a ratio of exactly $1/20$, but explicitly selects Mindlin theory to quantify the shear correction and investigate thin-plate numerical behavior.

### 06 | When does the linear model become inadequate?

**Core answer.** A nonlinear model is needed when the strain–displacement relationship, the configuration used for equilibrium, or the material stress–strain relationship can no longer be adequately linearized. Membrane action may become important when plate deflection is no longer small relative to thickness. Near buckling, stability analysis may be needed even before displacements become large.

**Follow-up.** Does small strain imply small rotation? No. Large rotation with small strain can still require geometric nonlinearity.

## B. Constitutive relations, energy, and weak form

### 07 | Why is bending rigidity proportional to thickness cubed?

**At the whiteboard.** Substitute $\boldsymbol\varepsilon=-z\boldsymbol\kappa$ into the through-thickness energy integral:

$$
\mathbf D_b=\int_{-t/2}^{t/2}z^2\mathbf C_{ps}\,dz
=\frac{t^3}{12}\mathbf C_{ps},\qquad
D=\frac{Et^3}{12(1-\nu^2)}.
$$

**Physical meaning.** Material farther from the midsurface experiences greater strain and acts through a larger moment arm. In contrast, the shear rigidity $S=k_sGt$ is proportional to the first power of thickness.

### 08 | Why use plane-stress elasticity, and what is the shear correction factor?

**Core answer.** This reduced elastic plate model uses the in-plane plane-stress relation to construct its bending law, neglecting the leading influence of thickness-normal stress on that relation. It does not provide a complete three-dimensional stress recovery. Mindlin kinematics give constant transverse shear strain through the thickness, which does not match the actual distribution; $k_s$ corrects the shear energy.

**Follow-up.** Is $5/6$ always valid? No. It is a common energy-equivalent value for a homogeneous rectangular section and should not be transferred unconditionally to arbitrary layered materials or cross-sections.

### 09 | How do you derive the finite element equations from total potential energy?

$$
\Pi=\frac12\int_\Omega\boldsymbol\kappa^T\mathbf D_b\boldsymbol\kappa\,dA
+\frac12\int_\Omega\boldsymbol\gamma^T\mathbf D_s\boldsymbol\gamma\,dA
-\int_\Omega qw\,dA.
$$

**Core answer.** Within the displacement space satisfying the essential boundary conditions, require $\delta\Pi=0$ for every admissible virtual displacement. Substitute $\boldsymbol\kappa=\mathbf B_b\mathbf a_e$ and $\boldsymbol\gamma=\mathbf B_s\mathbf a_e$ to obtain the stiffness integrals and nodal equilibrium equations.

**Follow-up.** Why does the strain energy contain $1/2$ while the load potential does not? Linear-elastic strain energy is the integral of stress increasing linearly with strain. The potential of a prescribed dead load is $-\mathbf f^T\mathbf u$. The actual work accumulated while that load is increased linearly from zero is instead $\tfrac12\mathbf f^T\mathbf u$. These are different quantities.

### 10 | What is the difference between the strong and weak forms?

**Core answer.** The strong form requires the differential equations to hold pointwise. The weak form requires integral equilibrium against all admissible test functions. Integration by parts reduces the derivative order required of the unknown fields and reveals the natural boundary conditions.

**At the whiteboard.** With the moment tensor $\mathsf M=\left[\begin{smallmatrix}M_x&M_{xy}\\M_{xy}&M_y\end{smallmatrix}\right]$, the present sign convention gives

$$
-\nabla\!\cdot\mathbf Q=q,\qquad
\nabla\!\cdot\mathsf M+\mathbf Q=\mathbf0.
$$

**Follow-up.** Does the weak form weaken the physical requirements? For suitable function spaces and sufficiently smooth solutions, the two forms correspond. A numerical approximation enforces the weak form only within a finite test space.

### 11 | Why can Mindlin elements use $C^0$ interpolation, while Kirchhoff elements are more difficult?

**Core answer.** Mindlin energy contains only first derivatives of $w$ and $\theta$, so a standard conforming displacement discretization can use continuous $C^0$ shape functions. If rotations are eliminated in Kirchhoff theory, the energy contains second derivatives of $w$. A direct conforming displacement method then generally requires continuity of first derivatives across element boundaries.

**Follow-up.** How does DKQ address this? It constructs curvature using discrete Kirchhoff constraints while retaining corner-node degrees of freedom. It should not be described as simply taking exact second derivatives of an arbitrary bilinear Q4 deflection field.

### 12 | Is the stiffness matrix always symmetric and positive definite?

**Core answer.** In this conservative, linear-elastic discretization, stiffness is the Hessian of the energy and is therefore symmetric. Before constraints are applied, rigid-body modes generally make it only positive semidefinite. The reduced matrix becomes positive definite after rigid motions are removed, provided that no additional mechanisms remain and the discretization is stable.

**Follow-up.** What are the three rigid-body modes in plate bending? They can be represented by $w=c_0+c_1x+c_2y$, $\theta_x=c_1$, and $\theta_y=c_2$. Both curvature and transverse shear strain are zero.

## C. Discretization and numerical methods

### 13 | What do shape functions and the $\mathbf B$ matrix do?

**Core answer.** Shape functions interpolate nodal degrees of freedom into displacement fields inside an element. The $\mathbf B$ matrix maps those degrees of freedom to strains or curvatures. The former describes the field; the latter describes how it varies.

**At the whiteboard.** For node $i$, write

$$
\mathbf B_{b,i}=\begin{bmatrix}0&N_{i,x}&0\\0&0&N_{i,y}\\0&N_{i,y}&N_{i,x}\end{bmatrix},\quad
\mathbf B_{s,i}=\begin{bmatrix}N_{i,x}&-N_i&0\\N_{i,y}&0&-N_i\end{bmatrix}.
$$

### 14 | What does the Jacobian do?

**Core answer.** It connects natural and physical coordinates, transforming both derivatives and area:

$$
\nabla_xN_i=\mathbf J^{-T}\nabla_\xi N_i,\qquad dA=\det\mathbf J\,d\xi d\eta.
$$

Here, the rows of $\mathbf J$ correspond to physical coordinates and its columns to natural coordinates. With another storage convention, the transpose locations must be adjusted accordingly.

**Common mistake.** Assuming that $\det\mathbf J>0$ proves good mesh quality. It is only one necessary condition. Severe distortion can still damage accuracy and conditioning, and checking finitely many sampling points is not a general proof of geometric validity everywhere.

### 15 | Why use Gauss quadrature? Is a two-by-two rule always exact?

**Core answer.** Quadrature approximates element energy and load integrals using a small number of evaluation points. A two-point one-dimensional Gauss rule integrates polynomials of degree three or less exactly; its tensor-product extension applies this property in each coordinate.

**Follow-up.** What about distorted Q4 elements? The inverse mapping can make the integrand non-polynomial. The sinusoidal load used later is also non-polynomial, so two-by-two load quadrature introduces a small error.

### 16 | Why can consistent nodal loads not always be distributed equally?

**Core answer.** Equality of external virtual work gives $\mathbf f_e=\int\mathbf N_w^Tq\,dA$. Equal transverse forces at the four nodes occur only in special cases, such as a regular rectangular element under a uniform surface load.

**Follow-up.** Does pressure alone directly produce nodal rotation loads? Not for the Mindlin interpolation and external work term $\int qw\,dA$ used here. Distributed couples or different kinematic interpolation require their own virtual-work derivation.

### 17 | What exactly is shear locking?

**Core answer.** A thin plate requires $\boldsymbol\theta\approx\nabla w$. With ordinary low-order interpolation, these fields are not sufficiently compatible, leaving spurious shear strain. Relative to bending, the shear term is weighted increasingly strongly as $L^2/t^2$ grows. The model becomes artificially stiff and underpredicts deflection.

**At the whiteboard.** Compare $SL^2/D\sim(L/t)^2$.

**Common mistake.** Saying that the absolute shear rigidity increases as the plate becomes thinner. In fact, $S\sim t$ decreases too; what matters is its ratio to $D\sim t^3$.

### 18 | How does MITC4 reduce locking?

**Core answer.** It samples covariant shear components at specified tying points and interpolates an assumed shear field inside the element, improving the discrete compatibility of the thin-plate constraint. This project samples $\gamma_\xi$ at the upper and lower edge midpoints and $\gamma_\eta$ at the left and right edge midpoints, then transforms back to physical components.

**Follow-up.** Is this the same as reduced integration? No. MITC changes the shear field being integrated; one-point integration changes the sampling rule. A mixed-method projection interpretation of MITC4 is also described in the [official GetFEM documentation](https://getfem.org/userdoc/model_Mindlin_plate.html).

### 19 | Why not always use reduced integration?

**Core answer.** Underintegration can allow non-rigid deformation patterns to have zero strain energy at every sampling point, producing additional zero-energy modes. Reducing locking is not enough: stability, eigenvalues, and sensitivity to distortion also need checking.

**Common mistake.** Declaring a method correct merely because its deflection becomes larger and closer to an analytical value.

## D. Assembly, boundaries, recovery, and research judgment

### 20 | How are elements assembled into a global system?

**Core answer.** Adjacent elements share the degrees of freedom of a global node. Their virtual-work contributions are added using the same global numbering:

$$
\mathbf K=\sum_e\mathbf A_e^T\mathbf K_e\mathbf A_e,\qquad
\mathbf f=\sum_e\mathbf A_e^T\mathbf f_e.
$$

**Follow-up.** Why addition? Displacements are compatible at a shared node, while internal force contributions from all adjoining elements together balance the external nodal force.

### 21 | How do essential, natural, hard simply supported, and soft simply supported boundaries differ?

**Core answer.** Essential conditions directly restrict the unknown fields and their virtual variations. Natural conditions specify forces or moments conjugate to those variations after integration by parts. The Mindlin hard simply supported condition used here is $w=0,\theta_t=0$. The unrestricted normal rotation $\theta_n$ has the natural condition $M_{nn}=0$.

Fixing only $w$, leaving both rotations free, and applying no edge moments gives a different condition: $\mathsf M\mathbf n=0$. It must not automatically be compared with the hard simply supported sinusoidal solution. A clamped boundary instead imposes $w=\theta_x=\theta_y=0$.

### 22 | How are nonzero prescribed displacements handled, and how are reactions obtained?

$$
\mathbf K_{ff}\mathbf u_f=\mathbf f_f-\mathbf K_{fc}\mathbf u_c,
\qquad \mathbf r=\mathbf K\mathbf u-\mathbf f.
$$

**Core answer.** Prescribed displacements contribute to the right-hand side through the coupling stiffness. It is not sufficient to delete the corresponding rows and columns while ignoring that contribution. After reconstructing the full displacement vector, residuals at constrained degrees of freedom give the support reactions; residuals at free degrees of freedom should be near zero.

**Follow-up.** What about an inclined edge? Transform the rotations into local normal and tangential coordinates before imposing the relevant constraints.

### 23 | How do displacements become moments and stresses, and how do you avoid sign errors?

**Core answer.** Follow $\mathbf a_e\rightarrow\boldsymbol\kappa,\boldsymbol\gamma\rightarrow\mathbf M,\mathbf Q\rightarrow\boldsymbol\sigma(z)$. Here $\mathbf M=\mathbf D_b\boldsymbol\kappa$, while $\boldsymbol\varepsilon=-z\boldsymbol\kappa$ gives

$$
\boldsymbol\sigma(z)=-z\mathbf C_{ps}\boldsymbol\kappa
=-\frac{12z}{t^3}\mathbf M.
$$

Thus, the stress at $z=t/2$ is $-6\mathbf M/t^2$. Under this convention, $\mathbf M=-\int z\boldsymbol\sigma\,dz$. Do not mix these signs with a textbook using the opposite definition of the moment resultant.

### 24 | Does a smoother stress contour mean greater accuracy?

**Core answer.** No. Stress depends on displacement derivatives and is usually more sensitive. Extrapolation and nodal averaging produce smoother displays but may conceal differences between elements. Retain the original Gauss-point values, identify the recovery location, and check convergence of the quantity of interest.

**Follow-up.** Why can shear force be zero at the point of maximum deflection? Symmetry can make the center shear force zero. The locations of maximum moment and deflection need not coincide with maximum shear force.

### 25 | Do a small residual and balanced reactions prove an accurate answer?

**Core answer.** They show that the assembled discrete equations have been solved consistently and satisfy the corresponding global equilibrium. Incorrect material properties, incorrect boundaries, or locking elements can still pass these checks.

**Connection to the example.** The fully integrated $4\times4$ mesh later has a relative residual of approximately $10^{-14}$, but underpredicts center deflection by about $78\%$. Analytical benchmarks and mesh convergence remain necessary.

### 26 | How would you build a credible verification procedure?

**Core answer.** Check units and geometry; shape-function consistency; rigid-body modes; patch tests for reproducible fields; matrix symmetry; residuals and global equilibrium; energy; analytical benchmarks; mesh and thickness studies; and sensitivity to distortion.

**Follow-up.** What is the difference between verification and validation? Verification asks whether the chosen equations are solved correctly. Validation asks whether the chosen physical model describes reality and needs experimental or other physical evidence. The worked example provides mathematical and numerical verification.

### 27 | How would you design a convincing convergence study?

**Core answer.** Keep the physical problem, boundaries, and error definition fixed while varying mesh size. Compare the same response at the same location across multiple meshes. Before claiming a particular convergence order, establish that the results are in the asymptotic regime.

**Follow-up.** Are mesh convergence and freedom from locking the same? No. Mesh convergence refines the mesh at fixed thickness. A locking study also varies thickness on a fixed or controlled mesh and checks whether errors deteriorate as $t/L\to0$. A single fixed-thickness convergence table cannot establish uniform locking-free behavior.

### 28 | How could this project lead to a PhD research question?

**Core answer.** First identify the thicknesses, distortions, materials, or coupled conditions for which the current method loses accuracy, then propose a testable improvement. For example: can a particular mixed discretization maintain stability and a uniform error bound when both thickness-to-span ratio and mesh distortion vary?

**Follow-up.** What mathematics must be added for nonlinear plate and shell analysis? Updated kinematics, an internal-force residual, a consistent tangent, and iterative equilibrium solution are needed. Material nonlinearity also requires evolving material states; buckling requires an appropriate geometric stiffness and stability problem. These are possible extensions, not capabilities already present in this linear plate core.

# Part II. One worked example connecting the full calculation

## 1. Problem statement

Consider a homogeneous, isotropic, constant-thickness rectangular plate with midsurface $0\le x\le a$, $0\le y\le b$. All four edges are **hard simply supported in the Mindlin sense**: $w=0$ and tangential rotation $\theta_t=0$. Normal rotation is free, with no applied normal bending moment, so $M_{nn}=0$.

Apply the following sinusoidal surface load, positive in the same direction as $w$:

$$q(x,y)=q_0\sin\frac{\pi x}{a}\sin\frac{\pi y}{b}.$$

Use these numerical values:

$$
a=b=2.0\ \mathrm m,\quad t=0.10\ \mathrm m,\quad
E=30\times10^9\ \mathrm{Pa},\quad \nu=0.30,
$$

$$k_s=\frac56,\qquad q_0=10\,000\ \mathrm{N/m^2}.$$

Use small-deformation, linear-elastic Reissner–Mindlin theory. Neglect membrane action and any loading additional to the prescribed $q$; detailed three-dimensional thickness-normal stress is outside the model. The value $q_0$ is the peak load intensity, not its spatial average.

**Tasks:**

1. Derive curvature, shear strain, constitutive relations, and total potential energy from the displacement assumption.
2. Choose three amplitudes satisfying the boundary conditions, derive a three-equation linear system, and obtain its analytical solution.
3. Calculate center deflection, rotation amplitudes, center moments, and normal stresses on the upper and lower surfaces.
4. Find the shear force at an edge midpoint and check global force equilibrium, energy, and the thin-plate limit.
5. Explain how Q4 finite elements solve the same problem, including element calculations, assembly, and numerical convergence.

**Why is this example broadly useful?** Dimensions, thickness, material properties, and load amplitude can all be varied. Under linear conditions, more general loads can be expanded into sinusoidal terms and their responses superposed. The energy and discretization reasoning transfers to other structural problems; this particular analytical formula remains restricted to the rectangular geometry, stated supports, and constant-coefficient plate model.

## 2. From displacement to strain: describe the deformation first

Start with

$$u_x=-z\theta_x(x,y),\quad u_y=-z\theta_y(x,y),\quad u_z=w(x,y)$$

Differentiate term by term:

$$\varepsilon_x=\frac{\partial u_x}{\partial x}=-z\theta_{x,x},\qquad
\varepsilon_y=-z\theta_{y,y},$$

$$\gamma_{xy}=\frac{\partial u_x}{\partial y}+\frac{\partial u_y}{\partial x}
=-z(\theta_{x,y}+\theta_{y,x}),$$

$$\gamma_{xz}=\frac{\partial u_x}{\partial z}+\frac{\partial w}{\partial x}
=w_{,x}-\theta_x,\qquad
\gamma_{yz}=w_{,y}-\theta_y.$$

The in-plane strain vector is therefore $-z\boldsymbol\kappa$, and the transverse shear strain is $\boldsymbol\gamma=\nabla w-\boldsymbol\theta$. In this pure plate-bending model, the midsurface has no in-plane extension degrees of freedom. Bending strains on opposite sides of it have opposite signs.

## 3. From strain to rigidity: how the material resists deformation

The plane-stress constitutive matrix is

$$
\mathbf C_{ps}=\frac{E}{1-\nu^2}
\begin{bmatrix}1&\nu&0\\\nu&1&0\\0&0&(1-\nu)/2\end{bmatrix}.
$$

Integrate through the thickness, using

$$\int_{-t/2}^{t/2}z^2\,dz=\left[\frac{z^3}{3}\right]_{-t/2}^{t/2}=\frac{t^3}{12},$$

to obtain

$$\mathbf D_b=D\begin{bmatrix}1&\nu&0\\\nu&1&0\\0&0&(1-\nu)/2\end{bmatrix},\quad
D=\frac{Et^3}{12(1-\nu^2)},\quad \mathbf D_s=S\mathbf I,\quad S=k_sGt.$$

Numerically,

$$G=\frac{30\times10^9}{2(1+0.3)}=1.153846154\times10^{10}\ \mathrm{Pa},$$

$$D=\frac{30\times10^9(0.1)^3}{12(1-0.3^2)}
=2.747252747\times10^6\ \mathrm{N\,m},$$

$$S=\frac56(1.153846154\times10^{10})(0.1)
=9.615384615\times10^8\ \mathrm{N/m}.$$

The work-conjugate resultants are $\mathbf M=\mathbf D_b\boldsymbol\kappa$ and $\mathbf Q=S\boldsymbol\gamma$. A bending moment resultant is a moment per unit length. Even when its numerical unit is written as N, it should be understood as N·m/m.

## 4. Weak form: where the equilibrium equations come from

Combine the strain energy and load potential:

$$\Pi=U_b+U_s-\int_\Omega qw\,dA.$$

Requiring $\delta\Pi=0$ gives

$$
\int_\Omega\delta\boldsymbol\kappa^T\mathbf M\,dA
+\int_\Omega(\nabla\delta w-\delta\boldsymbol\theta)^T\mathbf Q\,dA
=\int_\Omega q\delta w\,dA.
$$

Integrate by parts the terms containing $\nabla\delta w$ and derivatives of the virtual rotations. Their interior coefficients are

$$-\nabla\cdot\mathbf Q-q=0,\qquad-\nabla\cdot\mathsf M-\mathbf Q=\mathbf0.$$

The boundary contribution is

$$\int_{\partial\Omega}\left[(\mathbf Q\cdot\mathbf n)\delta w
+(\mathsf M\mathbf n)\cdot\delta\boldsymbol\theta\right]ds.$$

Hard simple support imposes $\delta w=\delta\theta_t=0$, but leaves $\delta\theta_n$ arbitrary. With no applied normal bending moment, the natural condition is consequently $M_{nn}=0$. Since tangential rotation is constrained, a corresponding twisting-moment reaction may exist; do not additionally impose $M_{nt}=0$.

## 5. Choose trial functions: reduce the continuum to three unknowns

Define

$$\alpha=\frac\pi a,\quad\beta=\frac\pi b,\quad k^2=\alpha^2+\beta^2.$$

Introduce three independent amplitudes $W,A,B$:

$$w=W\sin\alpha x\sin\beta y,$$

$$\theta_x=A\cos\alpha x\sin\beta y,\qquad
\theta_y=B\sin\alpha x\cos\beta y.$$

The amplitude $W$ has units of m; $A$ and $B$ are dimensionless rotation amplitudes.

**Check each edge.** At $x=0,a$, both $w$ and $\theta_y$ vanish. At $y=0,b$, both $w$ and $\theta_x$ vanish. These are exactly the required zero tangential rotations. The expressions for $M_x$ and $M_y$ derived below both contain $\sin\alpha x\sin\beta y$, so the normal bending moment also vanishes on the corresponding edges.

The curvature components are

$$\kappa_x=-\alpha A\sin\alpha x\sin\beta y,\qquad
\kappa_y=-\beta B\sin\alpha x\sin\beta y,$$

$$\kappa_{xy}=(\beta A+\alpha B)\cos\alpha x\cos\beta y.$$

The shear strains are

$$\gamma_{xz}=(\alpha W-A)\cos\alpha x\sin\beta y,$$

$$\gamma_{yz}=(\beta W-B)\sin\alpha x\cos\beta y.$$

We have not set $A=\alpha W$: doing so at this stage would prematurely remove the Mindlin shear deformation.

## 6. Evaluate the area integrals: obtain the three-amplitude energy

Using $\sin^2s=(1-\cos2s)/2$ gives

$$\int_0^a\sin^2\alpha x\,dx=\frac a2,\quad
\int_0^a\cos^2\alpha x\,dx=\frac a2.$$

The same holds in the $y$ direction. The relevant squared products therefore have area integral $H=ab/4$. With $c=(1-\nu)/2$,

$$U_b=\frac{HD}{2}\left[\alpha^2A^2+\beta^2B^2
+2\nu\alpha\beta AB+c(\beta A+\alpha B)^2\right],$$

$$U_s=\frac{HS}{2}\left[(\alpha W-A)^2+(\beta W-B)^2\right],$$

$$\int_\Omega qw\,dA=Hq_0W.$$

Physically, the first expression penalizes spatial changes in rotation, the second penalizes mismatch between the deflection gradient and rotation, and the last describes the applied load acting on the chosen deflection mode through its generalized force.

## 7. Differentiate with respect to the amplitudes: write every equation

Set $\partial\Pi/\partial W=0$ and cancel the nonzero factor $H$:

$$S\left[k^2W-\alpha A-\beta B\right]=q_0.$$

Set $\partial\Pi/\partial A=0$:

$$D\left[(\alpha^2+c\beta^2)A+(\nu+c)\alpha\beta B\right]+S(A-\alpha W)=0.$$

Set $\partial\Pi/\partial B=0$:

$$D\left[(\nu+c)\alpha\beta A+(\beta^2+c\alpha^2)B\right]+S(B-\beta W)=0.$$

In matrix form,

$$
\begin{bmatrix}
Sk^2&-S\alpha&-S\beta\\
-S\alpha&S+D(\alpha^2+c\beta^2)&D(\nu+c)\alpha\beta\\
-S\beta&D(\nu+c)\alpha\beta&S+D(\beta^2+c\alpha^2)
\end{bmatrix}
\begin{bmatrix}W\\A\\B\end{bmatrix}
=\begin{bmatrix}q_0\\0\\0\end{bmatrix}.
$$

This is a Ritz system using three global trigonometric trial functions. Like finite elements, it follows from stationary potential energy. Finite elements instead use local, piecewise shape functions, allowing more general geometries and loads.

## 8. Solve the equations: separate bending and shear deflections

Set $A=\alpha C$ and $B=\beta C$, where $C$ has units of m. Check this form by substituting it into both rotation equations; do not assume that $C=W$.

Since $\nu+2c=1$, the bending term in the first rotation equation becomes

$$D\alpha\left[\alpha^2+(c+\nu+c)\beta^2\right]C=D\alpha k^2C.$$

The second rotation equation simplifies in the same way. Both reduce to

$$S(W-C)=Dk^2C.$$

The deflection equation becomes

$$Sk^2(W-C)=q_0.$$

Substitute the former into the latter:

$$Dk^4C=q_0\quad\Rightarrow\quad C=\frac{q_0}{Dk^4}.$$

Back-substitution gives

$$\boxed{W=\frac{q_0}{Dk^4}+\frac{q_0}{Sk^2}=W_b+W_s},$$

$$\boxed{A=\alpha W_b,\qquad B=\beta W_b}.$$

Here, $W_b=C$ is the bending contribution and $W_s=W-C$ is the additional deflection due to transverse shear. Both have units of m.

**Why is this more than a Ritz approximation?** The resulting fields satisfy all boundary conditions. Furthermore,

$$\mathbf Q=\frac{q_0}{k^2}
\begin{bmatrix}\alpha\cos\alpha x\sin\beta y\\\beta\sin\alpha x\cos\beta y\end{bmatrix}$$

directly gives $-\nabla\cdot\mathbf Q=q$. Also,

$$\nabla\cdot\mathsf M=-Dk^2C\begin{bmatrix}\alpha\cos\alpha x\sin\beta y\\\beta\sin\alpha x\cos\beta y\end{bmatrix}=-\mathbf Q.$$

Thus the fields satisfy the Mindlin strong form pointwise throughout the domain. They are the analytical solution of this plate-model boundary-value problem.

## 9. Insert the numbers: deflection and rotations

$$\alpha=\beta=\frac\pi2=1.570796327\ \mathrm{m^{-1}},\quad
k^2=4.934802201\ \mathrm{m^{-2}},\quad k^4=24.35227276\ \mathrm{m^{-4}}.$$

$$W_b=\frac{10\,000}{(2.747252747\times10^6)(24.35227276)}
=1.494727016\times10^{-4}\ \mathrm m,$$

$$W_s=\frac{10\,000}{(9.615384615\times10^8)(4.934802201)}
=2.107480620\times10^{-6}\ \mathrm m.$$

At the center $(a/2,b/2)$, the product of the sine factors equals one, so

$$\boxed{w_c=W=0.151580182\ \mathrm{mm}}.$$

The rotation amplitudes are

$$\boxed{A=B=\frac\pi2(1.494727016\times10^{-4})
=2.347911707\times10^{-4}\ \mathrm{rad}}.$$

Both rotations at the center are zero because the relevant cosine factors vanish there. The values above are amplitudes of the rotation fields, not center rotations.

As a consistency check, $w_c/t=0.0015158$, and the rotations are also small. These numerical results are consistent with the small-deformation assumption used in this example.

## 10. Moments and surface stresses: keep track of the signs

At the center,

$$\kappa_x=\kappa_y=-\alpha^2C=-3.688091085\times10^{-4}\ \mathrm{m^{-1}},\qquad\kappa_{xy}=0.$$

The constitutive equations give

$$M_x=D(\kappa_x+\nu\kappa_y),\quad M_y=D(\nu\kappa_x+\kappa_y),\quad M_{xy}=Dc\kappa_{xy}.$$

Therefore,

$$\boxed{M_{x,c}=M_{y,c}=-1317.175387\ \mathrm{N\,m/m},\qquad M_{xy,c}=0}.$$

On the surface $z=+t/2=+0.05\ \mathrm m$,

$$\sigma_{x,+}=\sigma_{y,+}=-\frac6{t^2}M_{x,c}
=-\frac6{0.1^2}(-1317.175387)
=790305.232\ \mathrm{Pa}.$$

Thus,

$$\boxed{\sigma_{x,+}=\sigma_{y,+}=+0.790305\ \mathrm{MPa}},$$

$$\boxed{\sigma_{x,-}=\sigma_{y,-}=-0.790305\ \mathrm{MPa}}.$$

Tension is positive. Here, “upper” and “lower” strictly mean $z=+t/2$ and $z=-t/2$. Reversing the load direction reverses the corresponding displacement, moment, and stress signs. The stated signs follow directly from $\boldsymbol\sigma=-z\mathbf C_{ps}\boldsymbol\kappa$.

## 11. Shear forces and reactions: check global equilibrium

Using $\mathbf Q=S\boldsymbol\gamma$, at the left-edge midpoint $(0,b/2)$,

$$Q_x=\frac{q_0\alpha}{k^2}=3183.098862\ \mathrm{N/m},\qquad Q_y=0.$$

The outward normal on the left edge is $\mathbf n=(-1,0)$. The transverse boundary action on the plate, $\mathbf Q\cdot\mathbf n$, is therefore negative and opposes the positive surface load.

The total applied load is

$$F=\int_0^a\int_0^bq_0\sin\alpha x\sin\beta y\,dy\,dx
=q_0\frac2\alpha\frac2\beta
=\frac{4q_0ab}{\pi^2}=16211.389383\ \mathrm N.$$

Adding the reactions from all four edges gives

$$R=-\frac{4q_0\alpha}{\beta k^2}-\frac{4q_0\beta}{\alpha k^2}
=-\frac{4q_0(\alpha^2+\beta^2)}{\alpha\beta k^2}
=-\frac{4q_0}{\alpha\beta}=-F.$$

For the symmetric square plate, each edge carries a resultant reaction of $-4052.847346\ \mathrm N$. Multiplying the peak intensity $q_0$ by the full plate area would not give the correct total load for this problem.

To recover transverse shear stress, the project uses the parabolic assumption

$$\tau_{xz}(z)=\frac{3Q_x}{2t}\left(1-\frac{4z^2}{t^2}\right).$$

At the left-edge midpoint, its value on the thickness midplane is $0.0477465\ \mathrm{MPa}$, and it vanishes at both surfaces. This is an assumed distribution reconstructed from a shear resultant, not an exact description of the local three-dimensional stress at the support.

## 12. Energy and limiting behavior: two further checks

For $A=\alpha C$ and $B=\beta C$, the bracket in the bending energy reduces to $k^4C^2$, giving

$$U_b=\frac H2Dk^4C^2=0.747363508\ \mathrm J,$$

$$U_s=\frac H2Sk^2W_s^2=0.010537403\ \mathrm J.$$

The total strain energy is

$$U=0.757900911\ \mathrm J=\frac12Hq_0W.$$

Hence $2U=\int qw\,dA=1.515801822\ \mathrm J$. If the load is increased linearly from zero, the accumulated external work is $U$.

The ratio of the shear and bending deflection contributions is

$$\frac{W_s}{W_b}=\frac{Dk^2}{S}
=\frac{t^2k^2}{6k_s(1-\nu)}=0.014099435.$$

The Mindlin deflection exceeds the Kirchhoff deflection by approximately $1.410\%$. If Mindlin deflection is used as the denominator, neglecting shear underpredicts it by approximately $1.390\%$. Always specify the denominator of a percentage error.

At fixed in-plane geometry, this ratio tends to zero as $t\to0$, recovering the Kirchhoff limit. **This is a relative limit.** At fixed load, taking $t\to0$ increases the actual deflection and eventually violates the small-deformation assumptions. In a thickness study, scaling the load as $q_0\propto t^4$ keeps bending-dominated $w/t$ small, allowing comparison of normalized linear responses.

## 13. Solve the same problem with Q4 finite elements

### 13.1 Nodes, shape functions, and mapping

Start with a regular $4\times4$ mesh: 16 elements, 25 nodes, and 75 degrees of freedom. Each square element has side length $h=0.5\ \mathrm m$, with nodes ordered counterclockwise. The lower-left element has coordinates $(0,0)$, $(0.5,0)$, $(0.5,0.5)$, and $(0,0.5)$.

The natural-coordinate shape functions are

$$N_1=\tfrac14(1-\xi)(1-\eta),\quad N_2=\tfrac14(1+\xi)(1-\eta),$$

$$N_3=\tfrac14(1+\xi)(1+\eta),\quad N_4=\tfrac14(1-\xi)(1+\eta).$$

Interpolate $w=\sum N_iw_i$, $\theta_x=\sum N_i\theta_{xi}$, and $\theta_y=\sum N_i\theta_{yi}$, using the nodal degree-of-freedom order $[w_i,\theta_{xi},\theta_{yi}]$.

For this regular element,

$$x=\frac h2(1+\xi),\quad y=\frac h2(1+\eta),\quad
\mathbf J=\begin{bmatrix}h/2&0\\0&h/2\end{bmatrix},\quad\det\mathbf J=0.0625\ \mathrm{m^2}.$$

At its center, $N_i=1/4$, and

$$[N_{i,x}]_{i=1}^4=[-1,1,1,-1]\ \mathrm{m^{-1}},$$

$$[N_{i,y}]_{i=1}^4=[-1,-1,1,1]\ \mathrm{m^{-1}}.$$

For example, the first-node blocks at the element center are

$$\mathbf B_{b,1}=\begin{bmatrix}0&-1&0\\0&0&-1\\0&-1&-1\end{bmatrix},\quad
\mathbf B_{s,1}=\begin{bmatrix}-1&-1/4&0\\-1&0&-1/4\end{bmatrix}.$$

The derivative entries and shape-function entries have different units. Interpret each together with the deflection or rotation that it multiplies.

### 13.2 Integrating stiffness: calculate an individual coefficient

The ordinary fully integrated formulation uses

$$\mathbf K_e=\sum_{g=1}^4\left(\mathbf B_{b,g}^T\mathbf D_b\mathbf B_{b,g}
+\mathbf B_{s,g}^TS\mathbf B_{s,g}\right)\det\mathbf J_g,$$

where $(\xi_g,\eta_g)=(\pm1/\sqrt3,\pm1/\sqrt3)$ and each quadrature weight is one. MITC4 replaces $\mathbf B_s$ in the second term with the reconstructed $\widetilde{\mathbf B}_s$.

For the fully integrated $w_1$ diagonal coefficient, the bending contribution is zero because the bending strain matrix has no $w$ entries:

$$K_{w_1w_1}=S\int_{\Omega_e}(N_{1,x}^2+N_{1,y}^2)\,dA.$$

Since $N_{1,x}=-(1-\eta)/(2h)$ and $N_{1,y}=-(1-\xi)/(2h)$,

$$K_{w_1w_1}=\frac S{16}\int_{-1}^1\int_{-1}^1
\left[(1-\eta)^2+(1-\xi)^2\right]d\xi d\eta.$$

For each squared term, $\int_{-1}^1(1-s)^2ds=8/3$. The combined double integral is therefore $32/3$, and

$$K_{w_1w_1}=\frac{2S}{3}=6.410256410\times10^8\ \mathrm{N/m}.$$

The fully integrated $\theta_{x1}$ diagonal coefficient is

$$K_{\theta_{x1}\theta_{x1}}=D\left(\frac13+\frac c3\right)+S\frac{h^2}{9}
=2.794566545\times10^7\ \mathrm{N\,m}.$$

For the same coefficient, the project's MITC4 implementation gives $2.126831502\times10^7\ \mathrm{N\,m}$. The difference comes from the shear field. Neither an individual diagonal coefficient nor arbitrarily removing the shear term is sufficient to establish global accuracy.

### 13.3 The MITC4 reconstruction

First transform physical shear strain to covariant components: $\boldsymbol\gamma_{cov}=\mathbf J^T\boldsymbol\gamma$. Using the upper and lower edge midpoints for one component and the left and right midpoints for the other, construct

$$\widetilde\gamma_\xi(\xi,\eta)=\frac{1-\eta}{2}\gamma_\xi(0,-1)+\frac{1+\eta}{2}\gamma_\xi(0,1),$$

$$\widetilde\gamma_\eta(\xi,\eta)=\frac{1-\xi}{2}\gamma_\eta(-1,0)+\frac{1+\xi}{2}\gamma_\eta(1,0).$$

At the current quadrature point, transform back using $\widetilde{\boldsymbol\gamma}=\mathbf J^{-T}\widetilde{\boldsymbol\gamma}_{cov}$, from which $\widetilde{\mathbf B}_s$ follows. This retains the selected low-order shear structure and reduces incompatible discrete constraints.

### 13.4 Consistent loading: the four actual nodal force values

At the same four Gauss points, evaluate

$$f_{w_i,e}\approx\sum_{g=1}^4N_i(\xi_g,\eta_g)q(x_g,y_g)\det\mathbf J_g.$$

For the lower-left element, this gives

$$[f_{w_1},f_{w_2},f_{w_3},f_{w_4}]
=[40.001710,\ 77.918896,\ 151.777371,\ 77.918896]\ \mathrm N.$$

All eight rotational load entries are zero. Node 3, which is closest to the plate center, receives the largest contribution because the sinusoidal load increases toward the center.

### 13.5 Assembly, supports, and solution

Add the element stiffness and load contributions at shared global nodes. Using $i,j=0,\ldots,4$ as the node indices in the $x$ and $y$ directions, let $n=5j+i$. Its degree-of-freedom indices are $(3n,3n+1,3n+2)$. The first element uses nodes $(0,1,6,5)$, and its stiffness entries are accumulated at the corresponding global positions. These indices are zero-based.

There are 16 boundary nodes whose $w$ values are fixed. Constrain $\theta_y$ on the left and right edges and $\theta_x$ on the lower and upper edges, giving 20 rotational constraints. Of the 75 degrees of freedom, 36 are constrained and 39 remain free.

All prescribed values are zero, so solve $\mathbf K_{ff}\mathbf u_f=\mathbf f_f$. The center node is $n=12$, and $u_{36}$ is its deflection. MITC4 gives $0.146214651\ \mathrm{mm}$. This value results from solving the assembled system of 39 equations; it is not obtained by assigning analytical displacements to the nodes.

Recover the reactions from $\mathbf r=\mathbf K\mathbf u-\mathbf f$. To calculate stresses, evaluate $\mathbf B_b\mathbf a_e$ and the constitutive law at specified element points. Do not label the analytical center stress as a recovered finite element stress for this mesh.

## 14. Numerical verification: analytical and project results

The following results were recomputed with this project. Every mesh explicitly uses the **Mindlin branch**, hard simple supports, and the same sinusoidal loading. Both shear formulations use the project's two-by-two load quadrature. Define percentage error as $(w_h/W-1)\times100\%$.

| Elements per direction | Full integration: center deflection, mm | Error | MITC4: center deflection, mm | Error |
|---|---:|---:|---:|---:|
| 2 | 0.00988890 | −93.476% | 0.12578850 | −17.015% |
| 4 | 0.03278439 | −78.372% | 0.14621465 | −3.540% |
| 8 | 0.07948889 | −47.560% | 0.15029349 | −0.849% |
| 16 | 0.12355963 | −18.486% | 0.15126174 | −0.210% |
| Analytical | 0.15158018 | 0 | 0.15158018 | 0 |

**Interpretation.** MITC4 center deflection approaches the analytical value as the mesh is refined. Ordinary full integration is substantially too stiff at this thickness-to-span ratio. Refinement improves it, but at greater computational cost. This is evidence for regular square meshes, one thickness, and one response quantity. It does not establish identical accuracy for distorted meshes or every stress component.

Additional checks for the $16\times16$ MITC4 mesh are:

| Quantity | Recomputed value | Interpretation |
|---|---:|---|
| Discrete total load | 16211.378214 N | Two-by-two quadrature |
| Sum of transverse reactions | −16211.378214 N | Balances the discrete load |
| Analytical total load | 16211.389383 N | About 0.01117 N above the discrete value |
| Relative free-residual norm | Approximately $1.12\times10^{-12}$ | Normalized by the full load-vector norm |
| Discrete strain energy | 0.751464614 J | From $\tfrac12\mathbf u^T\mathbf K\mathbf u$ |
| Analytical strain energy | 0.757900911 J | A global check distinct from center deflection |

The energy identity $\mathbf u^T\mathbf K\mathbf u=\mathbf u^T\mathbf f$ and reaction equilibrium both passed the numerical checks. These are necessary internal consistency checks; comparison with the analytical solution remains essential.

# Part III. Follow-up questions using the same example

| Changed condition or question | How to respond |
|---|---|
| Double the load | In this linear model, deflections, rotations, moments, shear forces, and stresses double; strain energy increases fourfold. |
| Double $E$, with everything else fixed | Both $D$ and $S$ double. Deflections and rotations halve; moments and shear forces in this load-controlled analytical problem remain unchanged. |
| Double the thickness | $W_b$ becomes $1/8$ of its previous value and $W_s$ becomes $1/2$. In this problem, moments remain unchanged and surface bending stresses become $1/4$. |
| Replace the sinusoidal load with a uniform load | Expand the constant load into a double sine series. Do not insert the uniform intensity directly into the single-mode formula. |
| Change the edges to clamped | The original rotation trial functions no longer satisfy the boundary conditions. Choose admissible functions or solve again with appropriate finite element constraints. |
| Constrain deflection only at the four corners | This is not a plate simply supported along all four edges. Its result cannot be compared directly with this solution. |
| Does zero center shear mean the plate carries no load there? | No. Center curvature and bending moment are nonzero; shear force is related to the spatial variation of moment. |
| The residual is small, but the deflection is very small | Check the theory, load, units, boundaries, and locking. Changing solver tolerances alone is not a sufficient response. |
| Extend the problem to an arbitrary quadrilateral | Separation of variables generally no longer applies. Finite element geometry mapping, distortion studies, and error control become particularly important. |

For a uniform load $q_u$, orthogonality gives the double sine-series coefficients:

$$q_{mn}=\frac4{ab}\int_0^a\int_0^bq_u\sin\frac{m\pi x}{a}\sin\frac{n\pi y}{b}\,dy\,dx.$$

For odd $m,n$, $q_{mn}=16q_u/(mn\pi^2)$; otherwise the coefficient is zero. For each term, set $\alpha_m=m\pi/a$, $\beta_n=n\pi/b$, and $k_{mn}^2=\alpha_m^2+\beta_n^2$, and use

$$W_{mn}=\frac{q_{mn}}{Dk_{mn}^4}+\frac{q_{mn}}{Sk_{mn}^2}.$$

Superpose the modal fields and check convergence with respect to series truncation. This is another reason for choosing the sinusoidal example: it serves as a building block for more general loading on rectangular, simply supported plates.

## A practical rehearsal sequence

**Round 1: explain it aloud.** Without looking at the equations, explain all nine stages of the calculation. Prioritize Questions 01, 03, 04, 09, 17, 21, 25, and 28.

**Round 2: derive it at the whiteboard.** Independently derive the strains, the cubic thickness dependence, the three-amplitude energy equations, and $W=W_b+W_s$. At every step, explain what is being calculated and which assumption is being used.

**Round 3: challenge the result.** After changing the load, thickness, or supports, which conclusions still apply? What does your verification demonstrate, and what remains untested?

A possible closing explanation of your project is:

> I begin with the plate kinematics to define curvature and transverse shear strain. Elastic constitutive relations and through-thickness integration then give the bending and shear energies. After discretization, assembly, and application of boundary conditions, I solve for nodal deflections and rotations and recover internal forces and stresses. For thin plates, I pay particular attention to shear locking. Using the same simply supported plate under sinusoidal loading, I compare finite element results with an analytical solution and check reaction equilibrium, energy, and mesh convergence. This allows me to distinguish solving the discrete equations correctly from obtaining a sufficiently accurate discretization.

# Sources and scope of verification

The formulas, implementation, and numerical example were checked against the current working project on 10 September 2026. This English edition preserves those calculations. Earlier study notes were used to identify the existing revision structure; historical test counts are not presented as a new full regression-test result.

| Mathematical topic | Project reference |
|---|---|
| Displacement, curvature, and shear strain | `mindlin_plate/kinematics.py` |
| Material rigidity and stress signs | `mindlin_plate/material.py` |
| Q4 mapping, strain matrices, and MITC4 | `mindlin_plate/q4.py` |
| Stiffness and consistent load integration | `mindlin_plate/element.py` |
| Assembly, hard simple supports, and reactions | `mindlin_plate/assembly.py` |
| Automatic theory selection | `mindlin_plate/theory.py` |
| Existing revision framework | The project's mathematical overview and interview Q&A notes in `docs/` |

External theoretical cross-reference: the [official GetFEM description of the Mindlin–Reissner model](https://getfem.org/userdoc/model_Mindlin_plate.html), used to cross-check the weak-form and MITC shear-field concepts. The amplitude derivation, numerical substitution, and convergence table in this guide were independently prepared and checked locally.

No solver implementation was changed. The numerical checks cover the three-amplitude system, two shear formulations on regular meshes, load and reaction equilibrium, energy, and center deflection. They do not establish performance for every thickness, distorted mesh, or three-dimensional effect.
