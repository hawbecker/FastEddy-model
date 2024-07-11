==================================================
Transport and dispersion over idealized topography
==================================================

Background
----------

This is an example of transport and dispersion of a scalar in the presence of idealized topography. The idealized terrain is based on the Witch of Agnesi, and a passive tracer is released on the lee side of the hill.

Input parameters
----------------

* Number of grid points: :math:`[N_x,N_y,N_z]=[504,498,90]`
* Isotropic grid spacings in the horizontal directions: :math:`[\Delta x,\Delta y]=[4,4]` m, the minimum vertical grid at the surface is :math:`\Delta z=3.6` m and stretched with verticalDeformFactor :math:`=0.26`
* Domain size: :math:`[2.16 \times 1.99 \times 1.44]` km
* Model time step: :math:`0.01` s
* Advection scheme: 5th-order upwind
* Time scheme: 3rd-order Runge Kutta
* Geostrophic wind: :math:`[U_g,V_g]=[8,0]` m/s
* Latitude: :math:`54.0^{\circ}` N
* Surface potential temperature: :math:`300` K
* Potential temperature profile:

.. math::
  \partial{\theta}/\partial z =
    \begin{cases}
      0 & \text{if $z$ $\le$ 500 m}\\
      0.003 & \text{if $z$ > 500 m}
    \end{cases} 

* Rayleigh damping layer: uppermost :math:`600` m of the domain
* Initial perturbations: :math:`\pm 0.25` K 
* Depth of perturbations: :math:`375` m
* Top boundary condition: free slip
* Lateral boundary conditions: periodic
* Time period: :math:`1` h

Execute FastEddy
----------------

Note that this example requires customization of the initial condition file. The first step consists in runnig for 1 timestep to create the *FE_DISPERSION.0* file. After that, exectue the Jupyer notebook provided in **/tutorial/notebooks/Dispersion_Prep1.ipynb** to create the topography file *./Topography_504x498.dat*. Run again for 1 timestep but pointing to the generated binary file containing the terrain information in the parameter file :code:`topoPos`. This will create a new initial file *FE_DISPERSION.0* that contains the toography and with the terrain-following vertical grid correspondingly adjusted. Then, exectute the Jupyer notebook provided in **/tutorial/notebooks/Dispersion_Prep2.ipynb** to modify the surface roughness distribution over some portion of the hill. Execute the Jupyer notebook provided in **/tutorial/notebooks/Dispersion_Prep3.ipynb** to create the input file required for the source specification. Now, point to the modified initial condition file (:code:`inPath`,:code:`inFile`) and the generated passive tracer source file (:code:`rcAuxScFile`, release starts :math:`45` min into the simulation) and run FastEddy for the :math:`1` h of the simulation by changing :code:`frqOutput`, :code:`Nt`, and :code:`NtBatch` back to their original values.

See :ref:`run_fasteddy` for instructions on how to build and run FastEddy on NSF NCAR's High Performance Computing machines.
