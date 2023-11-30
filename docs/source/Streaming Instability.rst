.. _Streaming_Instability:

Streaming Instability
===========

The streaming instability is an efficient mechanism for forming planetesimals out of small (mm-cm) sized dust grains, and has been used to explain how planetesimal growth can overcome bouncing/fragmentation barriers as well as the radial drift barrier.

.. figure:: _static/barriers.png
    :align: center
    :class: with-shadow with-border
    :width: 1200px

    Figure 1: The barriers to planetesimal growth which the streaming instability can overcome.

The strong clumping that arises from the streaming instability is caused by the differential motion in the azimuthal direction between the dust and the gas, which is moving at keplerian speeds while the gas is supported by pressure and in turn moves slower. The figure below demonstrates this dynamic, in which the headwind felt the dust results in drag, reducing its angular momentum and causing it to drift inward and drag the gas with it (step 1). As the gas is deflected by the Coriolis force (step 2), it accelerates the dust which reduces the dust drift effciency as it begins orbiting at sup-Keplerian velocities (step 3). This has the effect of deflecting the dust outwards and in turn yields an dust clump overdensity which can quickly grow to reach critical Roche densities (step 4). When this happens, the gravitational attraction of the pebble clump is enough to overcome the tidal shear of the disk and pebble accretion takes place.


.. figure:: _static/streaming_instability_schematic.png
    :align: center
    :class: with-shadow with-border
    :width: 1200px

    Figure 2: Illustration of the dust-gas aerodynamcis that give rise to the streaming instability. In a cylindrical shear flow with a radially decreasing pressure gradient, gas and dust exhibit differential rotational speeds, orbiting at sub-Keplerian and Keplerian speeds, respectively. The interaction between gas and dust, particularly the backreaction of dust on gas, influences the angular momentum of solids, leading to the concentration of particles and the potential formation of self-gravitating clumps that eventually give rise to planetesimals.



Simulations
===========

For this work we analyzed streaming instability simulations, both with and without self-grabity. For our streaming instability study without self-gravity, we use archival data from a single-species shearing box simulation conducted and published by Yang & Johansen 2014 using the Pencil Code, a high-order non- conservative finite-difference code for astrophysics fluid dynamics. The simulation ran for a duration of 100 orbital periods and was configured with 17 million superparticles and a spatial resolution of 256 grid cells in each dimension, where :math:`L_x` = :math:`L_y` = :math:`L_z` = 1.6H. The simulation was conducted using a Stokes number of St = 0.314 with a pressure gradient parameter of :math:`\Pi` = 0.05 and an initial solid-to-gas ratio of Z = 0.02


.. only:: html

   .. figure:: _static/sim_without_sg.gif

      Figure 3: Streaming instability simulation without self-gravity.

.. figure:: _static/si_simulation_no_sg.png
    :align: center
    :class: with-shadow with-border
    :width: 1200px

    Figure 4: Azimuthally averaged dust column density (left) and maximum particle density (right) as the simulation progresses in time.

Disk Model 
===========


.. figure:: _static/NewDisk_Model_.png
    :align: center
    :class: with-shadow with-border
    :width: 1200px

    Figure 5: Protoplanetary disk model.

