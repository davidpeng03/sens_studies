.. _PythiaInterface:

Pythia Simulation Interface
===========================

The Pythia Simulation Interface provides a structured way to configure and simulate particles, particularly heavy neutral leptons (HNL) and Dark Photons. This interface allows users to set up simulations via custom configurations and utilize pre-defined or new particle configuration classes.

The main functionalities of this interface include:

- **Configurable particle types**: Users can register and manage various particle configurations.
- **Simulation execution**: Based on these configurations, the interface sets up and executes Pythia simulations.
- **Argument parsing**: The interface supports flexible configuration input via command-line arguments using `argparse`.

Classes Overview
----------------

ParticleConfigFactory
---------------------

.. autoclass:: ParticleConfigFactory
   :members:
   :undoc-members:
   :show-inheritance:

This factory class manages and retrieves configurations for different particle types.

### Example Usage:

.. code-block:: python

    from Pythia_conf import ParticleConfigFactory

    # Retrieve HNL configuration
    hnl_config = ParticleConfigFactory.get_particle_config("HNL", params)

PythiaSimulation
----------------

.. autoclass:: PythiaSimulation
   :members:
   :undoc-members:
   :show-inheritance:

The `PythiaSimulation` class sets up and runs a Pythia simulation based on particle-specific configurations.

### Example Usage:

.. code-block:: python

    simulation = PythiaSimulation(hnl_config)
    simulation.setup_simulation()

Command-Line Arguments
-----------------------

This interface provides several command-line arguments for flexible configuration. These options allow users to specify the particle type, mass, couplings, and other settings related to the simulation.

Here is a breakdown of the available arguments:

.. code-block:: bash

    usage: Pythia Simulation for HNL Particles
           [-h] [--model MODEL] [--particle PARTICLE] [--mass MASS]
           [--coupling COUPLING [COUPLING ...]] [--process PROCESS]
           [--may_decay MAY_DECAY] [--epsilon EPSILON]
           [--MesonMother MESONMOTHER]

    optional arguments:
      -h, --help            Show this help message and exit
      --model MODEL         Particle model (default: "HNL")
      --particle PARTICLE   Particle name (default: "N1")
      --mass MASS           Mass of the HNL particle (default: 1.0)
      --coupling COUPLING [COUPLING ...]
                            Three couplings for the HNL particle (default: 
                            [0.447e-9, 7.15e-9, 1.88e-9])
      --process PROCESS     Process selection for the simulation (default: "c")
      --may_decay MAY_DECAY
                            True or False, whether we consider particle decays 
                            (default: False)
      --epsilon EPSILON     Epsilon mixing value for DarkPhoton (default: 8e-08)
      --MesonMother MESONMOTHER
                            Choose DarkPhoton production meson source (default: True)

### Example Usage:

.. code-block:: bash

    python pythia_simulation.py --model HNL --particle N1 --mass 1.0 --coupling 0.447e-9 7.15e-9 1.88e-9 --process c --may_decay False

Particle Configuration Registration
------------------------------------

Particle configurations can be registered dynamically using the `register_config` decorator. This allows users to easily extend the supported particles by creating new configuration classes.

### Example:

.. code-block:: python

    from Pythia_conf import register_config

    @register_config("NewParticle")
    class NewParticleConfig:
        # Custom configuration logic here
        pass

Automatic Documentation with Sphinx
-----------------------------------

Automatic documentation : 

.. automodule:: Pythia_conf
   :members:
   :undoc-members:
   :show-inheritance:

