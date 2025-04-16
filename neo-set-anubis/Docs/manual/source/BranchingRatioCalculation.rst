.. _BRcalculator:

Branching Ratio Calculation
===========================

The Branching Ratio (BR) calculation is a crucial part of determining the phenomenology of long-lived particles (LLP). This manual explains how to use the BR calculation interface, allowing you to compute different types of branching ratios for any model that has UFO as input.

The BR interface supports two main types of calculations for LLPs:

- **Production BR**: This is the branching ratio corresponding to the production of a given particle from a decaying mother particle.
- **Decay BR**: This refers to the branching ratio for the decay of the particle into specific final states.

Both types of branching ratios can be calculated using Python-based functions. The interface is designed to allow flexibility in the methods used for the computation, with options for pure Python calculations or input-output files to handle large data sets or precomputed results.

Command-Line Arguments
-----------------------

This interface provides several command-line arguments for flexible configuration. These options allow users to specify the particle type, mass, couplings, and other settings related to the calculation.

Here is a breakdown of the available arguments:

.. include:: BR_calculator.rst

### Example Usage:

.. code-block:: bash

    python Pipeline/BR_calculator.py --model HNL --params 1 1 1 --mass 1 --calc_type DecayTot


Setting up the BR Interface
---------------------------

To use the interface, first instantiate the `BRInterface` class:

.. code-block:: python

    testbr = BRInterface()

Once initialized, the interface requires you to set a calculation method and a model before performing any computations.

Choosing a Calculation Method
------------------------------

The interface supports two calculation strategies:

- **Python-based calculations**: This is the recommended method for flexibility and custom analysis.
- **File-based calculations**: For situations where results have already been precomputed and stored in files.

To select a calculation method, use:

.. code-block:: python

    testbr.set_calculation_method(method_name, pythoncalculation_prod_in_file=False, pythoncalculation_decay_in_file=False)

Here, `method_name` can be either `"Python"` or `"File"`. If using Python-based calculations, you can further customize whether the production or decay BR computations should rely on precomputed data files using the `pythoncalculation_prod_in_file` and `pythoncalculation_decay_in_file` flags.

For example:

.. code-block:: python

    testbr.set_calculation_method("Python", pythoncalculation_prod_in_file=True, pythoncalculation_decay_in_file=False)

Defining the Model
------------------

Next, set the model you wish to use for your calculations. Models should be specified by their name as a string, corresponding to the UFO input you have prepared.

.. code-block:: python

    testbr.set_model("HNL")

In this example, the model is set to `"HNL"` for heavy neutral leptons.

Setting Parameters and Masses
-----------------------------

Once the model is set, you must define the parameters and masses relevant to your computation. Parameters are typically coupling constants, mixing angles, or other model-specific values. They can be passed as a dictionary:

.. code-block:: python

    testbr.set_params({"Ve": 1, "Vmu": 1, "Vta": 1})

Masses for the particles in your model should also be provided as a dictionary:

.. code-block:: python

    testbr.set_masses({"N1": 1})

This example sets the mass of particle `N1` to 1 GeV.

Performing Branching Ratio Calculations
---------------------------------------

With the calculation method, model, parameters, and masses set, you can now compute the branching ratios.

### Decay Total Calculation

To calculate the total decay width (`DecayTot`) of a particle, use the following command:

.. code-block:: python

    testbr.calculate("DecayTot", "N1")

Here, `"DecayTot"` indicates the type of calculation, and `"N1"` is the particle for which the decay width is being computed.

### Decay Branching Ratio

To compute the branching ratio for a specific decay channel, you can use:

.. code-block:: python

    testbr.calculate("BR", "N1", channel=(211, 13))

In this example, the decay branching ratio of particle `N1` into a pion (`211`) and a muon (`13`) is calculated. The `channel` argument specifies the final state particles.

### Production Branching Ratio

To compute the production branching ratio from a mother particle, use the `ProdBR` calculation:

.. code-block:: python

    testbr.calculate("ProdBR", "N1", mother_particle=24)

This computes the branching ratio for producing `N1` from the decay of the W boson (`24`).

Error Handling and Customization
--------------------------------

The interface is designed to handle invalid inputs and calculations gracefully. If you attempt to perform a calculation without setting a model or calculation method, an appropriate error will be raised. Additionally, you can adjust individual parameters using:

.. code-block:: python

    testbr.set_one_param("Ve", 0.8)

This will update the parameter `"Ve"` to `0.8` while leaving other parameters unchanged.

Conclusion
----------

The `BRInterface` class provides a flexible and customizable way to compute branching ratios for various models. By offering different strategies for performing the calculations and allowing for detailed control over input parameters, this interface is a powerful tool for LLP studies.

