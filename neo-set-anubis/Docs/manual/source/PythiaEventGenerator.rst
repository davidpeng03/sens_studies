.. _PythiaEventGenerator:

Pythia Event Generator
======================

This script facilitates the generation of particle physics events using the Pythia simulation framework. It provides flexibility for handling multiple configuration files and outputting results in `.lhe` and `.hepmc` formats. The script also supports custom file naming conventions, including timestamping, and allows specification of input/output directories.

Overview
--------

The script performs the following key tasks:

- **Event Generation**: Uses Pythia to generate events based on configuration files.
- **Customizable Output**: Allows users to define output filenames with optional timestamp inclusion.
- **Directory Management**: Ensures necessary directories for `.lhe` and `.hepmc` files exist before generating events.
- **Flexible Input**: Accepts a single `.cmnd` configuration file, a list of files, or an entire directory of `.cmnd` files.

Functions
---------

create_pythia_generator
-----------------------

.. autofunction:: create_pythia_generator

This function initializes and returns a Pythia generator with the specified parameters.

- **Arguments**:

    - `config_file`: Path to the Pythia configuration `.cmnd` file.
    - `lhe_output`: Path for the output `.lhe` file.
    - `hepmc_output`: Path for the output `.hepmc` file.
    - `num_events`: Number of events to generate.

- **Returns**: A configured Pythia generator instance.

process_file
------------

.. autofunction:: process_file

Generates events for a given Pythia configuration file and writes the output to `.lhe` and `.hepmc` files.

- **Arguments**:

    - `config_file`: Path to the Pythia configuration `.cmnd` file.
    - `output_lhe_dir`: Directory for the output `.lhe` file.
    - `output_hepmc_dir`: Directory for the output `.hepmc` file.
    - `num_events`: Number of events to generate.
    - `suffix`: Suffix to append to the output filenames.
    - `include_time`: Boolean flag to determine whether to include a timestamp in the filenames.

ensure_directories
------------------

.. autofunction:: ensure_directories

Ensures that the specified directories exist, creating them if necessary.

- **Arguments**:

    - `base_dir`: Base directory where the subdirectories should be created.
    - `sub_dirs`: List of subdirectories to ensure exist within `base_dir`.

- **Returns**: List of full paths to the subdirectories.

Command-Line Arguments
----------------------

The script accepts several command-line arguments to configure the input, output, and behavior of the event generation process.

.. code-block:: bash

    usage: Pythia Event Generator [-h] [-n NUM_EVENTS] [-i INPUT] [-d INPUT_DIR]
                                  [-o OUTPUT_DIR] [-s SUFFIX] [-t]

    optional arguments:
      -h, --help            Show this help message and exit
      -n NUM_EVENTS, --num_events NUM_EVENTS
                            Number of events to generate (default: 2000)
      -i INPUT, --input INPUT
                            Input .cmnd file, list of files, or directory containing .cmnd files
      -d INPUT_DIR, --input_dir INPUT_DIR
                            Input directory containing .cmnd files (default: db/Temp/Pythia/cmnd)
      -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                            Base output directory for generated files (default: db/Temp/Pythia)
      -s SUFFIX, --suffix SUFFIX
                            Suffix for the output files (default: output)
      -t, --time            Include timestamp in output filenames

### Example Usage:

1. **Generate events from a single configuration file:**

    .. code-block:: bash

        python pythia_event_generator.py --input config_file.cmnd --num_events 5000

2. **Generate events from multiple `.cmnd` files listed in the argument:**

    .. code-block:: bash

        python pythia_event_generator.py --input "config1.cmnd,config2.cmnd" --output_dir /path/to/output

3. **Generate events from all `.cmnd` files in a directory, adding a timestamp to the output files:**

    .. code-block:: bash

        python pythia_event_generator.py --input /path/to/config_dir --time

Directory Management
--------------------

The script creates subdirectories within the base output directory for storing `.lhe` and `.hepmc` files. These directories are automatically created if they do not exist, ensuring a structured organization for the output files.

### Example:

If the base output directory is `db/Temp/Pythia`, the following subdirectories will be created:

- `db/Temp/Pythia/lhe`
- `db/Temp/Pythia/hepmc`

Conclusion
----------

The Pythia Event Generator script is a powerful tool for particle physics simulations, allowing users to customize input and output handling through a flexible command-line interface. By supporting multiple input formats and structured output, this script simplifies the process of generating Pythia events for different particle configurations.
