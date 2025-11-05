Installation
============

Prerequisites
-------------

* Python 3.8 or higher
* pip (Python package installer)

Installation Methods
-------------------

From PyPI (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pip install leishdomains

From Source
~~~~~~~~~~~

1. Clone the repository:

   .. code-block:: bash

      git clone https://github.com/your-org/leish-domains.git
      cd leish-domains

2. Install in development mode:

   .. code-block:: bash

      pip install -e .

3. Install with development dependencies:

   .. code-block:: bash

      pip install -e ".[dev]"

Dependencies
------------

LeishDomains requires the following packages:

* biopython >= 1.79
* numpy >= 1.24.0
* matplotlib >= 3.7.0
* scipy >= 1.10.0
* xlwt >= 1.3.0

Optional dependencies for development:

* pytest >= 7.0.0
* black >= 22.0.0
* flake8 >= 5.0.0
* mypy >= 1.0.0

Verification
------------

To verify the installation, run:

.. code-block:: bash

   leishdomains --help

This should display the help message and available commands.
