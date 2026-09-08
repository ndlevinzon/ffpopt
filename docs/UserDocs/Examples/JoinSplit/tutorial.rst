.. _JoinSplit-tutorial:


Other ffpopt Json File operations
=================================

This tutorial will show how split and join json files.


Learning Objectives
-------------------


Tutorial
--------
1. User Options
~~~~~~~~~~~~~~~

 .. code-block::

    [user@comp] ffpopt-JsonJoin.py --help
    usage: ffpopt-JsonJoin.py [-h] [--strict] --out OUT inpjson [inpjson ...]

    Split systems within a json file into separate files (or extract a molecule by name)

    positional arguments:
      inpjson            One or more json files to join

    options:
      -h, --help         show this help message and exit
      --strict           Do not rename molecules
      --out OUT, -o OUT  Output json file

      
 .. code-block::
    
    [user@code] ffpopt-JsonSplit.py --help
    usage: ffpopt-JsonSplit.py [-h] [--name NAME] inp

    Split systems within a json file into separate files (or extract a molecule by name)

    positional arguments:
      inp          Input json file

    options:
      -h, --help   show this help message and exit
      --name NAME  Extract structure with a given name

      

2. Split a json file
~~~~~~~~~~~~~~~~~~~~
Suppose foo.json contains 3 structures with names "s000", "s001", and "s002".

.. code-block::

   [user@comp] ffpopt-JsonSplit.py foo.json
   [user@comp] ls *json
   foo.json s000.json s001.json s002.json

   
3. Extract a json file for a specific structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Suppose foo.json contains 3 structures with names "s000", "s001", and "s002".

.. code-block::

   [user@comp] ffpopt-JsonSplit.py --name=s001 foo.json
   [user@comp] ls *json
   foo.json s001.json

   
4. Join several json files
~~~~~~~~~~~~~~~~~~~~~~~~~~
Suppose you have s000.json, s001.json, and s002.json.


.. code-block::

   [user@comp] ffpopt-JsonJoin.py --strict --out=foo.json s000.json s001.json s002.json




