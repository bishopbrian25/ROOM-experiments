# ROOM-experiments
This is a PSAMM implementation of Regulatory on/off minimization of metabolic flux changes (ROOM).

## Important Files

* fba_test.py - Deletes each gene one by one and return the percent biomass of the wild type.
* room_test.py - Deletes each gene one by one and return the percent biomass of the wild type.
* room.py - ROOM implementation from the paper.
* gene_deletion.py - Set up and execution of all the expirements.

### Prerequisites

Updated version of [PSAMM](https://github.com/zhanglab/psamm) and a LP solver (CPLEX).

## Authors

* **Brian Bishop** - ROOM implementation

## Acknowledgments

* [Jon Lund Steffensen](https://github.com/jonls) - Initial code for the expirements
* [Zhang Lab](https://github.com/zhanglab/psamm) - PSAMM implementation
