erg_loc
=======

A numerical toolbox to distinguish between ergodic and localized phases in quantum many-body systems.

Prerequisites: Python 3.7.6+, QuSpin 0.3.6+

Static Systems
--------------

The ``programs`` along with their corresponding ``tools``:

* **inst_ham** -- instantaneous Hamiltonian
	* ener -- disorder-averaged energy spectrum
	* ener_spac -- disorder-concatenated energy spacings
	* ent -- disorder-averaged entanglement entropy spectrum
* **W_flow** -- disorder evolution
	* ener_W_flow -- disorder-averaged energy spectrum against disorder amplitude
	* r_W_flow -- disorder-averaged first moment of the energy spacings against disorder amplitude
	* ent_W_flow -- disorder-averaged entanglement entropy spectrum against disorder amplitude
* **t_flow** -- time evolution
	* ener_t_flow -- disorder-averaged energy spectrum against time
	* ent_t_flow -- disorder-averaged entanglement entropy spectrum against time

The ``plot`` scripts:

* **ener_spec** -- energy spectrum

.. image:: figures/ener_spec/heisenberg/ener_spec_heisenberg_L_8_obc_J_1_1_1_W_0.5_comparison.png
	:align: center
	:width: 80%

* **ent_arc** -- entanglement entropy arc

.. image:: figures/ent_arc/heisenberg/ent_arc_heisenberg_L_12_obc_J_1_1_1_W_0.5_comparison.png
	:align: center
	:width: 80%

* **ener_stat** -- energy spacing statistics

.. image:: figures/ener_stat/heisenberg/ener_stat_heisenberg_L_8_Nup_4_pauli_0_obc_dis_10000_J_1_1_1_W_0.5_comparison.png
	:align: center
	:width: 80%

* **ener_stat_W_flow** -- energy spacing statistics disorder evolution

.. image:: figures/ener_stat_W_flow/heisenberg/ener_stat_W_flow_heisenberg_L_8_Nup_4_pauli_0_obc_dis_10000_J_1_1_1_W_0.5_12.5_24.png
	:align: center
	:width: 80%

* **ent_W_flow** -- entanglement entropy disorder evolution

Driven Systems
--------------

Coming soon...

Example Command
---------------

``python code/inst_ham.py -mod heisenberg -thr 1 -L 8 -bc o -dis 1 -W 0.5``

Cluster Usage
-------------

Programs to install:

* (htop)[https://htop.dev/]
* (parallel)[https://www.gnu.org/software/parallel/]

Modules to load:

