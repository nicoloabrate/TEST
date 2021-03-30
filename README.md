## TEST (Transport Equation Solver in Torino)
This is a suite of functions and procedures to solve 1D neutron transport equation in simple yet physically significant cases. TEST can be profitably used both as a test bench for different computational methods and a tool for approximated neutron transport calculations.

## Basic features

1. Multi-layer slab geometry;
2. Multi-group energy approximation;
3. PN and diffusion angular approximations;
4. Higher-order harmonics for different eigenvalue problems (criticality, collision, time...);
5. Perturbative calculations (GPT method);

## Installation and pre-requisites
To install the `master` version of the code, type:
        ```bash
	git clone https://nicolo_abrate@bitbucket.org/nicolo_abrate/test.git
        cd TEST
        python setup.py install 
	```

To install the `development` version of the code, type:
        ```bash
	git clone https://nicolo_abrate@bitbucket.org/nicolo_abrate/test.git
        cd TEST
        git checkout devel
        python setup.py install 
	```

To ensure the correct code installation, type
        ```bash
	python setup.py test
	```
