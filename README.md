Code to simulate the mathematical model of the controlled circulation written by A Kennard, Z Ahmad, K C Peskin and K Ong.

The Matlab code is deprecated. Code to simulate the controlled model of the steady state circulation can be found in steady_state. Code to simulate the ODE model is currently under development in its own branch.

This code requires python 3 and the packages listed in the requirements.txt file.   

To run the code:

```
git clone https://github.com/alexoblomov/circulation_with_gravity
cd circulation with gravity
```

create a python virtual environment (tested on ubuntu 18.04.6)

```
pyenv virtualenv 3.9.0 circulation_model
pyenv activate circulation_model
pip install -r requirements.txt

```
checkout the version of the model you want to run:
steady state model is on the main branch in the steady-state folder

```
cd ../steady_state/
python varyG_conventional_units.py 
```

to run the dynamic model (work in progress)
```
git checkout ode-model
cd dynamic_model
python dynamic_model.py
```
