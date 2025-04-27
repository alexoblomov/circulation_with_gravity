Code to simulate the mathematical model of the controlled circulation written by A Kennard, Z Ahmad, K C Peskin and K Ong.

Code to simulate the controlled model of the steady state circulation can be found in steady_state. Code to simulate the ODE model is currently under development in its own branch.

This code requires python 3 and the packages listed in the requirements.txt file.   

To run the code:

```
git clone https://github.com/alexoblomov/circulation_with_gravity
cd circulation_with_gravity
```

create a python virtual environment (tested on ubuntu 18.04.6)

```
pyenv virtualenv 3.9.0 circulation_model
pyenv activate circulation_model
```

Install the requirements

```
pip install -r requirements.txt

```
checkout the version of the model you want to run:
e.g. steady state model is on the main branch in the steady-state folder

```
cd ../steady_state/
python varyG_conventional_units.py 
```
for plots used in the paper, checkout the steady_state_plots branch

to run the dynamic model (work in progress)
```
git checkout ode-model
cd dynamic_model
python dynamic_model.py
```



# HOW TO SETUP IN WINDOWS

1. Clone the repository into Git Bash:
```
git clone https://github.com/alexoblomov/circulation_with_gravity
cd circulation with gravity
```

2. Install python into windows - go to COMMAND PROMPT in Start Menu (this will download and install python)
```
python
```

3. Create a python virtual environment using venv:
```
python -m venv circulation_model
```

4. Activate it with the same command:
```
python -m venv circulation_model
```

5. Install requirements:
```
pip install -r requirements.txt
```

6. Now go to the steady state directory
```
cd steady_state
```

7. Run the program
```
python varyG_conventional_units.py 
```

8. View the results (graphical *.png files in the same directory)
```
dir
```

