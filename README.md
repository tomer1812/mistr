# MISTR: Heterogeneous Treatment Effect in Time-to-Event Outcomes

This repository contains the implementation of the procedure and the analyses presented in the paper:

**"Heterogeneous Treatment Effect in Time-to-Event Outcomes: Harnessing Censored Data with Recursively Imputed Trees"**  
Tomer Meir, Uri Shalit, and Malka Gorfine, 42nd International Conference of Machine Learning (ICML) 2025.

## Repository Structure

The repository is organized as follows:

- `src/`: Contains implementations of the MISTR and MISTR-IV procedures
- `simulations/`: Contains the code to replicate the simulation studies from the paper
- `use_cases/`: Contains the code for the use case analyses presented in the paper

## Environment Setup

The procedure requires the "min_observed_leaf" parameter, which sets the minimum number of observations with observed events in each leaf. Thus, **it depends on custom versions of scikit-survival and scikit-learn.**

1. Create a new environment with both Python and R

2. run pip install -r requirements.txt

3. run Rscript requirements.R

4. git clone https://github.com/tomer1812/scikit-survival.git

5. git clone https://github.com/tomer1812/scikit-learn.git

6. <u> In local git/scikit-survival dir: </u>

git submodule update --init

python setup.py build_ext --inplace

pip install -e . --verbose --no-build-isolation --config-settings editable-verbose=true

pip uninstall scikit-learn

7. <u> In local git/scikit-learn dir </u>

python setup.py build_ext --inplace

pip install -e . --verbose --no-build-isolation --config-settings editable-verbose=true


## Running Examples Using Command Line

After installing and activating the environment, navigate to the mistr project directory and run: 

bash run_examples.sh