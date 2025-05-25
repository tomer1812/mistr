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

1. Create a new environment with both Python and R. For example, when using anaconda, type:
```bash
conda create -n mistr-env -c conda-forge python=3.10 r-base=4.2.3
```

2. Activate the new environment. For example, when using anaconda, type:
```
conda activate mistr-env
```

3. clone mistr repository:

```
git clone https://github.com/tomer1812/mistr.git
```

4. Install Python dependencies. In mistr project directory, type:
```bash
pip install -r requirements.txt
```

5. Install R dependencies. In mistr project directory, type:
```
Rscript requirements.R
```

6. Clone the custom version of scikit-survival:
```
git clone https://github.com/tomer1812/scikit-survival.git
```

7. <u> Install the custom version of scikit-survival. </u>

In local scikit-survival dir type: 
```bash
python setup.py build_ext --inplace

pip install -e . --verbose --no-build-isolation --config-settings editable-verbose=true

pip uninstall scikit-learn
```

8. Clone the custom version of scikit-survival:
```bash
git clone https://github.com/tomer1812/scikit-learn.git
```

9. <u> Install the custom version of scikit-learn. </u>

In local scikit-learn dir type: 
```bash
python setup.py build_ext --inplace

pip install -e . --verbose --no-build-isolation --config-settings editable-verbose=true
```

## Running Examples Using Command Line

After installing and activating the environment, navigate to the mistr project directory and run: 
```bash
bash run_examples.sh
```