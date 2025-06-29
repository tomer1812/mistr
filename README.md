# MISTR: Heterogeneous Treatment Effect in Time-to-Event Outcomes

This repository contains the implementation of the procedure and the analyses presented in the paper:

**"Heterogeneous Treatment Effect in Time-to-Event Outcomes: Harnessing Censored Data with Recursively Imputed Trees"**  
Tomer Meir, Uri Shalit, and Malka Gorfine, 42nd International Conference of Machine Learning (ICML) 2025.

![MISTR](figure1.png)  
Figure 1. Our goal is to estimate the heterogeneous treatment effect, defined as the expected difference in survival times (or their transformation) with and without treatment, conditional on a set of covariates. We propose a multiple-imputation-based estimator that effectively leverages censored observations, outperforms existing methods, and is applicable in settings with instrumental variable adjustment for unobserved confounders.

<br><br>
Set up the environment and execute run_examples.sh for quick end-to-end usage examples. See detailed instructions below.

## Repository Structure

The repository is organized as follows:

- `src/`: Contains implementations of the MISTR and MISTR-IV procedures
- `simulations/`: Contains the code of the simulation studies from the paper
- `use_cases/`: Contains the code for the use case analyses presented in the paper

## Environment Setup

The procedure requires the "min_observed_leaf" parameter, which sets the minimum number of observations with observed events in each leaf. Thus, **it depends on custom versions of scikit-survival and scikit-learn.**

1. Create a new environment with both Python and R. For example, when using anaconda, type:
```bash
conda create -n ENV-NAME -c conda-forge python=3.10 r-base=4.4.3
```

2. Activate the new environment. For example, when using anaconda, type:
```
conda activate ENV-NAME
```

3. Install xcode
```bash
xcode-select --install
```

4. clone mistr repository:

```
git clone https://github.com/tomer1812/mistr.git
```

5. Clone the custom version of scikit-survival:
```
git clone https://github.com/tomer1812/scikit-survival.git
```

6. Clone the custom version of scikit-learn:
```bash
git clone https://github.com/tomer1812/scikit-learn.git
```

7. Install dependencies:
```bash
pip install ecos joblib numexpr numpy osqp pandas scipy ninja scikit-learn packaging Cython lifelines tableone matplotlib meson-python 
```

8. <u> Install the custom version of scikit-survival. </u>

In local scikit-survival dir type: 
```bash
git submodule update --init

python setup.py build_ext --inplace

pip install -e . --verbose --no-build-isolation --config-settings editable-verbose=true
```

9. <u> Make sure scikit-survival did not install the external scikit-learn. </u>

In local scikit-learn dir type: 
```bash
pip uninstall scikit-learn

python setup.py build_ext --inplace

pip install -e . --verbose --no-build-isolation --config-settings editable-verbose=true
```

10. Install Python dependencies. In mistr project directory, type:
```bash
pip install -r requirements.txt
```

11. Install R dependencies. In mistr project directory, type:
```bash
Rscript requirements.R
```

## Running Examples Using Command Line

After installing and activating the environment, navigate to the mistr project directory and run: 
```bash
bash run_examples.sh
```

## Results

After the `run_examples.sh` script completes successfully, examine the `example_data/output` directory.

For each data type, there should be a file ending with `_test_with_htes.csv` or `_test_with_htes_iv.csv`.  
This file includes the following key columns:

- `X.i` with i=1, ..., 5: feature columns.
- `W`: treatment assignment status.
- `Y`: last observed time.
- `D`: observed event indicator.
- **`cate`**: the Monte Carlo approximation of the true heterogeneous treatment effect (HTE)
- **`MEAN_HTE`**: the HTE estimates produced by MISTR
- **`TOTAL_VAR_MISTR`** or **`TOTAL_VAR_MISTRIV`**: the estimated variance of the MISTR\MISTR-IV estimates, respectively.
