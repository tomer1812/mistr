import matplotlib

from simulations.constants import Tmax_dict, p_dict

matplotlib.use('agg')
from src.qrist import QRIST
#from src.plots_and_tables import create_simple1_tableone
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
from time import time
import logging
import os
import numpy as np


def example(datatype='type8', p=5, Tmax=7, seed=1):
    n_patients = 5000
    n_imputations = 1000
    n_htes = 200
    n_jobs = 3
    q_recursion_steps = 3

    base_dir = "./example_data/"

    data_dir = os.path.join(base_dir, f"data")
    output_dir = os.path.join(base_dir, f"output/{datatype}")

    features = [f'X.{i}' for i in range(1, p + 1)]

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    np.random.seed(seed)

    data_file = os.path.join(data_dir, f"{datatype}/n_{n_patients}_seed_{seed}.csv")
    patients_df = pd.read_csv(data_file, index_col=0)

    adms_censoring_idx = patients_df[patients_df['Y'] > Tmax].index
    patients_df.loc[adms_censoring_idx, 'Y'] = Tmax
    patients_df.loc[adms_censoring_idx, 'D'] = 0

    trees_kwargs = {'splitter': 'random', 'min_observed_leaf': 15}

    ci = QRIST(Tmax=Tmax, event_col='D', duration_col='Y', treatment_col='W')
    ci.fit(train_df=patients_df, features=features, q_recursion_steps=q_recursion_steps,
           n_imputations=n_imputations, treatment_col='W', n_jobs=n_jobs,
           output_dir=output_dir,
           validate_min_samples_leaf=False, trees_kwargs=trees_kwargs)

    imputed_df = ci.impute_censored(patients_df, features=features, treatment_col=ci.treatment_col,
                                    n_samples=n_htes)

    imputed_idx = patients_df[(patients_df[ci.event_col] == 0) &
                              (patients_df[ci.duration_col] < ci.Tmax)].index
    not_imputed_idx = patients_df[~patients_df.index.isin(imputed_idx)].index

    for ii in range(n_htes):
        imputed_df[f"Di{ii}"] = imputed_df[ci.event_col]
        imputed_df.loc[imputed_idx, f"Di{ii}"] = 1
        isTmax = imputed_df[imputed_df[f"Ti{ii}"] == ci.Tmax].index
        imputedTmax = isTmax.intersection(imputed_idx)
        imputed_df.loc[imputedTmax, f"Di{ii}"] = 0

    imputed_df.loc[not_imputed_idx, [f'Ti{ii}' for ii in range(n_htes)]] = np.repeat(
        np.expand_dims(imputed_df.loc[not_imputed_idx, ci.duration_col].values, axis=1), n_htes, axis=1)

    imputed_df.loc[patients_df.index].to_csv(os.path.join(output_dir, f'n_{n_patients}_seed_{seed}_imputed.csv'))


if __name__ == "__main__":
    example(datatype='type8', p=5, Tmax=7, seed=1)
    example(datatype='type200', p=3, Tmax=9, seed=1)
