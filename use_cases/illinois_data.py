import numpy as np
import scipy
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
from src.qrist import QRIST
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
from lifelines.utils import restricted_mean_survival_time
from lifelines.fitters.kaplan_meier_fitter import KaplanMeierFitter
from time import time
import logging
import os


def main(seed=1):

    group = 'hie'

    data_file = f"example_data/{group}_data.csv"

    output_dir = f'example_data/illinois/'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Configure the logger with a custom time format for messages
    logging.basicConfig(filename=os.path.join(output_dir, 'log.txt'),
                        level=logging.DEBUG,
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        filemode='w')

    Tmax = 26 * 7
    h = 25 * 7
    duration_col = 'Y'
    event_col = 'D'
    treatment_col = 'W'
    iv_col = 'Z'


    n_imputations = 2500
    n_htes = 200
    q_recursion_steps = 3
    features = ['age', 'male', 'white', 'black']
    n_jobs = 25

    np.random.seed(seed)

    patients_df = pd.read_csv(data_file, index_col=0)

    logger = logging.getLogger('MyLogger')
    logger.info(f"Tmax Illinois data: {Tmax}")

    after_Tmax_idx = patients_df[patients_df[duration_col] > Tmax].index
    patients_df.loc[after_Tmax_idx, duration_col] = Tmax
    patients_df.loc[after_Tmax_idx, event_col] = 0

    patients_df = patients_df[features + [treatment_col, duration_col, event_col, iv_col]]

    for min_observed_leaf in [12]:
        trees_kwargs = {'splitter': 'random', 'min_observed_leaf': min_observed_leaf}

        start = time()
        ci = QRIST(Tmax=Tmax, event_col=event_col, duration_col=duration_col, treatment_col=treatment_col)
        ci.fit(train_df=patients_df, features=features, q_recursion_steps=q_recursion_steps,
                       n_imputations=n_imputations, treatment_col=treatment_col, n_jobs=n_jobs,
                       output_dir=output_dir,
                       validate_min_samples_leaf=False, trees_kwargs=trees_kwargs)
        end = time()
        logger.info(f'Full recursive imputation took {int(end - start)} seconds.')

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

        imputed_df = pd.concat([imputed_df.loc[patients_df.index], patients_df.loc[patients_df.index, 'Z']], axis=1)
        imputed_df.to_csv(os.path.join(output_dir, f'{group}_imputed_mol_{min_observed_leaf}_nimp_{n_imputations}_rep_{seed}.csv'))

        del ci


if __name__ == '__main__':
    reps = 1
    for seed in range(1, reps+1):
        print(f'Starting seed {seed}')
        s0 = time()
        main(seed=seed)
        s1 = time()
        print(f"Seed {seed} took {int(s1-s0)} seconds.")
        print('********************************************************************s')
