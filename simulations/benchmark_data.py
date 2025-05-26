import matplotlib

from simulations.constants import Tmax_dict, p_dict

matplotlib.use('agg')
from src.qrist import QRIST
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
from time import time
import logging
import os
import numpy as np


def main_class():
    n_patients = 5000
    n_test = 5000
    start_reps = 1
    reps = 100
    n_imputations = 1000
    n_htes = 200
    n_jobs = 25
    q_recursion_steps = 3

    base_dir = "examples_data/"

    types_list = [f'type18']
    data_dir = os.path.join(base_dir, f"data")

    for seed in range(start_reps, reps + 1):

        seed_start = time()

        for datatype in types_list:

            print(f"Starting {datatype} of seed {seed}")
            datatype_start = time()

            output_dir = os.path.join(base_dir, f"output/{datatype}")

            # Configure the logger with a custom time format for messages
            logging.basicConfig(filename=os.path.join(output_dir, 'log.txt'),
                                level=logging.DEBUG,
                                format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                                datefmt='%Y-%m-%d %H:%M:%S',
                                filemode='w')

            # Create a logger object
            p = p_dict[datatype]
            features = [f'X.{i}' for i in range(1, p + 1)]

            logger = logging.getLogger('MyLogger')
            logger.info(f"Tmax[{datatype}]: {Tmax_dict[datatype]}")
            logger.info(f"p_dict[{datatype}]: {p_dict[datatype]}")

            if not os.path.exists(output_dir):
                os.mkdir(output_dir)

            np.random.seed(seed)

            data_file = os.path.join(base_dir, f"data/{datatype}/n_{n_patients}_seed_{seed}.csv")
            patients_df = pd.read_csv(data_file, index_col=0)

            # Change column names for mimic data
            if datatype in [f'type140', f'type141', f'type142', f'type143', f'type144']:
                patients_df.columns = [f'X.{rr}' for rr in range(1, p+1)] + list(patients_df.columns[p:])

            Tmax = Tmax_dict[datatype]

            adms_censoring_idx = patients_df[patients_df['Y'] > Tmax].index
            patients_df.loc[adms_censoring_idx, 'Y'] = Tmax
            patients_df.loc[adms_censoring_idx, 'D'] = 0

            trees_kwargs = {'splitter': 'random', 'min_observed_leaf': 15}

            start = time()
            ci = QRIST(Tmax=Tmax, event_col='D', duration_col='Y', treatment_col='W')
            ci.fit(train_df=patients_df, features=features, q_recursion_steps=q_recursion_steps,
                   n_imputations=n_imputations, treatment_col='W', n_jobs=n_jobs,
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

            imputed_df.loc[patients_df.index].to_csv(os.path.join(output_dir,
                                                                  f'n_{n_patients}_seed_{seed}_imputed.csv'))
            datatype_end = time()
            logger.info(f'Datatype {datatype} of seed {seed} took {int(datatype_end - datatype_start)} seconds.')

        seed_end = time()
        logger.info(f'Full seed took {int(seed_end - seed_start)} seconds.')
        logger.info(f'***********************************************************************************')


if __name__ == "__main__":
    main_class()
