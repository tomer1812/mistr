import numpy as np
import scipy
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
from src.qrist import QRIST
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
from lifelines.fitters.kaplan_meier_fitter import KaplanMeierFitter
from time import time
import logging
import os


def main(seed=1):
    base_dir = 'example_data/'
    data_file = os.path.join(base_dir, "ACTG175.csv")

    time_unit = 28

    output_dir = os.path.join(base_dir, 'hiv')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    continuous_covariates = ['age', 'wtkg', 'karnof', 'cd40', 'cd80']
    binary_covariates = ['gender', 'homo', 'race', 'symptom', 'drugs', 'hemo', "str2"]
    unobserved_covariate = 'z30'

    # Configure the logger with a custom time format for messages
    logging.basicConfig(filename=os.path.join(output_dir, 'log.txt'),
                        level=logging.DEBUG,
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        filemode='w')

    Tmax = 868 // time_unit  # (31 months)
    h = 840 // time_unit  # (30 months)
    duration_col = 'U'
    event_col = 'Delta'
    treatment_col = 'W'

    colors_list = ['tab:blue', 'tab:green', 'darkorange', 'purple']
    arms_names_list = ["ZDV", "ZDV+ddI", "ZDV+Zal", "ddI"]

    arms_list = [[0, 1], [0, 2], [0, 3]]
    fig, axes = plt.subplots(3, 3, figsize=(12, 10))

    for idp, add_perc in enumerate([0.6]):

        for ida, arms in enumerate(arms_list):

            n_imputations = 2500
            n_htes = 200
            q_recursion_steps = 3
            features = continuous_covariates + binary_covariates
            n_jobs = 25

            int_seed = idp*10000 + ida*1000 + n_imputations + seed
            np.random.seed(int_seed)

            ax = axes[idp, ida]
            # ax = axes
            ax.set_xlabel('T [days]', fontsize=12)
            ax.set_ylabel('P(t > T)', fontsize=12)
            ax.axvline(Tmax, label='Tmax', color='k', ls='--')
            ax.axvline(h, label='h', color='tab:green', ls='--')

            patients_df = pd.read_csv(data_file, index_col=0).set_index('pidnum')

            # treatment arm (0=zidovudine, 1=zidovudine and didanosine, 2=zidovudine and zalcitabine, 3=didanosine).
            patients_df = patients_df[patients_df['arms'].isin(arms)]
            patients_df['arms'].replace({arms[0]: 0, arms[1]: 1}, inplace=True)
            patients_df['months'] = patients_df['days'] // time_unit

            patients_df.rename(columns={'months': duration_col, 'arms': treatment_col, 'cens': event_col}, inplace=True)
            patients_df['preanti'] = (patients_df['preanti'] > 0).values.astype(int)

            if add_perc > 0:
                print(len(patients_df[(patients_df[event_col] == 1) & (patients_df[treatment_col] == 0)]),
                      len(patients_df[(patients_df[event_col] == 1) & (patients_df[treatment_col] == 1)]))

                for idx in patients_df.index:
                    unobs_censoring_perc = 0.25 * patients_df.loc[idx, unobserved_covariate]
                    if (patients_df.loc[idx, duration_col] > 2) and (np.random.random() < (add_perc+unobs_censoring_perc)):
                        obs_event_time = np.min([Tmax // 5, patients_df.loc[idx, duration_col]])
                        sampled_censoring = np.random.randint(1, obs_event_time - 1)
                        patients_df.loc[idx, duration_col] = sampled_censoring
                        patients_df.loc[idx, event_col] = 0

                print(len(patients_df[(patients_df[event_col] == 1) & (patients_df[treatment_col] == 0)]),
                      len(patients_df[(patients_df[event_col] == 1) & (patients_df[treatment_col] == 1)]))

            censoring_perc = 100*(1 - len(patients_df[(patients_df[event_col] == 1) & (patients_df[duration_col] < h)]) / len(patients_df))
            print(100*(len(patients_df[(patients_df[event_col] == 0) & (patients_df[duration_col] < h)]) / len(patients_df)))

            kmf = KaplanMeierFitter()

            untreated_group = patients_df[patients_df[treatment_col] == 0]
            kmf.fit(untreated_group[duration_col], event_observed=untreated_group[event_col], label=arms_names_list[arms[0]])
            ax.step(kmf.survival_function_.index, np.squeeze(kmf.survival_function_.values),
                    color=colors_list[arms[0]], label=arms_names_list[arms[0]])

            # Plotting for group B on the same axes
            treated_group = patients_df[patients_df[treatment_col] == 1]
            kmf.fit(treated_group[duration_col], event_observed=treated_group[event_col],
                    label=arms_names_list[arms[1]])
            ax.step(kmf.survival_function_.index, np.squeeze(kmf.survival_function_.values), color=colors_list[arms[1]],
                    label=arms_names_list[arms[1]])

            ax.set_ylim([0.61, 1.04])
            ax.set_xlim([-40, 1240])
            ax.set_title(f'HIV {ida+1} KM, {censoring_perc:.1f}% censoring at h', fontsize=14)

            # Create a logger object
            logger = logging.getLogger('MyLogger')
            logger.info(f"Tmax HIV data: {Tmax}")

            after_Tmax_idx = patients_df[patients_df[duration_col] > Tmax].index
            patients_df.loc[after_Tmax_idx, duration_col] = Tmax
            patients_df.loc[after_Tmax_idx, event_col] = 0

            features = continuous_covariates + binary_covariates
            patients_df = patients_df[features + [treatment_col, duration_col, event_col]]

            for min_observed_leaf in [3]:

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

                imputed_df.loc[patients_df.index].to_csv(os.path.join(output_dir, f'hiv_imputed_p{int(100*add_perc)}_mol_{min_observed_leaf}_nimp_{n_imputations}_arms_{arms[0]}{arms[1]}_rep_{seed}.csv'))
                del ci
            ax.legend(loc='lower left', fontsize=9)

    fig.tight_layout()
    # fig.savefig(os.path.join(output_dir, f'kaplan_meier_fig_arms_rep_{seed}.png'), dpi=300)


if __name__ == '__main__':
    reps = 10
    for seed in range(1, reps+1):
        print(f'Starting seed {seed}')
        s0 = time()
        main(seed=seed)
        s1 = time()
        print(f"Seed {seed} took {int(s1-s0)} seconds.")
        print('********************************************************************s')
