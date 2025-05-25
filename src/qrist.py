import pandas as pd
from sksurv.tree import SurvivalTree
import numpy as np
from joblib import Parallel, delayed
from time import time
import itertools
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('agg')
import os
import logging
import sys


logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

DEFAULT_LOGGER = logging.getLogger('MyLogger')


SURVIVAL_TREES_KWARGS = dict(splitter='random')


class QRIST(object):

    def __init__(self, Tmax, duration_col='U', event_col='Delta', treatment_col='W', logger=DEFAULT_LOGGER):
        self.imputation_trees = []
        self.features = None
        self.Tmax = Tmax
        self.duration_col = duration_col
        self.event_col = event_col
        self.treatment_col = treatment_col
        self.trees_kwargs = None
        self.times = []
        self.logger = logger
        self.low_memory = False
        self.dtype = np.float64
        self.min_observed_per_leaf = None
        self.validate_min_samples_leaf = False

    def fit(self, train_df, q_recursion_steps, n_imputations,
            features, treatment_col='W', n_jobs=-1, trees_kwargs=SURVIVAL_TREES_KWARGS,
            output_dir="", low_memory=True, validate_min_samples_leaf=False, min_observed_per_leaf=8):

        self.low_memory = low_memory
        self.features = features
        self.treatment_col = treatment_col
        if self.low_memory:
            self.dtype = np.float32
            train_df = self._reduce_memory_footprint(train_df)
            self.Tmax = np.float32(self.Tmax)

        original_train_df = train_df.copy()

        self.min_observed_per_leaf = min_observed_per_leaf
        self.validate_min_samples_leaf = validate_min_samples_leaf

        self.logger.info(f'Fitting imputation trees...')
        start_fit = time()
        self.fit_imputation_trees(train_df, n_estimators=n_imputations,
                                  features=features, treatment_col=treatment_col,
                                  event_col=self.event_col, duration_col=self.duration_col,
                                  trees_kwargs=trees_kwargs,
                                  n_jobs=n_jobs)
        end_fit = time()
        self.logger.info(f'Fitting imputation trees took {int(end_fit - start_fit)} seconds.')

        imputation_candidates = original_train_df[original_train_df[self.event_col] == 0]
        imputation_candidates = imputation_candidates[imputation_candidates[self.duration_col] < self.Tmax]
        self.times = self._update_times(train_df, n_imputations)

        for q_step in range(1, q_recursion_steps+1):
            self.logger.info(f'Recursion step {q_step}')
            if len(imputation_candidates) == 0:
                self.logger.warning(f'No loss of followup censoring - imputation stage is skipped')
                break
            self.logger.info(f'Loss of followup censoring in fold: {len(imputation_candidates)}')

            train_df = original_train_df.copy()
            start = time()
            self.logger.info('Timing imputations start.')

            train_df = self.impute_censored(train_df, features=features, treatment_col=treatment_col,
                                            n_samples=n_imputations, n_jobs=n_jobs)
            end = time()

            self.logger.info(f'Imputation took {int(end - start)} seconds.')

            not_imputed_idx = train_df[train_df['Ti0'].isnull()].index
            imputed_idx = train_df[train_df['Ti0'].notnull()].index
            for ii in range(n_imputations):
                train_df[f"Di{ii}"] = train_df[self.event_col]
                if len(imputed_idx) > 0:
                    train_df.loc[imputed_idx, f"Di{ii}"] = 1
                    isTmax = train_df[train_df[f"Ti{ii}"] == self.Tmax].index
                    imputedTmax = isTmax.intersection(imputed_idx)
                    train_df.loc[imputedTmax, f"Di{ii}"] = 0

            train_df.loc[not_imputed_idx, [f'Ti{ii}' for ii in range(n_imputations)]] = np.repeat(
                np.expand_dims(train_df.loc[not_imputed_idx, self.duration_col].values, axis=1), n_imputations, axis=1)

            self.times = self._update_times(train_df, n_imputations)

            self.logger.info(f'Fitting tree by tree')
            start_fit = time()
            self.trees_kwargs = trees_kwargs
            self.imputation_trees = []
            parallel = Parallel(n_jobs=n_jobs, verbose=1)
            self.imputation_trees = parallel(
                delayed(self.fit_imputation_tree)(train_df[features + [treatment_col, f"Di{i}", f"Ti{i}"]],
                                                  features, treatment_col, f"Di{i}", f"Ti{i}", trees_kwargs)
                                                  for i in range(n_imputations))
            end_fit = time()
            self.logger.info(f'Fitting imputation trees took {int(end_fit - start_fit)} seconds.')

        del train_df, imputation_candidates, original_train_df
        return

    def _validate_observed_in_nodes(self, train_df, event_col, n_jobs=-1):
        parallel = Parallel(n_jobs=n_jobs, verbose=1)
        X = train_df[self.features + [self.treatment_col]].values.astype(self.dtype)
        D = train_df[event_col].values
        min_observed_in_nodes = parallel(
            delayed(self._get_min_observed_in_node)(self.imputation_trees[i].tree_, X, D)
            for i in range(len(self.imputation_trees)))
        return np.min(min_observed_in_nodes)

    @staticmethod
    def _get_min_observed_in_node(tree, X, D):
        nodes = tree.apply(X)
        observed_in_nodes = []
        samples_in_nodes = []
        for node in np.unique(nodes):
            observed_in_nodes.append(D[nodes == node].sum())
            samples_in_nodes.append(sum(nodes == node))
        return np.min(observed_in_nodes)

    @staticmethod
    def _reduce_memory_footprint(df):
        for col in df.select_dtypes(include=['float64']):
            df[col] = df[col].astype('float32')
        for col in df.select_dtypes(include=['int64']):
            df[col] = df[col].astype('int32')
        return df

    def _update_times(self, train_df, n_imputations):
        time_arr = train_df[train_df[self.event_col] == 1][self.duration_col].values
        if f"Ti0" in train_df.columns:
            for i in range(n_imputations):
                time_arr = np.append(time_arr, train_df[train_df[f"Di{i}"] == 1][f"Ti{i}"].values)
        final = np.array(sorted(np.unique(list(np.unique(time_arr)) + [0, self.Tmax])), dtype=self.dtype)
        return final

    def fit_imputation_trees(self, df, features, treatment_col, event_col, duration_col, n_estimators,
                             trees_kwargs=SURVIVAL_TREES_KWARGS, n_jobs=-1):
        self.trees_kwargs = trees_kwargs
        self.imputation_trees = []
        required_df = df[features + [treatment_col, event_col, duration_col]]
        parallel = Parallel(n_jobs=n_jobs, verbose=1)
        self.imputation_trees = parallel(delayed(self.fit_imputation_tree)(required_df, features, treatment_col,
                                                                           event_col, duration_col, trees_kwargs)
                                         for i in range(n_estimators))

    def fit_imputation_tree(self, df, features, treatment_col, event_col, duration_col, trees_kwargs):
        X = df[features + [treatment_col]]
        y = df[[event_col, duration_col]].copy()
        y.columns = ['event', 'time']
        y['event'] = y['event'].astype(bool, copy=False)
        y['time'] = y['time'].astype(self.dtype, copy=False)
        tree = SurvivalTree(**trees_kwargs)
        tree.fit(X, y.to_records(index=False))
        if self.validate_min_samples_leaf:
            min_observed = self._get_min_observed_in_node(tree, X.astype(np.float32).values, df[event_col].values)
            self.logger.info(f'Min observed: {int(min_observed)}.')
        return tree

    def impute_censored(self, df, features, treatment_col, n_samples=10, n_jobs=-1):
        lof_censoring = (df[self.event_col] == 0) & (df[self.duration_col] < self.Tmax)
        sampled_times = self.sample_event_times(df[lof_censoring], features, treatment_col, n_samples=n_samples,
                                                n_jobs=n_jobs)
        return pd.concat([df[~lof_censoring], sampled_times])

    def sample_event_times(self, df, features, treatment_col, n_samples=10, n_jobs=-1):
        if len(df) == 0:
            for i in range(n_samples):
                df[f'Ti{i}'] = np.nan
            return df
        _features = features + [treatment_col]
        surv_funcs = self.predict_tree_by_tree(df, _features)

        censoring_times = df[self.duration_col].values
        censoring_indices = np.searchsorted(self.times, censoring_times, side='right') - 1

        sampled_times_array = np.zeros((surv_funcs.shape[0], n_samples), dtype=self.dtype)

        for i in range(surv_funcs.shape[0]):
            factor = surv_funcs[i, censoring_indices[i]]
            conditional_surv = surv_funcs[i, censoring_indices[i]:] / factor

            random_samples = np.random.uniform(0, 1, n_samples).astype(self.dtype)

            conditional_surv_matrix = np.tile(conditional_surv, (n_samples, 1)).T
            comparison_matrix = conditional_surv_matrix >= random_samples

            indices_first_false = np.argmax(~comparison_matrix, axis=0)
            sampled_times_array[i, :] = self.times[censoring_indices[i]:][indices_first_false]

        sampled_times_df = pd.DataFrame(sampled_times_array,
                                        columns=[f'Ti{s}' for s in range(n_samples)],
                                        index=df.index)

        full_df = pd.concat([df, sampled_times_df], axis=1)
        return full_df

    def predict_tree_by_tree(self, df, features):
        df_features = df[features]
        surv_funcs_generator = (self.apply_per_tree(tree, df_features, features) for tree in self.imputation_trees)
        mean_surv_funcs = np.mean(np.array(list(surv_funcs_generator)), axis=0, dtype=self.dtype)
        return mean_surv_funcs

    def apply_per_tree(self, tree, df, features):
        surv_df = pd.DataFrame(tree.predict_survival_function(df[features], return_array=True),
                               columns=tree.unique_times_)

        if 0 not in tree.unique_times_:
            surv_df.insert(0, 0, pd.Series(np.ones(len(df))))

        if self.Tmax not in tree.unique_times_:
            surv_df.insert(len(surv_df.columns), self.Tmax, pd.Series(np.zeros(len(df))))
        else:
            surv_df.loc[:, self.Tmax] = 0

        surv_df = surv_df.reindex(columns=self.times, method='ffill')

        if self.low_memory:
            return surv_df.values.astype(self.dtype)
        return surv_df.values
