import numpy as np
import pandas as pd
import glob
import os

import NaiveDE


def aggregate_data(df, annotations):
    X = df.groupby(['cellID'])['abs_X', 'abs_Y'].mean()
    X = X.sort_index()
    exp = pd.DataFrame(df.groupby(['barcode_id', 'cellID'])['barcode_id'].count())
    exp.columns = ['count']

    # add gene names
    tmp = exp.reset_index(['barcode_id'])
    barcode_ids = tmp['barcode_id'] - 1
    RNA_names = annotations.iloc[barcode_ids]
    exp['gene'] = RNA_names.values
    exp = exp.reset_index(['barcode_id']).drop(['barcode_id'], axis=1)
    exp = exp.set_index(['gene'], append=True)

    # reshape to matrix format and filling missing values
    exp = exp.unstack('gene')
    exp = exp.fillna(0)
    exp.columns = exp.columns.droplevel()

    return {'X': X, 'exp': exp}


# def process_mer_fish(exp_file, annotation_file, data_dir):
def process_mer_fish(exp_file, annotation_file, d):

    # read data
    df = pd.read_csv(exp_file, index_col=0)
    annotations = pd.read_csv(annotation_file)['gene']

    # aggregate data
    tmp = aggregate_data(df, annotations)
    X, exp = tmp['X'], tmp['exp']

    # filter practically unobserved genes
    exp = exp.T[exp.sum(0) >= 3].T

    # Get total counts per cell
    tot = pd.DataFrame(exp.sum(1))
    tot.columns = ['total_count']

    # Convert data to log-scale, and account for depth
    dfm = NaiveDE.stabilize(exp.T).T
    res = NaiveDE.regress_out(tot, dfm.T, 'np.log(total_count)').T

    # Add total_count as pseudogene for reference
    res['log_total_count'] = np.log(tot['total_count'])

    res.to_csv(d+'/expressions.txt', sep=' ', header=True, index=False)
    X.to_csv(d+'/positions.txt', sep=',', header=False, index=False)

if __name__ == '__main__':
    data_dir = '/Users/damienarnol1/Documents/local/pro/PhD/other_data/mer-FISH/data/big_image/'
    save_dir = '/Users/damienarnol1/Documents/local/pro/PhD/other_data/mer-FISH/Vprocessing/'
    annotation_file = '/Users/damienarnol1/Documents/local/pro/PhD/other_data/mer-FISH/data//small_image/codebook_no_header.csv'

    all_files = glob.glob(data_dir+'/*')

    i = 0
    for f in all_files:
        d = save_dir + '/' + 'image_' +str(i) + '/'
        print d
        os.makedirs(d)

        process_mer_fish(f, annotation_file, d)

        i+=1
