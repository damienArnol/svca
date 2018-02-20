import pandas as pd
import numpy as np
import scipy.io as sio
import glob
import os
# from PIL import Image
from tifffile import imread

import NaiveDE


def get_centroids(mask):
    cell_ix = np.sort(np.unique(mask))
    cell_ix = cell_ix[1:]

    X = np.zeros([len(cell_ix), 2])

    for cell_i in range(len(cell_ix)):

        tmp = mask == cell_ix[cell_i]

        # get x position
        tmp2 = tmp.sum(0)
        ix = np.array([i for i in range(len(tmp2)) if tmp2[i] !=0])
        X[cell_i, 0] =  ix.mean()

        # get y position
        tmp2 = tmp.sum(1)
        ix = np.array([i for i in range(len(tmp2)) if tmp2[i] !=0])
        X[cell_i, 1] =  ix.mean()

    return pd.DataFrame(X)

def aggregate_image(tiff_dir, mask_file):
    tiff_files = glob.glob(tiff_dir+'/*')

    # grep protein name
    tmp = [f.split('/')[-1] for f in tiff_files]
    prots = [s.split('(')[0] for s in tmp]

    # reading mask file
    try:
        mask = sio.loadmat(mask_file)['Image']
    except:
        tmp_m = imread(mask_file)
        mask = np.array(tmp_m)

    mask = mask[2:, :]  # removing first two rows

    # getting centroids
    X = get_centroids(mask)

    # flatten mask for following steps
    mask = mask.flatten()

    # reading and aggregating tiff files
    dfs = []
    for i in range(len(tiff_files)):
        prot = prots[i]
        tf = tiff_files[i]

        im = imread(tf)
        imarray = np.array(im)
        imarray = imarray[2:, :]
        imarray = imarray.flatten()

        tmp = pd.DataFrame()
        tmp['exp'] = imarray
        tmp['cell'] = mask
        tmp['protein'] = prot

        # remove non-cell data
        tmp = tmp[tmp['cell'] != 0]

        dfs.append(tmp)

    aggregated_data = pd.concat(dfs)
    # TODO check whats best
    # exp = aggregated_data.groupby(['cell', 'protein']).sum()
    exp = aggregated_data.groupby(['cell', 'protein']).median()
    exp = exp.unstack('protein')
    exp.columns = exp.columns.droplevel()

    return {'exp': exp, 'X':X}

def filter_out(exp):
    # to_drop = ['Akt',
            #    'AMPKalpha',
            #    'Bad',
            #    'p53',
            #    'p38',
            #    'CD31',
            #    'EGFR',
            #    'SHP2',
            #    'CD45',
            #    'Ru1',
            #    'Ru2',
            #    'Ru3',
            #    'Ru4',
            #    'Ru5',
            #    'Ru6',
            #    'Ru7']
    # return exp.drop(to_drop, axis=1)
    return exp

def process_image(tiff_dir, mask_file, d):

    # read and aggregate data
    data = aggregate_image(tiff_dir, mask_file)
    X, exp = data['X'], data['exp']
    exp = filter_out(exp)

    # Get total counts per cell
    tot = pd.DataFrame(exp.sum(1))
    tot.columns = ['total_count']

    # remove cells with total count bellow 3
    X = X[tot.values > 3]
    exp = exp[tot.values > 3]
    tot = tot[tot.values > 3]

    # Convert data to log-scale, and account for depth
    dfm = NaiveDE.stabilize(exp.T).T

    res = NaiveDE.regress_out(tot, dfm.T, 'np.log(total_count)').T

    # Add total_count as pseudogene for reference
    # res['log_total_count'] = np.log(tot['total_count'])

    res.to_csv(d+'/expressions.txt', sep=' ', header=True, index=False)
    X.to_csv(d+'/positions.txt', sep=',', header=False, index=False)

def process_all_images(images, masks, d):
    image_files = np.sort(glob.glob(images+'/*'))
    masks_files = np.sort(glob.glob(masks+'/*'))

    assert len(image_files) == len(masks_files), 'number of image != number of masks'

    ms = [m.split('/')[-1] for m in masks_files]
    ms =[m.split('_')[0] for m in ms]

    im = [i.split('/')[-1] for i in image_files]
    im =[i.split('_')[0] for i in im]

    for i in range(len(image_files)):
        assert im[i] == ms[i], 'mask image mismatch'

        di = d + '/' +im[i]
        os.makedirs(di)

        tmp = glob.glob(image_files[i]+'/*')
        if 'CellProfiler' in tmp:
            im_dir = [t for t in tmp if t.split('/')[-1] != 'CellProfiler'][0]
        else:
            im_dir = image_files[i]

        process_image(im_dir, masks_files[i], di)



if __name__ == '__main__':
    images_dir = '/Users/damienarnol1/Documents/local/pro/PhD/spatial/data/normal_tissue/Data_6_NormalTissueSamples/'
    masks_dir ='/Users/damienarnol1/Documents/local/pro/PhD/spatial/data/normal_tissue/masks'
    d = '/Users/damienarnol1/Documents/local/pro/PhD/spatial/data/normal_tissue/normal_processed/'


    process_all_images(images_dir, masks_dir, d)
