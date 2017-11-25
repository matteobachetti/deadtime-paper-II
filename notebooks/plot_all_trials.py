from pandas import read_csv
import matplotlib.pyplot as plt
import numpy as np

def plot_scatter_ratios(results, channel, toler_d_clean=0.01,
                            filter={}):
    plt.figure()
    plt.title(channel)
    plt.axhline(0.0, ls='-')
    plt.axvline(1.0, ls='-')

    good_clean = np.abs(results['delta clean']) < toler_d_clean
    minrate, maxrate = np.min(results['incident_rate']), np.max(results['incident_rate'])
    colors = results['incident_rate'] / maxrate
    for segsize in list(set(results['segment size'][good_clean])):
        good_segsize = results['segment size'] == segsize
        good = good_segsize & good_clean
        marker = np.int(np.log2(segsize)//2)
        plt.scatter(results['channel_ratio'][good], results[channel][good],
                    c=colors[good],
                    s=10 * results['length'][good] / np.min(results['length'][good]),
                    marker=(marker, 0), cmap='magma')
        plt.xlabel('Ratio of flux in the two channels')
        plt.ylabel('Relative change of measured rms'
                   ' ({} - clean)/clean'.format(channel))
    plt.colorbar()

def plot_all(file, toler_d_clean=0.01):
    results = read_csv(file)

    plot_scatter_ratios(results, 'delta clean', toler_d_clean=toler_d_clean)
    plot_scatter_ratios(results, 'CS', toler_d_clean=toler_d_clean)
    plot_scatter_ratios(results, 'PDS 1', toler_d_clean=toler_d_clean)
    plot_scatter_ratios(results, 'PDStot', toler_d_clean=toler_d_clean)

def pairplot(file, toler_d_clean=0.5):
    from seaborn import pairplot
    results = read_csv(file)
    good_clean = np.abs(results['delta clean']) < toler_d_clean
    meanrate = results['incident_rate'] // 100 * 100 + 50
    results["rate"] = meanrate
    g = pairplot(results[good_clean], hue="rate",
             plot_kws={"s": 5},
                 vars=['channel_ratio', #'delta clean',
                   'CS', 'PDS 1', 'PDStot'],size=2,
             palette="GnBu_d")
    g.fig.subplots_adjust(right=0.9, left=0.1)

def scatter_matrix(file, toler_d_clean=0.5):
    from pandas.plotting import scatter_matrix
    results = read_csv(file)
    good_clean = np.abs(results['delta clean']) < toler_d_clean

    scatter_matrix(results[good_clean], alpha=0.2, figsize=(6, 6),
                   diagonal='kde')


def stats(file, toler_d_clean=0.01):
    results = read_csv(file)
    good_clean = np.abs(results['delta clean']) < toler_d_clean
    for data in ['CS', 'PDStot', 'PDS 1']:
        print(data + ':')
        filtered = results[data][good_clean] * 100
        print('Mean _relative_ shift '
              '(({} - clean)/clean) = '
              '({:.1f} +- {:.1f})%'.format(data, np.mean(filtered),
                                           np.std(filtered)))


if __name__ == '__main__':
    import sys
    import time
    file = sys.argv[1]
    toler_d_clean = 0.05
    plt.ion()
    while True:
#        plot_all(file, toler_d_clean=toler_d_clean)
        pairplot(file, toler_d_clean=toler_d_clean)
#        stats(file, toler_d_clean=toler_d_clean)
#        scatter_matrix(file, toler_d_clean=toler_d_clean)
        plt.pause(200)
        plt.close('all')
