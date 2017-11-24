from pandas import read_csv
import matplotlib.pyplot as plt
import numpy as np

def plot_scatter_ratios(results, channel, tolerance_delta_clean=0.05,
                            filter={}):
    plt.figure()
    plt.title(channel)
    plt.axhline(0.0, ls='-')
    plt.axvline(1.0, ls='-')

    good_clean = np.abs(results['delta clean']) < tolerance_delta_clean
    for segsize in list(set(results['segment size'][good_clean])):
        good_segsize = results['segment size'] == segsize
        good = good_segsize & good_clean
        marker = np.int(np.log2(segsize)//2)
        plt.scatter(results['channel_ratio'][good], results[channel][good],
                    c=results['incident_rate'][good],
                    s=10 * results['length'][good] / np.min(results['length'][good]),
                    marker=(marker, 0))
    plt.colorbar()

def plot_all(file):
    results = read_csv(file)

    plot_scatter_ratios(results, 'delta clean')
    plot_scatter_ratios(results, 'CS')
    plot_scatter_ratios(results, 'PDS 1')
    plot_scatter_ratios(results, 'PDStot')

def pairplot(file, tolerance_delta_clean=0.01):
    import seaborn as sns
    results = read_csv(file)
    good_clean = np.abs(results['delta clean']) < tolerance_delta_clean

    sns.pairplot(results[good_clean], hue="smoothing length",
                 plot_kws={"s": 10},
                 vars=['channel_ratio', 'delta clean', 'CS', 'PDS 1', 'PDStot'])

def stats(file, tolerance_delta_clean=0.01):
    import seaborn as sns
    results = read_csv(file)
    good_clean = np.abs(results['delta clean']) < tolerance_delta_clean
    for data in ['CS', 'PDStot', 'PDS 1']:
        print(data + ':')
        filtered = results[data][good_clean] * 100
        print('Mean shift = ({:.1f} +- {:.1f})%'.format(np.mean(filtered),
                                                     np.std(filtered)))


if __name__ == '__main__':
    import sys
    file = sys.argv[1]
    plot_all(file)
    pairplot(file)
    stats(file)
    plt.show()
