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


def pairplot(file, toler_d_clean=0.5, filters={'channel_ratio': [0.8, 1.2]}):
    from seaborn import pairplot
    results = read_csv(file)
    good_clean = np.abs(results['delta clean']) < toler_d_clean
    for key, val in filters.items():
        good_clean = good_clean & (results[key] > val[0])&(results[key] < val[1])
    meanrate = results['incident_rate'] // 100 * 100 + 50
    results["rate"] = meanrate
    results["median_ratio"] = results['median'] / results['clean_median']
    vars_unfilt = ['channel_ratio', 'frequency', 'Frac_rms',
            'CS', 'PDS 1', 'PDStot', 'slope', 'median_ratio']
    vars = []
    for v in vars_unfilt:
        if v in results:
            if np.max(results[v][good_clean]) > np.min(results[v][good_clean]):
                vars.append(v)
            
    g = pairplot(results[good_clean], hue="rate",
                 plot_kws={"s": 5},
                 vars=vars, size=1.7,
                 palette="GnBu_d", dropna=True, diag_kind='kde')
    g.fig.subplots_adjust(right=0.9, left=0.1)
    return g.fig


def pairplot_gain(file, toler_d_clean=0.01):
    from seaborn import pairplot
    results = read_csv(file)
    good_clean = np.isclose(results['pds_m'], results['cs_m'], rtol=0.05)&\
        np.isclose(results['cs_q'], 0, atol=2*np.std(results['cs_q']))&\
        np.isclose(results['pds_q'], np.median(results['pds_q']), atol=2*np.std(results['pds_q']))&\
        (results['slope'] == results['slope'])

#   good_clean = np.abs(results['delta clean']) < toler_d_clean
#   good_clean = np.isclose(results['pds_m'], results['cs_m'], rtol=0.02)
    meanrate = results['incident_rate'] // 100 * 100 + 50
    results["rate"] = meanrate
    g = pairplot(results[good_clean], hue="rate",
             plot_kws={"s": 5},
                 vars=['incident_rate', 'rms', 'pds_m',
                       'pds_q', 'cs_m', 'cs_q', 'slope', 'intercept'], size=1.7,
             palette="GnBu_d")
    g.fig.subplots_adjust(right=0.9, left=0.1)
    return g.fig


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

def line_to_zero(x, a=0):
    return a * x


def square(x, a=0):
    return a * x ** 2

def line(x, m, q):
    return m*x + q
    
def offset(x, q):
    return q


def fit_incident_vs_delta(file, toler_d_clean=0.01, 
                          filters={'channel_ratio': [0.8, 1.2]}):
    from scipy.optimize import curve_fit
    from matplotlib.pyplot import cm
    from matplotlib.gridspec import GridSpec
    import matplotlib as mpl

    results = read_csv(file)
    good_clean = np.abs(results['delta clean']) < toler_d_clean
    for key, val in filters.items():
        good_clean = good_clean & (results[key] > val[0])&(results[key] < val[1])
 
    fig = plt.figure(figsize=(8, 8))

    plt.title('rms error vs rate ratio correlation')
    gs = GridSpec(2, 2, width_ratios=[9.5, 0.5])
    ax0 = plt.subplot(gs[0, 0])
    ax1 = plt.subplot(gs[1, 0])
    ax2 = plt.subplot(gs[:, 1])
    nsub = 9
    frac_rmss = np.linspace(0.05,
                            np.max(results['Frac_rms']), nsub)
    colors=iter(cm.magma(np.linspace(0,1,nsub)))
    coeffs = []
    for frac_rms_min, frac_rms_max in zip(frac_rmss[:-1], frac_rmss[1:]):
        good_rms = (results['Frac_rms'] >= frac_rms_min)&(results['Frac_rms'] < frac_rms_max)
        good = good_rms & good_clean
        x = results['incident_rate'][good] + 1# Equivalent to Rmeas/R = (dR/R + R)/R
        y = results['CS'][good] + 1
        color = next(colors)
        ax0.scatter(x, y, c=color)
        try:
            par, pcov = curve_fit(line_to_zero, x, y - 1)
            x = np.linspace(0, max(x), 50)
            ax0.plot(x, x * par[0] + 1, color=color)
            coeffs.append(par[0])
        except:
            coeffs.append(np.nan)
         
    coeffs = np.array(coeffs)
    ax0.set_xlabel('Incident count rate')
    ax0.set_ylabel('rms relative error')
    ax0.set_xlim([0, np.max(results['incident_rate'])])

    colors=cm.magma(np.linspace(0,1,nsub))

    frac_rmss = np.mean(list(zip(frac_rmss[:-1], frac_rmss[1:])), axis=1)
    ax1.scatter(frac_rmss, coeffs, c=colors)
    good = coeffs == coeffs
    par, pcov = curve_fit(square, frac_rmss[good], coeffs[good], p0=[0.0])
    x = np.linspace(0, np.max(results['Frac_rms']), 50)
    ax1.plot(x, x**2 * par[0])
    ax1.set_xlabel('Frac. rms')
    ax1.set_ylabel('Slope of rms error vs rate ratio correlation')
    ax1.set_xlim([0, np.max(results['Frac_rms'])])
    ax1.set_ylim([np.min([np.min(coeffs), 0]), np.max(coeffs)])
    print(par[0])
    cmap = mpl.cm.magma
    norm = mpl.colors.Normalize(vmin=np.min(results['Frac_rms']),
                                vmax=np.max(results['Frac_rms']))

    cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm)
    cb.set_label('Fractional rms')


def fit_pds_quantities_vs_rms(file, toler_d_clean=0.01, 
                         filters={'channel_ratio': [0.8, 1.2]}, quantity='slope', offset=1, fix=True,
                         slope_fit_kind='square'):
    from scipy.optimize import curve_fit
    from matplotlib.pyplot import cm
    from matplotlib.gridspec import GridSpec
    import matplotlib as mpl

    results = read_csv(file)
    good_clean = (results['slope'] == results['slope'])&\
        (np.abs(results['delta clean']) < toler_d_clean)
    for key, val in filters.items():
        good_clean = good_clean & (results[key] > val[0])&(results[key] < val[1])
 
    fig = plt.figure(figsize=(8, 8))

    plt.title('rms error vs rate ratio correlation')
    gs = GridSpec(3, 2, width_ratios=[9.5, 0.5])
    ax0 = plt.subplot(gs[0, 0])
    ax1 = plt.subplot(gs[1, 0])
    ax3 = plt.subplot(gs[2, 0])
    ax2 = plt.subplot(gs[:, 1])
    nsub = 9
    frac_rmss = np.linspace(np.min(results['Frac_rms'][good_clean]),
                            np.max(results['Frac_rms']), nsub + 1)
    colors=iter(cm.viridis(np.linspace(0,1,nsub)))
    ms = []
    qs = []
    for frac_rms_min, frac_rms_max in zip(frac_rmss[:-1], frac_rmss[1:]):

        good_rms = (results['Frac_rms'] >= frac_rms_min)&(results['Frac_rms'] < frac_rms_max)
        good = good_rms & good_clean

        x = results['incident_rate'][good]
        y = results[quantity][good]
        color = next(colors)
        ax0.scatter(x, y, c=color, s=10)
        try:
            x_plot = np.linspace(0, max(x), 50)
            if fix:
                par, pcov = curve_fit(line_to_zero, x, y - offset)
                ax0.plot(x_plot, x_plot * par[0] + offset, color=color)
            else:
                par, pcov = curve_fit(line, x, y, p0=[0, offset])
                ax0.plot(x_plot, x_plot * par[0] + par[1], color=color)
            ms.append(par[0])
            if fix:
                qs.append(offset)
            else:
                qs.append(par[1])
        except Exception as e:
            print(e)
            ms.append(np.nan)
            qs.append(np.nan)
         
    ms = np.array(ms)
    qs = np.array(qs)
    
    ax0.set_xlabel('Incident rate')
    ax0.set_ylabel(quantity)

    colors=cm.viridis(np.linspace(0,1,nsub))

    frac_rmss = np.mean(list(zip(frac_rmss[:-1], frac_rmss[1:])), axis=1)
    good = ms == ms
    ax1.scatter(frac_rmss[good], ms[good], c=colors[good], cmap='viridis')
    ax1.set_ylim([np.min(ms), np.max(ms)])
 
    x = np.linspace(0, np.max(results['Frac_rms']), 50)

    if slope_fit_kind == 'square':
        par, pcov = curve_fit(square, frac_rmss[good], ms[good], p0=[0.0])
        ax1.plot(x, x**2 * par[0])
    else:
        par, pcov = curve_fit(line, frac_rmss[good], ms[good], p0=[0.0, 0.0])
        ax1.plot(x, x * par[0] + par[1])
        
    ax1.set_xlabel('Frac_rms')
    ax1.set_ylabel('Slope of '+ quantity)
    
    retval = par[0]
    ax3.scatter(frac_rmss[good], qs[good], c=colors[good], cmap='viridis')
    ax3.set_ylim([np.min(qs), np.max(qs)])
    par, pcov = curve_fit(line, frac_rmss[good], qs[good])

    x = np.linspace(0, np.max(results['Frac_rms']), 50)

    ax3.plot(x, x * par[0] + par[1])
    ax3.set_xlabel('Frac_rms')
    ax3.set_ylabel('Intercept of '+ quantity)

    cmap = mpl.cm.viridis
    norm = mpl.colors.Normalize(vmin=np.min(results['Frac_rms']),
                                vmax=np.max(results['Frac_rms']))

    cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm)
    cb.set_label('Fractional rms')
    return retval

if __name__ == '__main__':
    import sys
    import time
    from hendrics.base import r_in, r_det
    import warnings
    warnings.filterwarnings('ignore')
    infile = sys.argv[1]

    toler_d_clean = 0.01
#    plt.ion()
    post_filters = {'channel_ratio': [0.7, 1.5]}
    for f in np.logspace(1, np.log10(500), 5):
        filters ={'frequency': [f, f + 30]}
        table = read_csv(infile)
        table['ratio'] = table['median'] / table['clean_median']
        table['inc_to_det'] = table['incident_rate'] / r_det(2.5e-3, table['incident_rate'])
        table['median_corr'] = table['median'] / table['inc_to_det']
        table['ratio_corr'] = table['median_corr'] / table['clean_median']
        table.to_csv('bubu.csv')
        file = 'bubu.csv'

    #        plot_all(file, toler_d_clean=toler_d_clean)
#        pairplot(file, toler_d_clean=toler_d_clean, filters=filters)
    #       pairplot_gain(file, toler_d_clean=toler_d_clean)
    #        stats(file, toler_d_clean=toler_d_clean)
    #        scatter_matrix(file, toler_d_clean=toler_d_clean)
    #       fit_incident_vs_delta(file, toler_d_clean=toler_d_clean)
        slope = fit_pds_quantities_vs_rms(file, toler_d_clean=0.01,
                                          filters={**filters, **post_filters},
                                          quantity='CS', offset=0)
        print(f + 15, slope)
    #        fit_pds_quantities_vs_rms(file, toler_d_clean=0.01,
    #                             filters={**filters, **post_filters},
    #                             quantity='median', offset=0, fix=False)
    #        fit_pds_quantities_vs_rms(file, toler_d_clean=0.01,
    #                             filters={**filters, **post_filters},
    #                             quantity='clean_median', offset=0, fix=False)
    #        fit_pds_quantities_vs_rms(file, toler_d_clean=0.01,
    #                                  filters={**filters, **post_filters},
    #                                  quantity='ratio', offset=1, fix=True,
    #                                  slope_fit_kind='linear')
    #        fit_pds_quantities_vs_rms(file, toler_d_clean=0.01,
    #                                  filters={**filters, **post_filters},
    #                                  quantity='inc_to_det', offset=1, fix=True,
    #                                  slope_fit_kind='linear')
    #        fit_pds_quantities_vs_rms(file, toler_d_clean=0.01,
    #                                  filters={**filters, **post_filters},
    #                                  quantity='ratio_corr', offset=1, fix=True,
    #                                  slope_fit_kind='linear')

        plt.show()
