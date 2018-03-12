import os
import numpy
import pickle
import pandas
import matplotlib.pyplot as plt

def iqr(data):
    q75 = numpy.percentile(data, 75, axis=0)
    q25 = numpy.percentile(data, 25, axis=0)
    return q75 - q25

def get_interp(t_meas, x_meas, t_interp):
    filt = (pandas.isnull(x_meas) == False)
    t_valid = t_meas[filt]
    x_valid = x_meas[filt]
    x_interp = numpy.interp(t_interp, t_valid, x_valid)
    return x_interp

def interpolate_cells(df, groups, t_interp):
    grouped = df.groupby(groups)
    df_interp = pandas.DataFrame(columns=(groups + ['Values']))
    for conditions, group in grouped:
        t_meas = group['Readout time [min]']
        x_meas = group['FoxO3a Ratio [log10 C/N]']
        x_interp = get_interp(t_meas, x_meas, t_interp)
        # Make a row for the new data frame
        row = {}
        for condition, value in zip(groups, conditions):
            row[condition] = value
        row['Values'] = x_interp
        # Add the row to the new data frame
        df_interp = df_interp.append(row, ignore_index=True)
    return df_interp

def calculate_iqrs(df, groupd):
    grouped = df.groupby(groups)
    df_iqr = pandas.DataFrame(columns=(groups + ['IQRs']))
    for conditions, group in grouped:
        values = numpy.vstack(group['Values'].as_matrix())
        iqrs = iqr(values)
        row = {}
        for condition, value in zip(groups, conditions):
            row[condition] = value
        row['IQRs'] = iqrs
        df_iqr = df_iqr.append(row, ignore_index=True)
    return df_iqr

def merge_iqr_pulse(df_iqr, df_pulse_frac, groups):
    df_iqr['Fraction of pulsing'] = \
            pandas.Series(numpy.random.randn(len(df_iqr))*numpy.inf)
    iqr_grouped = df_iqr.groupby(groups)
    pulse_grouped = df_pulse_frac.groupby(groups)
    for conditions, group in iqr_grouped:
        pfrac = pulse_grouped.get_group(conditions)['Fraction of pulsing']
        pfrac = pfrac.as_matrix()[0]
        filt = True
        for c, v in zip(groups, conditions):
            filt &= (df_iqr[c] == v)
        df_iqr.loc[filt, 'Fraction of pulsing'] = pfrac
    return df_iqr

def plot_pulsing_iqr(df, groups, ligand):
    df = df[df['Drug'] == 'None']
    df = df[df['Ligand'] == ligand]
    df = df[df['Ligand dose [ng/ml]'] >= 5]
    plt.figure()
    plt.ion()
    markers = {'CI-1040': 's', 'MK-2206': '^', 'None': 'o'}
    colors = {'BTC': 'm', 'EGF': 'r', 'EPR': 'k', 'FGF': 'g',
              'HGF': 'b', 'HRG': 'c', 'IGF': 'y'}
    grouped = df.groupby(groups)
    all_iqrs = []
    all_pulse = []
    for conditions, group in grouped:
        iqr = group['IQRs'].as_matrix()[0]
        pulse_frac = group['Fraction of pulsing'].as_matrix()[0]
        all_iqrs.append(iqr)
        all_pulse.append(pulse_frac)
    all_iqrs = numpy.array(all_iqrs)
    all_pulse = numpy.array(all_pulse)

    times = [1, 11, 12, 13, 14, 15]

    for t in times:
        label = '%d min' % (interp_times[t])
        if t == 1:
            plt.plot(all_pulse, all_iqrs[:, t], '--', label=label)
        else:
            plt.plot(all_pulse, all_iqrs[:, t], 'o-', label=label)

    plt.legend()

    plt.xlim([0, 0.8])
    plt.ylim([0.004, 0.016])
    plt.xlabel('Fraction of pulsing cells')
    plt.ylabel('IQR of FOXO3a C/N')
    plt.legend(numpoints=1, loc='upper left')
    plt.show()

if __name__ == '__main__':
    single_cell_fname = '../rawdata/130722_SCdyn.csv'
    pulsing_fname = '../rawdata/130722_Pav.csv'
    df_pulse_frac = pandas.DataFrame.from_csv(pulsing_fname, index_col=None)

    groups = ['Cell line', 'Ligand', 'Ligand dose [ng/ml]',
                  'Drug', 'Drug dose [nM]']
    groups_single_cell = groups + ['Cell ID']
    interp_times = [0, 5, 10, 15, 20, 30, 45, 60, 90, 120, 180, 300,
                    360, 420, 480, 540, 600]
    rerun = False
    if rerun or not os.path.exists('single_intp.pkl'):
        df_single_cell = pandas.DataFrame.from_csv(single_cell_fname,
                                                   index_col=None)
        df_single_intp = interpolate_cells(df_single_cell, groups_single_cell,
                                           interp_times)
        with open('single_intp.pkl', 'wb') as fh:
            pickle.dump(df_single_intp, fh)
    else:
        with open('single_intp.pkl', 'rb') as fh:
            df_single_intp = pickle.load(fh)

    df_iqr = calculate_iqrs(df_single_intp, groups)

    df_merged = merge_iqr_pulse(df_iqr, df_pulse_frac, groups)
    df_merged.to_csv('iqr_pulse_merged.csv', index=False)
    plot_groups = ['Ligand', 'Drug', 'Ligand dose [ng/ml]']
    plot_pulsing_iqr(df_merged, plot_groups, 'EGF')
    plot_pulsing_iqr(df_merged, plot_groups, 'BTC')
