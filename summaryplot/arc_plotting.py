# ============================================================================ #
#  This script is intended to by called by                                     #
#  spt_dq_system/summaryplot/cache_archive_data.py, which provides a list of   #
#  pickle files that contain data points to be plotted, the location to        #
#  save the figures, and so on.                                                #
#                                                                              #
# ============================================================================ #


import matplotlib
matplotlib.use('agg')
matplotlib.rcParams['mathtext.default'] = 'regular'

import pickle
import datetime
import numpy
import os

from spt3g      import std_processing
from matplotlib import pyplot



def get_xlims(start_obs_id, end_obs_id):

    diff   = end_obs_id - start_obs_id
    margin = diff * 0.025
    xliml  = start_obs_id - margin
    xlimr  = end_obs_id   + margin

    return xliml, xlimr



def get_xticks(start_obs_id, end_obs_id):

    xticks_major  = []
    xticks_minor  = []
    total_seconds = end_obs_id - start_obs_id

    current_tick = start_obs_id
    end_obs_id  += 1
    while current_tick <= end_obs_id:
        xticks_major.append(current_tick)
        if total_seconds <= 2.5 * 24 * 3600:
            for hour in [1, 2, 3, 4, 5]:
                tick_minor = current_tick + hour * 3600
                if tick_minor <= end_obs_id:
                    xticks_minor.append(tick_minor)
            current_tick += 6 * 3600
        elif total_seconds <= 5.0 * 24 * 3600:
            for hour in [3, 6, 9]:
                tick_minor = current_tick + hour * 3600
                if tick_minor <= end_obs_id:
                    xticks_minor.append(tick_minor)
            current_tick += 12 * 3600
        elif total_seconds <= 10 * 24 * 3600:
            for hour in [6, 12, 18]:
                tick_minor = current_tick + hour * 3600
                if tick_minor <= end_obs_id:
                    xticks_minor.append(tick_minor)
            current_tick += 24 * 3600
        elif total_seconds <= 40 * 24 * 3600:
            for day in [1, 2, 3, 4, 5, 6]:
                tick_minor = current_tick + day * 24 * 3600
                if tick_minor <= end_obs_id:
                    xticks_minor.append(tick_minor)
            current_tick += 7 * 24 * 3600

    return xticks_major, xticks_minor



def get_xtick_labels(ticks):

    labels = []
    for tick in ticks:
        time_str = str(std_processing.obsid_to_g3time(tick))
        time_obj = datetime.datetime.strptime(
                       time_str, '%d-%b-%Y:%H:%M:%S.000000000')
        new_str  = time_obj.strftime('%m/%d\n%H:00')
        labels.append(new_str)

    return labels



def run_weather(input_files=None,
                output_dir=None,
                start_obs_id=None,
                end_obs_id=None,
                logger=None):

    logger.info('')
    logger.info('--------------------------')
    logger.info(' These are the arguments: ')
    logger.info('--------------------------')
    logger.info('')
    logger.info('Start time as in observation ID: %d', start_obs_id)
    logger.info('End   time as in observation ID: %d', end_obs_id)
    logger.info('Output directory:')
    logger.info(' %s', output_dir)
    logger.info('Input files:')
    for counter_f, f in enumerate(input_files, 1):
        logger.info(' file %04d: %s', counter_f, f)
    logger.info('')


    timestreams = {'tipper_tau'    : [],
                   'tipper_tatm'   : [],
                   'tipper_oid'    : [],
                   'weather_tatm'  : [],
                   'weather_tint'  : [],
                   'weather_humi'  : [],
                   'weather_pres'  : [],
                   'weather_wdir'  : [],
                   'weather_wmag'  : [],
                   'weather_oid'   : [],
                   'cb_rx_4k_head' : [],
                   'cb_rx_4k_sqstr': [],
                   'cb_rx_4k_plafa': [],
                   'cb_rx_50k_head': [],
                   'cb_rx_50k_harn': [],
                   'cb_ox_4k_head' : [],
                   'cb_ox_4k_lens' : [],
                   'cb_ox_4k_lyot' : [],
                   'cb_ox_50k_head': [],
                   'cb_ox_wbp_far' : [],
                   'cb_ox_50k_tube': [],
                   'cb_oid'        : [],
                   'scu_tcryoshell': [],
                   'scu_ticecrate' : [],
                   'scu_tsecondary': [],
                   'scu_tbackerack': [],
                   'scu_tcabin'    : [],
                   'scu_tundercb'  : [],
                   'scu_oid'       : []}

    logger.info('')
    logger.info('Gathering the data points to plot...')
    for counter_f, f in enumerate(input_files, 1):
        with open(f, 'rb') as fobj:
            d = pickle.load(fobj)
            for k, v in d.items():
                timestreams[k].append(v)
        logger.info(' Got the data from %s (file %04d).', f, counter_f)
    logger.info('Done.')


    logger.info('')
    logger.info('Converting some units...')
    conversion_methods = {'tipper_tau'    : lambda a: a,
                          'tipper_tatm'   : lambda a: a,
                          'weather_tatm'  : lambda a: a + 273.15,
                          'weather_tint'  : lambda a: a + 273.15,
                          'weather_humi'  : lambda a: a,
                          'weather_pres'  : lambda a: a * 9.8692e-4,
                          'weather_wdir'  : lambda a: a,
                          'weather_wmag'  : lambda a: a * 1.9438,
                          'cb_rx_4k_head' : lambda a: a,
                          'cb_rx_4k_sqstr': lambda a: a,
                          'cb_rx_4k_plafa': lambda a: a,
                          'cb_rx_50k_head': lambda a: a,
                          'cb_rx_50k_harn': lambda a: a,
                          'cb_ox_4k_head' : lambda a: a,
                          'cb_ox_4k_lens' : lambda a: a,
                          'cb_ox_4k_lyot' : lambda a: a,
                          'cb_ox_50k_head': lambda a: a,
                          'cb_ox_wbp_far' : lambda a: a,
                          'cb_ox_50k_tube': lambda a: a,
                          'scu_tcryoshell': lambda a: a,
                          'scu_ticecrate' : lambda a: a,
                          'scu_tsecondary': lambda a: a,
                          'scu_tbackerack': lambda a: a,
                          'scu_tcabin'    : lambda a: a,
                          'scu_tundercb'  : lambda a: a}

    all_okay_indices = {}
    for k, ts in timestreams.items():
        if '_oid' in k:
            a = numpy.asarray(ts)
            okay_indices = numpy.where((start_obs_id<=a)&(a<=end_obs_id))
            timestreams[k] = a[okay_indices]
            all_okay_indices[k.split('_')[0]] = okay_indices

    for k, ts in timestreams.items():
        if '_oid' not in k:
            okay_indices   = all_okay_indices[k.split('_')[0]]
            converted_ts   = conversion_methods[k](numpy.asarray(ts))
            timestreams[k] = converted_ts[okay_indices]
    logger.info('Done.')


    logger.info('')
    logger.info('Deciding on the appearance of each plot...')
    xliml, xlimr = get_xlims(start_obs_id, end_obs_id)
    xticks_major, xticks_minor = get_xticks(start_obs_id, end_obs_id)
    xtick_labels = get_xtick_labels(xticks_major)

    plot_groups = {'tipper_tau'  : ['tipper_tau'  ],
                   'tipper_tatm' : ['tipper_tatm' ],
                   'weather_tatm': ['weather_tatm'],
                   'weather_humi': ['weather_humi'],
                   'weather_pres': ['weather_pres'],
                   'weather_wdir': ['weather_wdir'],
                   'weather_wmag': ['weather_wmag'],
                   'cabin_temps' : ['scu_tcryoshell',
                                    'scu_ticecrate' ,
                                    'scu_tsecondary',
                                    'scu_tbackerack',
                                    'scu_tcabin'    ,
                                    'scu_tundercb'  ],
                   'fourk_temps' : ['cb_rx_4k_head' ,
                                    'cb_rx_4k_sqstr',
                                    'cb_rx_4k_plafa',
                                    'cb_ox_4k_head' ,
                                    'cb_ox_4k_lens' ,
                                    'cb_ox_4k_lyot' ],
                   'fiftyk_temps': ['cb_rx_50k_head',
                                    'cb_rx_50k_harn',
                                    'cb_ox_50k_head',
                                    'cb_ox_wbp_far' ,
                                    'cb_ox_50k_tube']}

    legends = {'tipper_tau'    : '',
               'tipper_tatm'   : '',
               'weather_tatm'  : '',
               'weather_humi'  : '',
               'weather_pres'  : '',
               'weather_wdir'  : '',
               'weather_wmag'  : '',
               'cb_rx_4k_head' : 'Receiver 4K head',
               'cb_rx_4k_sqstr': 'Receiver 4K SQUID strap',
               'cb_rx_4k_plafa': 'Receiver 4K plate far side',
               'cb_ox_4k_head' : 'Optics 4K head',
               'cb_ox_4k_lens' : 'Optics 4K lens tab',
               'cb_ox_4k_lyot' : 'Optics 4K Lyot',
               'cb_rx_50k_head': 'Receiver 50K head',
               'cb_rx_50k_harn': 'Receiver 50K harness',
               'cb_ox_50k_head': 'Optics 50K head',
               'cb_ox_wbp_far' : 'WBP far side',
               'cb_ox_50k_tube': 'Optics 50K tube',
               'scu_tcryoshell': 'Cryostat shell',
               'scu_ticecrate' : 'Above icecrate 015',
               'scu_tsecondary': 'Behind secondary',
               'scu_tbackerack': 'Back of e-rack',
               'scu_tcabin'    : 'Near cabin thermostat',
               'scu_tundercb'  : 'Under cryoboard crate'}

    yticks_major = {'tipper_tau'  : numpy.linspace(1.0, 5.0, 5),
                    'tipper_tatm' : numpy.linspace(195, 255, 7),
                    'weather_tatm': numpy.linspace(195, 255, 7),
                    'weather_humi': numpy.linspace(0.10, 0.80, 6),
                    'weather_pres': numpy.linspace(0.60, 0.80, 5),
                    'weather_wdir': numpy.linspace(0, 360, 7),
                    'weather_wmag': numpy.linspace(0, 40, 5),
                    'cabin_temps' : numpy.linspace(-20, 20, 5),
                    'fourk_temps' : numpy.linspace(2.0, 6.0, 5),
                    'fiftyk_temps': numpy.linspace(30, 70, 5)}

    ylims = {}
    for quantity, ticklist in yticks_major.items():
        vmax   = numpy.max(ticklist)
        vmin   = numpy.min(ticklist)
        margin = (vmax - vmin) * 0.05
        ylims[quantity] = (vmin-margin, vmax+margin)

    yticks_minor = {}
    for quantity, ticklist_major in yticks_major.items():
        delta = ticklist_major[1] - ticklist_major[0]
        vmax  = numpy.max(ticklist_major)
        vmin  = numpy.min(ticklist_major)
        ntck  = len(ticklist_major)
        ticklist_minor = numpy.linspace(vmin, vmax, 4*(ntck-1)+1)
        ticklist_minor = [t for t in ticklist_minor if t not in ticklist_major]
        yticks_minor[quantity] = ticklist_minor

    ylabels = {'tipper_tau'  : '{:s} [{:s}]'.format(
                                   r'$\tau$', r'${airmass}^{-1}$'),
               'tipper_tatm' : '{:s} [K]'.format(r'$T_{atm}$'),
               'weather_tatm': 'Air temperature [K]',
               'weather_humi': 'Relative humidity [%]',
               'weather_pres': 'Atmospheric pressure [atm]',
               'weather_wdir': 'Wind direction [degree]',
               'weather_wmag': 'Wind speed [Knot]',
               'cabin_temps' : 'Temperature [{:s}]'.format(r'${}^{\circ} C$'),
               'fourk_temps' : 'Temperature [K]',
               'fiftyk_temps': 'Temperature [K]'}

    titles = {'tipper_tau'  : '{:s} ({:s} {:s}) {:s}\n{:s}'.format(
                                  'Atmospheric opacity',
                                  'zenith optical depth',
                                  r'$\tau$',
                                  '@ 350 micron',
                                  'obtained from the tipper measurements'),
              'tipper_tatm' : '{:s} ({:s})\n{:s}'.format(
                                  'Effective atmospheric brightness',
                                  r'$T_{atm}$',
                                  'obtained from the tipper measurements'),
              'weather_tatm': 'Air temperature around the weather station',
              'weather_humi': 'Relative humidity around the weather station',
              'weather_pres': 'Atmospheric pressure around the weather station',
              'weather_wdir': 'Azimuth from which wind is blowing',
              'weather_wmag': 'Wind speed measured by the weather station',
              'cabin_temps' : 'Temperature sensed by cabin thermometers',
              'fourk_temps' : 'Temperature of some parts nominally at 4 [K]',
              'fiftyk_temps': 'Temperature of some parts nominally at 50 [K]'}
    logger.info('Done.')


    logger.info('')
    logger.info('Starting to plot data points...')
    for group, quantities in plot_groups.items():
        logger.info(' Plotting data related to %s...', group)

        figure = pyplot.figure(figsize=(12, 8))
        plot   = figure.add_subplot(111)

        not_many_outliers   = True
        planned_ylim_bottom = ylims[group][0]
        planned_ylim_top    = ylims[group][1]

        for quantity in quantities:
            data_type   = quantity.split('_')[0]
            time_stamps = timestreams[data_type+'_oid']
            time_series = timestreams[quantity]
 
            quarter_N = len(time_stamps) // 4
            if len(numpy.where(time_series>planned_ylim_top)[0]) > quarter_N:
                not_many_outliers = False
            if len(numpy.where(time_series<planned_ylim_bottom)[0]) > quarter_N:
                not_many_outliers = False

            plot.plot(time_stamps, time_series, label=legends[quantity],
                      linestyle='solid', linewidth=1.5, alpha=0.50)

        plot.set_xlim(left=xliml, right=xlimr)

        plot.tick_params(which='both', direction='in', pad=5.0,
                         labelsize=20, labelcolor='black')
        plot.set_xticks(xticks_major, minor=False)
        plot.set_xticks(xticks_minor, minor=True)
        plot.set_xticklabels(xtick_labels, rotation=-45)

        if not_many_outliers:
            plot.set_ylim(bottom=planned_ylim_bottom, top=planned_ylim_top)
            plot.set_yticks(yticks_major[group], minor=False)
            plot.set_yticks(yticks_minor[group], minor=True)

        plot.grid(which='major', linestyle='dotted', linewidth=1.0)
        plot.grid(which='minor', linestyle='dotted', linewidth=0.5)

        if len(quantities) > 1:
            plot.legend(loc='upper right', fontsize=18, framealpha=0.1)

        plot.set_xlabel('\nTime (UTC)', fontsize=22)
        plot.set_ylabel(ylabels[group]+'\n', fontsize=22)

        if not_many_outliers:
            full_title = titles[group]
        else:
            full_title  = titles[group]
            full_title += '\n(Usual y axis limits not used '
            full_title += 'because of the many outliers!)'
        plot.set_title(full_title+'\n', fontsize=24)

        figure.tight_layout(pad=1.0)
        figure.savefig(os.path.join(output_dir, group+'.png'))
        pyplot.close(figure)

    logger.info('')
    logger.info('All plots have been made!')
    logger.info('')



def run_fridge_cycles(input_file=None, output_dir=None, logger=None):

    logger.info('')
    logger.info('Loading the data...')

    with open(input_file, 'rb') as fobj:
        data = pickle.load(fobj)

    oids = data['OID']
    beg  = std_processing.time_to_obsid(data['INFO']['start'].\
               strftime('20%y%m%d_%H%M%S'))
    end  = std_processing.time_to_obsid(data['INFO']['stop'].\
               strftime('20%y%m%d_%H%M%S'))
    end  = (end  - beg) / 60.0
    mins = (oids - beg) / 60.0

    common_title  = 'Temperature of various thermometers '
    common_title += 'during cycle_tune({})\n'.format(data['INFO']['args'])
    common_title += '(From {} till {})'.format(
                        data['INFO']['start'].strftime('20%y%m%d_%H%M%S'),
                        data['INFO']['stop'].strftime('20%y%m%d_%H%M%S'))

    logger.info('Making the first plot...')

    figure = pyplot.figure(figsize=(12, 6))
    plot   = figure.add_subplot(1, 1, 1)

    for thermometer in ['UCHEAD',  'ICHEAD', 'HE4HEAD', 'HE4FB',
                        'HE4PUMP', 'ICPUMP', 'UCPUMP',
                        'HE4SW',   'ICSW',   'UCSW']:
        plot.plot(mins, data[thermometer], label=thermometer,
                  linewidth=1.5, alpha=0.75)

    plot.set_xlim(left=-15, right=375)
    plot.set_ylim(bottom=0.2, top=70)
    plot.set_yscale('log')
    xticks_major = numpy.arange(7) * 60
    xticks_minor = (numpy.arange(39) - 1) * 10
    xticks_minor = xticks_minor[xticks_minor%60!=0]
    plot.set_xticks(xticks_major, minor=False)
    plot.set_xticks(xticks_minor, minor=True)
    plot.tick_params(which='both', direction='in', pad=6.0,
                     labelsize=12, labelcolor='black')
    plot.grid(axis='both', which='major', linestyle='dotted', linewidth=1.0)
    plot.grid(axis='both', which='minor', linestyle='dotted', linewidth=0.5)
    plot.set_xlabel('\nTime since the beginning of the cycle [minutes]',
                    fontsize=14)
    plot.set_ylabel('Temperature [K]\n', fontsize=14)
    plot.axvline(end, color='black', linestyle='dashed', linewidth=1.25,
                 label='END')
    plot.legend(loc='upper right', fontsize=13)
    sub_title = '  -  Plot 1: He10 fridge parts, entire period\n'
    plot.set_title(common_title+sub_title, fontsize=15)
    figure.tight_layout()
    figure.savefig('{}/he10_full.png'.format(output_dir),
                   bbox_inches='tight')
    pyplot.close(figure)


    logger.info('Making the second plot...')

    figure = pyplot.figure(figsize=(12, 6))
    plot   = figure.add_subplot(1, 1, 1)

    for thermometer in ['UCHEAD',  'ICHEAD', 'HE4HEAD', 'HE4FB',
                        'HE4PUMP', 'ICPUMP', 'UCPUMP',
                        'HE4SW',   'ICSW',   'UCSW']:
        plot.plot(mins, data[thermometer], label=thermometer,
                  linewidth=1.5, alpha=0.75)

    plot.set_xlim(left=-15, right=175)
    plot.set_ylim(bottom=0.2, top=70)
    plot.set_yscale('log')
    xticks_major = numpy.arange(9) * 20
    xticks_minor = (numpy.arange(36) - 1) * 5
    xticks_minor = xticks_minor[xticks_minor%20!=0]
    plot.set_xticks(xticks_major, minor=False)
    plot.set_xticks(xticks_minor, minor=True)
    plot.tick_params(which='both', direction='in', pad=6.0,
                         labelsize=12, labelcolor='black')
    plot.grid(axis='both', which='major', linestyle='dotted', linewidth=1.0)
    plot.grid(axis='both', which='minor', linestyle='dotted', linewidth=0.5)
    plot.set_xlabel('\nTime since the beginning of the cycle [minutes]',
                    fontsize=14)
    plot.set_ylabel('Temperature [K]\n', fontsize=14)
    plot.legend(loc='upper right', fontsize=13)
    sub_title = '  -  Plot 2: He10 fridge parts, 1st half\n'
    plot.set_title(common_title+sub_title, fontsize=15)
    figure.tight_layout()
    figure.savefig('{}/he10_half.png'.format(output_dir),
                   bbox_inches='tight')


    logger.info('Making the third plot...')

    figure = pyplot.figure(figsize=(12, 6))
    plot   = figure.add_subplot(1, 1, 1)

    for thermometer in ['UCHEAD', 'ICHEAD', 'UC_STAGE', 'IC_STAGE', 'LC_TOWER']:
        plot.plot(mins, data[thermometer], label=thermometer,
                  linewidth=1.5, alpha=0.75)

    plot.set_xlim(left=-15, right=285)
    plot.set_ylim(bottom=0.2, top=3.1)
    xticks_major = numpy.arange(10) * 30
    xticks_minor = (numpy.arange(31) - 1) * 10
    xticks_minor = xticks_minor[xticks_minor%30!=0]
    plot.set_xticks(xticks_major, minor=False)
    plot.set_xticks(xticks_minor, minor=True)
    yticks_major = (numpy.arange(10) + 1) * 0.3
    yticks_minor = (numpy.arange(30)) * 0.1 + 0.2
    yticks_minor = yticks_minor[yticks_minor%0.3!=0.0]
    plot.set_yticks(yticks_major, minor=False)
    plot.set_yticks(yticks_minor, minor=True)
    plot.tick_params(which='both', direction='in', pad=6.0,
                         labelsize=12, labelcolor='black')
    plot.grid(axis='both', which='major', linestyle='dotted', linewidth=1.0)
    plot.grid(axis='both', which='minor', linestyle='dotted', linewidth=0.5)
    plot.set_xlabel('\nTime since the beginning of the cycle [minutes]',
                    fontsize=14)
    plot.set_ylabel('Temperature [K]\n', fontsize=14)
    plot.legend(loc='upper right', fontsize=13)
    sub_title = '  -  Plot 3: Low temperature parts\n'
    plot.set_title(common_title+sub_title, fontsize=15)
    figure.tight_layout()
    figure.savefig('{}/low_temps.png'.format(output_dir),
                   bbox_inches='tight')
    pyplot.show(plot)

    logger.info('Done.')
    logger.info('')

