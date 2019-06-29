This project is a website for data quality monitoring of SPT-3G in the north and at pole.

# Installation
1.) In order to run the data quality page locally, the only dependencies you will need that are not installed automatically in a later step are node.js and its package manager npm. The data quality system has been tested with node.js versions >=7.10.1. These packages are already installed on anal at pole and amundsen/scott in the north. For work on other machines, you can install node.js from its website: [https://nodejs.org/en/](https://nodejs.org/en/)

2.) The python components of the data quality system require python3. The easiest way to enforce the use of python3 on SPT machines is to use the clustertools system install at pole or on scott/amundsen in the north. For instructions on installing clustertools on other machines, refer to the clustertools github page: [https://github.com/SouthPoleTelescope/clustertools](https://github.com/SouthPoleTelescope/clustertools)

3.) After cloning the repository, run the setup script to get server- and client-side dependencies, and to build a local version of `spt3g_software`. Select `pole` or `north` depending on the install location.
```
./setup.sh [pole/north]
```
The rationale behind building a local version of `spt3g_software` is that user versions of the repo are subject to change and we would like the copy used by the server to remain unchanged unless deliberately updated.

4.) Create a `config.yaml` file. In practice, the easiest way to do this is to copy one of two template files, `config.yaml.template_north` and `config.yaml.template_pole`, intended as examples for running and pole and in the north, and to modify it to have paths that suit your needs. The fields have the following meanings:

* `scanify_db_path`: Path to SQLite version of scanify database.
* `auxtransfer_db_path`: Path to SQLite version of aux files transfer database.
* `autoproc_db_path`: Path to SQLite version of autoprocessing status database.
* `python_location`: Path to version of python to use for plot generation. Must be the same version that `spt3g_software` is compiled against.
* `key_file`: `.key` file used for password-protected login. See step 4 for instructions to create this. Leaving this field and `cert_file` blank will default to no password protection.
* `cert_file`: `.pem` certificate file used for password-protected login. See step 4 for instructions to create this. Leaving this field and `key_file` blank will default to no password protection.
* `plot_dir`: Directory where plots are saved. Must have write access to this directory.
* `static_plot_dir`: Directory where "static" plots used for summary page are saved. Must have write access to this directory.
* `calib_data_dir`: Location of autoprocessed calibration data.
* `bolo_data_dir`: Location of scanified data.
* `port`: Port from which website will be accessible. This must be unused; note that production versions of the data quality system run on port 3000 at pole and in the north, and a development version is often running on port 3001, so you may want to choose a port different from these.
* `min_time_summary_plots`: Date string of format `YYYYMMDD` denoting time after which to generate summary plots.
* `min_time_maps`: Date string of format `YYYYMMDD` denoting time after which to generate map coadds and plots..

5.) Launch the server. If you are certain that you have the correct environment set up, this can be done simply with:
```
node db_server.js
```
To ensure that the proper environment is being used, it is better to launch the server with the following:
```
# north
/bin/bash -c 'unset SPT3G_SOFTWARE_PATH; unset SPT3G_SOFTWARE_BUILD_PATH; unset SPT3G_BUILD_ROOT; unset PATH; unset LD_LIBRARY_PATH; unset PYTHONPATH; eval `/cvmfs/spt.opensciencegrid.org/py3-v3/setup.sh` /PATH/TO/YOUR/INSTALL/spt_dq_system/spt3g_software/build/env-shell.sh /PATH/TO/YOUR/INSTALL/spt_dq_system/db_server.js'

# pole
/bin/bash -c 'unset SPT3G_SOFTWARE_PATH; unset SPT3G_SOFTWARE_BUILD_PATH; unset SPT3G_BUILD_ROOT; unset PATH; unset LD_LIBRARY_PATH; unset PYTHONPATH; eval `/software/clustertools/py3-v3/setup.sh` /PATH/TO/YOUR/INSTALL/spt_dq_system/spt3g_software/build/env-shell.sh /PATH/TO/YOUR/INSTALL/spt_dq_system/db_server.js'
```
To avoid conflicts with environment variables set by other copies of `spt3g_software`, this invocation unsets all variables used by other copies, and runs using the environment of the version built inside the `spt_dq_server` directory.

6.) Visit the page. If your port is 3002, and you are using password protection and scott, for example, then go to:
```
http://scott.grid.uchicago.edu:3002/
```


# Documentation
For a thorough description of how to use and extend the website, please visit the wiki page at https://pole.uchicago.edu/spt3g/index.php/SPT_Data_Quality

