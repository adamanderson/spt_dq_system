This project is a website for data quality monitoring of SPT-3G in the north and at pole.

# Installation
1.) In order to run the data quality page locally, the only dependencies you will need that are not installed automatically in a later step are node.js and its package manager npm. The data quality system has been tested with node.js versions >=7.10.1. These packages are already installed on anal at pole and amundsen/scott in the north. For work on other machines, you can install node.js from its website: [https://nodejs.org/en/](https://nodejs.org/en/)

2.) The python components of the data quality system require python3. The easiest way to enforce the use of python3 on SPT machines is to use run the clustertools environment setup script:
```
eval `/cvmfs/spt.opensciencegrid.org/py3-v1/setup.sh`
```

For instructions on installing clustertools on other machines, refer to the clustertools github page: [https://github.com/SouthPoleTelescope/clustertools](https://github.com/SouthPoleTelescope/clustertools)

3.) After cloning the repository, install the javascript dependencies via npm:
```
cd spt_dq_system/
npm install
```

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
* `min_time_static_plots`: Date string of format `YYYYMMDD` at which to start static plot generation. 

5.) If you do not want password-protected login, then skip this step. Otherwise, run the following in a terminal:

```
openssl genrsa -des3 -passout pass:x -out server.pass.key 2048                
openssl rsa -passin pass:x -in server.pass.key -out key.pem                   
rm server.pass.key                                                            
openssl req -new -key key.pem -out server.csr                                 
openssl x509 -req -sha256 -days 365 -in server.csr -signkey key.pem -out cert.pem                                                                             
rm server.csr
```

Then create a hash for the password, with the commands below. Note that the login username is fixed to be 'spt'. Change 'your_password' below to the password you want to use.

```
echo "var bcrypt = require('bcryptjs');var salt = bcrypt.genSaltSync(10);var hash = bcrypt.hashSync('your_password', salt);console.log(hash);" > temp.js
node temp.js > hash
rm temp.js
```

6.) Add `spt_dq_system` to your `PYTHONPATH`:
```
export PYTHONPATH=$PYTHONPATH:/path/to/spt_dq_system/
```

7.) Launch the server with:
```
node db_server.js
```

8.) Visit the page. If your port is 3002, and you are using password protection and scott, for example, then go to:
```
https://scott.grid.uchicago.edu:3002/
```

Without password-protection, go to:
```
http://scott.grid.uchicago.edu:3002/
```


# Documentation
For a thorough description of how to use and extend the website, please visit the wiki page at https://pole.uchicago.edu/spt3g/index.php/SPT_Data_Quality

