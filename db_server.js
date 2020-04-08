#!/usr/bin/env node

var express = require('express');
var sqlite3 = require('sqlite3');
var assert = require('assert');
var bodyParser = require('body-parser');
var squel = require("squel");
var execFile = require('child_process').execFile;
var exec = require('child_process').exec;
var fs = require('fs');
var readline = require('readline');
var util = require('util');
var auth = require('express-basic-auth');
var moment = require('moment');
var SSE = require('express-sse');
var sse = new SSE();
var yaml = require('yamljs');
var serveIndex = require('serve-index')
var path = require('path');

// TODO error handling

// load configuration info from yaml file
config = yaml.load('config.yaml')

// Hack environment variables if needed, in case they aren't pulled from the
// parent bash shell. Seems to be due to a bug in node.js as far as I can tell,
// although it would be nice to understand this better. It only seems to be a
// problem on spt3g-dq.grid.uchicago.edu.
if(config.ld_library_path != null) {
	process.env.LD_LIBRARY_PATH = config.ld_library_path
}
if(config.pythonpath != null) {
	process.env.PYTHONPATH = config.pythonpath
}

// setup express
var app = express();
app.use(bodyParser.json());
app.use(bodyParser.urlencoded({extended: true}));
// server side events for passing messages
app.get('/sse', sse.init);
// counter to separate messages by plot request
var sseid = 0;

// password authentication
app.use(auth({
    users: {'spt': 'sptjamz'},
    challenge: true
}));

// database filenames
db_files = {scanify: config.scanify_db_path,
	        aux_transfer: config.auxtransfer_db_path,
	        autoproc: config.autoproc_db_path,
            schedule: config.schedule_db_path}

// setup home page
app.use(express.static('public'));
app.get('/index.html', function (req, res) {
  res.sendFile( path.join(__dirname, "index.html") );
})
app.get('/', function (req, res) {
  res.sendFile( path.join(__dirname, "index.html") );
})

// needed to load js and css files
app.use('/js', express.static( path.join(__dirname, '/js') ));
app.use('/css', express.static( path.join(__dirname, '/css') ));

// file browser for log files
if(config.site == 'pole') {
  app.use('/scanify_logs', express.static(config.scanify_log_path),
  		  serveIndex(config.scanify_log_path, {'icons': true}));
}
app.use('/autoproc_logs', express.static(config.autoproc_log_path),
		serveIndex(config.autoproc_log_path, {'icons': true}));

// create directory in /tmp to store plots
// TODO: don't cache every plot. Remove old plots from time to time
// NOTE: not a high priority because right now ~1.2TB free in /tmp
if (!fs.existsSync(config.plot_dir)) {
      fs.mkdirSync(config.plot_dir);
}
app.use('/img', express.static(config.plot_dir));


// serve directory for static plots
if (!fs.existsSync(config.static_plot_dir)) {
    fs.mkdirSync(config.static_plot_dir);
}
app.use('/staticimg', express.static(config.static_plot_dir));
app.use('/readmeimg', express.static(path.join(__dirname, 'READMEfigs')));


// get subdirectory information for previous plots
app.get('/staticdirs', function(req, res) {
	// just do an ls on the plots directory to figure out the other
	// subdirectories of plots
    var interval_path = path.join(config.static_plot_dir,
                                  req.query.subdirectory,
                                  req.query.interval);

    fs.readdir(interval_path, function(err, filelist) {
        // treat file not found like an empty file
        if(err && err.code == "ENOENT") {
            res.send([]);
            return;
        }

        var dirlist = [];
        // get the list of directories that contain plots
	    for (jfile=0; jfile<filelist.length; jfile++) {
            var test_path = path.join(interval_path,
                                     filelist[jfile]);

		    if(fs.readdirSync(test_path).length)
			    dirlist.push( path.join(req.query.subdirectory,  req.query.interval, filelist[jfile]) );
	    }
	    res.send(dirlist);
    });
});

// get time that calibration plots were lastmodified
app.get('/lastmodified_calibration', function(req, res) {
	var stat_dir = path.join(config.static_plot_dir, '/plots/last_n/last_03/');

    fs.readdir(stat_dir, function(err, filelist) {
        lastupdate_time = check_lastmodified(stat_dir, filelist);
	    res.send({time:lastupdate_time});
    });
});

// get time that winter fields related plots  were lastmodified
app.get('/lastmodified_maps_winter', function(req, res) {
	var stat_dir = path.join(config.coadds_figs_dir, '/last_n/last_07/');

    fs.readdir(stat_dir, function(err, filelist) {
        lastupdate_time = check_lastmodified(stat_dir, filelist);
	    res.send({time:lastupdate_time});
    });
});

// get time that summer fields related plots  were lastmodified
app.get('/lastmodified_maps_summer', function(req, res) {
    var stat_dir = path.join(config.coadds_figs_summer_dir, '/last_n/last_07/');

    fs.readdir(stat_dir, function(err, filelist) {
        lastupdate_time = check_lastmodified(stat_dir, filelist);
	    res.send({time:lastupdate_time});
    });
});

// get time that weather plots were lastmodified
app.get('/lastmodified_weather', function(req, res) {
	var stat_dir = path.join(config.arc_figs_dir, '/last_n/last_03/');

    fs.readdir(stat_dir, function(err, filelist) {
        lastupdate_time = check_lastmodified(stat_dir, filelist);
	    res.send({time:lastupdate_time});
    });
});


// check the last-modified property for a list of files and return the latest one
function check_lastmodified(stat_dir, filelist) {
    var lastupdate_time;
    for(var jfile=0; jfile < filelist.length; jfile++) {
		file_info = fs.statSync( path.join(stat_dir, filelist[jfile]) );
		if(typeof lastupdate_time !== 'undefined') {
			if(file_info.mtime > lastupdate_time)
				lastupdate_time = file_info.mtime;
		}
		else
			lastupdate_time = file_info.mtime;
	}
    return lastupdate_time;
}



// logs messages along with a timestamp. keeps writesteam open
var ws = fs.createWriteStream('./db.log', {'flags': 'a'});
ws.on('error', function (err) {
  console.error('Error logging message');
  console.error(err);
});
function log(msg) {
    var d = new Date();
    ws.write(d.toLocaleString() + '\t' + msg + '\n');
}


app.get('/dbpage', function(req, res) {
	// open the database
	db = new sqlite3.Database(db_files[req.query.dbname],
				  sqlite3.OPEN_READONLY);

	// If the autoprocessing database was requested, then merge with the
	// scanify database so that we can get the date of the observations.
	// Note that we must wrap this in a 'serialize' call because otherwise
	// subsequent queries will execute before the scanify database is
	// attached.
	if(req.query.dbname == "autoproc") {
	    db2 = new sqlite3.Database(db_files["scanify"],
				       sqlite3.OPEN_READONLY);
	    db.serialize(function() {
		    db.run("ATTACH \""+db2.filename+"\" as scanify");
		});
	}


	// get data from the database
	query = squel.select().from(req.query.dbname);
	parseSearch(query, req.query, req.query.dbname);
	db.all(query.toParam()['text'],
	       query.toParam()['values'],
	       function(err, rows) {
			   if(req.query.dbname == 'scanify') {
				   // append log file information to database results
				   for(var jrow in rows) {
					   rows[jrow]['log_file'] = path.join('scanify_logs', rows[jrow]['source'], rows[jrow]['observation']);
				   }
				   
				   // interpret transfer_fullrate and transfer_downsampled
				   // integer flags as booleans for simpler parsing by the 
				   // transfer status table
				   for(var jrow in rows) {
					   if(rows[jrow]['transfer_fullrate'] > 0)
						   rows[jrow]['transfer_fullrate'] = true;
					   else
						   rows[jrow]['transfer_fullrate'] = false;
					   if(rows[jrow]['transfer_downsampled'] > 0)
						   rows[jrow]['transfer_downsampled'] = true;
					   else
						   rows[jrow]['transfer_downsampled'] = false;
				   }
			   }
			   if(req.query.dbname == 'autoproc') {
				   // append log file information to database results
				   for(var jrow in rows) {
					   rows[jrow]['log_file'] = path.join('autoproc_logs/', rows[jrow]['source'], rows[jrow]['observation']);
				   }
			   }
			   res.send(rows);
	       });

	// close the database
	db.close();
    });


// page for displaying plots/data
app.get('/display.html', function(req, res) {
  res.sendFile( path.join(__dirname, "display.html") );
});

// summary page
app.get('/summary.html', function(req, res) {
  res.sendFile( path.join(__dirname, "summary.html") );
});
//monthly summary
app.get('/monthly.html', function(req, res) {
  res.sendFile( path.join(__dirname, "monthly.html") );
});

// request new sseid
app.get('/sseid', function(req, res) {
  res.json(sseid);
  sseid++;
});

// used to make plotting requests. User requests a specific plot(s) and
// the server executes a plotting script that saves the plot in the img
// directory. The request times out after 20s in case of a bad script or
// the server is being slow.
app.get('/data_req', function (req, res) {
  // save id to send messages to
  var id = req.query['sseid'];
  options = {'timeout':60000};
  tab = req.query['table'];
  if (tab == 'scanify' || tab =='autoproc') {
    obs = req.query['observation'];
    source = req.query['source'];
  }
  else if (tab == 'aux') {
    filename = req.query['filename'];
  }
  plot_type = req.query['plot_type'];

  log(util.inspect(req.query));
  // execute python plotting script. Safe because user input
  // is passed as arguments to the python script and the python
  // script handles the arguments safely.
  var args;
  if (tab == 'scanify' || tab == 'autoproc')
    args = ['-B', './plot/_plot.py',
	    source,
	    config.plot_dir,
	    config.calib_data_dir,
	    config.bolo_data_dir,
	    obs].concat(
        plot_type.split(' '));
  else if (tab == 'aux')
    args = ['-B', './plot/_plot.py', 'aux', filename].concat(
        plot_type.split(' '));
  var err = null;
  var child = execFile(config.python_location, args, options);

  child.stdout.on('data', function(data) {
    // sometimes stdout combines messages so split them up and send them
    // individually
    var msgs = data.split('\n');

    // remove empty strings
    msgs = msgs.filter(entry => entry.trim() != '');

    for (var i = 0; i < msgs.length; i++)
      sse.send(msgs[i], 'out' + id);
  });

  child.stderr.on('data', function(data) {
    // log error messages
    log(data)
    sse.send(data, 'err' + id);
  });

  // send a blank response
  res.send("0");
  return;
})

// get all available plot types, removing the driver file and .py
app.get('/plot_list', function(req, res) {
  if (req.query['type'] == '')
    req.query['type'] = 'any';
  var json = JSON.parse(fs.readFileSync('./plot/plot_config.json', 'utf8'));
  // check if plot_config has an entry for this source. Otherwise return 'any'
  if (json[req.query['tab']][req.query['type']] == null)
    req.query['type'] = 'any';
  res.json(json[req.query['tab']][req.query['type']]);
});

// get list of available sources
app.get('/sourcelist', function(req, res) {
	db = new sqlite3.Database(config.scanify_db_path,
				  sqlite3.OPEN_READONLY);
	query = squel.select()
	    .from('scanify')
	    .field('source')
	    .distinct();
	db.all(query.toString(), function(err, rows) {
		res.send(rows);
		    });
    });

// get the site location: north/pole
app.get('/site', function(req, res) {
	res.send(config.site)
});

function parseSearch(query, searchJSON, tab) {
    // special option to just get the single most recent obserrvation
    if (searchJSON.search.hasOwnProperty('date') &&
	searchJSON.search['date'] == 'latest') {
	query.where(tab + '.date = ?', squel.select()
		    .field("MAX(" + tab + '.date)')
		    .from(tab)
		    .where(tab+'.source == ?', searchJSON.search['source']));
	return query;
    }

    for(var column in searchJSON.search) {
	// special handling for ranges of dates
	if(column == 'date' || column == 'modified' || column == 'sch_start' || column == 'sch_stop') {
	    // otherwise assume that some date range is specified
	    var min_time = searchJSON.search[column]['min'];
	    var max_time = searchJSON.search[column]['max'];
	    if(!min_time)
		min_time = "2017-01-01";
	    if(!max_time)
		max_time = moment().format('YYYY-MM-DD');
	    // special handling to deal with joined 'autoproc' and 'scanify' databases
	    if(tab == 'autoproc' && column == 'date')
		query.where("date(scanify."+column+") BETWEEN date(?) AND date(?)", min_time, max_time);
	    else
		query.where("date("+tab+'.'+column+") BETWEEN date(?) AND date(?)", min_time, max_time);
	}
	else {	    
	    // handle other columns that span ranges
	    if(searchJSON.search[column].hasOwnProperty('min') && searchJSON.search[column]['min']!='')
		query.where(tab+'.'+column + ' >= ?', searchJSON.search[column]['min']);
	    if(searchJSON.search[column].hasOwnProperty('max') && searchJSON.search[column]['max']!='')
		query.where(tab+'.'+column + ' <= ?', searchJSON.search[column]['max']);
	    
	    // special handling for filenames or other things with wildcards
	    if(column == 'filename' && searchJSON.search['filename'])
		query.where("filename LIKE ?", '%' + searchJSON.search['filename'] + '%');
	
	    // if neither min nor max is present assume we are looking for an exact
	    // match on this column
	    if(searchJSON.search[column].hasOwnProperty('min')==false && 
	       searchJSON.search[column].hasOwnProperty('max')==false &&
	       searchJSON.search[column] != '' && 
	       column != 'filename')
		query.where(tab+'.'+column + " == ?", searchJSON.search[column]);
	}
    }

    // joining commands
    if(tab == 'autoproc') {
	query.field('autoproc.*').field('scanify.date');
     	query.join('scanify', null, squel.expr().and("autoproc.observation == scanify.observation"));
    }

    // sorting commands
    var sort = searchJSON['sort'];
    var sort_dir = searchJSON['sort_dir'];
    if (sort) {
	if (sort_dir == 'desc')
	    query.order(sort, false);
	else
	    query.order(sort, true);
    }
    else if (tab == 'scanify' || tab == 'aux')
	query.order('scanify.date', false);
    else if (tab == "autoproc")
	query.order('autoproc.modified', false);
    
    return query;
}

app.listen(parseInt(config.port))


// static timeseries plots update
is_update_running = false;
is_winter_map_update_running = false;
is_summer_map_update_running = false;
is_arc_update_running = false;
var child;
var child_winter_maps;
var child_summer_maps;
var child_arcs;

function updateSummaryPlots() {
    args = ['-B',
            'update_summary.py',
            'summarystats',
            config.static_plot_dir,
            config.calib_data_dir,
            config.bolo_data_dir,
            config.min_time_summary_plots]

    if(config.no_data_update) {
        args.push('--no-data-update')
    }

    if(is_update_running == false) {
        is_update_running = true;
        // update data skims
        child = execFile(config.python_location,
                         args,
                         {maxBuffer: 1024*1024*8},
                         function(err) {
                             console.log(err);
                             console.log('Finished updating data skims and plots.');
                             is_update_running = false;
                         });
        console.log('Updating summary plots...');
    }
    else {
        console.log('Summary plot updater already running, so not spawning again!');
    }
}

function updateWinterMapPlots() {
    args = ['-B',
            'update_summary.py',
            'maps',
            'winter',
            config.maps_data_dir,
            config.coadds_data_dir,
            config.coadds_figs_dir,
            config.coadds_logs_dir,
            config.calib_data_dir,
            config.bolo_data_dir,
            config.min_time_maps]

    if(is_winter_map_update_running == false) {
        is_winter_map_update_running = true;
        // update data skims
        child = execFile(config.python_location,
                         args, {maxBuffer: 1024*1024*8},
                         function(err) {
                             console.log(err);
                             console.log('Finished winter map coadds and plots.');
                             is_winter_map_update_running = false;
                         });
        console.log('Updating winter maps...');
    }
    else {
        console.log('Winter map updater already running, so not spawning again!');
    }
}

function updateSummerMapPlots() {
    args = ['-B',
            'update_summary.py',
            'maps',
            'summer',
            config.maps_data_dir,
            config.coadds_data_dir,
            config.coadds_figs_dir,
            config.coadds_logs_dir,
            config.calib_data_dir,
            config.bolo_data_dir,
            config.min_time_maps]

    if(is_summer_map_update_running == false) {
        is_summer_map_update_running = true;
        // update data skims
        child = execFile(config.python_location,
                         args,
                         {maxBuffer: 1024*1024*8},
                         function(err) {
                             console.log(err);
                             console.log('Finished summer map coadds and plots.');
                             is_summer_map_update_running = false;
                         });
         console.log('Updating summer maps...');
    }
    else {
        console.log('Summer map updater already running, so not spawning again!');
    }
}

function updateArcPlots() {
    args = ['-B',
            'update_summary.py',
            'arcs',
            config.arc_data_dir,
            config.arc_pkls_dir,
            config.arc_figs_dir,
            config.arc_logs_dir,
            config.min_time_arcs]

    if(is_arc_update_running == false) {
        is_arc_update_running = true;
        // update arc file skims
        child = execFile(config.python_location,
                         args,
                         {maxBuffer: 1024*1024*8},
                         function(err) {
                             console.log(err);
                             console.log('Finished updating arc file skims and plots');
                             is_arc_update_running = false;
                         });
         console.log('Updating arc file skims and plots...');
    }
    else {
        console.log('Arc file updater already running, so not spawning again!');
    }
}

// update both types of plots in parallel every 10 minutes
setInterval(updateSummaryPlots, 600000);
if(config.site == 'pole') { // only update map and arc plots at pole
    setInterval(updateWinterMapPlots, 600000);
    setInterval(updateSummerMapPlots, 600000);
    setInterval(updateArcPlots,      3600000);
}
