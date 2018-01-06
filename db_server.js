var express = require('express');
var sqlite3 = require('sqlite3');
var assert = require('assert');
var bodyParser = require('body-parser');
var exphbs = require('express-handlebars');
var squel = require("squel");
var execFile = require('child_process').execFile;
var fs = require('fs');
var readline = require('readline');
var util = require('util');
var auth = require('express-basic-auth');
var bcrypt = require('bcryptjs');
var https = require('https');
var moment = require('moment');
var SSE = require('express-sse');
var sse = new SSE();
var yaml = require('yamljs');

// TODO error handling

// load configuration info from yaml file
config = yaml.load('config.yaml')

// setup express
var app = express();
app.use(bodyParser.json());
app.use(bodyParser.urlencoded({extended: true}));
// server side events for passing messages
app.get('/sse', sse.init);
// counter to separate messages by plot request
var sseid = 0;

// determines if username and password are correct
function authorizer(username, password, cb) {
  if (username != 'spt')
    return cb(null, false);
  // get hash
  filePath = __dirname + '/hash';
  fs.readFile(filePath, {encoding: 'utf-8'}, function(err, hash){
    if (!err) {
      // check password vs hash
      bcrypt.compare(password, hash.trim(), function(err, res) {
        return cb(null, res);
      });
    } else {
      log(err);
    }
  });
}

// database filenames
db_files = {transfer: config.transfer_db_path,
	    aux_transfer: config.auxtransfer_db_path,
	    autoproc: config.autoproc_db_path}

// response if wrong credentials
function getUnauthorizedResponse(req) {
  return 'Credentials rejected';
}

// password protect the website if certificate is given in config.yaml
if(config.key_file && config.cert_file) {
  app.use(auth({
    authorizer: authorizer,
    authorizeAsync: true,
    challenge: true,
    unauthorizedResponse: getUnauthorizedResponse,
  }))
}

// setup home page
app.use(express.static('public'));
app.get('/index.html', function (req, res) {
  res.sendFile( __dirname + "/" + "index.html" );
})
app.get('/', function (req, res) {
  res.sendFile( __dirname + "/" + "index.html" );
})

// handlebars for templating
app.engine('handlebars', exphbs({defaultLayout: 'main'}));
app.set('view engine', 'handlebars');

// needed to load js and css files
app.use('/js',express.static(__dirname + '/js'));
app.use('/css',express.static(__dirname + '/css'));

// create directory in /tmp to store plots
// TODO: don't cache every plot. Remove old plots from time to time
// NOTE: not a high priority because right now ~1.2TB free in /tmp
if (!fs.existsSync(config.plot_dir)){
      fs.mkdirSync(config.plot_dir);
}
app.use('/img', express.static(config.plot_dir));

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
	db = new sqlite3.Database(db_files[req.query.dbname]);

	// get data from the database
	query = squel.select().from(req.query.dbname);
	parseSearch(query, req.query, req.query.dbname);
	db.all(query.toParam()['text'],
	       query.toParam()['values'],
	       function(err, rows) {
		   res.send(rows);
	       });
	
	// close the database
	db.close();
    });




// page for displaying plots/data
app.get('/display.html', function(req, res) {
  res.sendFile( __dirname + "/" + "display.html" );
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
  options = {'timeout':20000};
  tab = req.query['table'];
  if (tab == 'transfer') {
    obs = req.query['observation'];
    source = req.query['source'];
  }
  else if (tab == 'aux') {
    filename = req.query['filename'];
  }
  else if (tab == 'autoproc') {
      filename = req.query['filename'];
  }
  plot_type = req.query['plot_type'];
  func_val = req.query['func'];

  log(util.inspect(req.query));
  // execute python plotting script. Safe because user input
  // is passed as arguments to the python script and the python
  // script handles the arguments safely.
  var python = '/cvmfs/spt.opensciencegrid.org/py3-v1/RHEL_6_x86_64/bin/python';
  var args;
  if (func_val == 'individual' && tab == 'transfer')
    args = ['-B', './plot/_plot.py',
	    func_val,
	    source,
	    config.plot_dir,
	    config.calib_data_dir,
	    config.bolo_data_dir,
	    obs].concat(
        plot_type.split(' '));
  else if (func_val == 'timeseries' && tab == 'transfer')
    args = ['-B', './plot/_plot.py',
	    func_val,
	    source,
	    config.plot_dir,
	    config.calib_data_dir,
	    config.bolo_data_dir,
	    plot_type].concat(
        obs.split(' '));
  else if (func_val == 'individual' && tab == 'aux')
    args = ['-B', './plot/_plot.py', func_val, 'aux', filename].concat(
        plot_type.split(' '));
  else if (func_val == 'timeseries' && tab == 'aux')
    args = ['-B', './plot/_plot.py', func_val, 'aux', plot_type].concat(
        filename.split(' '));
  var err = null;
  var child = execFile(python, args, options);

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
  if (json[req.query['tab']][req.query['func']][req.query['type']] == null)
    req.query['type'] = 'any';
  res.json(json[req.query['tab']][req.query['func']][req.query['type']]);
});

// get list of available sources
app.get('/sourcelist', function(req, res) {
	db = new sqlite3.Database(config.transfer_db_path);
	query = squel.select()
	    .from('transfer')
	    .field('source')
	    .distinct();
	db.all(query.toString(), function(err, rows) {
		res.send(rows);
		    });
    });

// turns search into a sql query
function parseSearch(query, searchJSON, tab) {
  if (tab == 'transfer') {
    if(searchJSON['source']) {
      query.where("source == ?", searchJSON['source']);
    }

    // user could specify min obs, max obs or both
    if(searchJSON['observation']['min'] && searchJSON['observation']['max']) {
      query.where('observation >= ?', searchJSON['observation']['min'])
          .where('observation <= ?', searchJSON['observation']['max']);
    }
    else if (searchJSON['observation']['min']) {
      query.where('observation >= ?', searchJSON['observation']['min']);
    }
    else if (searchJSON['observation']['max']) {
      query.where('observation <= ?', searchJSON['observation']['max']);
    }

    // user could specify min date, max date, or both
    var min_time = searchJSON['date']['min'];
    var max_time = searchJSON['date']['max'];
    if(min_time && max_time) {
	query.where("date(date) BETWEEN date(?) AND date(?)", min_time, max_time);
    }
    else if (min_time) {
	// set max time to current date
	max_time = moment().format('YYYY-MM-DD'); 
	query.where("date(date) BETWEEN date(?) AND date(?)", min_time, max_time);
    }
    else if (max_time) {
	// set min time before any observations
	min_time = "2000-01-01";
	query.where("date(date) BETWEEN date(?) AND date(?)", min_time, max_time);
    }
  }
  else if (tab == 'aux') {
    if (searchJSON['filename']) {
      query.where("filename LIKE ?", '%' + searchJSON['filename'] + '%');
    }
    if (searchJSON['type']) {
      query.where("type == ?", searchJSON['type']);
    }

    // user could specify min date, max date, or both
    var min_time = searchJSON['date']['min'];
    var max_time = searchJSON['date']['max'];
    if(min_time && max_time) {
	query.where("date(date) BETWEEN date(?) AND date(?)", min_time, max_time);
    }
    else if (min_time) {
	// set max time to current date
	max_time = moment().format('YYYY-MM-DD'); 
	query.where("date(date) BETWEEN date(?) AND date(?)", min_time, max_time);
    }
    else if (max_time) {
	// set min time before any observations
	min_time = "2000-01-01";
	query.where("date(date) BETWEEN date(?) AND date(?)", min_time, max_time);
    }
  }
  else if (tab == 'autoproc') {
      // user could specify min date, max date, or both
    var min_time = searchJSON['modified']['min'];
    var max_time = searchJSON['modified']['max'];
    if(min_time && max_time) {
	query.where("date(modified) BETWEEN date(?) AND date(?)", min_time, max_time);
    }
    else if (min_time) {
	// set max time to current date
	max_time = moment().format('YYYY-MM-DD'); 
	query.where("date(modified) BETWEEN date(?) AND date(?)", min_time, max_time);
    }
    else if (max_time) {
	// set min time before any observations
	min_time = "2000-01-01";
	query.where("date(modified) BETWEEN date(?) AND date(?)", min_time, max_time);
    }

  }

  var sort = searchJSON['sort'];
  var sort_dir = searchJSON['sort_dir'];
  if (sort) {
    if (sort_dir == 'desc')
      query.order(sort, false);
    else
      query.order(sort, true);
  }
  else if (tab == 'transfer' || tab == 'aux')
    query.order('date', false);
  else if (tab == "autoproc")
      query.order('modified', false);

  return query;
}

// use https if certificate info is given in config.yaml
if(config.key_file && config.cert_file) {
  var options = {
    key: fs.readFileSync(config.key_file),
    cert: fs.readFileSync(config.cert_file)
  };

  // run server
  log('Listening on port ' + config.port);
  https.createServer(options, app).listen(parseInt(config.port));
}
else {
    app.listen(parseInt(config.port))
};