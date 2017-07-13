var express = require('express');
// TODO disable verbose for production
var sqlite3 = require('sqlite3').verbose();
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

// start the server
var app = express();
app.use(bodyParser.json());
app.use(bodyParser.urlencoded({extended: true}));

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

function getUnauthorizedResponse(req) {
  return 'Credentials rejected';
}

app.use(auth({
  authorizer: authorizer,
  authorizeAsync: true,
  challenge: true,
  unauthorizedResponse: getUnauthorizedResponse,
}))

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

// open the database
db = new sqlite3.Database('/spt/data/transfer_database/transfer.db');

// needed to load chosen
app.use('/js',express.static(__dirname + '/js'));
app.use('/css',express.static(__dirname + '/css'));

// create directory in /tmp to store plots
// TODO: don't cache every plot. Remove old plots from time to time
// NOTE: not a high priority because right now ~1.2TB free in /tmp
var plot_dir = '/tmp/spt_dq/';
if (!fs.existsSync(plot_dir)){
      fs.mkdirSync(plot_dir);
}

// logs messages along with a timestamp. keeps writesteam open
var ws = fs.createWriteStream('./db.log', {'flags': 'a'});
ws.on('error', function (err) {
  console.error('Error logging message:');
  console.error(msg);
  console.error(err);
});
function log(msg) {
    var d = new Date();
    ws.write(d.toLocaleString() + '\t' + msg + '\n');
}

app.use('/img', express.static(plot_dir));

// TODO error handling here
app.get('/page', function(req, res) {
  // get all data from the database
  query = squel.select()
                .from('transfer');
  parseSearch(query, req.query);
  db.all(query.toString(), function(err, rows) {
    var max_pages = Math.ceil(Object.keys(rows).length / req.query['size']);
    if (max_pages == 0)
      max_pages = 1;

    query = query.offset((req.query['page']-1)*req.query['size'])
                .limit(req.query['size'])
                .order('observation', false);
    log(query.toString());
    db.all(query.toString(), function(err, rows) {
      assert.equal(null, err);
      var data = {'last_page': max_pages, 'data': rows};
      res.send(data);
    });
  });

});

// page for displaying plots/data
app.get('/display.html', function(req, res) {
  res.sendFile( __dirname + "/" + "display.html" );
});

// page for displaying plots/data
app.get('/ps.html', function(req, res) {
  res.sendFile( __dirname + "/" + "ps.html" );
})


// request data, request contains observation to make a plot of
// the python script saves a plot and this gets read and sent
// back
app.get('/data_req', function (req, res) {
    options = {'timeout':20000};
    //path = req.query['path'];
    obs = req.query['observation'];
    source = req.query['source'];
    plot_type = req.query['plot_type'];

    log(util.inspect(req.query));
    // execute python plotting script. Safe because user input
    // is passed as arguments to the python script and the python
    // script handles the arguments safely.
    var python = '/cvmfs/spt.opensciencegrid.org/py3-v1/RHEL_6_x86_64/bin/python';
    var args = ['-B', './plot/_plot.py', source, obs].concat(plot_type.split(' '));
    var err = null;
    execFile(python, args, options, (error, stdout, stderr) => {
      if (error) {
        err = stdout;
        log('exec python error: ' + error.toString());
      }
      res.json(err);
    });
    return;
})

// get all available plot types, removing the driver file and .py
app.get('/plot_list', function(req, res) {
  if (req.query['type'] == '')
    req.query['type'] = 'any';
  var json = JSON.parse(fs.readFileSync('./plot/plot_config.json', 'utf8'));
  res.json(json[req.query['type']]);
})

function parseSearch(query, searchJSON) {
  if(searchJSON['source']) {
    query.where('source == \'' + searchJSON['source'] + '\'');
  }

  if(searchJSON['observation']['min'] && searchJSON['observation']['max']) {
    query.where('observation >= ' + searchJSON['observation']['min'])
        .where('observation <= ' + searchJSON['observation']['max']);
  }
  else if (searchJSON['observation']['min']) {
    query.where('observation >= ' + searchJSON['observation']['min']);
  }
  else if (searchJSON['observation']['max']) {
    query.where('observation <= ' + searchJSON['observation']['max']);
  }

  var min_time = searchJSON['date']['min'];
  var max_time = searchJSON['date']['max'];
  if(min_time && max_time) {
    query.where("date(date) BETWEEN date('" + min_time +
        "') AND date('" + max_time + "')");
  }
  else if (min_time) {
    // set max time to current date
    max_time = moment().format('YYYY-MM-DD'); 
    query.where("date(date) BETWEEN date('" + min_time +
        "') AND date('" + max_time + "')");
  }
  else if (max_time) {
    // set min time before any observations
    min_time = "2000-01-01";
    query.where("date(date) BETWEEN date('" + min_time +
        "') AND date('" + max_time + "')");
  }

  return query;
}

var options = {
  key: fs.readFileSync('key.pem'),
  cert: fs.readFileSync('cert.pem')
};

log('Listening on port 3000');
https.createServer(options, app).listen(3000);
