var express = require('express');
var sqlite3 = require('sqlite3').verbose();
var assert = require('assert');
var bodyParser = require('body-parser');
var exphbs = require('express-handlebars');
var squel = require("squel");
var moment = require("moment");
var execFile = require('child_process').execFile
var fs = require('fs')
var readline = require('readline');

// start the server
var app = express()
app.use(bodyParser.json());
app.use(bodyParser.urlencoded({extended: true}));

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
db = new sqlite3.Database('master_database.db');

// needed to load chosen
app.use('/js',express.static(__dirname + '/js'));
app.use('/css',express.static(__dirname + '/css'));

// create directory in /tmp to store plots
// TODO: be nice and clean up /tmp every now
// and then so it doesn't fill up
var plot_dir = '/tmp/spt_dq/';
if (!fs.existsSync(plot_dir)){
      fs.mkdirSync(plot_dir);
}

app.use('/img', express.static(plot_dir));

app.get('/page', function(req, res) {
  // get all data from the database
  query = squel.select()
                .from('master_database')
  parseSearch(query, req.query);
  var max_pages;
  db.all(query.toString(), function(err, rows) {
    max_pages = Math.ceil(Object.keys(rows).length / req.query['size']);
  });

  query = query.offset((req.query['page']-1)*req.query['size'])
              .limit(req.query['size'])
              .order('observation', false);
  console.log(query.toString());
  db.all(query.toString(), function(err, rows) {
    assert.equal(null, err);
    var data = {'last_page': max_pages, 'data': rows};
	  res.send(data);
  });
});

// page for displaying plots/data
app.get('/display.html', function(req, res) {
  res.sendFile( __dirname + "/" + "display.html" );
})


// request data, request contains observation to make a plot of
// the python script saves a plot and this gets read and sent
// back
app.get('/data_req', function (req, res) {
    options = {'timeout':10000};
    path = req.query['path'];
    obs = req.query['observation'];
    source = req.query['source'];
    plot_type = req.query['plot_type'];
    console.log('Request for ' + source + ' ' + obs);

    // setup return value
    var err = {num:null, msg:null};

    if (fs.existsSync(path + 'nominal_online_cal.g3')) { 
      console.log(req.query);
      console.log('python -B ./plot/_plot.py ' + ' ' + path + ' ' + source + ' ' + obs + ' ' + plot_type);
      // execute python plotting script. Safe because user input
      // is passed as arguments to the python script and the python
      // script handles the arguments safely.
      var python = '/cvmfs/spt.opensciencegrid.org/py3-v1/RHEL_6_x86_64/bin/python';
      var args = ['-B', './plot/_plot.py', path, source, obs].concat(plot_type.split(' '))
      execFile(python, args, options, (error, stdout, stderr) => {
        if (error) {
          err = {num:error.code, msg:stderr, type:stdout};
          if (err.num == 3) {
            err.msg = 'Error making plot ' + err.type + '. Data does not exist in file.';
          }
          if (err.num == 4) {
            err.msg = 'Error. Could not find observation file.';
          }
          console.error(`exec error: ${error}`);
        }
        res.json(err);
      });
      return;
    }
    console.error('File does not exist');
    err.num = 1;
    err.msg = 'Error. Data file not found.';
    res.json(err);
})

// get all available plot types, removing the driver file and .py
app.get('/plot_list', function(req, res) {
  fs.readdir('./plot/', function(err, items) {
    index = items.indexOf('plot.py');
    items.splice(index, 1);
    for (var i = 0; i < items.length; i++) {
      items[i] = items[i].slice(0, -3);
    }
    res.json(items);
  });
})

function parseSearch(query, searchJSON) {
  if(searchJSON['source']) {
    query.where('source == \'' + searchJSON['source'] + '\'')
  }
  if(searchJSON['observation']['min'] && searchJSON['observation']['max']) {
    query.where('observation > ' + searchJSON['observation']['min'])
        .where('observation < ' + searchJSON['observation']['max'])
  }
  if(searchJSON['date']['min'] && searchJSON['date']['max']) {
    query.where('start_time > ' + moment(searchJSON['date']['min'], 'DD/MM/YYYY').unix()*1e8)
        .where('start_time < ' + moment(searchJSON['date']['max'], 'DD/MM/YYYY').unix()*1e8)
  }
  return query
}

app.listen(3000, function() {
  console.log('Listening on port 3000');
});
