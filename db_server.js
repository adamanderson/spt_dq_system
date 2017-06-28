var express = require('express');
var sqlite3 = require('sqlite3').verbose();
var assert = require('assert');
var bodyParser = require('body-parser');
var exphbs = require('express-handlebars');
var squel = require("squel");
var moment = require("moment");
var spawnSync = require('child_process').spawnSync
var exec = require('child_process').exec
var fs = require('fs')

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

app.post('/series', function (req, res) {
  res.render('series')
})

// get all data from the database
app.get('/all', function(req, res) {
  query = squel.select()
  console.log(query.toString())
  /*
  db.all(query.toString(), function(err, rows) {
  	assert.equal(null, err);
  	res.send(rows);
  });
  */
});

app.get('/page', function(req, res) {
  // get all data from the database
  query = squel.select()
                .from('master_database')
  parseSearch(query, req.query);
  var max_pages;
  /*
  db.all(query.toString(), function(err, rows) {
    max_pages = Math.ceil(Object.keys(rows).length / req.query['size']);
  });
  */

  query = query.offset((req.query['page']-1)*req.query['size'])
              .limit(req.query['size'])
              .order('observation', false);
  console.log(query.toString());
  /*
  db.all(query.toString(), function(err, rows) {
    assert.equal(null, err);
    var data = {'last_page': max_pages, 'data': rows}
	  res.send(data);
  });
  */
});

// used to populate list of available observations
// only gets offline calibration dates currently
app.get('/date_list', function (req, res) {
    console.log("date list request");
    var path = "/spt/data/bolodata/downsampled/RCW38-pixelraster/";
    fs.readdir(path, function(err, items) {
        res.json(items);
    });
})

// request data, request contains observation to make a plot of
// the python script saves a plot and this gets read and sent
// back
app.get('/data_req', function (req, res) {
    console.log('Request for date: ' + req._parsedUrl.query);
    options = {'timeout':10000};
    //TODO sanitize input (sanitized for now but will have to be better
    // once we get input for more files)
    input = req._parsedUrl.query;
    if (fs.existsSync('/spt/data/bolodata/downsampled/RCW38-pixelraster/'
                + input + '/nominal_online_cal.g3')) { 
        exec('python test.py ' + input , options, (error, stdout, stderr) => {
          if (error) {
            console.error(`exec error: ${error}`);
            res.json('error');
            return;
          }
          var img = fs.readFileSync('./' + stdout).toString('base64');
          res.writeHead(200, {'Content-Type': 'image/png' });
          res.end(img, 'binary');
        });
        return;
    }
    console.error('File does not exist');
    res.json('error');
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
