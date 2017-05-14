var express = require('express');
var sqlite3 = require('sqlite3').verbose();
var assert = require('assert');
var bodyParser = require('body-parser');
var exphbs = require('express-handlebars');
var squel = require("squel");
var moment = require("moment");

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
  db.all(query.toString(), function(err, rows) {
  	assert.equal(null, err);
  	res.send(rows);
  });
});

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
  console.log(query.toString())
  db.all(query.toString(), function(err, rows) {
    assert.equal(null, err);
    var data = {'last_page': max_pages, 'data': rows}
	  res.send(data);
  });
});

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
