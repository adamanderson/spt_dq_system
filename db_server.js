var express = require('express');
var sqlite3 = require('sqlite3').verbose();
var assert = require('assert');
var bodyParser = require('body-parser');
var exphbs = require('express-handlebars');
var squel = require("squel");

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

app.get('/all', function(req, res) {
    // get all data from the database
    db.all('SELECT * FROM master_database;', function(err, rows) {
	assert.equal(null, err);
	console.log('querying all rows from observation file table...')
	res.send(rows);
    });
});

app.post('/search', function(req, res) {
    // build the query
    query = squel.select()
                .from('master_database')
    for(var condition in req.body) {
      if(req.body[condition] && condition == 'source') {
        query.where('source == \'' + req.body[condition] + '\'')
      }
      if(condition == 'observation' || condition == 'date') {
        if(req.body[condition]['min'] && req.body[condition]['min']) {
          query.where(condition + ' > ' + req.body[condition]['min'])
              .where(condition + ' < ' + req.body[condition]['max'])
        }
      }
    }

    // do the query
    console.log(query.toString())
    db.all(query.toString(), function(err, rows) {
	       assert.equal(null, err);
	       res.send(rows);
	   });
});

// app.all('/', function(req, res) {
//   res.redirect("/index.html");
// });

app.listen(3000, function() {
  console.log('Listening on port 3000');
});
