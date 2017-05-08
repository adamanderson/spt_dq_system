var express = require('express');
var sqlite3 = require('sqlite3').verbose();
var assert = require('assert');
var bodyParser = require('body-parser');
var exphbs = require('express-handlebars')

// start the server
var app = express()
app.use(bodyParser.json());
app.use(bodyParser.urlencoded({extended: true}));

app.use(express.static('public'));
app.get('/index.html', function (req, res) {
   res.sendFile( __dirname + "/" + "index.html" );
})

// handlebars for templating
app.engine('handlebars', exphbs({defaultLayout: 'main'}));
app.set('view engine', 'handlebars');

// open the database
db = new sqlite3.Database('obsfiles.db');

app.get('/series.html', function (req, res) {
  res.render('series')
})

app.get('/all', function(req, res) {
    // get all data from the database
    db.all('SELECT * FROM obsfiles;', function(err, rows) {
	assert.equal(null, err);
	console.log('querying all rows from observation file table...')
	res.send(rows);
    });
});

app.post('/search', function(req, res) {
    // do a search
    console.log(req.body.type)
    db.all('SELECT * FROM obsfiles WHERE type = \'' + req.body.type + '\'',
	   function(err, rows) {
	       assert.equal(null, err);
	       console.log('Getting rows satisfying: type = '+req.body.type);
	       res.send(rows);
	   });
});

app.all('/', function(req, res) {
  res.redirect("http://127.0.0.1:3000/index.html");
});

app.listen(3000, function() {
  console.log('Listening on port 3000');
});
