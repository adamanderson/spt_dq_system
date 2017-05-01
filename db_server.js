var express = require('express');
var sqlite3 = require('sqlite3').verbose();
var assert = require('assert');
var bodyParser = require('body-parser');

// open the database
db = new sqlite3.Database('obsfiles.db');

// start the server
var app = express()
app.use(bodyParser.json());
app.use(bodyParser.urlencoded({extended: true}));

app.use(express.static('public'));
app.get('/index.html', function (req, res) {
   res.sendFile( __dirname + "/" + "index.html" );
})

app.get('/all', function(req, res) {
    // get all data from the database
    db.all('SELECT * FROM obsfiles;', function(err, rows) {
	assert.equal(null, err);
	console.log('querying all rows from observation file table...')
	res.send(rows);
    });
});

app.listen(3000, function() {
  console.log('Listening on port 3000');
});
