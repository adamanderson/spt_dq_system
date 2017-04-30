var express = require('express');
var sqlite3 = require('sqlite3').verbose();
var assert = require('assert');

// open the database
db = new sqlite3.Database('obsfiles.db');

// start the server
var app = express()

app.use(express.static('public'));
app.get('/index.html', function (req, res) {
   res.sendFile( __dirname + "/" + "index.html" );
})

app.get('/', function(req, res) {
  res.send('hello world');
});


app.get('/all', function(req, res) {
    // query the database with an example query
    db.all('SELECT * FROM obsfiles;', function(err, rows) {
	    assert.equal(null, err);
	    console.log(rows);
	    res.send(rows);
    });
});

app.listen(3000, function() {
  console.log('Listening on port 3000');
});
