var express = require('express')
var fs = require('fs')
// var sql = require('sql.js');
var MongoClient = require('mongodb').MongoClient
var assert = require('assert');

// var filebuffer = fs.readFileSync('filelist.db')
// var database = new SQL.Database(filebuffer);

var url = 'mongodb://localhost:27017/myproject';

var app = express()

app.use(express.static('public'));
app.get('/index.html', function (req, res) {
   res.sendFile( __dirname + "/" + "index.html" );
})

app.get('/', function(req, res) {
  res.send('hello world');
});


app.get('/search', function(req, res) {
  // query the database with an example query
  // var contents = database.exec('SELECT * FROM filetable;');
  //   console.log('test');
  //   console.log(contents);

  // Use connect method to connect to the server
  MongoClient.connect(url, function(err, db) {
    assert.equal(null, err);
    console.log("Connected successfully to server");

    //insertDocuments(db, function() {
    findDocuments(db, function(contents) {
      res.send(contents);
    //  });
    });
  });


});


var findDocuments = function(db, callback) {
  // Get the documents collection
  var collection = db.collection('documents');
  // Find some documents
  collection.find({'a': 3}).toArray(function(err, docs) {
    assert.equal(err, null);
    console.log("Found the following records");
    console.log(docs);
    callback(docs);
  });
}


app.listen(3000, function() {
  console.log('Listening on port 3000');

});

// var server = http.createServer(function(req, res) {
//     res.writeHead(200);
//
//     // query the database with an example query
//     var contents = database.exec('SELECT * FROM filetable;')
//     console.log(contents.values)
//
//     res.end(contents.values);
//     console.log('Returned a message!')
// });
// server.listen(8080);
