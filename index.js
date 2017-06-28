// index.js
const path = require('path')  
const express = require('express')  
const exphbs = require('express-handlebars')
const jquery = require('jquery')
const spawnSync = require('child_process').spawnSync
const fs = require('fs')

const app = express()

// spawns a python process and handles input / output
// waits for process to finish before returning
function run_py(file) {
    options = {"input":JSON.stringify(file)};
    py = spawnSync('python', ['test.py'], options);
    if (py.status != 0)
        console.log(py.stderr.toString());
    return py.stdout;
}

app.engine('.hbs', exphbs({  
  defaultLayout: 'index',
  extname: '.html',
  layoutsDir: path.join(__dirname, 'views/layouts')
}))
app.set('view engine', '.hbs')  
app.set('views', path.join(__dirname, 'views'))

// home page
app.get('/', function (req, res) {
  res.sendFile('./index.html', { root: __dirname })
})

// request data, request contains observation to make a plot of
// the python script saves a plot and this gets read and sent
// back
app.get('/data_req', function (req, res) {
    console.log("hello request");
    console.log(req._parsedUrl.query);
    // make plot
    var plot = run_py(req._parsedUrl.query);
    // read image file made by plot and convert to base64
    var img = fs.readFileSync('./' + plot).toString("base64");
    // send base64 image to client
    res.writeHead(200, {'Content-Type': 'image/png' });
    res.end(img, 'binary');
})

// used to populate list of available observations
// only gets offline calibration dates currently
app.get('/date_list', function (req, res) {
    console.log("date list request");
    var path = "/spt/data/bolodata/downsampled/RCW38-pixelraster/";
    fs.readdir(path, function(err, items) {
        res.json(items);
    });
})

app.use((err, request, response, next) => {
  // log the error, for now just console.log
  console.log(err)
  response.status(500).send('Something broke!')
})

app.listen('3000')
