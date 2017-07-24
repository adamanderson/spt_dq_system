// redraw determines if the selected tabulator table needs to be redrawn
var redraw = false;
// switch tabs
function open_tab(evt, tab) {
  // remove tabs from display
  var tabcontent = document.getElementsByClassName("tabcontent");
  for (var i = 0; i < tabcontent.length; i++)
    tabcontent[i].style.display = "none";

  // remove active from tablinks
  var tablinks = document.getElementsByClassName("tablinks");
  for (var i = 0; i < tablinks.length; i++)
    tablinks[i].className = tablinks[i].className.replace(" active", "");

  // add new tab to display and active
  document.getElementById(tab).style.display = "block";
  evt.currentTarget.className += " active";

  if (redraw) {
    // redraw tabulator because it needs to have a display window when
    // created
    $("#" + tab).tabulator("redraw", true);
    redraw = false;
  }
}


function make_t_table(select) {
  //create Tabulator on DOM element with id "transfer_table"
  $("#transfer_table").tabulator({
    //pagination:"remote",  // use this for normal pagination
    ajaxURL:"/tpage", //set the ajax URL
    ajaxParams: {source: $("#obstype-search").val(),
                 date:  {min: $("#date-from").val(),
                         max: $("#date-to").val()},
                 observation:    {min: $("#obsid-from").val(),
                         max: $("#obsid-to").val()}},
    ajaxConfig:'GET',
    ajaxSorting: true,
    //paginationSize:40,
    selectable: select,
    index:"_id",
    height:"400px", // set height of table (optional)
    fitColumns:true, //fit columns to width of table (optional)
    columns:[ //Define Table Columns
      {title:"Observation ID", field:"observation", sorter:"number"},
      {title:"Source", field:"source"},
      {title:"Fullrate status", field:"status_fullrate"},
      {title:"Downsampled status", field:"status_downsampled"},
      {title:"Fullrate transfer", field:"transfer_fullrate",
          formatter:"tickCross"},
      {title:"Downsampled transfer", field:"transfer_downsampled",
          formatter:"tickCross"},
      {title:"Date", field:"date", sorter:"date",
          sorterParams:{format:"YYYY-MM-DD hh:mm:ssZZ"}}
    ]
  });
}

function make_aux_table() {
// create Tabulator on DOM element with id "example-table"
  $("#aux_table").tabulator({
    //pagination:"remote",  // use this for normal pagination
    ajaxURL:"/apage", //set the ajax URL
    ajaxParams: {date:  {min: $("#date-from").val(),
                         max: $("#date-to").val()}},
    ajaxConfig:'GET',
    ajaxSorting: true,
    //paginationSize:40,
    index:"_id",
    height:"400px", // set height of table (optional)
    fitColumns:true, //fit columns to width of table (optional)
    columns:[ //Define Table Columns
      {title:"Filename", field:"filename"},
      {title:"Type", field:"type"},
      {title:"Size", field:"size", sorter:"number"},
      {title:"Status", field:"status"},
      {title:"Modified", field:"modified", sorter:"date",
          sorterParams:{format:"YYYY-MM-DD hh:mm:ssZZ"}},
      {title:"Date", field:"date", sorter:"date",
          sorterParams:{format:"YYYY-MM-DD hh:mm:ssZZ"}}
    ]
  });
}

// get plotting mode and build the table
var func_val = $('input[name="func"]:checked').val();
if (func_val == "timeseries")
  make_t_table(true);
else if (func_val == "individual")
  make_t_table(4);

make_aux_table();

// wrapper for two search functions
function search() {
  var tab = document.getElementsByClassName("tablinks active")[0].id;
  if (tab == 'aux')
    asearch();
  else if (tab == 'transfer')
    tsearch();
}

// sends search request and rebuilds table and plot types
function tsearch() {
  // clear selections in table
  $("#transfer_table").tabulator("deselectRow");
  // package up the search fields into a json to send to server
  querydata = {source: $("#obstype-search").val(),
               date:  {min: $("#date-from").val(),
                       max: $("#date-to").val()},
               observation:    {min: $("#obsid-from").val(),
                       max: $("#obsid-to").val()}}
  $("#transfer_table").tabulator("setData", "/tpage", querydata);
  plot_list();
};

function asearch() {
  // package up the search fields into a json to send to server
  querydata = {date:  {min: $("#date-from").val(),
                       max: $("#date-to").val()}}
  $("#aux_table").tabulator("setData", "/apage", querydata);
};

// sends plot request to server and opens images in new window
function plot() {
  // check if transfer or aux
  var tab = document.getElementsByClassName("tablinks active")[0];
  if (tab.id == "aux") {
    alert("There are no plotting capabilities for the aux database.");
    return;
  }
  // get selected observations
  var rows = $("#transfer_table").tabulator("getSelectedData");
  if (rows.length == 0) {
    alert("Please select at least one observation.");
    return;
  }

  // get plot type
  var selected_values = [];    
  $("#plot_select :selected").each(function(){
    selected_values.push($(this).val()); 
  });
  if (selected_values.length == 0) {
    alert("Please select at least one plot type.");
    return;
  }

  // get plotting mode
  var func_val = $('input[name="func"]:checked').val();

  // items holds the image info needed for photoswipe
  var items = [];
  // deferred holds the jquery requests so we can know when they are return
  var deferred = [];
  
  if (func_val == "individual") {
    // loop over each observation
    $.each(rows, function(i, obsdata) {
      // combine all requested plot types into a string
      obsdata['plot_type'] = selected_values.join(' ');

      obsdata['func'] = func_val;
      // request the plot
      deferred.push($.get("data_req", obsdata, function(err, status) {
        if (err != null) {
          alert(err);
          return;
        }
        // put the image info in the list
        for (var i = 0; i < selected_values.length; i++) {
          src = 'img/' + obsdata['source'] + obsdata['observation'] +
              selected_values[i] + '.png';
          items.push({src: src, w: 0, h: 0, obs: obsdata['observation']});
        }
      }));
    });
  } else if (func_val == "timeseries") {
    // don't allow "any" source. Only one source allowed for timeseries mode
    if ($("#obstype-search").val() == '') {
      alert('Please search for a specific source in timeseries mode.');
      return;
    }
    // loop over each plot type
    $.each(selected_values, function(i, type) {
      // combine requested observations into a string
      var obs = [];
      for (var j = 0; j < rows.length; j++)
        obs.push(rows[j]['observation']);

      var obsdata = {plot_type: type, observation: obs.join(' '),
          source: rows[0]['source'], func: func_val};
      // request plots
      deferred.push($.get("data_req", obsdata, function(err, status) {
        if (err != null) {
          alert(err);
          return;
        }
        // put the image info in the list
        for (var i = 0; i < selected_values.length; i++) {
          src = 'img/' + obsdata['source'] + obsdata['observation'] +
              selected_values[i] + '.png';
          items.push({src: src, w: 0, h: 0, obs: obsdata['observation']});
        }
      }));
    });
  }
  // after all plots have finished
  $.when.apply(null, deferred).then( function() {
    // do nothing if all images generated errors
    if (items.length == 0)
      return;
    // open display window and load images
    var w = window.open('display.html', '_blank');
    // need to open image viewer after page has loaded
    w.onload = function() {
      // sort by observation
      items.sort(compare);
      w.start_image(items);
    }
  });
}

// used to sort image items
function compare(a, b) {
  if (a.obs > b.obs)
    return -1;
  else if (a.obs < b.obs)
    return 1;
  else
    return 0;
}


