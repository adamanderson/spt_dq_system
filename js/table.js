
function make_table(select) {
  //create Tabulator on DOM element with id "example-table"
  $("#example-table").tabulator({
    pagination:"remote",  // use this for normal pagination
    ajaxURL:"/page", //set the ajax URL
    ajaxParams: {source: $("#obstype-search").val(),
                 date:  {min: $("#date-from").val(),
                         max: $("#date-to").val()},
                 observation:    {min: $("#obsid-from").val(),
                         max: $("#obsid-to").val()}},
    ajaxConfig:'GET',
    paginationSize:40,
    selectable: select,
    index:"_id",
    height:"400px", // set height of table (optional)
    fitColumns:true, //fit columns to width of table (optional)
    columns:[ //Define Table Columns
      {title:"Observation ID", field:"observation", sorter:"number"},
      {title:"Source", field:"source"},
      {title:"Fullrate status", field:"status_fullrate"},
      {title:"Downsampled status", field:"status_downsampled"},
      {title:"Fullrate transfer", field:"transfer_fullrate", formatter:"tickCross"},
      {title:"Downsampled transfer", field:"transfer_downsampled", formatter:"tickCross"},
      {title:"Date", field:"date", sorter:"date"}
    ]
  });
}

var func_val = $('input[name="func"]:checked').val();
if (func_val == "timeseries")
  make_table(true);
else if (func_val == "individual")
  make_table(4);

function search() {
  // clear selections in table
  $("#example-table").tabulator("deselectRow");
  // package up the search fields into a json to send to server
  querydata = {page:1,
               size:20,
               source: $("#obstype-search").val(),
               date:  {min: $("#date-from").val(),
                       max: $("#date-to").val()},
               observation:    {min: $("#obsid-from").val(),
                       max: $("#obsid-to").val()}}
  $("#example-table").tabulator("setData", "/page", querydata);
  plot_list();
};

function plot() {
  // get selected observations
  var rows = $("#example-table").tabulator("getSelectedData");
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

  var func_val = $('input[name="func"]:checked').val();

  // items holds the image info needed for photoswipe
  var items = [];
  // deferred holds the jquery requests so we know when they are done
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
        for (var i = 0; i < selected_values.length; i++) {
          src = 'img/' + obsdata['source'] + obsdata['observation'] + selected_values[i] + '.png';
          items.push({src: src, w: 0, h: 0, obs: obsdata['observation']});
        }
      }));
    });
  } else if (func_val == "timeseries") {
    if ($("#obstype-search").val() == '') {
      alert('Please search for a specific source in timeseries mode.');
      return;
    }
    // loop over each plot
    $.each(selected_values, function(i, type) {
      // combine requested observations into a string
      var obs = [];
      for (var j = 0; j < rows.length; j++)
        obs.push(rows[j]['observation']);

      var obsdata = {plot_type: type, observation: obs.join(' '), source: rows[0]['source'], func: func_val};
      // request plots
      deferred.push($.get("data_req", obsdata, function(err, status) {
        if (err != null) {
          alert(err);
          return;
        }
        for (var i = 0; i < selected_values.length; i++) {
          src = 'img/' + obsdata['source'] + obsdata['observation'] + selected_values[i] + '.png';
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
    // need to add images after page has loaded
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


