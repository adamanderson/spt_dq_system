// redraw determines if the selected tabulator table needs to be redrawn
var redraw = false;
// switch tabs
function open_tab(evt, tab) {
  // deselect all from current tab
  $("#" + tab).tabulator("deselectRow");
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
  // remake plot_list
  plot_list();

  // hide/show obs id,source and type,filename search fields
  if (tab == 'transfer_table') {
    document.getElementById("obsid-row").style.display = "table-row";
    document.getElementById("source-row").style.display = "table-row";
    document.getElementById("type-row").style.display = "none";
    document.getElementById("file-row").style.display = "none";
    document.getElementById("modified-row").style.display = "none";
  }
  else if (tab == 'aux_table') {
    document.getElementById("obsid-row").style.display = "none";
    document.getElementById("source-row").style.display = "none";
    document.getElementById("type-row").style.display = "table-row";
    document.getElementById("file-row").style.display = "table-row";
    document.getElementById("modified-row").style.display = "none";
  }
  else if (tab == 'autoproc_table') {
    document.getElementById("obsid-row").style.display = "table-row";
    document.getElementById("source-row").style.display = "table-row";
    document.getElementById("type-row").style.display = "none";
    document.getElementById("file-row").style.display = "none";
    document.getElementById("modified-row").style.display = "table-row";
  }
}


function make_table(type, select) {
    //create Tabulator on DOM element with id "transfer_table"
    $("#transfer_table").tabulator({
	    //pagination:"remote",  // use this for normal pagination
	    ajaxURL:"/dbpage", 
		ajaxParams: {search: {source: $("#obstype-search").val(),
			date:  {min: $("#date-from").val(),
			    max: $("#date-to").val()},
			observation:    {min: $("#obsid-from").val(),
			    max: $("#obsid-to").val()}},
		    dbname: "transfer"},
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


function make_t_table(select) {
  //create Tabulator on DOM element with id "transfer_table"
  $("#transfer_table").tabulator({
    //pagination:"remote",  // use this for normal pagination
    ajaxURL:"/dbpage", //set the ajax URL
	      ajaxParams: {search: {source: $("#obstype-search").val(),
		      date:  {min: $("#date-from").val(),
			  max: $("#date-to").val()},
		      observation:    {min: $("#obsid-from").val(),
			  max: $("#obsid-to").val()}},
		  dbname: "transfer"},
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

function make_aux_table(select) {
// create Tabulator on DOM element with id "example-table"
  $("#aux_table").tabulator({
    //pagination:"remote",  // use this for normal pagination
    ajaxURL:"/dbpage", //set the ajax URL
	      ajaxParams: {search: {date:  {min: $("#date-from").val(),
			  max: $("#date-to").val()},
		      filename: $("#file-search").val(),
		      type: $("#auxtype-search").val()},
		  dbname: "aux_transfer"},
    ajaxConfig:'GET',
    ajaxSorting: true,
    //paginationSize:40,
    selectable: select,
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

function make_autoproc_table(select) {
// create Tabulator on DOM element with id "example-table"
  $("#autoproc_table").tabulator({
    //pagination:"remote",  // use this for normal pagination
    ajaxURL:"/dbpage", //set the ajax URL
    ajaxParams: {modified:  {min: $("#date-from").val(),
		      max: $("#date-to").val()},
                 observation:    {min: $("#obsid-from").val(),
		      max: $("#obsid-to").val()},
		  source: $("#obstype-search").val(),
		  dbname: "autoproc"},
    ajaxConfig:'GET',
    ajaxSorting: true,
    //paginationSize:40,
    selectable: select,
    index:"_id",
    height:"400px", // set height of table (optional)
    fitColumns:true, //fit columns to width of table (optional)
    columns:[ //Define Table Columns
      {title:"Source", field:"source"},
      {title:"Observation", field:"observation"},
      {title:"Status", field:"status"},
      {title:"Modified", field:"modified", sorter:"date",
          sorterParams:{format:"YYYY-MM-DD hh:mm:ssZZ"}},
    ]
  });
}


// get plotting mode and build the table
var func_val = $('input[name="func"]:checked').val();
if (func_val == "timeseries") {
  make_t_table(true);
  make_aux_table(true);
  make_autoproc_table(true);
}
else if (func_val == "individual") {
  make_t_table(10);
  make_aux_table(10);
  make_autoproc_table(10);
}


// wrapper for two search functions
function search() {
  var tab = document.getElementsByClassName("tablinks active")[0].id;
  if (tab == 'aux')
    asearch();
  else if (tab == 'transfer')
    tsearch();
  else if (tab == 'autoproc')
      autoprocsearch();
}

// sends search request and rebuilds table and plot types
function tsearch() {
  // clear selections in table
  $("#transfer_table").tabulator("deselectRow");
  // package up the search fields into a json to send to server
  querydata = {search: {source: $("#obstype-search").val(),
			date:  {min: $("#date-from").val(),
				max: $("#date-to").val()},
			observation:    {min: $("#obsid-from").val(),
					 max: $("#obsid-to").val()}},
	       dbname: "transfer"}
  $("#transfer_table").tabulator("setData", "/dbpage", querydata);
  plot_list();
};

function asearch() {
    // package up the search fields into a json to send to server
    querydata = {search: {date:  {min: $("#date-from").val(),
				  max: $("#date-to").val()},
			  filename: $("#file-search").val(),
			  type: $("#auxtype-search").val()},
		 dbname: "aux_transfer"};
    $("#aux_table").tabulator("setData", "/dbpage", querydata);
    plot_list();
};

function autoprocsearch() {
    // package up the search fields into a json to send to server
    querydata = {search: {modified:  {min: $("#date-from").val(),
				      max: $("#date-to").val()},
			  observation:    {min: $("#obsid-from").val(),
					   max: $("#obsid-to").val()}},
			  dbname: "autoproc"};
		 $("#autoproc_table").tabulator("setData", "/dbpage", querydata);
		 plot_list();
};


// create sse listener
var es = new EventSource("/sse");

// sends plot request to server and opens images in new window
function plot() {
  // holds selected observations
  var rows;
  // check if transfer or aux
  var tab = document.getElementsByClassName("tablinks active")[0];
  if (tab.id == "aux") {
    rows = $("#aux_table").tabulator("getSelectedData");
  }
  else if (tab.id == "transfer") {
    rows = $("#transfer_table").tabulator("getSelectedData");
  }

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

  // open new window immediately so it isn't considered a pop-up
  var display_win = window.open('display.html', '_blank');

  // items holds the image info needed for photoswipe
  var items = [];

  // counts the number of images to be created
  var image_ctr = 0;

  var sseid;
  $.get("sseid", function(id, status) {
    sseid = id;
    // counts the images that have finished being created
    var plot_ctr = 0;

    es.addEventListener('out' + sseid, function (event) {
      // events include opening and closing quotations so slice them out

      if (event.data.slice(1, 4) == 'msg')
        display_win.show(event.data.slice(4, -1));

      else if (event.data.split('fln').length == 2) {
	  path =  event.data.split('fln')[1].slice(0,-1);
	  items.push({src: 'img/' + path, w: 0, h: 0});
      }

      else if (event.data.slice(1, 4) == 'plt') {
	  display_win.load_progress();
        plot_ctr++;

        if (image_ctr == plot_ctr) {
          // done making plots
          items.sort(compare);
          display_win.start_image(items);
        }
      }

    });

    es.addEventListener('err' + sseid, function (event) {
      display_win.show_err(event.data);
    });
  });

  // make the plots and everything else once display window has opened
  // so that we can track plotting progress and send images to the new window.
  display_win.onload = function () {

  // get plotting mode
  var func_val = $('input[name="func"]:checked').val();

  // set the number of images and variables to keep track of creation
  if (func_val == "individual")
    display_win.images_to_load = selected_values.length * rows.length;
  else if (func_val == "timeseries")
    display_win.images_to_load = selected_values.length;
  display_win.loaded_count = 0;
  display_win.loading_progress = 0;

  if (func_val == "individual") {
    // loop over each observation
    $.each(rows, function(i, obsdata) {
      // combine all requested plot types into a string
      obsdata['plot_type'] = selected_values.join(' ');

      obsdata['table'] = tab.id;

      obsdata['func'] = func_val;

      obsdata['sseid'] = sseid;

      // request the plot
      $.get("data_req", obsdata, function(data, status) {
        image_ctr += selected_values.length;
      });
    });
  } else if (func_val == "timeseries") {
    // don't allow "any" source. Only one source allowed for timeseries mode
    if (tab.id == 'transfer') {
      if ($("#obstype-search").val() == '') {
        alert('Please search for a specific source in timeseries mode.');
        return;
      }
    }
    // loop over each plot type
    $.each(selected_values, function(i, type) {
      var obsdata;
      if (tab.id == 'transfer') {
        // combine requested observations into a string
        var obs = [];
        for (var j = 0; j < rows.length; j++)
          obs.push(rows[j]['observation']);

        obsdata = {plot_type: type, observation: obs.join(' '),
            source: rows[0]['source'], func: func_val, table: tab.id,
            sseid: sseid};
      } else if (tab.id == 'aux') {
        // combine requested observations into a string
        var obs = [];
        for (var j = 0; j < rows.length; j++)
          obs.push(rows[j]['filename']);

        obsdata = {plot_type: type, filename: obs.join(' '),
            date: rows[0]['date'], func: func_val, table: tab.id,
            sseid: sseid};
      }
      // request plots
      $.get("data_req", obsdata, function(data, status) {
        image_ctr++;
      });
    });
  }
  }
}

// used to sort image items
function compare(a, b) {
  return a.src.localeCompare(b.src);
}


