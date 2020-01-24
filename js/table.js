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
  if (tab == 'scanify_table') {
    document.getElementById("obsid-row").style.display = "table-row";
    document.getElementById("source-row").style.display = "table-row";
    document.getElementById("type-row").style.display = "none";
    document.getElementById("file-row").style.display = "none";
    document.getElementById("modified-row").style.display = "none";
    document.getElementById("plot-type").style.display = "block";
    document.getElementById("plot-button").style.display = "inline";
  }
  else if (tab == 'schedule_table') {
    document.getElementById("obsid-row").style.display = "none";
    document.getElementById("source-row").style.display = "none";
    document.getElementById("type-row").style.display = "none";
    document.getElementById("file-row").style.display = "none";
    document.getElementById("modified-row").style.display = "table-row";
    document.getElementById("plot-type").style.display = "none";
    document.getElementById("plot-button").style.display = "none";
  }
  else if (tab == 'aux_table') {
    document.getElementById("obsid-row").style.display = "none";
    document.getElementById("source-row").style.display = "none";
    document.getElementById("type-row").style.display = "table-row";
    document.getElementById("file-row").style.display = "table-row";
    document.getElementById("modified-row").style.display = "none";
    document.getElementById("plot-type").style.display = "none";
    document.getElementById("plot-button").style.display = "none";
  }
  else if (tab == 'autoproc_table') {
    document.getElementById("obsid-row").style.display = "table-row";
    document.getElementById("source-row").style.display = "table-row";
    document.getElementById("type-row").style.display = "none";
    document.getElementById("file-row").style.display = "none";
    document.getElementById("modified-row").style.display = "table-row";
    document.getElementById("plot-type").style.display = "block";
    document.getElementById("plot-button").style.display = "inline";
  }
}


function make_t_table(select, site) {
  // create Tabulator on DOM element with id "transfer_table"
  var site, obsid_col;

  jQuery.ajax({
	url: '/site',
    success: function(data) {
	  site = data;
	},
    async: false
  });

  if(site == 'north') {
	obsid_col = {title:"Observation ID", field:"observation", sorter:"number"};
  }
  else if(site == 'pole') {
    obsid_col = {title:"Observation ID", field:"observation", sorter:"number",
    			 formatter:"link",
				 formatterParams:{labelField:"observation", target:"_blank",
								  url:function(cell) {
									  return window.location.href + cell.getData().log_file;
								  }}};
  }

  $("#scanify_table").tabulator({
    ajaxURL:"/dbpage", // set the ajax URL
	      ajaxParams: {search: {source: $("#obstype-search").val(),
		      date:  {min: $("#date-from").val(),
			  max: $("#date-to").val()},
		      observation:    {min: $("#obsid-from").val(),
			  max: $("#obsid-to").val()}},
		  dbname: "scanify"},
    ajaxConfig:'GET',
    selectable: select,
    index:"_id",
    height:"400px", // set height of table (optional)
	layout:"fitColumns",
    columns:[ // define table columns
	  obsid_col,
      {title:"Source", field:"source"},
	  {title:"Scanify status", field:"status_scanify"},
      {title:"Transfer status<br>(fullrate)", field:"status_fullrate"},
      {title:"Transfer status<br>(downsampled)", field:"status_downsampled"},
      {title:"Marked for transfer<br>(fullrate)", field:"transfer_fullrate",
          formatter:"tickCross"},
      {title:"Marked for transfer<br>(downsampled)", field:"transfer_downsampled",
          formatter:"tickCross"},
      {title:"Date (UTC)", field:"date", sorter:"date",
          sorterParams:{format:"YYYY-MM-DD hh:mm:ssZZ"}},
    ],
	initialSort:[{column:"date", dir:"desc"}]
  });
}

function make_aux_table(select) {
// create Tabulator on DOM element with id "example-table"
  $("#aux_table").tabulator({
    ajaxURL:"/dbpage", // set the ajax URL
	      ajaxParams: {search: {date:  {min: $("#date-from").val(),
			  max: $("#date-to").val()},
		      filename: $("#file-search").val(),
		      type: $("#auxtype-search").val()},
		  dbname: "aux_transfer"},
    ajaxConfig:'GET',
    selectable: select,
    index:"_id",
    height:"400px", // set height of table (optional)
	layout:"fitColumns",
    columns:[ // define table columns
      {title:"Filename", field:"filename"},
      {title:"Type", field:"type"},
      {title:"Size", field:"size", sorter:"number"},
      {title:"Status", field:"status"},
      {title:"Modified (UTC)", field:"modified", sorter:"date",
          sorterParams:{format:"YYYY-MM-DD hh:mm:ssZZ"}},
      {title:"Date (UTC)", field:"date", sorter:"date",
          sorterParams:{format:"YYYY-MM-DD hh:mm:ssZZ"}}
	],
	initialSort:[{column:"date", dir:"desc"}]
  });
}

function make_autoproc_table(select) {
// create Tabulator on DOM element with id "example-table"
  $("#autoproc_table").tabulator({
    ajaxURL:"/dbpage", // set the ajax URL
	      ajaxParams: {search: {modified:  {min: $("#date-from").val(),
			  max: $("#date-to").val()},
		      date:  {min: $("#date-from").val(),
			  max: $("#date-to").val()},
		      observation:    {min: $("#obsid-from").val(),
			  max: $("#obsid-to").val()},
		      source: $("#obstype-search").val()},
		  dbname: "autoproc"},
    ajaxConfig:'GET',
    selectable: select,
    index:"_id",
    height:"400px", // set height of table (optional)
	layout:"fitColumns",
    columns:[ // define table columns
	  {title:"Source", field:"source", sorter:"number",
	   formatter:function(cell, formatterParams, onRendered) {
		   logFileName = cell.getData().log_file
		   if(logFileName.indexOf("calframe") == -1)
			   return "<a href=" + window.location.href + logFileName + ">" + cell.getData().source + "</a>";
		   else
		       return cell.getData().source;
	   }},
      {title:"Observation", field:"observation"},
      {title:"Status", field:"status"},
      {title:"Modified (UTC)", field:"modified", sorter:"date",
          sorterParams:{format:"YYYY-MM-DD hh:mm:ssZZ"}},
      {title:"Observation date (UTC)", field:"date", sorter:"date",
          sorterParams:{format:"YYYY-MM-DD hh:mm:ssZZ"}},
    ],
	initialSort:[{column:"date", dir:"desc"}]
  });
}

function make_schedule_table(select) {
// create Tabulator on DOM element with id "example-table"
  $("#schedule_table").tabulator({
    ajaxURL:"/dbpage", // set the ajax URL
	      ajaxParams: {search: {modified:    {min: $("#modified-from").val(),
		                                      max: $("#modified-to").val()},
		                        sch_start:   {min: $("#date-from").val(),
		                                      max: $("#date-to").val()},
		                        name: $("#obstype-search").val()},
		               dbname: "sched_table"},
    ajaxConfig:'GET',
    selectable: select,
    index:"_id",
    height:"400px", // set height of table (optional)
	layout:"fitColumns",
    columns:[ // define table columns
	    {title:"Name", field:"name", sorter:"number"},
        {title:"Arguments", field:"args"},
        {title:"Aborted", field:"aborted"},
        {title:"Start", field:"sch_start", sorter:"date",
         sorterParams:{format:"YYYY-MM-DD hh:mm:ssZZ"}},
        {title:"Stop", field:"sch_stop", sorter:"date",
         sorterParams:{format:"YYYY-MM-DD hh:mm:ssZZ"}},
        {title:"Modified (UTC)", field:"modified", sorter:"date",
         sorterParams:{format:"YYYY-MM-DD hh:mm:ssZZ"}},
    ],
	initialSort:[{column:"modified", dir:"desc"}]
  });
}


// get plotting mode and build the table
make_t_table(10);
make_schedule_table(10);
make_aux_table(10);
make_autoproc_table(10);


// wrapper for two search functions
function search() {
    var tab = document.getElementsByClassName("tablinks active")[0].id;
    if (tab == 'aux')
        asearch();
    else if (tab == 'scanify')
        tsearch();
    else if (tab == 'schedule')
        schedule_search();
    else if (tab == 'autoproc')
        autoprocsearch();
}

// sends search request and rebuilds table and plot types
function tsearch() {
    // clear selections in table
    $("#scanify_table").tabulator("deselectRow");
    // package up the search fields into a json to send to server
    querydata = {search: {source:      $("#obstype-search").val(),
			              date:        {min: $("#date-from").val(),
				                        max: $("#date-to").val()},
			              observation: {min: $("#obsid-from").val(),
					                    max: $("#obsid-to").val()}},
	             dbname: "scanify"};
    $("#scanify_table").tabulator("setData", "/dbpage", querydata);
    plot_list();
};

function asearch() {
    // package up the search fields into a json to send to server
    querydata = {search: {date:     {min: $("#date-from").val(),
				                     max: $("#date-to").val()},
			              filename: $("#file-search").val(),
			              type:     $("#auxtype-search").val()},
		         dbname: "aux_transfer"};
    $("#aux_table").tabulator("setData", "/dbpage", querydata);
    plot_list();
};

function autoprocsearch() {
    // package up the search fields into a json to send to server
    querydata = {search: {date:        {min: $("#date-from").val(),
				                        max: $("#date-to").val()},
			              modified:    {min: $("#modified-from").val(),
				                        max: $("#modified-to").val()},
			              observation: {min: $("#obsid-from").val(),
					                    max: $("#obsid-to").val()},
			              source: $("#obstype-search").val()},
			     dbname: "autoproc"};
	$("#autoproc_table").tabulator("setData", "/dbpage", querydata);
	plot_list();
};

function schedule_search() {
    querydata = {search: {modified:  {min: $("#modified-from").val(),
				                      max: $("#modified-to").val()},
                          sch_start: {min: $("#date-from").val(),
                                      max: $("#date-to").val()}},
		         dbname: "sched_table"};
    $("#schedule_table").tabulator("setData", "/dbpage", querydata);
    plot_list();
};


// create sse listener
var es = new EventSource("/sse");

// sends plot request to server and opens images in new window
function plot() {
  // holds selected observations
  var rows;
  // check if scanify or aux
  var tab = document.getElementsByClassName("tablinks active")[0];
  if (tab.id == "aux") {
    rows = $("#aux_table").tabulator("getSelectedData");
  }
  else if (tab.id == "scanify") {
    rows = $("#scanify_table").tabulator("getSelectedData");
  }
    else if (tab.id == "autoproc") {
    rows = $("#autoproc_table").tabulator("getSelectedData");
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
  display_win.images_to_load = selected_values.length * rows.length;
  display_win.loaded_count = 0;
  display_win.loading_progress = 0;

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
}
}

// used to sort image items
function compare(a, b) {
  return a.src.localeCompare(b.src);
}


