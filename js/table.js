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
  selectable: 4,
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


function search() {
  // package up the search fields into a json to send to server
  querydata = {page:1,
               size:20,
               source: $("#obstype-search").val(),
               date:  {min: $("#date-from").val(),
                       max: $("#date-to").val()},
               observation:    {min: $("#obsid-from").val(),
                       max: $("#obsid-to").val()}}
  $("#example-table").tabulator("setData", "/page", querydata);
  // load available plot types
  $("#plot_select").empty();
  var type = {"type": $("#obstype-search").val()};
  $.get("plot_list", type, function( plots ) {
    for (var i = 0; i < plots.length; i++) {
      var option = document.createElement("option");
      option.setAttribute("value", plots[i]);
      option.text = plots[i];
      $("#plot_select").append(option);
    }
    // let chosen know the options have changed
    $("#plot_select").trigger("chosen:updated");
  });
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

  // items holds the image info needed for photoswipe
  var items = [];
  // deferred holds the jquery requests so we know when they are done
  var deferred = [];
  // loop over each observation
  $.each(rows, function(i, obsdata) {
    // combine all requested plot types into a string
    obsdata['plot_type'] = selected_values.join(' ');
    // request the plot
    deferred.push($.get("data_req", obsdata, function(err, status) {
      if (err != null) {
        alert(err)
        return;
      }
      for (var i = 0; i < selected_values.length; i++) {
        src = 'img/' + obsdata['source'] + obsdata['observation'] + selected_values[i] + '.png';
        items.push({src: src, w: 0, h: 0, obs: obsdata['observation']});
      }
    }));
  });
  // after all plots have finished
  $.when.apply(null, deferred).then( function() {
    // open display window and load images
    var w = window.open('display.html', '_blank');
    // need to add images after page has loaded
    w.onload = function() {
      // sort by observation
      items.sort(compare)
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


