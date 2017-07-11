//create Tabulator on DOM element with id "example-table"
$("#example-table").tabulator({
  pagination:"remote",  // use this for normal pagination
  // progressiveRender:"remote",  // progressive rendering with AJAX requests; cool!!
  ajaxURL:"/page", //set the ajax URL
  ajaxParams: {source: $("#obstype-search").val(),
               date:  {min: $("#date-from").val(),
                       max: $("#date-to").val()},
               observation:    {min: $("#obsid-from").val(),
                       max: $("#obsid-to").val()}},
  ajaxConfig:'GET',
  paginationSize:40,
  index:"_id",
  height:"400px", // set height of table (optional)
  fitColumns:true, //fit columns to width of table (optional)
  columns:[ //Define Table Columns
    {title:"Observation ID", field:"observation", sorter:"number"},
    {title:"Source", field:"source"},
    //{title:"Path on scott", field:"path"},
    {title:"Fullrate status", field:"status_fullrate"},
    {title:"Downsampled status", field:"status_downsampled"},
    {title:"Fullrate transfer", field:"transfer_fullrate", formatter:"tickCross"},
    {title:"Downsampled transfer", field:"transfer_downsampled", formatter:"tickCross"},
    {title:"Date", field:"date", sorter:"date"}
    //{title:"Time", field:"time", sorter:"time"}
  ],
  //trigger an alert message when the row is clicked
  rowClick:function(e, row, id, obs){
    obsdata = row.getData()
    // get plot type
    var selected_values = [];    
    $("#plot_select :selected").each(function(){
      selected_values.push($(this).val()); 
    });
    if (selected_values.length == 0) {
        alert("Please select at least one plot type.");
        return;
    }
    obsdata['plot_type'] = selected_values.join(' ');
    // request the plot
    $.get("data_req", obsdata, function(err, status) {
      if (err != null) {
        alert(err)
        return;
      }
      // open display window and load images
      var w = window.open('display.html', '_blank');
      // need to add images after page has loaded
      var items = [
        {src: 'img/' + obsdata['source'] + obsdata['observation'] + selected_values[0] + '.png',}
                ];
      w.onload = function() {
        var items = [];
        for (var i = 0; i < selected_values.length; i++) {
          var img = new Image();
          img.src = 'img/' + obsdata['source'] + obsdata['observation'] + selected_values[i] + '.png';
          items.push({src: img.src, w: img.width, h: img.height});
        }
        w.start_image(items);
      };
    });
  },
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
