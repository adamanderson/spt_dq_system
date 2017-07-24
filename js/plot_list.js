// function to populate/update plot_list
function plot_list() {
  $("#plot_select").empty();
  // type holds obs source and plotting mode
  var type = {"type": $("#obstype-search").val(),
      'func': $('input[name="func"]:checked').val()};
  // load available plot types and add to multiselect
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
}

// need to delay creating until page is loaded so display is active
function make_plot_list() {
  // make obs select chosen
  $(".chosen-select").chosen({disable_search: true});

  // create the drop down menu of plot types
  var plot_div = document.getElementById("plot-type");
  var plot_select = document.createElement('select');
  plot_select.setAttribute("id", "plot_select");
  plot_select.setAttribute("class", "chosen-multiselect");
  plot_select.setAttribute("multiple", "");
  plot_div.appendChild(plot_select);

  // set multiselect options
  $(".chosen-multiselect").chosen({
    width: "90%",
    placeholder_text_multiple: "Select Plot Type",
  });

  //create initial plot list
  plot_list();
}
