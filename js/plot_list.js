// make obs select chosen
$(".chosen-select").chosen({disable_search: true});

// populate the drop down menu of plot types
var plot_div = document.getElementById("plot-type");
var plot_select = document.createElement('select');
plot_select.setAttribute("id", "plot_select");
plot_select.setAttribute("class", "chosen-multiselect");
plot_select.setAttribute("multiple", "");
plot_div.appendChild(plot_select);

// note: do not have to decode "plots" which comes in json            
var type = {'type': $("#obstype-search").val()};
$.get("plot_list", type, function( plots ) {
  for (var i = 0; i < plots.length; i++) {
    var option = document.createElement("option");
    option.setAttribute("value", plots[i]);
    option.text = plots[i];
    plot_select.appendChild(option);
  }
  $(".chosen-multiselect").chosen({
    width: "90%",
    placeholder_text_multiple: "Select Plot Type",
  });
});
