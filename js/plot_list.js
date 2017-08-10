// function to populate/update plot_list
function plot_list() {
  // type holds obs source and plotting mode
  var type = {"type": $("#obstype-search").val(),
      'func': $('input[name="func"]:checked').val(),
      'tab': document.getElementsByClassName("tablinks active")[0].id};
  if (document.getElementsByClassName("tablinks active")[0].id == 'aux')
    type['type'] = $("#auxtype-search").val();
  // load available plot types and add to multiselect
  $.get("plot_list", type, function( plots ) {
    // empty options
    $("#plot_select").empty();
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
  $(".chosen-select").chosen({width: '150px', disable_search: true});

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


// Constructs the list of sources from all unique sources in the transfer DB
function make_source_list() {
    $.get("sourcelist", function(sourceData) {
	    sourceList = [];
	    for(jsource = 0; jsource < sourceData.length; jsource++) {
		sourceList.push(sourceData[jsource].source);
	    };
	    sourceList.sort();

	    for(var jsource = 0; jsource < sourceList.length; jsource++) {
		var el = document.createElement("option");
		el.setAttribute('value', sourceList[jsource]);
		el.text = sourceList[jsource];
		$('#obstype-search').append(el);
	    };
	    $('#obstype-search').trigger('chosen:updated');
	});
}
