function open_table(table) {
    $('#recent').hide();
    $('#lastobs').hide();
    $('#'+table).show();

    if(table == "lastobs") {
	load_db_for_latestobs('calibrator');
	load_latestobs_plots();
    }
}

// create sse listener
var es = new EventSource("/sse");

// loads the database information for the latest observation of a given type
function load_db_for_latestobs(sourcename) {
    querydata = {search: {date: "latest",
                          source: sourcename},
                 dbname: "transfer"};
    $.get('/dbpage', querydata, function(data, status) {
	    console.log(data);
	});
}

// loads the plots for the latest set of observations of a given type
function load_latestobs_plots() {
    // toy request for demonstration purposes
    datareq = [{observation: '32537722',
		source: 'calibrator',
		table: 'transfer',
		plot_type: 'CalHistogram',
		func:'individual'}];
    
    var items = [];
    var nTotalPlots = 0;
    var sseid;

    // gets the IDs for each plot to be produced and keeps track of the number
    // that have been created
    $.get("sseid", function(id, status) {
	    sseid = id;
	    // counts the images that have finished being created         
	    var plot_ctr = 0;

	    es.addEventListener('out' + sseid, function (event) {
		    if (event.data.split('fln').length == 2) {
			path =  event.data.split('fln')[1].slice(0,-1);
			items.push({src: 'img/' + path, w: 0, h: 0});
		    }
		    else if (event.data.slice(1, 4) == 'plt') {
			plot_ctr++;
			
			if (nTotalPlots == plot_ctr) {
			    $.each(items, function(i, img) {
				    console.log(img);
				    $('#calSN').attr('src', img.src);
				});
			}
		    }
		});
	    es.addEventListener('err' + sseid, function (event) {
		    display_win.show_err(event.data);
		});

	    // loop over each observation
	    $.each(datareq, function(i, obsdata) {
		    obsdata['sseid'] = sseid;
		    
		    // request the plots
		    $.get("data_req", obsdata, function(data, status) {
			    nTotalPlots += datareq.length;
			});
		});
	});
}
