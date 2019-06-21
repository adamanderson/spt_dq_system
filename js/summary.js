function get_old_timeseries_plotdirs(interval) {
    // clear contents of the datalist div before we rebuild
    $('#datalist_'+interval+'_stats').empty()

    // now rebuild the div
    $.get('/oldstaticdirs', {subdirectory:"plots", interval:interval}, function(data, status) {
	    data.reverse();
	    for (jdir=0; jdir<data.length; jdir++) {
			datestring = data[jdir].split('/')[2];
			$('#datalist_'+interval+'_stats').append("<button class='btn' onclick=\"javascript:set_variable('weekdir', '" + data[jdir] + "')\";>" + datestring + "</button>");
	    }
	});
}

function get_old_map_plotdirs(interval) {
    // clear contents of the datalist div before we rebuild
    $('#datalist_'+interval+'_maps').empty()

    // now rebuild the div
    $.get('/oldstaticdirs', {subdirectory:"maps/figures", interval:interval}, function(data, status) {
	    data.reverse();
	    for (jdir=0; jdir<data.length; jdir++) {
			datestring = data[jdir].split('/')[3];
			$('#datalist_'+interval+'_maps').append("<button class='btn' onclick=\"javascript:set_variable('mapweekdir', '" + data[jdir] + "')\";>" + datestring + "</button>");
	    }
	});
}


// create sse listener
var es = new EventSource("/sse");

// loads the database information for recent data of a given type
function load_db_for_latestobs(sources, fcallback_on_output, rangetype) {
    $.each(sources, function(sourcename, plots) {
	    var plot_mode;
	    if (rangetype == 'lastobs') {
		querydata = {search: {date: "latest",
				      source: sourcename},
			     dbname: "scanify"};
		plot_mode = 'individual';
	    }
	    else if (rangetype == 'lastweek') {
		time_lastweek = moment().add(-14, 'days');
		time_now = moment();
		querydata = {search: {date: {min: time_lastweek.format('YYYY-MM-DD'),
					     max: time_now.format('YYYY-MM-DD')},
				      source: sourcename},
			     dbname: "scanify"};
		plot_mode = 'timeseries';
	    }
	    $.get('/dbpage', querydata, function(data, status) {
		    observation_list = [];
		    $.each(data, function(jobs, obsdata) {
			    observation_list.push(obsdata.observation.toString());
			});
		    console.log(data);
		    datareq = [];
		    $.each(plots, function(i, plot) {
			    datareq.push({observation: observation_list.join(' '),
					source: data[0].source,
					table: 'scanify',
					plot_type: plot,
					func: plot_mode});
			});
		    console.log(datareq);
		    fcallback_on_output(datareq);
		});
	});
}


// loads the plots for the latest set of observations of a given type
function load_plots(datareq) {
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
				    // this is a horrible kludge; solve with
				    // better message passing
				    for (iplot = 0; iplot < datareq.length; iplot++) {
					if (img.src.indexOf(datareq[iplot].plot_type) != -1) {
					    $('#'+datareq[iplot].plot_type).attr('src', img.src);
					}
				    }
				});
			}
		    }
		});
	    es.addEventListener('err' + sseid, function (event) {
		    console.log(event.data);
		});

	    // loop over each observation
	    $.each(datareq, function(i, obsdata) {
		    obsdata['sseid'] = sseid;
		    
		    // request the plots
		    $.get("data_req", obsdata, function(data, status) {
			    nTotalPlots += 1
			});
		});
	});
}
