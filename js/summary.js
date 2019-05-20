function get_timeseries_plotdirs(interval) {
    // clear contents of the datalist div before we rebuild
    $('#datalist_'+interval).empty()

    // now rebuild the div
    $.get('/oldstaticdirs', {interval:interval}, function(data, status) {
	    data.reverse();
	    for (jdir=0; jdir<data.length; jdir++) {
			datestring = data[jdir].split('/')[2];
			if (interval == 'last_3days') {
				$('#datalist_'+interval).append("<button class='btn' onclick=\"javascript:set_variable('weekdir', '" + data[jdir] + "')\";>last 3 days</button>");
			}
			else if (datestring != 'current') {
				$('#datalist_'+interval).append("<button class='btn' onclick=\"javascript:set_variable('weekdir', '" + data[jdir] + "')\";>" + datestring + "</button>");
			}
	    }
	});
}
