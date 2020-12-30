// check if cookies are defined, otherwise set them to defaults
var default_cookies = {'wafer': 'all',
                       'weekdir': 'plots/last_n/last_07/',
                       'mapweekdir': 'maps/figures_winter/last_n/last_07/',
                       'mapweekdirsummer': 'maps/figures_summer/last_n/last_07/',
                       'mapweekdirsummerb': 'maps/figures_summerb/last_n/last_07/',
                       'arcdir': 'arcs/figs/last_n/last_07/',
                       'cycledir': 'arcs/figs/cycles/newest'};
for (var cookie_name in default_cookies) {
    if (Cookies.get(cookie_name) == undefined) 
        Cookies.set(cookie_name, default_cookies[cookie_name], {expires: 1});
}

// get cookie values
var wafer = Cookies.get('wafer');
var weekdir = {calibration: Cookies.get('weekdir'),
               winter: Cookies.get('mapweekdir'),
               summer: Cookies.get('mapweekdirsummer'),
               summerb: Cookies.get('mapweekdirsummerb'),
               weather: Cookies.get('arcdir')};
var cycledir = Cookies.get('cycledir');

// function for getting the templates used to fill out the map tabs
Handlebars.getTemplate = function(name) {
    $.ajax({
        url : 'templates/' + name + '.handlebars',
        success : function(data) {
            if (Handlebars.templates === undefined) {
                Handlebars.templates = {};
            }
            Handlebars.templates[name] = Handlebars.compile(data);
        },
        async : false
    });
    return Handlebars.templates[name];
};

// tab names
var map_tab_names = ['winter', 'summer', 'summerb'];
for(var jtab = 0; jtab < map_tab_names.length; jtab++)
{
    var compiled_map_template = Handlebars.getTemplate('maptab');
    var compiled_time_selector_template = Handlebars.getTemplate('time_selector');
    var compiled_map_pointing_template = Handlebars.getTemplate('map_pointing_' + map_tab_names[jtab]);

}
var compiled_calibration_template = Handlebars.getTemplate('calibration_tab');
var compiled_weather_template = Handlebars.getTemplate('weather_tab');

/**
 * Updates the source of all images defined on the summary page.
 */
function update_figs() {
    // Compile the templates for the tabs
    for(var jtab = 0; jtab < map_tab_names.length; jtab++)
    {
        var context = { 'map_tab_name' : map_tab_names[jtab],
                        'map_week_dir': weekdir[map_tab_names[jtab]],
                        'yearly': weekdir[map_tab_names[jtab]].includes('yearly') };

        var html = compiled_map_template(context);
        $("#figs_maps" + map_tab_names[jtab]).html(html);

        var html = compiled_map_pointing_template(context);
        $('#maps_' + map_tab_names[jtab] + '_pointing').html(html);
    }

    var context = { 'week_dir': weekdir.calibration,
                    'wafer': wafer,
                    'yearly': weekdir.calibration.includes('yearly') };
    var html = compiled_calibration_template(context);
    $("#figs_calibration").html(html);

    var context = { 'arc_dir': weekdir.weather};
    var html = compiled_weather_template(context);
    $("#figs_weather").html(html);

    document["cycle_he10_full"].src = 'staticimg/'+cycledir+'/he10_full.png';
    document["cycle_he10_half"].src = 'staticimg/'+cycledir+'/he10_half.png';
    document["cycle_low_temps"].src = 'staticimg/'+cycledir+'/low_temps.png';
}


function update_lastmodified() {
    $.get('/lastmodified_calibration', function(data) {
        $("#lastmodified_calibration").replaceWith('Plots last modified: ' + data.time + ' (UTC)');
    });
    for(var jtab = 0; jtab < map_tab_names.length; jtab++)
    {
        $.get('/lastmodified_maps', {field: map_tab_names[jtab]}, function(data) {
            $("#lastmodified_maps").replaceWith('Plots last modified: ' + data.time + ' (UTC)');
        });
    }
    $.get('/lastmodified_weather', function(data) {
        $("#lastmodified_weather").replaceWith('Plots last modified: ' + data.time + ' (UTC)');
    });
}

/**
 * Sets a global variable to records its value as a cookie for retrieval later.
 */
function set_variable(variable, newVal)
{
  variable = newVal;
  Cookies.set(variable, newVal, { expires: 1 });
  update_figs();
}

/**
 * Builds the buttons on the summary page that select different time intervals
 * of data to display. These are constructed by appending the DOM directly.
 * Note also that this function also initializes the jQuery UI "controlgroup"
 * and binds the click event to it after appending all the buttons.
 * 
 * Arguments:
 *  interval : time interval to traverse, 'weekly' or 'monthly'
 *  subdirectory : subdirectory that we should traverse to get dated data
 *  tab : tab in which to add buttons, 'summary' or 'maps'
 */
function add_date_buttons(interval, subdirectory, tab)
{
    // now rebuild the div
    $.get('/staticdirs', {subdirectory:subdirectory, interval:interval},
          function(data, status) {
              div_id = '#datalist_'+tab+'_'+interval;
             
              data.reverse();
              for (jdir=0; jdir<data.length; jdir++)
              {
                  if (tab == 'calibration')
                      datestring = data[jdir].split('/')[2];
                  else if (tab == 'winter')
                      datestring = data[jdir].split('/')[3];
                  else if (tab == 'summer')
                      datestring = data[jdir].split('/')[3];
                  else if (tab == 'summerb')
                      datestring = data[jdir].split('/')[3];
                  else if (tab == 'weather')
                      datestring = data[jdir].split('/')[3];
                  else if (tab == 'fridge')
                      datestring = data[jdir].split('/')[3];

                  $(div_id).append("<input type='radio' id='dates-"+tab+"-"+datestring+"' name='dates' value='"+data[jdir]+"'>\n"+  "<label for='dates-"+tab+"-"+datestring+"'>"+datestring+"</label>");

                  // can use jquery <select ...></select> with this line instead of radio buttons
                  // $(div_id).append("<option value='"+data[jdir]+"'>"+datestring+"</option>");
              }
              $(div_id).controlgroup();
              $("[id^=dates-"+tab+"-]").click(function(event) {
                weekdir[tab] = event.target.value;
                // Cookies.set(variable, newVal, { expires: 1 });
                update_figs();
                  //set_variable(boundVar, event.target.value);
              });
          });
}

// Page initialization
$( document ).ready(function()
{
    // Compile the templates for the tabs
    for(var jtab = 0; jtab < map_tab_names.length; jtab++)
    {
        var context = { 'map_tab_name' : map_tab_names[jtab],
                    'map_week_dir': weekdir[map_tab_names[jtab]],
                    'yearly': weekdir[map_tab_names[jtab]].includes('yearly'), 
                    'recent_selector': true, 
                    'yearly_selector': true, 
                    'monthly_selector': true, 
                    'weekly_selector': true};

        var html = compiled_time_selector_template(context);
        $("#time_selector_maps" + map_tab_names[jtab]).html(html);

        var html = compiled_map_template(context);
        $("#figs_maps" + map_tab_names[jtab]).html(html);

        var html = compiled_map_pointing_template(context);
        $('#maps_' + map_tab_names[jtab] + '_pointing').append(html)
    }

    var context = { 'map_tab_name': 'calibration',
                    'week_dir': weekdir.calibration,
                    'wafer': wafer,
                    'yearly': weekdir.calibration.includes('yearly'), 
                    'recent_selector': true, 
                    'yearly_selector': true, 
                    'monthly_selector': true, 
                    'weekly_selector': true};
    var html = compiled_time_selector_template(context);
    $("#time_selector_calibration").html(html);
    var html = compiled_calibration_template(context);
    $("#figs_calibration").html(html);

    var context = { 'map_tab_name': 'weather',
                    'arc_dir': weekdir.weather, 
                    'recent_selector': true, 
                    'yearly_selector': false, 
                    'monthly_selector': true, 
                    'weekly_selector': true};
    var html = compiled_time_selector_template(context);
    $("#time_selector_weather").html(html);
    var html = compiled_weather_template(context);
    $("#figs_weather").html(html);
    
    // Initialize jQuery UI elements and make dynamic modifications to the DOM
    $("#tabs").tabs();
    $("#waferlist").controlgroup();
    add_date_buttons('last_n', 'plots', 'calibration');
    add_date_buttons('yearly', 'plots', 'calibration');
    add_date_buttons('monthly', 'plots', 'calibration');
    add_date_buttons('weekly', 'plots', 'calibration');

    add_date_buttons('last_n', 'maps/figures', 'winter');
    add_date_buttons('yearly', 'maps/figures', 'winter');
    add_date_buttons('monthly', 'maps/figures', 'winter');
    add_date_buttons('weekly', 'maps/figures', 'winter');

    add_date_buttons('last_n', 'maps/figures_summer', 'summer');
    add_date_buttons('yearly', 'maps/figures_summer', 'summer');
    add_date_buttons('monthly', 'maps/figures_summer', 'summer');
    add_date_buttons('weekly', 'maps/figures_summer', 'summer');

    add_date_buttons('last_n', 'maps/figures_summerb', 'summerb');
    add_date_buttons('yearly', 'maps/figures_summerb', 'summerb');
    add_date_buttons('monthly', 'maps/figures_summerb', 'summerb');
    add_date_buttons('weekly', 'maps/figures_summerb', 'summerb');

    add_date_buttons('last_n', 'arcs/figs', 'weather');
    add_date_buttons('monthly', 'arcs/figs', 'weather');
    add_date_buttons('weekly', 'arcs/figs', 'weather');
    
    add_date_buttons('cycles', 'arcs/figs', 'fridge', 'cycledir');

    // Bind the click event to the wafer buttons
    $("[id^=wafers-]").click(function(event) {
        wafer = event.target.value;
        update_figs();
    });

    // Bind the click event to the tabs so that we can keep track of the active tab
    $("#tabs").click(function(event) {
        Cookies.set('activeTab', $("#tabs").tabs("option", "active"), {expires: 1});
    });

    if(Cookies.get('activeTab') !== undefined)
        $("#tabs").tabs("option", "active", Cookies.get('activeTab'))
    else {
        Cookies.set('activeTab', $("#tabs").tabs("option", "active"), {expires: 1});
    }


    // Update all the figures. Need to call this on load because the values of
    // wafer and weekdir might be pulled from a cookie, which probably differs
    // from the default values.
    update_figs();
    update_lastmodified();
});
