// check if cookies are defined, otherwise set them to defaults
var default_cookies = {'wafer': 'all',
                       'timedir_calibration': 'plots/last_n/last_07/',
                       'timedir_winter': 'maps/figures_winter/last_n/last_07/',
                       'timedir_summer': 'maps/figures_summer/last_n/last_07/',
                       'timedir_summerb': 'maps/figures_summerb/last_n/last_07/',
                       'arcdir': 'arcs/figs/last_n/last_07/',
                       'cycledir': 'arcs/figs/cycles/newest'};
for (var cookie_name in default_cookies) {
    if (Cookies.get(cookie_name) == undefined) 
        Cookies.set(cookie_name, default_cookies[cookie_name], {expires: 1});
}

// get cookie values
var wafer = Cookies.get('wafer');
var timedir = {calibration: Cookies.get('timedir_calibration'),
               winter: Cookies.get('timedir_winter'),
               summer: Cookies.get('timedir_summer'),
               summerb: Cookies.get('timedir_summerb'),
               weather: Cookies.get('arcdir'),
               fridgecycle: Cookies.get('cycledir')};

/** 
 * Retrieves a precompiled template, or if it does not exist,
 * looks for a template in the `templates` directory and compiles it.
 * @param {string} name Name of template to search for
 * @return {Handlebars.templates} 
 */
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

// Compile the templates, at least one per tab
var map_tab_names = ['winter', 'summer', 'summerb'];
for(var jtab = 0; jtab < map_tab_names.length; jtab++)
{
    // Loop over map tabs. Each tab has a template for figures common to all
    // map tabs, pointing figures that differ between tabs, and the time 
    // selector.
    var compiled_map_template = Handlebars.getTemplate('maptab');
    var compiled_time_selector_template = Handlebars.getTemplate('time_selector');
    var compiled_map_pointing_template = Handlebars.getTemplate('map_pointing_' + map_tab_names[jtab]);

}
var compiled_calibration_template = Handlebars.getTemplate('calibration_tab');
var compiled_weather_template = Handlebars.getTemplate('weather_tab');
var compiled_fridgecycle_template = Handlebars.getTemplate('fridgecycle_tab');

/**
 * Updates the content in a tab by refilling HTML templates that comprise the
 * tab.
 * @param {string} name Name of tab to update
 * @param {bool} init Do initialization steps that should occur only on the
 * first time the lab is loaded
 */
function update_tab(name, init) {
    // prototype context
    var context = { wafer: wafer,
                    map_tab_name: name,
                    week_dir: timedir[name],
                    yearly: timedir[name].includes('yearly'), 
                    recent_selector: true, 
                    yearly_selector: true, 
                    monthly_selector: true, 
                    weekly_selector: true,
                    cycles_selector: true};

    // customizations for each tab
    if(name == 'calibration')
    {
        context.weekly_selector = false;
        context.cycles_selector = false;

        var html = compiled_calibration_template(context);
        $("#figs_calibration").html(html);
    }
    else if(name == 'winter' || name == 'summer' || name == 'summerb')
    {
        context.weekly_selector = false;
        context.cycles_selector = false;

        var html = compiled_map_template(context);
        $("#figs_maps" + name).html(html);
        var html = compiled_map_pointing_template(context);
        $('#maps_' + name + '_pointing').append(html)
    }
    else if(name == 'weather')
    {
        context.yearly_selector = false;
        context.cycles_selector = false;

        var html = compiled_weather_template(context);
        $("#figs_weather").html(html);
    }
    else if(name == 'fridgecycle')
    {
        context.recent_selector = false;
        context.yearly_selector = false;
        context.monthly_selector = false;
        context.weekly_selector = false;
        var html = compiled_fridgecycle_template(context);
        $("#figs_fridgecycle").html(html);
        
    }

    if(init)
    {
        var html = compiled_time_selector_template(context);
        $("#time_selector_" + name).html(html);

        // Initialize jQuery UI elements and make dynamic modifications to the DOM
        if(name == 'calibration')
        {
            add_date_buttons('last_n', 'plots', 'calibration');
            add_date_buttons('yearly', 'plots', 'calibration');
            add_date_buttons('monthly', 'plots', 'calibration');
            add_date_buttons('weekly', 'plots', 'calibration');
        }
        else if(name == 'winter')
        {
            add_date_buttons('last_n', 'maps/figures', 'winter');
            add_date_buttons('yearly', 'maps/figures', 'winter');
            add_date_buttons('monthly', 'maps/figures', 'winter');
            add_date_buttons('weekly', 'maps/figures', 'winter');
        }
        else if(name == 'summer')
        {
            add_date_buttons('last_n', 'maps/figures_summer', 'summer');
            add_date_buttons('yearly', 'maps/figures_summer', 'summer');
            add_date_buttons('monthly', 'maps/figures_summer', 'summer');
            add_date_buttons('weekly', 'maps/figures_summer', 'summer');
        }
        else if(name == 'summerb')
        {
            add_date_buttons('last_n', 'maps/figures_summerb', 'summerb');
            add_date_buttons('yearly', 'maps/figures_summerb', 'summerb');
            add_date_buttons('monthly', 'maps/figures_summerb', 'summerb');
            add_date_buttons('weekly', 'maps/figures_summerb', 'summerb');
        }
        else if(name == 'weather')
        {
            add_date_buttons('last_n', 'arcs/figs', 'weather');
            add_date_buttons('monthly', 'arcs/figs', 'weather');
            add_date_buttons('weekly', 'arcs/figs', 'weather');
        }
        else if(name == 'fridgecycle')
        {
            add_date_buttons('cycles', 'arcs/figs', 'fridgecycle');
        }
    }
}

/**
 * Updates the "last modified" field at the top of each tab. This could
 * probably be combined directly with update_tab.
 */
function update_lastmodified() {
    $.get('/lastmodified_calibration', function(data) {
        $("#lastmodified_calibration").replaceWith('Plots last modified: ' + data.time + ' (UTC)');
    });
    $.each(map_tab_names, function(jtab, tab_name)
    {
        $.get('/lastmodified_maps', {field: tab_name}, function(data) {
            $("#lastmodified_maps_" + tab_name).replaceWith('Plots last modified: ' + data.time + ' (UTC)');
        });
    });
    $.get('/lastmodified_weather', function(data) {
        $("#lastmodified_weather").replaceWith('Plots last modified: ' + data.time + ' (UTC)');
    });
}

/**
 * Builds the buttons on the summary page that select different time intervals
 * of data to display. These are constructed by appending the DOM directly.
 * Note also that this function also initializes the jQuery UI "controlgroup"
 * and binds the click event to it after appending all the buttons.
 * @param {string} interval Time interval to traverse, 'weekly' or 'monthly'
 * @param {string} interval Subdirectory that we should traverse to get dated data
 * @param {string} interval Tab in which to add buttons, 'summary' or 'maps'
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
                  else if (tab == 'fridgecycle')
                      datestring = data[jdir].split('/')[3];

                  $(div_id).append("<input type='radio' id='dates-"+tab+"-"+datestring+"' name='dates' value='"+data[jdir]+"'>\n"+  "<label for='dates-"+tab+"-"+datestring+"'>"+datestring+"</label>");
              }
              $(div_id).controlgroup();
              $("[id^=dates-"+tab+"-]").click(function(event) {
                timedir[tab] = event.target.value;
                update_tab(tab, false);
              });
          });
}


// Page initialization
$( document ).ready(function()
{
    // initialize the tabs in jquery-ui with a function that updates the tab
    // when clicked
    $("#tabs").tabs({
        activate: function( event, ui )
        {
            var tab_name = $("#tabs .ui-state-active a").attr('href');
            tab_name = tab_name.replace('#tabs-','');
            update_tab(tab_name, false);
            update_lastmodified();

            Cookies.set('activeTab', $("#tabs").tabs("option", "active"), {expires: 1});
        }
    });

    // Bind the click event to the wafer buttons
    $("#waferlist").controlgroup();
    $("[id^=wafers-]").click(function(event) {
        wafer = event.target.value;
        update_tab('calibration', false);
    });

    // load all the tab data
    $.each(timedir, function(tab_name, val) 
    {
        update_tab(tab_name, true);
        update_lastmodified();
    });
});