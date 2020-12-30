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
// var weekdir = Cookies.get('weekdir');
var weekdir = {summary: Cookies.get('weekdir'),
               winter: Cookies.get('mapweekdir'),
               summer: Cookies.get('mapweekdirsummer'),
               summerb: Cookies.get('mapweekdirsummerb')};
var arcdir = Cookies.get('arcdir');
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
    var context = { 'map_tab_name' : map_tab_names[jtab],
                    'map_week_dir': weekdir[map_tab_names[jtab]]};
    var compiled_map_template = Handlebars.getTemplate('maptab');
    var compiled_time_selector_template = Handlebars.getTemplate('time_selector');
    var compiled_map_pointing_template = Handlebars.getTemplate('map_pointing_' + map_tab_names[jtab]);
}

/**
 * Updates the source of all images defined on the summary page.
 */
function update_figs() {
    // Compile the templates for the tabs
    console.log(weekdir);
    for(var jtab = 0; jtab < map_tab_names.length; jtab++)
    {
        var context = { 'map_tab_name' : map_tab_names[jtab],
                        'map_week_dir': weekdir[map_tab_names[jtab]],
                        'yearly': weekdir[map_tab_names[jtab]].includes('yearly')};

        var html = compiled_map_template(context);
        $("#figs_maps" + map_tab_names[jtab]).html(html)

        var html = compiled_map_pointing_template(context);
        $('#maps_' + map_tab_names[jtab] + '_pointing').html(html)
    }

    document["ts_median_cal_sn_highel"].src        = 'staticimg/'+weekdir+'/median_cal_sn_4Hz_highel_'+wafer+'.png';
    document["ts_median_cal_response_highel"].src  = 'staticimg/'+weekdir+'/median_cal_response_4Hz_highel_'+wafer+'.png';
    document["ts_alive_bolos_cal_highel"].src      = 'staticimg/'+weekdir+'/alive_bolos_cal_4Hz_highel_'+wafer+'.png';
      
    document["ts_median_cal_sn_lowel"].src         = 'staticimg/'+weekdir+'/median_cal_sn_4Hz_lowel_'+wafer+'.png';
    document["ts_median_cal_response_lowel"].src   = 'staticimg/'+weekdir+'/median_cal_response_4Hz_lowel_'+wafer+'.png';
    document["ts_alive_bolos_cal_lowel"].src       = 'staticimg/'+weekdir+'/alive_bolos_cal_4Hz_lowel_'+wafer+'.png';
    
    document["ts_median_elnod_slope"].src          = 'staticimg/'+weekdir+'/median_elnod_response_'+wafer+'.png';
    document["ts_median_elnod_sn"].src             = 'staticimg/'+weekdir+'/median_elnod_sn_'+wafer+'.png';
    document["ts_median_elnod_opacity"].src        = 'staticimg/'+weekdir+'/median_elnod_opacity_'+wafer+'.png';
    document["ts_median_elnod_iq"].src             = 'staticimg/'+weekdir+'/median_elnod_iq_phase_'+wafer+'.png';
    document["ts_alive_bolos_elnod"].src           = 'staticimg/'+weekdir+'/alive_bolos_elnod_'+wafer+'.png';
    
    document["ts_rcw38_sky_transmission"].src      = 'staticimg/'+weekdir+'/rcw38_sky_transmission_'+wafer+'.png';
    document["ts_median_rcw38_fluxcal"].src        = 'staticimg/'+weekdir+'/median_rcw38_fluxcal_'+wafer+'.png';
    document["ts_median_rcw38_intflux"].src        = 'staticimg/'+weekdir+'/median_rcw38_intflux_'+wafer+'.png';
    document["ts_alive_bolos_rcw38"].src           = 'staticimg/'+weekdir+'/alive_bolos_rcw38_'+wafer+'.png';
    
    document["ts_mat5a_sky_transmission"].src      = 'staticimg/'+weekdir+'/mat5a_sky_transmission_'+wafer+'.png';
    document["ts_median_mat5a_fluxcal"].src        = 'staticimg/'+weekdir+'/median_mat5a_fluxcal_'+wafer+'.png';
    document["ts_median_mat5a_intflux"].src        = 'staticimg/'+weekdir+'/median_mat5a_intflux_'+wafer+'.png';
    document["ts_alive_bolos_mat5a"].src           = 'staticimg/'+weekdir+'/alive_bolos_mat5a_'+wafer+'.png';

    document["ts_w28a2_sky_transmission"].src      = 'staticimg/'+weekdir+'/w28a2_sky_transmission_'+wafer+'.png';
    document["ts_median_w28a2_fluxcal"].src        = 'staticimg/'+weekdir+'/median_w28a2_fluxcal_'+wafer+'.png';
    document["ts_median_w28a2_intflux"].src        = 'staticimg/'+weekdir+'/median_w28a2_intflux_'+wafer+'.png';
    document["ts_alive_bolos_w28a2"].src           = 'staticimg/'+weekdir+'/alive_bolos_w28a2_'+wafer+'.png';

    document["ts_iras17258_sky_transmission"].src  = 'staticimg/'+weekdir+'/iras17258_sky_transmission_'+wafer+'.png';
    document["ts_median_iras17258_fluxcal"].src    = 'staticimg/'+weekdir+'/median_iras17258_fluxcal_'+wafer+'.png';
    document["ts_median_iras17258_intflux"].src    = 'staticimg/'+weekdir+'/median_iras17258_intflux_'+wafer+'.png';
    document["ts_alive_bolos_iras17258"].src       = 'staticimg/'+weekdir+'/alive_bolos_iras17258_'+wafer+'.png';

    document["focus_fwhm_vs_bench_90"].src         = 'staticimg/'+weekdir+'/focus_fwhm_vs_bench_90.png';
    document["focus_fwhm_vs_bench_150"].src        = 'staticimg/'+weekdir+'/focus_fwhm_vs_bench_150.png';
    document["focus_fwhm_vs_bench_220"].src        = 'staticimg/'+weekdir+'/focus_fwhm_vs_bench_220.png';

    document["focus_ellip_vs_bench_90"].src        = 'staticimg/'+weekdir+'/focus_ellip_vs_bench_90.png';
    document["focus_ellip_vs_bench_150"].src       = 'staticimg/'+weekdir+'/focus_ellip_vs_bench_150.png';
    document["focus_ellip_vs_bench_220"].src       = 'staticimg/'+weekdir+'/focus_ellip_vs_bench_220.png';

    document["focus_min_fwhm"].src                 = 'staticimg/'+weekdir+'/focus_FWHM.png';
    document["focus_min_fwhm_bench"].src           = 'staticimg/'+weekdir+'/focus_BenchF.png';
    document["focus_min_ellip"].src                = 'staticimg/'+weekdir+'/focus_Ellipticity.png';
    document["focus_min_ellip_bench"].src          = 'staticimg/'+weekdir+'/focus_BenchE.png';

    document["ts_median_net_01Hz_to_05Hz"].src     = 'staticimg/'+weekdir+'/median_NET_0.1Hz_to_0.5Hz_'+wafer+'.png';
    document["ts_median_net_1Hz_to_2Hz"].src       = 'staticimg/'+weekdir+'/median_NET_1.0Hz_to_2.0Hz_'+wafer+'.png';
    document["ts_median_net_3Hz_to_5Hz"].src       = 'staticimg/'+weekdir+'/median_NET_3.0Hz_to_5.0Hz_'+wafer+'.png';
    document["ts_median_net_10Hz_to_15Hz"].src     = 'staticimg/'+weekdir+'/median_NET_10.0Hz_to_15.0Hz_'+wafer+'.png';
    
    document["ts_median_nep_01Hz_to_05Hz"].src     = 'staticimg/'+weekdir+'/median_NEP_0.1Hz_to_0.5Hz_'+wafer+'.png';
    document["ts_median_nep_1Hz_to_2Hz"].src       = 'staticimg/'+weekdir+'/median_NEP_1.0Hz_to_2.0Hz_'+wafer+'.png';
    document["ts_median_nep_3Hz_to_5Hz"].src       = 'staticimg/'+weekdir+'/median_NEP_3.0Hz_to_5.0Hz_'+wafer+'.png';
    document["ts_median_nep_10Hz_to_15Hz"].src     = 'staticimg/'+weekdir+'/median_NEP_10.0Hz_to_15.0Hz_'+wafer+'.png';
    
    document["ts_median_nei_01Hz_to_05Hz"].src     = 'staticimg/'+weekdir+'/median_NEI_0.1Hz_to_0.5Hz_'+wafer+'.png';
    document["ts_median_nei_1Hz_to_2Hz"].src       = 'staticimg/'+weekdir+'/median_NEI_1.0Hz_to_2.0Hz_'+wafer+'.png';
    document["ts_median_nei_3Hz_to_5Hz"].src       = 'staticimg/'+weekdir+'/median_NEI_3.0Hz_to_5.0Hz_'+wafer+'.png';
    document["ts_median_nei_10Hz_to_15Hz"].src     = 'staticimg/'+weekdir+'/median_NEI_10.0Hz_to_15.0Hz_'+wafer+'.png';

    document["ts_number_of_lines"].src             = 'staticimg/'+weekdir+'/number_of_lines_found_'+wafer+'.png';
    


    document["tipper_tau"].src   = 'staticimg/'+arcdir+'/tipper_tau.png';
    document["tipper_tatm"].src  = 'staticimg/'+arcdir+'/tipper_tatm.png';
    document["elnod_tau_cp"].src = 'staticimg/'+arcdir+'/median_elnod_opacity_all.png';

    document["weather_tatm"].src = 'staticimg/'+arcdir+'/weather_tatm.png';
    document["weather_humi"].src = 'staticimg/'+arcdir+'/weather_humi.png';
    document["weather_pres"].src = 'staticimg/'+arcdir+'/weather_pres.png';
    document["weather_wmag"].src = 'staticimg/'+arcdir+'/weather_wmag.png';
    document["weather_wdir"].src = 'staticimg/'+arcdir+'/weather_wdir.png';

    document["cabin_temps"].src  = 'staticimg/'+arcdir+'/cabin_temps.png';
    document["4K_temps"].src     = 'staticimg/'+arcdir+'/fourk_temps.png';
    document["50K_temps"].src    = 'staticimg/'+arcdir+'/fiftyk_temps.png';

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
  //window[variable] = newVal;
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
function add_date_buttons(interval, subdirectory, tab, boundVar)
{
    // now rebuild the div
    $.get('/staticdirs', {subdirectory:subdirectory, interval:interval},
          function(data, status) {
              div_id = '#datalist_'+tab+'_'+interval;
             
              data.reverse();
              for (jdir=0; jdir<data.length; jdir++)
              {
                  if (tab == 'summary')
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
                        'map_week_dir': weekdir[map_tab_names[jtab]]};
        console.log(weekdir[map_tab_names[jtab]]);

        var html = compiled_time_selector_template(context);
        $("#time_selector_maps" + map_tab_names[jtab]).html(html)

        var html = compiled_map_template(context);
        $("#figs_maps" + map_tab_names[jtab]).html(html)

        var html = compiled_map_pointing_template(context);
        $('#maps_' + map_tab_names[jtab] + '_pointing').append(html)
    }
    
    // Initialize jQuery UI elements and make dynamic modifications to the DOM
    $("#tabs").tabs();
    $("#waferlist").controlgroup();
    add_date_buttons('last_n', 'plots', 'summary', 'weekdir');
    add_date_buttons('yearly', 'plots', 'summary', 'weekdir');
    add_date_buttons('monthly', 'plots', 'summary', 'weekdir');
    add_date_buttons('weekly', 'plots', 'summary', 'weekdir');

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

    add_date_buttons('last_n', 'arcs/figs', 'weather', 'arcdir');
    add_date_buttons('monthly', 'arcs/figs', 'weather', 'arcdir');
    add_date_buttons('weekly', 'arcs/figs', 'weather', 'arcdir');
    
    add_date_buttons('cycles', 'arcs/figs', 'fridge', 'cycledir');

    // Bind the click event to the wafer buttons
    $("[id^=wafers-]").click(function(event) {
        set_variable("wafer", event.target.value);
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
