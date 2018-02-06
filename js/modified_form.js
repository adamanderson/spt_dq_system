$( function() {
  var dFormat = "yy-mm-dd",
    from = $( "#modified-from" )
      .datepicker({
        defaultDate: "+1w",
        changeMonth: true,
        numberOfMonths: 1,
        dateFormat: dFormat
      })
      .on( "change", function() {
        to.datepicker( "option", "minDate", getDate( this ) );
      }),
    to = $( "#modified-to" ).datepicker({
      defaultDate: "+1w",
      changeMonth: true,
      numberOfMonths: 1,
      dateFormat: dFormat
    })
    .on( "change", function() {
      from.datepicker( "option", "maxDate", getDate( this ) );
    });

  function getDate( element ) {
    var date;
    try {
      date = $.datepicker.parseDate( dFormat, element.value );
    } catch( error ) {
      date = null;
    }

    return date;
  }
});

/* set the select menu for the date picker month to use chosen
 * and to be chosen when it is changed as well. */
function date_chosen() {
  var tag = document.getElementsByClassName('ui-datepicker-month');
  for (var i = 0; i < tag.length; i++) {
    tag[i].classList.add("chosen-select");
    $(".chosen-select").chosen({disable_search: true});
    tag[i].setAttribute("onchange", "date_chosen();");
  }
  var arrows = document.querySelectorAll('.ui-datepicker-next, .ui-datepicker-prev');
  for (var i = 0; i < arrows.length; i++) {
    arrows[i].setAttribute("onclick", "date_chosen();");
  }
}
