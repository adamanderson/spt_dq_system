$( function() {
  var dFormat = "yy-mm-dd",
    from = $( "#date-from" )
      .datepicker({
        defaultDate: "+1w",
        changeMonth: true,
        numberOfMonths: 1,
        dateFormat: dFormat
      })
      .on( "change", function() {
        to.datepicker( "option", "minDate", getDate( this ) );
      }),
    to = $( "#date-to" ).datepicker({
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

function date_chosen() {
  var tag = document.getElementsByClassName('ui-datepicker-month');
  for (var i = 0; i < tag.length; i++) {
    tag[i].classList.add("chosen-select");
    $(".chosen-select").chosen({disable_search: true});
  }
}
