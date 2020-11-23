export email=$"Subject: Data Quality Shift Reminder
Hi!

This is an automated reminder that your shift monitoring data quality from the north is beginning soon. Please consult the shift schedule:
https://docs.google.com/spreadsheets/d/1RgkZXJxeFPv3qQRP_VjeauRam6S_J_7txtMtfaIBcXM/edit#gid=0

The data quality summary page is located here:
http://spt3g-dq.grid.uchicago.edu/summary.html

And notes for your shift should be recorded on the wiki:
https://pole.uchicago.edu/spt3g/index.php/Northern_Offsite_Shift"

#echo $@
echo "$email" | /usr/sbin/sendmail -v $@
