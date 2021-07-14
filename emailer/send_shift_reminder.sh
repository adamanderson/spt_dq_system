export email=$"Subject: Data Quality Shift Reminder
Hi!

This is an automated reminder that your data quality monitoring shift is beginning soon. Please read through the updated instructions on the wiki here:
https://pole.uchicago.edu/spt3g/index.php/Northern_Offsite_Shift

In addition, please read through the 'Read Me' tab of the data quality summary page:
http://spt3g-dq.grid.uchicago.edu/summary.html

Shifters who are new to SPT-3G data are especially encouraged to contact Wei Quan (weiquan@uchicago.edu) or Adam Anderson (adama@fnal.gov) if they are unfamiliar with the four types of calibration observations that we take, the calibration quantities derived from them, the structure of our observing cadence, or the typical noise levels of individual observations and full-season coadds. Familiarity with these details is important for evaluating the data quality, and we (and others) are happy to help you learn.

Other useful links include:
Shift schedule: https://docs.google.com/spreadsheets/d/1RgkZXJxeFPv3qQRP_VjeauRam6S_J_7txtMtfaIBcXM/edit#gid=0
Data quality summary page: http://spt3g-dq.grid.uchicago.edu/summary.html

Thanks!"

#echo $@
echo "$email" | /usr/sbin/sendmail -v $@
