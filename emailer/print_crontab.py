import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import os

df_schedule = pd.read_csv('northern_shift_schedule.tsv', sep='\t')
df_emails = pd.read_csv('northern_shift_emails.tsv', sep='\t')

dict_emails = {}
for ind, row in df_emails.iterrows():
    dict_emails[row['Shifter List']] = row['Email']

for ind, row in df_schedule.iterrows():
    time_start = datetime(int(row['Start Date'].split('/')[2]),
                          int(row['Start Date'].split('/')[0]),
                          int(row['Start Date'].split('/')[1]))
    dt = timedelta(days=1)
    time_remind = time_start - dt
    cron_string = '0 12 {} {} *'.format(time_remind.day, time_remind.month)

    script = os.path.join(os.getcwd(), 'send_shift_reminder.sh')
    cron_string += ' {} {} {}'.format(script,
                                      dict_emails[row['Shifter1']],
                                      dict_emails[row['Shifter2']])
    print(cron_string)
