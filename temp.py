from package.utility import getGMST
from datetime import datetime, timedelta

t = datetime(year=2020, month=4, day=15, hour=18, minute=22, second=4)+ timedelta(seconds=3900)
print(getGMST(t))
