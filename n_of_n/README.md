Checking if using neighbors-of-neighbors for blob growth works better
on branch u/yoachim/larger_contig
which is a commit ahead of tickets/OPSIM-1180



branch u/yoachim/larger_contig:
Flushed 13996 observations from queue for being stale
Completed 182147 observations
ran in 20 min = 0.3 hours
Writing results to  baseline_v3.6_1yrs.db

branch tickets/OPSIM-1180:
Flushed 15547 observations from queue for being stale
Completed 183367 observations
ran in 19 min = 0.3 hours
Writing results to  baseline_v3.6_1yrs.db

branch u/yoachim/blob_order:
Flushed 14100 observations from queue for being stale
Completed 183348 observations
ran in 19 min = 0.3 hours


------

OK, so looks like the alt,az thing isn't really needed. 


u/yoachim/blob_order:
Flushed 97929 observations from queue for being stale
Completed 2021920 observations
ran in 233 min = 3.9 hours
Writing results to  baseline_v3.6_10yrs.db


branch tickets/OPSIM-1180:
Flushed 99846 observations from queue for being stale
Completed 2022005 observations
ran in 220 min = 3.7 hours
Writing results to  baseline_v3.6_10yrs.db

branch u/yoachim/larger_contig:
Flushed 101689 observations from queue for being stale
Completed 2004304 observations
ran in 349 min = 5.8 hours
Writing results to  baseline_v3.6_10yrs.db


So, it looks like trying to increase the contiguous area results in a 
1% decrease in survey efficiency, and has a bit of a performance hit as well.
