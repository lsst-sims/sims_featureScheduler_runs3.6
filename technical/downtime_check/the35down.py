
import numpy as np
from astropy.time import Time, TimeDelta
import rubin_scheduler.site_models as site_models


def tma_movement(percent=70):
    """Get a dictionary of parameters to pass to `setup_telescope`
     defining altitude and azimuth speed, acceleration, and jerk
     in terms of 'percent' of total performance.

     Parameters
     ----------
     percent : `float`, optional
        Default performance for the scheduler simulations for operations
        has been 70% (70, default).
        Expected performance at the start of comcam on-sky
        science operations is about 10%.

    Returns
    -------
    tma : `dict` {`str`: `float`}
        A dictionary which can be passed as kwargs to
        KinematicModel.setup_telescope(**tma).
    """
    # See https://confluence.lsstcorp.org/display/LSSTCOM/TMA+Motion+Settings
    # Expected performance at end of comcam on-sky is probably 10%
    if percent > 125:
        percent = 125
        print("Cannot exceed 125 percent, by requirements.")
    tma = {}
    scale = percent / 100.0
    tma["azimuth_maxspeed"] = np.min([10.0 * scale, 7.0])
    tma["azimuth_accel"] = 10.0 * scale
    tma["azimuth_jerk"] = np.max([1.0, 40.0 * scale])
    tma["altitude_maxspeed"] = 5.0 * scale
    tma["altitude_accel"] = 5.0 * scale
    tma["altitude_jerk"] = np.max([1.0, 20.0 * scale])
    tma["settle_time"] = 3.0
    return tma


def rotator_movement(percent=100):
    """Get a dictionary of parameters to pass to `setup_camera`
     defining rotator max speed, acceleration and jerk,
     in terms of 'percent' of total performance.

     Parameters
     ----------
     percent : `float`, optional
        Default performance for the scheduler simulations for operations
        has been 100% (100, default).
        Expected performance at the start of comcam on-sky
        science operations is approximately full performance.

    Returns
    -------
    rot : `dict` {`str`: `float`}
        A dictionary which can be passed as kwargs to
        KinematicModel.setup_camera(**rot).
    """
    # Kevin and Brian say these can run 100%
    # and are independent of TMA movement
    if percent > 125:
        percent = 125
        print("Cannot exceed 125 percent, by requirements.")
    rot = {}
    rot["maxspeed"] = 3.5 * percent / 100
    rot["accel"] = 1.0 * percent / 100
    rot["jerk"] = 4.0 * percent / 100
    return rot


class UnscheduledDowntimeDataYearOne:
    """Handle (and create) the unscheduled downtime information.

    Parameters
    ----------
    sunsets : `np.ndarray`, (N,)
        Sunset information (in  mjd `float` format), for at least
        the N nights for which there should be random up and down times 
        within a night.
    sunrises : `np.ndarray`, (N,)
        Sunrise information (in  mjd `float` format), for at least
        the N nights for which there should be random up and down times 
        within a night.
    seed : `int`, optional
        The random seed for creating the random nights of unscheduled
        downtime. Default 43.
    """

    MINOR_EVENT = {"P": 0.0137, "length": 1, "level": "minor event"}
    INTERMEDIATE_EVENT = {"P": 0.00548, "length": 3, "level": "intermediate event"}
    MAJOR_EVENT = {"P": 0.00137, "length": 7, "level": "major event"}
    CATASTROPHIC_EVENT = {"P": 0.000274, "length": 14, "level": "catastrophic event"}

    def __init__(
        self,
        sunsets,
        sunrises,
        seed=43,
    ):
        self.seed = seed
        self.sunsets = sunsets
        self.sunrises = sunrises
        # Scheduled downtime data is a np.ndarray of start
        # / end / activity for each scheduled downtime.
        self.downtime = None
        self.make_data()

    def __call__(self):
        """Return the array of unscheduled downtimes.

        Parameters
        ----------
        time : `astropy.time.Time`
            Time in the simulation for which to find the current downtime.

        Returns
        -------
        downtime : `np.ndarray`
            The array of all unscheduled downtimes, with keys for
            'start', 'end', 'activity',  corresponding to
            `astropy.time.Time`, `astropy.time.Time`, and `str`.
        """
        return self.downtime

    def _downtime_status(self, time):
        """Look behind the scenes at the downtime status/next values"""
        next_start = self.downtime["start"].searchsorted(time, side="right")
        next_end = self.downtime["end"].searchsorted(time, side="right")
        if next_start > next_end:
            current = self.downtime[next_end]
        else:
            current = None
        future = self.downtime[next_start:]
        return current, future

    def make_data(self):
        """Configure the set of unscheduled downtimes.

        This function creates the unscheduled downtimes based on a set
        of probabilities of the downtime type occurance.

        The random downtime is calculated using the following
        probabilities:

        minor event : remainder of night and next day = 5/365 days
        e.g. power supply failure
        intermediate : 3 nights = 2/365 days e.g. repair filter
        mechanism, rotator, hexapod, or shutter
        major event : 7 nights = 1/2*365 days
        catastrophic event : 14 nights = 1/3650 days e.g. replace a raft
        """
        self.rng = np.random.default_rng(seed=self.seed)

        starts = []
        ends = []
        acts = []

        end_of_start = 380
        
        night_counted = np.zeros(len(self.sunsets))
        
        for night, (sunset, sunrise) in enumerate(zip(self.sunsets, self.sunrises)):
            prob = self.rng.random()
            hours_in_night = (sunrise - sunset) * 24.0
            if night_counted[night] == 1:
                continue
                
            if night < end_of_start:
                # Estimate a threshold probability of having some downtime - 
                # 50% at start, dropping until end_of_start, where it should be .. 5%?
                nightly_threshold = 0.5 * (1 - night / (end_of_start + 45))
                if prob <= nightly_threshold:
                    # Generate an estimate of how long the downtime should be
                    #prob_time = self.rng.uniform(low=1, high=hours_in_night, size=1)[0] 
                    prob_time = self.rng.gumbel(loc=1, scale=6, size=1)[0]
                    if prob_time >= hours_in_night:
                        prob_time = hours_in_night
                    if prob_time <= 1:
                        prob_time = 1.0
                    # And generate a starting time during the night for this event
                    tmax = hours_in_night - prob_time
                    if tmax <= 0:
                        starts.append(Time(sunset, format='mjd', scale='utc'))
                        ends.append(Time(sunrise, format='mjd', scale='utc'))
                        acts.append("Year1 Eng")
                    else:
                        offset = self.rng.uniform(low=sunset, high=sunset + tmax / 24.0)
                        starts.append(Time(offset, format='mjd', scale='utc'))
                        ends.append(Time(offset + prob_time/24.0, format='mjd', scale='utc'))
                        acts.append("Year1 Eng")
                night_counted[night] = 1
                continue
            # And also add the standard unscheduled downtime  
            start_time = Time(sunset, format='mjd', scale='utc')
            if prob < self.CATASTROPHIC_EVENT["P"]:
                starts.append(start_time)
                end_night = start_time + TimeDelta(self.CATASTROPHIC_EVENT["length"], format="jd")
                ends.append(end_night)
                acts.append(self.CATASTROPHIC_EVENT["level"])
                night_counted[night:night + self.CATASTROPHIC_EVENT["length"]] = 1
            elif prob < self.MAJOR_EVENT["P"]:                    
                starts.append(start_time)
                end_night = start_time + TimeDelta(self.MAJOR_EVENT["length"], format="jd")
                ends.append(end_night)
                acts.append(self.MAJOR_EVENT["level"])
                night_counted[night:night + self.MAJOR_EVENT["length"]] = 1
            elif prob < self.INTERMEDIATE_EVENT["P"]:
                starts.append(start_time)
                end_night = start_time + TimeDelta(self.INTERMEDIATE_EVENT["length"], format="jd")
                ends.append(end_night)
                acts.append(self.INTERMEDIATE_EVENT["level"])
                night_counted[night:night + self.INTERMEDIATE_EVENT["length"]] = 1
            elif prob < self.MINOR_EVENT["P"]:
                starts.append(start_time)
                end_night = start_time + TimeDelta(self.MINOR_EVENT["length"], format="jd")
                ends.append(end_night)
                acts.append(self.MINOR_EVENT["level"])
                night_counted[night:night + self.MINOR_EVENT["length"]] = 1
        
        self.downtime = np.array(
            list(zip(starts, ends, acts)),
            dtype=[("start", "O"), ("end", "O"), ("activity", "O")],
        )

    def total_downtime(self):
        """Return total downtime (in days).

        Returns
        -------
        total : `int`
            Total number of downtime days.
        """
        total = 0
        for td in self.downtime["end"] - self.downtime["start"]:
            total += td.jd
        return total


def new_downtimes(mjd_start=None):
    """return the array of new downtimes
    """
    almanac = site_models.Almanac(mjd_start=mjd_start)
    year1_sunsets = np.where((almanac.sunsets['night'] >= 0) & (almanac.sunsets['night'] < 366))
    sunset = almanac.sunsets[year1_sunsets]['sun_n12_setting']
    sunrise = almanac.sunsets[year1_sunsets]['sun_n12_rising']

    mjd_start_time = Time(mjd_start, format="mjd")
    sched_downtime_data = site_models.ScheduledDowntimeData(mjd_start_time)
    unsched_downtime = UnscheduledDowntimeDataYearOne(sunsets=sunset, sunrises=sunrise)
    unscheduled_downtimes = unsched_downtime()

    mjd_start_time = Time(mjd_start+365, format="mjd")
    reg_dt = site_models.UnscheduledDowntimeData(mjd_start_time)
    regular_downtimes = reg_dt()

    down_starts = []
    down_ends = []

    for dt in unscheduled_downtimes:
        down_starts.append(dt["start"].mjd)
        down_ends.append(dt["end"].mjd)

    for dt in regular_downtimes:
        down_starts.append(dt["start"].mjd)
        down_ends.append(dt["end"].mjd)

    for dt in sched_downtime_data():
        down_starts.append(dt["start"].mjd)
        down_ends.append(dt["end"].mjd)

    downtimes = np.array(
        list(zip(down_starts, down_ends)),
        dtype=list(zip(["start", "end"], [float, float])),
    )
    downtimes.sort(order="start")

    # Make sure there aren't any overlapping downtimes
    diff = downtimes["start"][1:] - downtimes["end"][0:-1]
    while np.min(diff) < 0:
        # Should be able to do this without a loop, but this works
        for i, dt in enumerate(downtimes[0:-1]):
            if downtimes["start"][i + 1] < dt["end"]:
                new_end = np.max([dt["end"], downtimes["end"][i + 1]])
                downtimes[i]["end"] = new_end
                downtimes[i + 1]["end"] = new_end

        good = np.where(downtimes["end"] - np.roll(downtimes["end"], 1) != 0)
        downtimes = downtimes[good]
        diff = downtimes["start"][1:] - downtimes["end"][0:-1]

    return downtimes

