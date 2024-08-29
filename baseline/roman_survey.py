import numpy as np
import rubin_scheduler.scheduler.basis_functions as bf
import rubin_scheduler.scheduler.detailers as detailers
from astropy.utils import iers
from rubin_scheduler.scheduler.surveys import (
    DeepDrillingSurvey,
)
from rubin_scheduler.utils import ra_dec2_hpid, _angular_separation

# So things don't fail on hyak
iers.conf.auto_download = False
# XXX--note this line probably shouldn't be in production
iers.conf.auto_max_age = None


class InSeasonBasisFunction(bf.BaseBasisFunction):
    """Only let a survey go if it is in a defined season"""

    def __init__(self, seasons=[]):
        super().__init__()
        self.seasons = seasons

    def check_feasibility(self, conditions):
        result = False
        for season_range in self.seasons:
            if np.min(season_range) <= conditions.mjd <= np.max(season_range):
                result = True
        return result


class AirmassPointRange(bf.BaseBasisFunction):
    """set an airmass limit for a single point"""

    def __init__(self, ra, dec, airmass_range=[1.05, 2.7], nside=32):
        super().__init__()
        self.hpid = ra_dec2_hpid(nside, ra, dec)
        self.airmass_range = airmass_range

    def check_feasibility(self, conditions):
        result = False
        airmass = conditions.airmass[self.hpid]
        if (np.min(self.airmass_range) <= airmass) & (
            airmass <= np.max(self.airmass_range)
        ):
            result = True
        return result


class MoonDistPointRange(bf.BaseBasisFunction):
    """set an airmass limit for a single point"""

    def __init__(self, ra, dec, moon_limit=15.0):
        super().__init__()
        self.ra = np.radians(ra)
        self.dec = np.radians(dec)
        self.moon_limit = np.radians(moon_limit)

    def check_feasibility(self, conditions):
        result = False
        moon_dist = _angular_separation(
            self.ra, self.dec, conditions.moon_ra, conditions.moon_dec
        )
        if moon_dist > self.moon_limit:
            result = True
        return result


def roman_info():
    """Manually enter some Roman RGES info."""
    # From the TVS Slack channel:
    # spring 2027, fall 2027, spring 2028, then fall 2030, spring 2031, fall 2031

    result = {}
    result["RA"] = 268.708
    result["dec"] = -28.975

    # Guessing these from the notebook in same dir.
    observing_season_mid_mjds = [61947.3, 62318.8, 62670.3, 63067.2, 63381.4, 63773.2]

    result["seasons_on"] = [[val - 32, val + 32] for val in observing_season_mid_mjds]

    result["seasons_off"] = []
    for i in range(len(result["seasons_on"]) - 1):
        result["seasons_off"].append(
            [result["seasons_on"][i][1] + 1, result["seasons_on"][i + 1][0] - 1]
        )

    return result


def gen_roman_on_season(
    nside=32,
    camera_ddf_rot_limit=75.0,
    exptime=30.0,
    nexp=2,
):

    field_info = roman_info()

    RA = field_info["RA"]
    dec = field_info["dec"]

    survey_name = "DD: RGES_onseason"

    # Add some feasability basis functions. Maybe just give it a set of nights where it can execute for now.
    basis_functions = []
    # These are crude hard limits. Nominally we would try to pre-schedule these
    # when they would be at the best airamss in the night.
    basis_functions.append(
        bf.HourAngleLimitBasisFunction(RA=RA, ha_limits=[[20, 24], [0, 4]])
    )
    basis_functions.append(bf.NotTwilightBasisFunction())
    # Force it to delay 30 minutes
    basis_functions.append(
        bf.ForceDelayBasisFunction(days_delay=30.0 / 24.0, survey_name=survey_name)
    )
    # Force it to be in a given observing season
    basis_functions.append(InSeasonBasisFunction(seasons=field_info["seasons_on"]))
    basis_functions.append(MoonDistPointRange(RA, dec))
    basis_functions.append(AirmassPointRange(RA, dec, nside=nside))

    # Add a dither detailer, so it dithers between each set of exposures I guess?
    details = []
    details.append(detailers.DitherDetailer(max_dither=0.5, seed=42, per_night=True))
    details.append(
        detailers.CameraRotDetailer(
            min_rot=-camera_ddf_rot_limit, max_rot=camera_ddf_rot_limit
        )
    )

    survey = DeepDrillingSurvey(
        basis_functions,
        RA=RA,
        dec=dec,
        sequence="giriz",
        nvis=[1, 1, 1, 1, 1],
        exptime=exptime,
        nexp=nexp,
        survey_name=survey_name,
        detailers=details,
    )
    return survey


def gen_roman_off_season(
    nside=32,
    camera_ddf_rot_limit=75.0,
    exptime=30.0,
    nexp=2,
):
    """Generate a ddf-like survey object to observe the roman field every ~3 days in the off-season"""

    field_info = roman_info()
    RA = field_info["RA"]
    dec = field_info["dec"]

    survey_name = "DD: RGES_offseason"

    # Add some feasability basis functions. Maybe just give it a set of nights where it can execute for now.
    basis_functions = []
    # These are crude hard limits. Nominally we would try to pre-schedule these
    # when they would be at the best airamss in the night.
    basis_functions.append(
        bf.HourAngleLimitBasisFunction(RA=RA, ha_limits=[[20, 24], [0, 4]])
    )
    basis_functions.append(bf.NotTwilightBasisFunction())
    # Force it to not go every day
    basis_functions.append(
        bf.ForceDelayBasisFunction(days_delay=3.0, survey_name=survey_name)
    )
    # Force it to be in a given observing season
    basis_functions.append(InSeasonBasisFunction(seasons=field_info["seasons_off"]))
    basis_functions.append(MoonDistPointRange(RA, dec))
    basis_functions.append(AirmassPointRange(RA, dec, nside=nside))

    # Add a dither detailer, so it dithers between each set of exposures I guess?
    details = []
    details.append(detailers.DitherDetailer(max_dither=0.5, seed=42, per_night=True))
    details.append(
        detailers.CameraRotDetailer(
            min_rot=-camera_ddf_rot_limit, max_rot=camera_ddf_rot_limit
        )
    )

    survey = DeepDrillingSurvey(
        basis_functions,
        RA=RA,
        dec=dec,
        sequence="griz",
        nvis=[1, 1, 1, 1],
        exptime=exptime,
        nexp=nexp,
        survey_name=survey_name,
        detailers=details,
    )
    return survey
