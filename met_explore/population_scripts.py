from met_explore.serializers import *
import numpy as np
import logging
import json

logger = logging.getLogger(__name__)

# Give the sample CSV file to populate the samples.
# KMcL: Working but need to condiser the filepath.

def populate_samples(sample_csv):

    sample_details = np.genfromtxt(sample_csv, delimiter=',', dtype=str)[2:]

    for sample in sample_details:
        sample_serializer = SampleSerializer(
            data={"name": sample[0], "group": sample[1], "life_stage": sample[2], "tissue": sample[3],
                  "mutant": sample[4]})
        if sample_serializer.is_valid():
            db_sample = sample_serializer.save()
            logger.info("sample saved ", db_sample.name)
        else:
            logger.error(sample_serializer.errors)


# This requires the input in the order taken from the construct_peak_df method/
# It requires all secondary_ids to be unique and reports any errors (throw?)
# To get the list and dict back from the json.dumps just use json.loads

def populate_peaks(peak_array):
    for peak in peak_array:
        print(peak)
        cmpd_id = json.dumps(peak[12])
        frank_annot = json.dumps(peak[13])

        peak_serializer = PeakSerializer(
            data={"psec_id": peak[1], "m_z": format(peak[2], '.9f'), "neutral_mass": format(peak[14], '.9f'),
                  "rt": peak[3], "polarity": peak[4],
                  "cmpd_name": peak[10], "cmpd_formula": peak[6], "cmpd_identifiers": cmpd_id, "identified": peak[8],
                  "frank_anno": frank_annot, "adduct": peak[7], "db": peak[11]})
        if peak_serializer.is_valid():
            db_peak = peak_serializer.save()
            logger.info("peak saved ", db_peak.psec_id)
        else:
            logger.error(peak_serializer.errors)