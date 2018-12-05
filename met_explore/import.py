from met_explore.serializers import *
import numpy as np
import logging

logger = logging.getLogger(__name__)


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

