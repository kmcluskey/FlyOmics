from django.db import models

class Sample(models.Model):
    """
    Model class defining an instance of an experimental Sample including the tissue and life-stage from which it came
    """

    name = models.CharField(max_length=250, blank=False)
    life_stage = models.CharField(max_length=250, blank=False)
    tissue = models.CharField(max_length=250)
    mutant = models.CharField(max_length=250, blank=True, null=True)

    def __unicode__(self):
        """
        Method to return a representation of the Sample
        """

        return "Sample" + self.sample.name


class Peak(models.Model):
    """
    Model class representing a basic peak including the compound as a simple string.
    """

    m_z = models.DecimalField(max_digits=20, decimal_places=10)
    rt = models.DecimalField(max_digits=20, decimal_places=10)
    polarity = models.CharField(max_length=8)
    anno_type = models.CharField(max_length=100)
    cmpd_name = models.CharField(max_length=600) # At this stage just a name for the metabolite
    cmpd_identifiers = models.CharField(max_length=600) # Any idenifiers we can associate with the peak

    def __unicode__(self):
        """
        Method to return a representation of the SamplePeak including the name of the compound
        :return: String:
        """

        return "Peak" + self.peak.id + "peak_compound " + self.peak.cmpd_name


class SamplePeak(models.Model):

    peak = models.ForeignKey(Peak, on_delete=models.CASCADE)
    sample = models.ForeignKey(Sample, on_delete=models.CASCADE)
    intensity = models.FloatField(null=True, blank=True)

    def __unicode__(self):

        """
        Method to return a representation of the SamplePeak including the name of the compound
        :return: String:
        """

        return "Sample" + self.sample.name + "peak" + self.peak.id + "peak_compound "+ self.peak.cmpd_name








