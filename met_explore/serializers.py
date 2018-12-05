from rest_framework import serializers
from met_explore.models import Sample, Peak, SamplePeak


class SampleSerializer(serializers.ModelSerializer):
    class Meta:
        model = Sample
        fields = ('name','life_stage', 'group','tissue','mutant')


class PeakSerializer(serializers.ModelSerializer):
    class Meta:
        model = Peak
        fields = ('pid','m_z','rt','polarity','cmpd_name', 'cmpd_formula','cmpd_identifiers','identified','frank_anno','adduct')


class SamplePeakSerializer(serializers.ModelSerializer):
    class Meta:
        model = SamplePeak
        fields = ('peak', 'sample','intensity')

