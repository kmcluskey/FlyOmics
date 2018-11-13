# Generated by Django 2.1.3 on 2018-11-12 11:52

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Peak',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('m_z', models.DecimalField(decimal_places=10, max_digits=20)),
                ('rt', models.DecimalField(decimal_places=10, max_digits=20)),
                ('polarity', models.CharField(max_length=8)),
                ('anno_type', models.CharField(max_length=100)),
                ('cmpd_name', models.CharField(max_length=600)),
                ('cmpd_identifiers', models.CharField(max_length=600)),
            ],
        ),
        migrations.CreateModel(
            name='Sample',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=250)),
                ('life_stage', models.CharField(max_length=250)),
                ('tissue', models.CharField(max_length=250)),
                ('mutant', models.CharField(blank=True, max_length=250, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='SamplePeak',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('intensity', models.FloatField(blank=True, null=True)),
                ('peak', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='met_explore.Peak')),
                ('sample', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='met_explore.Sample')),
            ],
        ),
    ]