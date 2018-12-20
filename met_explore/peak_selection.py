from IPython.display import display, HTML
from met_explore.models import *
from met_explore.serializers import *
from difflib import SequenceMatcher
from collections import OrderedDict

import pandas as pd
import numpy as np
import logging
import os
import json
import collections

logger = logging.getLogger(__name__)

PROTON = 1.007276
NAME_MATCH_SIG = 0.5

# A class to select the peaks required to use in further analysis

class PeakSelector(object):

    # The constructor just takes in a peak_json file from PiMP
    def __init__(self, peak_json_file):

        peak_details_df = pd.read_json(peak_json_file)

        # Filter on adduct types
        selected_adducts = (peak_details_df['adduct'] == 'M+H') | (peak_details_df['adduct'] == 'M-H')
        f_adducts = peak_details_df[selected_adducts].copy()

        # Select peaks that have been identified and/or has a FrAnk annotation associated with them
        with_annot = (f_adducts['frank_annot'].notnull()) | (f_adducts['identified'] == 'True')
        with_annots_df = f_adducts[with_annot].copy()

        # This is the df to be used to select peaks/compounds
        self.selected_df = self.add_neutral_masses(with_annots_df)
        self.unique_sec_ids = self.selected_df['sec_id'].unique()

        #This is the new df to store the chosen peaks.
        headers = list(self.selected_df.columns.values)
        self.final_df = pd.DataFrame(columns=headers)

        print ("in the beginning the final ID DF is ")
        display(self.final_df)

        logger.info("PeakSelector initialised with ", len(self.unique_sec_ids),"peaks/unique ids")


    def construct_peak_df(self):

        print ("Constructing the peak DF")
        unique_sec_ids = self.selected_df['sec_id'].unique()
        # unique_sec_ids = [3722]

        for sid in unique_sec_ids:
            # Collect a single sec_id into a DF
            sid_df = self.selected_df[self.selected_df.sec_id == sid]
            print("The single SID DF is")
            display(sid_df)
            # If the peak has an identified compound then keep that
            identified_df = sid_df[sid_df.identified == 'True']
            print("The identified df is: ")
            display(identified_df)
            new_row = None
            # If some of the rows have compounds that have identified=True
            if not identified_df.empty:

                # Check if there are more than one standard compounds for this sid
                standard_cmpds = sid_df[sid_df.db == 'stds_db']
                num_std_cmpds = standard_cmpds.shape[0]

                # If there is only one standard compound add this to the final DF and collect identifiers.
                if (num_std_cmpds == 1):
                    print("we have only one standard compound")
                    cmpd_id = sid_df[sid_df.db == "stds_db"]['cmpd_id'].values[0]
                    new_row = self.get_peak_by_cmpd_id(sid_df, cmpd_id)

                    # Here we have the senario that more that 1 standard compound has been identified and we
                # want to select a standard compound if possible
                if (num_std_cmpds > 1):
                    print("the number of standard compounds for sid is", sid, "is", num_std_cmpds)
                    new_row = self.select_standard_cmpd(sid_df, standard_cmpds)

                # If a new_row has been returned for this SID - add it to the final_df
                if new_row is not None:

                    print("we are adding the row for sid", sid)
                    display(pd.DataFrame(new_row).T)
                    self.final_df = self.final_df.append(new_row)

                    # If the new_row has not been determined for this SID
                else:

                    unique_cmpd_ids = sid_df['cmpd_id'].unique()

                    # For each unique compound id add a row to the final df, this will produce duplicates for later
                    for ucid in unique_cmpd_ids:
                        new_row = self.get_peak_by_cmpd_id(sid_df, ucid)
                        print("we are adding the row: for sid", sid)
                        display(pd.DataFrame(new_row).T)
                        self.final_df = self.final_df.append(new_row)

                        # Else nothing identified so look at the fragmentation data.
            else:
                # Get all the rows for this secondary ID
                print("nothing identified here so get best match FrAnk compound")
                new_row = self.select_on_frank(sid_df)
                print("we are adding the row: for sid", sid)
                display(pd.DataFrame(new_row).T)
                self.final_df = self.final_df.append(new_row)



                # final_df = check_for_duplicates(final_df, chosen_std_cmpds)

        print("This is the final df")
        display(self.final_df)

        print("There are", self.final_df['sec_id'].nunique(), "unique compounds out of", self.final_df.shape[0], "rows added")

        return self.final_df


    def select_standard_cmpd(self, sid_df, standard_cmpds):
        """
        Choose a standard compound based on the one which matches the name of the FrAnk annotation most closely.
        :param sid_df: a dataframe of the rows of the peaks for a single peak ID
        :param standard_cmpds: A dataframe for the peak containing only the standard compounds
        :return: A new row for the final df based on FrAnk or None (if FrAnk None)
        """

        print ("selecting standard compound")
        new_row = None
        display(standard_cmpds)
        name_match_dic = {}
        # For each of the standard compounds identified for the peak
        for i in standard_cmpds['identifier'].unique():
            pimp_cmpd_name = standard_cmpds[standard_cmpds.identifier == i]['compound'].iloc[0]
            annotation = standard_cmpds['frank_annot'].values[0]  # Get the value in the cell

            # If there is a FrAnK annotation get the best name match
            if (annotation is not None):

                frank_cmpd_name = annotation['frank_cmpd_name']
                m = SequenceMatcher(None, frank_cmpd_name, pimp_cmpd_name)
                name_match_dic[pimp_cmpd_name] = m.ratio()

                print(name_match_dic)
                max_value = max(name_match_dic.values())  # maximum value
                max_keys = [k for k, v in name_match_dic.items() if v == max_value]
                max_key = max_keys[0]
                print("max_key to grab row", max_key)

                new_row = standard_cmpds[standard_cmpds['compound'] == max_key].iloc[0]

                ucid = new_row['cmpd_id']
                cmpd_rows_df = sid_df[sid_df.cmpd_id == ucid]
                identifiers = self.get_all_identifiers(cmpd_rows_df)
                new_row.at['identifier'] = identifiers

        print("sending back a new row from select_standard_compound of type ", type(new_row))
        return new_row

    # Choose a compound bases on how closely it matches the name of the FrAnK annotation. If it is less than 50% return
    # the FrAnk details instead

    def select_on_frank(self, sid_df):
        """
        :param sid_df: a dataframe of the rows of the peaks for a single peak ID
        :return: A new row for the final df based on FrAnk (should not be None)
        """
        new_row = None
        name_match_dic = {}
        compound_names = sid_df['compound'].values
        single_annot = sid_df['frank_annot'].values[0]

        for pimp_cmpd_name in compound_names:
            # Find the best fit for to frank.
            frank_cmpd_name = single_annot['frank_cmpd_name']
            m = SequenceMatcher(None, frank_cmpd_name, pimp_cmpd_name)
            name_match_dic[pimp_cmpd_name] = m.ratio()

        if name_match_dic:
            print(name_match_dic)
            max_value = max(name_match_dic.values())  # maximum value
            max_keys = [k for k, v in name_match_dic.items() if v == max_value]
            max_key = max_keys[0]

            if (max_value >= 0.5):
                print("max_key to grab row", max_key)

                new_row = sid_df[sid_df['compound'] == max_key].iloc[0]
                ucid = new_row['cmpd_id']
                cmpd_rows_df = sid_df[sid_df.cmpd_id == ucid]
                identifiers = self.get_all_identifiers(cmpd_rows_df)
                new_row.at['identifier'] = identifiers

            # If the match is less than 50% just take the frank annotation instead.
            if (max_value < 0.5):
                print("max_value is < 0.5 (", max_value, ") and so taking frank annot")

                new_row = self.get_frank_annot(sid_df)

        print("The FrAnK probability score is", single_annot['probability'])

        return new_row

    def get_frank_annot(self, sid_df):
        '''
        :param sid_df: Dataframe containing all rows for a single peak (sid, secondary id)
        :return: A new row of the dataframe based on the FrAnk compound (instead of the PiMP one)
        '''
        new_row = None
        frank_annots = sid_df['frank_annot']
        single_annot = frank_annots.iloc[0]

        if any(single_annot):

            print("Collecting the compound info from", single_annot)

            new_row = sid_df.iloc[0]
            new_row.at['compound'] = single_annot['frank_cmpd_name']

            # Get all the frAnk identifiers for a single PiMP compound (could be one).
            identifiers = []
            identifier_keys = ['inchikey', 'cas_code', 'hmdb_id']
            for i in identifier_keys:
                identifiers.append(single_annot[i])
            new_row.at['identifier'] = identifiers

        return new_row


    def get_peak_by_cmpd_id(self, sid_df, ucid):
        """
        A method to return a row for an identified peak given the chosen compound id
        :param: A df containing a set of peaks with a single ID and the unique cmpd ID
        :returns: A new_row (peak) to be added to the final df.
        """
        # new_row = None  # Clear the new row at this stage
        cmpd_rows_df = sid_df[sid_df.cmpd_id == ucid]
        identifiers = self.get_all_identifiers(cmpd_rows_df)
        print("The returned identifiers are ", identifiers)
        cmpd_id = cmpd_rows_df['cmpd_id'] == ucid

        # Take the row with std_db just to remember this compound was identified.
        db = cmpd_rows_df['db'] == 'stds_db'
        new_row = cmpd_rows_df[cmpd_id & db].iloc[0]
        new_row.at['identifier'] = identifiers

        print('Returning new_row by cmpd_id')

        return new_row


    def get_all_identifiers(self, cmpd_rows_df):

        """ A method to return a list of identifiers for an identified
            Params: A df containing a unique 'pimp' compound
            returns: A list of identifiers relating to the compound
        """

        num_rows = cmpd_rows_df.shape[0]
        # Get all the identifiers for a single PiMP compound (could be one but want this as list).
        identifiers = []
        for i in range(0, (num_rows)):
            new_id = cmpd_rows_df.iloc[i]['identifier']
            identifiers.append(new_id)

            # Take one of the rows for this compound, add identifiers and save in the final df.

        return identifiers

    def add_neutral_masses(self, df):

        """ A method to add neutral masses to a DF
            Params: The df containing M+H or M-H adducts (no other adducts)
            Returns: A dataframe with the neutral masses added.
        """
        logger.info("Adding neutral masses to the peaks")

        masses = df['mass'].values
        adducts = df['adduct'].values

        neutral_masses = []
        joint_list = [masses, adducts]

        mass_adducts = list(zip(*joint_list))
        for ma in mass_adducts:
            mass = ma[0]
            adduct = ma[1]
            neutral_mass = self.get_neutral_mass(mass, adduct)
            neutral_masses.append(neutral_mass)
        print(type(neutral_masses))

        df['neutral_mass'] = np.asarray(neutral_masses)

        return df

    def get_neutral_mass(self, mass, adduct):

        """
        Small function to return a neutral mass given the m/z and the adduct.
        :param mass: m/z of the adduct
        :param adduct: type of adduct - only M+H, M-H currently accepted
        :return: The neutral mass np.float
        """
        if adduct == 'M+H':
            neutral_mass = mass - PROTON
        elif adduct == 'M-H':
            neutral_mass = mass + PROTON
        else:
            logger.warning("This is not the correct type of adduct and therefor skipping")
            return

        return neutral_mass