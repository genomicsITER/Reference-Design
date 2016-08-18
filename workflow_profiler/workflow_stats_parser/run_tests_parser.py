#!/usr/bin/env python

import unittest
from subprocess import call
from contextlib import contextmanager
import os
import shutil
from glob import glob
from hashlib import md5

###########################################################################
#Classes and functions that are used to test success and failures!

#Creating a folder for storing the automation test script outputs
def make_temp_output_folder(folder_name="test_output"):
    """
    PURPOSE: Creates a temporary output folder. If folder exists, delete it and create it.
    """
    if os.path.exists(folder_name):
        shutil.rmtree(folder_name)
        os.makedirs(folder_name)
    else: os.makedirs(folder_name)
    return folder_name

# Verifying the plots exist
def verify_plots_exist(folder):
    csv_files = glob(folder + '/*.csv')
    png_files = glob(folder + '/*.png')

    # Make sure there are some files
    if any((not csv_files, not png_files)):
        return False

    for output_file in csv_files + png_files:
        # Check filesize
        if os.path.getsize(output_file) <= 0:
            return False

    return True

def verify_plots_md5(input_folder, sample_folder):
    png_files = glob(input_folder + '/*.png')

    for output_file in png_files:
        with open(output_file, 'rb') as file_handle:
            output_md5 = md5(file_handle.read()).digest()

        bare_filename = os.path.basename(output_file)

        with open(os.path.join(sample_folder, bare_filename), 'rb') as file_handle:
            sample_md5 = md5(file_handle.read()).digest()

        if sample_md5 != output_md5:
            print("md5sum didn't match between the data generated and the sample output. Best to match the outputs manually\n")
            return False

    return True

###########################################################################
# Test for failures - Validating arguments mainly.
class TestFailures(unittest.TestCase):
    def test_missing_metric_argument(self):
        print("\nMissing Metric Argument.")
        out_folder = make_temp_output_folder()
        self.assertFalse(call("./workflow_stats_parser.py sample_multistage_input/ -o %s -N sample -p" % out_folder, shell=True) == 0)

    def test_input_folder_doesnt_exist(self):
        print("\nInput Folder doesn't exist.")
        out_folder = make_temp_output_folder()
        self.assertFalse(call("./workflow_stats_parser.py not_sample_multistage_input/ -o %s -N sample -isp" % out_folder, shell=True) == 0)

    def test_invalid_pipeline_specified(self):
        print("\nWrong pipeline name is provided.")
        out_folder = make_temp_output_folder()
        self.assertFalse(call("./workflow_stats_parser.py sample_multistage_input/ -o %s -N not_a_pipeline -isp" % out_folder, shell=True) == 0)

    def test_single_step_doesnt_exist(self):
        print("\nSingle Step doesn't exist.")
        out_folder = make_temp_output_folder()
        self.assertFalse(call("./workflow_stats_parser.py sample_multistage_input/ -o %s -N sample -S not_a_step_name -isp" % out_folder, shell=True) == 0)


###########################################################################
# Test for success - Match output data to sample output folders.
class TestRunsFinish(unittest.TestCase):
    def test_sample_multistage_with_all_metrics(self):
        out_folder = make_temp_output_folder("_test_multi_stage")
        cmd = "./workflow_stats_parser.py sample_multistage_input/ -N sample -isp -o '%s'" % out_folder
        print(("\nRunning: %s" % cmd))
        self.assertTrue(call(cmd, shell=True) == 0)
        self.assertTrue(verify_plots_exist(out_folder))
        self.assertTrue(verify_plots_md5(out_folder, "./sample_multistage_output"))

    def test_sample_single_stage_with_all_metrics(self):
        out_folder = make_temp_output_folder("_test_single_stage")
        cmd = "./workflow_stats_parser.py sample_onestage_input/ -N sample -S stage1 -isp -o '%s'" % out_folder
        print(("\nRunning: %s" % cmd))
        self.assertTrue(call(cmd, shell=True) == 0)
        self.assertTrue(verify_plots_exist(out_folder))
        self.assertTrue(verify_plots_md5(out_folder, "./sample_onestage_output"))


###########################################################################
def run_tests():
    unittest.main()

if __name__ == "__main__":
    run_tests()
