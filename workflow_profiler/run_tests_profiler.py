#!/usr/bin/env python

from glob import glob
import unittest
from subprocess import call
from contextlib import contextmanager
#from multiprocessing.connection import Client, Listener
import os
import re
import shutil
import threading
import signal
import sys
from hashlib import md5

##############################################################################
#Classes and functions that are used to test success and failures!

#Creating a folder for storing the automation test script outputs
def make_temp_output_folder():
    """
    PURPOSE: Creates a temporary output folder. If folder exists, delete it and create it.
    """
    folder_name = "./__test_folder"
    shutil.rmtree(folder_name, ignore_errors=True)
    os.makedirs(folder_name)

    return folder_name

# Verifying the workflow data against a baseline
def verify_workflow_data(folder):
    """
    PURPOSE: Verifies workflow output of test run by comparing MD5 hash sums of
             the *.bam and *.sam files to those of the baseline.
    INPUTS:  folder - path of test output folder
    OUTPUTS: True if all hashes compare; False if not
    """
    baseline_folder = "baseline_sim1M_pairs"
    bam_files = glob(baseline_folder + '/*.bam')
    sam_files = glob(baseline_folder + '/*.sam')
    bai_files = glob(baseline_folder + '/*.bai')

    #BAM files
    for output_file in bam_files:
        with open(output_file, 'rb') as file_handle:
            output_md5 = md5(file_handle.read()).digest()

        bare_filename = os.path.basename(output_file)

        with open(os.path.join(baseline_folder, bare_filename), 'rb') as file_handle:
            sample_md5 = md5(file_handle.read()).digest()

        if sample_md5 != output_md5:
            print("md5sum of BAM files didn't match between the data generated and the baseline. Best to match the outputs manually\n")
            return False

    #SAM files
    for output_file in sam_files:
        with open(output_file, 'rb') as file_handle:
            output_md5 = md5(file_handle.read()).digest()

        bare_filename = os.path.basename(output_file)

        with open(os.path.join(baseline_folder, bare_filename), 'rb') as file_handle:
            sample_md5 = md5(file_handle.read()).digest()

        if sample_md5 != output_md5:
            print("md5sum of SAM files didn't match between the data generated and the sample baseline. Best to match the outputs manually\n")
            return False

    #BAI files
    for output_file in bai_files:
        with open(output_file, 'rb') as file_handle:
            output_md5 = md5(file_handle.read()).digest()

        bare_filename = os.path.basename(output_file)

        with open(os.path.join(baseline_folder, bare_filename), 'rb') as file_handle:
            sample_md5 = md5(file_handle.read()).digest()

        if sample_md5 != output_md5:
            print("md5sum of BAI files didn't match between the data generated and the sample baseline. Best to match the outputs manually\n")
            return False

    return True

# Verifying the plots exist
def verify_plots_exist(folder):
    csv_files = glob(folder + '/sim1M*/post_processed_stats/*.csv')
    png_files = glob(folder + '/sim1M*/post_processed_stats/*.png')

    # Make sure there are some files
    if any((not csv_files, not png_files)):
        return False

    for output_file in csv_files + png_files:
        # Check filesize
        if os.path.getsize(output_file) <= 0:
            return False

    return True

# Verifying the plots md5sum against a baseline
def verify_plots_md5(folder):
    sample_folder = "baseline_sim1M_pairs/post_processed_stats"
    png_files = glob(folder + '/sim1M*/post_processed_stats/*.png')
    import pdb; pdb.set_trace()

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

# Test for kill scripts PIDs
def check_if_sar_or_iostat_running(username=None):
    if not username:
         username = os.popen("whoami").read()

    ps_output = os.popen("ps -u %s" % username).read().split('\n')
    for line in ps_output:
        if re.search(r'.*\ssar.*', line):
            return True
        if re.search(r'.*\siostat.*', line):
            return True
    return False

class MetricsMonitor(threading.Thread):
     def __init__(self):
         super(MetricsMonitor, self).__init__()
         self.sar_or_iostat_ran = False
         self.time_to_die = False

     def run(self):
         user = os.popen("whoami").read()
         while True:
             if check_if_sar_or_iostat_running(user):
                 self.sar_or_iostat_ran = True
             if self.time_to_die:
		return

##########################################################################
# Test for failures - Validating arguments mainly.
class TestFailures(unittest.TestCase):
    def test_01_workflow_script_non_existent(self):
        print("\nworkflow Script non-existent. Fail!")
        folder_name = make_temp_output_folder()
        self.assertFalse(call("./workflow_profiler.py data_collection.pl ohsu sim1M_pairs 16 /data/genomes/simulated/ %s/ -Ap" % folder_name, shell=True) != 0)

    def test_02_input_folder_doesnt_exist(self):
        print("\nInput Directory non-existent. Fail!")
        folder_name = make_temp_output_folder()
        self.assertFalse(call("./workflow_profiler.py data_collection_ohsu_dnapipeline.pl ohsu sim1M_pairs 16 /data/simulated/ %s/ -Ap" % folder_name, shell=True) != 0)

    def test_03_output_folder_doesnt_exist(self):
        print("\nOutput Directory non-existent. Fail!")
        self.assertFalse(call("./workflow_profiler.py data_collection_ohsu_dnapipeline.pl ohsu sim1M_pairs 16 /data/genomes/simulated/ folderdoesntexist/ -Ap", shell=True) != 0)

    def test_04_stats_non_existent(self):
        print("\nStats are not provided. Fail!")
        folder_name = make_temp_output_folder()
        self.assertFalse(call("./workflow_profiler.py data_collection_ohsu_dnapipeline.pl ohsu sim1M_pairs 16 /data/genomes/simulated/ %s/" % folder_name, shell=True) != 0)

#################################################################################
# Test for success - kill-scripts validation and run a full workflow.
class TestSuccesses(unittest.TestCase):
    def __init__(self, x):
        signal.signal(signal.SIGINT, self.signal_handler)
        super(TestSuccesses, self).__init__(x)
        self.metric_monitor = MetricsMonitor()

    def setUp(self):
        self.metric_monitor.start() # This actually causes the thread to run

    def tearDown(self):
        with open(os.devnull, "w") as fnull:
            call (["kill_scripts/kill-sar-script.sh"], stdout=fnull, stderr=fnull)
            call (["kill_scripts/kill-iostat-script.sh"], stdout=fnull, stderr=fnull)
        metric_seen = self.metric_monitor.sar_or_iostat_ran
        self._kill_children()
        print("Did we see that a metric was running at all? : " + str(metric_seen))
        self.assertTrue(metric_seen)


    def test_01_kill_scripts(self):
        print("\nTesting that metrics are killed")
        folder_name = make_temp_output_folder()
        folder_name = folder_name + "/test_kill_scripts"  # profiler validates the folder name, so don't want this to error out
        os.makedirs(folder_name)
        self.assertTrue(call("./workflow_profiler.py data_collection_1stage_dnapipeline.pl ohsu sim1M_pairs 16 /data/genomes/simulated/ %s -int 5 -pp 0 -A" % folder_name, shell=True) == 0)
        self.assertFalse(check_if_sar_or_iostat_running())

    def test_02_normal_run(self):
        print("\nTesting small workflow")
        folder_name = make_temp_output_folder()
        folder_name = folder_name + "/test_full"  # profiler validates the folder name, so don't want this to error out
        os.makedirs(folder_name)
        self.assertTrue(call("./workflow_profiler.py data_collection_4stages_dnapipeline.pl ohsu sim1M_pairs 16 /data/genomes/simulated/ %s -int 5 -Ap" % folder_name, shell=True) == 0)
        md5_workflow_verify = verify_workflow_data(folder_name)
        print("Did the workflow data verify against baseline? : " + str(md5_workflow_verify))
        made_plots = verify_plots_exist(folder_name)
        print("Did the plots get made? : " + str(made_plots))
        self.assertTrue(md5_workflow_verify)
        self.assertTrue(made_plots)

    def signal_handler(self, signal, frame):
        print('You pressed Ctrl+C!')
        self._kill_children()

    def _kill_children(self):
        self.metric_monitor.time_to_die = True
        self.metric_monitor.join()

#################################################################################
def run_tests():
    unittest.main()

if __name__ == "__main__":
    run_tests()
