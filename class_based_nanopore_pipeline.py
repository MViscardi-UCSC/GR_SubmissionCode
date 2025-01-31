"""
class_based_nanopore_pipeline.py
Marcus Viscardi,    March 24, 2023

This is a semi temporary file. The goal is to move the majority
of my pipeline functionality to a class based system. The main
reason for this is that there are many flags that should be carried
between the separate steps/functions of the pipeline.

I am not sure if this is the right approach to the problem at hand,
but it is the best way I can think of!

This will also provide an opportunity to fix some shape edges
that the current pipeline has, including:

 - Lack of comprehensive logging functionality...
    
    - A large hangup in the past has been getting this to
      work with the called shell functions, but I should be
      able to do with now with a better understanding of how
      this works.
 
 - Inability to have some steps function on their own
    
    - For example, the basecaller doesn't concatenate to
      a single fastq because the "concatenate" step is its
      own thing!
 
 - Ability to diagnose common errors that occur...
   
   - The problem here is that most (all) of my pathways to
     troubleshoot errors that come up are rattling around in
     my head. These should really be integrated into the code
     so that future users don't have to rely on me. At least,
     I can make an effort to solve the easy ones!
   
   - An easy example of this is an empty fastq being passed to
     early steps of pipeline. This error will not be caught until
     WAY down the line, and should be easy to identify and notify
     the user of.
 
 - Better path handling with pathlib. Half of the current errors
   that occur are due to missing files or lack of ownership/access.

Additionally, this will be an opportunity to update with some disparate
pieces of code I have around this repo.

 - Integrate new standards' functionality, drop old standards.
 
 - Integrate Nano3P functionality!!
 
 - Add some general QC steps to produce print-outs and plots.
"""
import subprocess
import sys
import traceback
from os import environ
from argparse import ArgumentParser
from pprint import pprint
from typing import List
from pathlib import Path
from glob import glob
import logging
import warnings
from collections import namedtuple
import time
from functools import wraps

from tqdm import tqdm

from nanoporePipelineCommon import find_newest_matching_file, get_dt, \
    gene_names_to_gene_ids, SamOrBamFile, FastqFile, gtf_to_df

import pandas as pd

# This is also in my .bashrc but that doesn't always seem to trigger,
#   ie. if I run the script from inside pycharm.
# If the HDF5 path isn't specified nanopolish freaks out, this solution is based on:
#   https://stackoverflow.com/questions/5971312/how-to-set-environment-variables-in-python
environ['HDF5_PLUGIN_PATH'] = '/usr/local/hdf5/lib/plugin'

GUPPY_PATH = Path("/opt/ont/guppy/bin/guppy_basecaller")
GUPPY_6_3_8_PATH = Path("/data16/marcus/scripts/ont-guppy/bin/guppy_basecaller")


# noinspection PyArgumentList
class NanoporePipeline:
    """
    This is a class to handle Marcus' nanopore pipeline. Including:
     - Parsing the command line arguments
     - Parsing the settings file
     - Running the requested steps
     - Logging the steps, their inputs, and their outputs
     - Handling errors
     - Maybe making some plots?
    """

    def __init__(self, skip_cli_dict: dict = None):
        """
        :param skip_cli_dict: This will allow the user to skip the command
                              line and instead pass a dictionary of arguments.
        """
        self.absolute_default_dict = {'printArgs': False,
                                      'nestedData': False,
                                      'regenerate': False,
                                      'altGenomeDirs': [],
                                      'threads': 20,
                                      'stepsToRun': "GMNCFPS",
                                      'sampleID': "sample1",
                                      'condition': "conditionA",
                                      'minimapParam': "-x splice -uf -k14",
                                      'guppyConfig': "rna_r9.4.1_70bps_hac.cfg",
                                      'tera3adapter': False,
                                      'tera5adapter': False,
                                      'extraGuppyOptions': False,
                                      'callWithJoshMethod': False,
                                      'freezeGuppyVersion6_3_8': False,
                                      'do_not_log': False,
                                      }

        self.arg_dict = dict()

        if skip_cli_dict:
            self.cli_arg_dict = skip_cli_dict
        else:
            self.cli_arg_dict = self._parse_args()
        self.settings_dict = self._parse_settings(**self.cli_arg_dict)

        self.arg_dict.update(self.absolute_default_dict)
        self.arg_dict.update(self.settings_dict)
        self.arg_dict.update(self.cli_arg_dict)

        # This reformatting helps to ensure that the types are correct:
        for key, arg in self.arg_dict.items():
            if self.arg_dict[key] == "True" or self.arg_dict[key] is True:
                self.arg_dict[key] = True
            elif self.arg_dict[key] == "False" or self.arg_dict[key] is False:
                self.arg_dict[key] = False
            elif key == "altGenomeDirs":
                self.arg_dict[key] = list(arg)
            elif key == "extraGuppyOptions":
                self.arg_dict[key] = arg.strip('"')
            else:
                try:
                    self.arg_dict[key] = int(arg)
                except ValueError or TypeError:
                    self.arg_dict[key] = str(arg)

        # Create attributes for all the arguments:
        self.settings_file = Path(self.arg_dict['settings_file'])
        self.output_dir = Path(self.arg_dict['outputDir'])
        self.genome_dir = Path(self.arg_dict['genomeDir'])
        self.alt_genome_dirs = [Path(alt_genome_dir) for alt_genome_dir in self.arg_dict['altGenomeDirs']]
        self.data_dir = Path(self.arg_dict['dataDir'])
        self.fast5_dir = self.data_dir / "fast5"
        self.threads = self.arg_dict['threads']
        self.guppy_config = self.arg_dict['guppyConfig']
        self.steps_to_run = self.arg_dict['stepsToRun']
        self.sample_id = self.arg_dict['sampleID']
        self.condition = self.arg_dict['condition']
        self.minimap_param = self.arg_dict['minimapParam']
        self.tera3adapter = self.arg_dict['tera3adapter']
        self.tera5adapter = self.arg_dict['tera5adapter']
        self.tera_trimming_ran = False
        self.do_not_log = self.arg_dict['do_not_log']

        self.extra_guppy_options = self.arg_dict['extraGuppyOptions']
        self.freeze_guppy_version_6_3_8 = self.arg_dict['freezeGuppyVersion6_3_8']
        if isinstance(self.extra_guppy_options, str):
            if not self.extra_guppy_options.endswith(" "):
                self.extra_guppy_options = self.extra_guppy_options + " "
        else:
            self.extra_guppy_options = ""

        if self.freeze_guppy_version_6_3_8 or "fast5_out" in self.extra_guppy_options:
            self.guppy_basecaller_path = GUPPY_6_3_8_PATH
            if "fast5_out" not in self.extra_guppy_options:
                self.extra_guppy_options += "--fast5_out "
        else:
            self.guppy_basecaller_path = GUPPY_PATH

        self.call_with_josh_method = self.arg_dict['callWithJoshMethod']
        self.regenerate = self.arg_dict['regenerate']
        self.drop_genes_with_hits_less_than = self.arg_dict['dropGeneWithHitsLessThan']

        if self.do_not_log:
            self.log_dir = Path(f"/tmp/{get_dt(extended_for_file=True)}_nanoporePipeline_logs")
        else:
            self.log_dir = self.output_dir / "logs"
        self.log_dir.mkdir(parents=True, exist_ok=True)
        self.run_setup_file = self.log_dir / f"{get_dt(extended_for_file=True)}_" \
                                             f"{self.sample_id}_{self.condition}.run.txt"
        self.run_setup_file.touch(exist_ok=True)
        with open(self.run_setup_file, 'w') as f:
            f.write(f"Run setup file for {self.sample_id} {self.condition}\n")
            f.write(f"Started at {get_dt(for_print=True)}\n")
            f.write(f"{'=' * 50}\n")
            f.write(f"{'Argument':>25}   {'Value'}\n")
            for key, value in self.arg_dict.items():
                f.write(f"{key:>25} = {str(value)}\n")
        self.log_file = self.log_dir / f"{get_dt(extended_for_file=True)}_" \
                                       f"{self.sample_id}_{self.condition}.log"
        self._setup_logging()

        # Define a named tuple to store expected output directory and explanation
        OutputStep = namedtuple('OutputStep', ['directory', 'description', 'function', 'step_order_number'])
        self.steps_dict = {"G": OutputStep("fastqs",
                                           "Guppy Basecalling",
                                           self.guppy_basecaller,
                                           1),
                           "A": OutputStep("cat_files",
                                           "Filtering Alt. Genomes (currently pretty slow)",
                                           self.filter_alt_genomes,
                                           2),
                           "T": OutputStep("cat_files",
                                           "Trimming 5TERA-Seq Adapters",
                                           self.trim_5TERA_adapters,
                                           3),
                           "M": OutputStep("cat_files",
                                           "Minimap2 and SamTools",
                                           self.minimap2_and_samtools,
                                           4),
                           "N": OutputStep("nanopolish",
                                           "Nanopolish Index and polyA Calling",
                                           self.nanopolish,
                                           5),
                           "F": OutputStep("featureCounts",
                                           "FeatureCounts",
                                           self.featureCounts,
                                           6),
                           "P": OutputStep("merge_files",
                                           "Merging Results w/ Pandas",
                                           self.merge_files,
                                           7),
                           "S": OutputStep("merge_files",
                                           "Mapping ENO2 Standards (version 2 assignment method)",
                                           self.assign_ENO2_standards,
                                           8),
                           "L": OutputStep("flair",
                                           "Calling Transcripts w/ Flair",
                                           self.flair_for_transcripts,
                                           9),
                           "X": OutputStep("extras",
                                           "Running random eXtra steps",
                                           self.extra_steps,
                                           10),
                           }
        # If the library was cDNA then we need to run tailfindr instead of nanopolish!:
        self.tailfindr_tag = "N" in self.steps_to_run.upper() and any([self.guppy_config.startswith("dna"),
                                                                       self.guppy_config.startswith("cdna")])
        if self.tailfindr_tag:
            self.steps_dict["N"] = OutputStep("tailfindr",
                                              "Tailfindr polyA Calling",
                                              self.tailfindr,
                                              5)
        # Only keep the steps that are in the steps_to_run list:
        self.steps_dict = {key: value for key, value in self.steps_dict.items()
                           if key.upper() in self.steps_to_run.upper()}
        # Create a list of the output directories and make them if they don't exist:
        self.sub_dirs = list(set([self.output_dir / value.directory for key, value in self.steps_dict.items()]))
        self.sub_dirs_dict = {str(directory.name): directory for directory in self.sub_dirs}
        self._setup_directories()
        self.fastq_dir = self.output_dir / "fastqs"
        self.cat_files_dir = self.output_dir / "cat_files"
        self.cat_files_dir.mkdir(exist_ok=True)  # I think we ALWAYS need the cat_files dir...
        self.tailfindr_dir = self.output_dir / "tailfindr"  # These are both here, but only one will be used.
        self.nanopolish_dir = self.output_dir / "nanopolish"  # ^^^
        self.featureCounts_dir = self.output_dir / "featureCounts"
        self.merge_files_dir = self.output_dir / "merge_files"
        self.flair_dir = self.output_dir / "flair"
        self.extras_dir = self.output_dir / "extras"

        # Add attributes for each step to keep track of whether it has been run:
        self.trim_5TERA_adapters_ran = False
        self.minimap_ran = False
        self.nanopolish_ran = False
        self.tailfindr_ran = False
        self.featureCounts_ran = False
        self.assign_standards_ran = False
        self.flair_ran = False

        # TODO: Make the below checks more robust:
        #       The idea would be to actually check for the output files/edits.
        # TODO: Should all of these checks, inclusing the regen stuff, happen
        #       here instead of in the individual functions?
        #       I think so, but I'm not sure. Revisit this.
        if "G" not in self.steps_to_run.upper():
            self.guppy_basecaller_ran = True
        if "T" not in self.steps_to_run.upper():
            self.trim_5TERA_adapters_ran = True
        if "M" not in self.steps_to_run.upper():
            self.minimap_ran = True
        if "S" not in self.steps_to_run.upper():
            self.assign_standards_ran = True
        if "L" not in self.steps_to_run.upper():
            self.flair_ran = True

    def _parse_args(self) -> dict:
        """Parse the arguments from the command line."""
        self.parser = ArgumentParser(description="A pipeline to handle the outputs "
                                                 "of nanopore runs!")
        # Required Arguments:
        self.parser.add_argument('settings_file', metavar='settings_file', type=str,
                                 help="A .txt file with inputs all in arg|value format. "
                                      "Any of the below arguments will also function if put "
                                      "into this file.")
        # Arguments that can be included in the settings file
        self.parser.add_argument('--outputDir', metavar='outputDir', type=str, default=None,
                                 help="Directory for the output of all resulting files.")
        self.parser.add_argument('--genomeDir', metavar='genomeDir', type=str, default=None,
                                 help="Path to genome directory.")
        self.parser.add_argument('--dataDir', metavar='dataDir', type=str, default=None,
                                 help="Path to sequencing data directory.")
        self.parser.add_argument('--threads', metavar='threads', type=int, default=None,
                                 help="Number of threads to be used by nanopolish and minimap2. [20]")
        self.parser.add_argument('--guppyConfig', metavar='guppyConfig', type=str,
                                 help="Configuration preset passed to the guppy_basecaller "
                                      "based on flowcell and kit used for run. Helpful "
                                      "table for picking a config @ https://denbi-nanopore-"
                                      "training-course.readthedocs.io/en/latest/basecalling/"
                                      "basecalling_1.html or just google 'guppy_basecalling'"
                                      "[rna_r9.4.1_70bps_hac.cfg]")
        self.parser.add_argument('--dropGeneWithHitsLessThan', metavar='dropGeneWithHitsLessThan',
                                 type=int, default=None,
                                 help="Minimum number of reads per gene to have in the "
                                      "compressedOnGenes outputs.")
        self.parser.add_argument('--altGenomeDirs', metavar='altGenomeDirs',
                                 nargs='*', type=List[str], default=None,
                                 help="Alternative genomes that can be used for filtering out reads "
                                      "that map to them")
        self.parser.add_argument('--stepsToRun', metavar='stepsToRun',
                                 type=str, default=None,
                                 help="Steps to run within the pipeline: (G)uppy basecalling, "
                                      "(T)rim TERA-Seq adapters, "
                                      "(M)inimap, (N)anopolish, (F)eature counts, (C)oncat files, "
                                      "merge with (P)andas, use f(L)air to ID transcripts, "
                                      "map ENO2 nanopore (S)tandards, and/or random e(X)tra "
                                      "steps (plotting). [GMNCFPS]")
        self.parser.add_argument('--sampleID', metavar='sampleID',
                                 type=int, default=None,
                                 help="sampleID to pass to FLAIR [sample1]")
        self.parser.add_argument('--condition', metavar='condition',
                                 type=int, default=None,
                                 help="condition to pass to FLAIR [conditionA]")
        self.parser.add_argument('--minimapParam', metavar='minimapParam',
                                 type=str, default=None,
                                 help="Arguments to pass to minimap2. If used on the command line, "
                                      "be sure to surround in double-quotes! [\"-x splice -uf -k14\"]")
        self.parser.add_argument('--tera3adapter', metavar='tera3adapter',
                                 type=int, default=None,
                                 help="Adapter to be trimmed for TERA3 [None]")
        self.parser.add_argument('--tera5adapter', metavar='tera5adapter',
                                 type=int, default=None,
                                 help="Adapter to be trimmed for 5TERA [None]")
        self.parser.add_argument('--extraGuppyOptions', metavar='extraGuppyOptions',
                                 type=int, default=None,
                                 help="String flags/options to be added to the guppy_basecaller call [None]")
        # Flag Arguments
        self.parser.add_argument('-p', '--printArgs', action='store_true',
                                 help="Boolean flag to show how arguments are overwritten/accepted")
        self.parser.add_argument('-n', '--nestedData', action='store_true',
                                 help="Boolean flag that will account for data coming out of gridIONs "
                                      "as these will produce large nested fast5 dictionaries.")
        self.parser.add_argument('-r', '--regenerate', action='store_true',
                                 help="Boolean flag to ignore previously produced files "
                                      "and generate all files anew")
        self.parser.add_argument('-j', '--callWithJoshMethod', action='store_true',
                                 help="Boolean flag to use Josh's read assignment method, rather"
                                      "than FeatureCounts")
        self.parser.add_argument('-F', '--freezeGuppyVersion6_3_8', action='store_true',
                                 help="A Boolean flag to force the use of guppy_basecaller v6.3.8."
                                      "This version was the last to offer --fast5_out support, which"
                                      "is required for tailfindr. Note that if '--fast5_out' is present"
                                      "in the '--extraGuppyOptions' string, this flag will automatically"
                                      "be activated!")
        self.parser.add_argument('--do_not_log', action='store_true',
                                 help="Boolean flag to turn off logging")
        self.parser.add_argument('-T', '--useTailfindr', action='store_true',
                                 help="Boolean flag to use tailfindr instead of nanopolish")

        # Spit out namespace object from argParse
        args = self.parser.parse_args()
        return {arg: vars(args)[arg] for arg in vars(args) if vars(args)[arg] is not None}

    def _parse_settings(self, settings_file, printArgs=False, **other_kwargs_from_cli) -> dict:
        """
        This function will parse a settings file and return a dictionary of the arguments
        :param settings_file: Path to settings file
        :param printArgs: Option to print the arguments as they are parsed
        :param other_kwargs_from_cli: Catch-all for other arguments passed from the command line
        :return: Dictionary of arguments that came from the settings file
        """
        # first, parse the settings file into _settingsDict
        settingsDict = {}

        with open(settings_file, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    line = line.strip()
                    if line != '':
                        line = line.split('|')
                        if len(line) == 2:
                            if line[0] == "altGenomeDirs":
                                settingsDict[line[0]] = line[1].split(",")
                            else:
                                settingsDict[line[0]] = line[1]
                        else:
                            raise NotImplementedError("\033[31;1m\nRemove pipes ('|') from settings "
                                                      "file arguments (or rewrite parser)\n\033[0m")
        if printArgs:
            print(f"\nSettings Arguments (file: '{settings_file}')")
            for key, arg in settingsDict.items():
                if not arg:
                    print(f"\t{key} = {arg} -> (Will not be passed)")
                else:
                    print(f"\t{key} = {arg}")
        settingsDict = {k: v for k, v in settingsDict.items() if v is not None and v != ''}
        return settingsDict

    def _setup_directories(self):
        """
        Will create the output directory and subdirectories
        """
        # Create the output directory
        if not self.output_dir.exists():
            self.output_dir.mkdir()
        # Create the subdirectories
        for sub_dir in self.sub_dirs:
            sub_dir.mkdir(exist_ok=True)

    def _setup_logging(self):
        """
        Will set up the logging file
        """
        # Create a custom logger
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)

        # Create handlers
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(logging.DEBUG)
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)

        # Create a formatter that includes the time, logging level, and message
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        console_handler.setFormatter(formatter)

        # Add the handlers to the logger
        self.logger.addHandler(file_handler)
        self.logger.addHandler(console_handler)

        self.logger.info(f"Initial set up complete.")
        self.logger.info(f"Check log file for more details: {self.log_file.parent.name}/{self.log_file.name}")
        self.logger.debug(f"Settings file: {self.settings_file}")
        self.logger.debug(f"Run parameters:  {self.run_setup_file}")

    def seconds_to_hms(self, seconds):
        """
        Converts seconds to hours:minutes:seconds
        """
        hours = int(seconds // 3600)
        minutes = int((seconds % 3600) // 60)
        seconds = seconds % 60
        return f"{hours:02}:{minutes:02}:{seconds:02.2f}"

    def pipeline_step_decorator(func):
        """
        This is a decorator that will be used to time each step of the pipeline.
        """
        # So this currently works w/ nanopore_run_obj down below,
        # but I'm not sure if it's the best way...
        # Or if it will even work in all cases!
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            start_time = None  # This is so we can check if we ever started the timer
            try:
                start_time = time.time()
                self.logger.debug(f"Starting {func.__name__}")
                result = func(self, *args, **kwargs)
                end_time = time.time()
                elapsed_seconds = end_time - start_time
                self.logger.debug(f"Finished {func.__name__}")
                self.logger.debug(f"Running  {func.__name__} took {self.seconds_to_hms(elapsed_seconds)}")
                return result
            except Exception as e:
                traceback_path = self.log_dir / f"{get_dt(extended_for_file=True)}_traceback.txt"
                error_msg = f"Error in {func.__name__}: {e.__repr__()}"
                self.logger.error(error_msg + f", find traceback in {traceback_path.parent.name}/{traceback_path.name}")
                with open(traceback_path, 'w') as f:
                    f.write(f"\n{error_msg}\n\n")
                    if start_time:
                        end_time = time.time()
                        elapsed_seconds = end_time - start_time
                        f.write(f"Error occurred at {get_dt(for_print=True)} after the function ran for "
                                f"{self.seconds_to_hms(elapsed_seconds)}.\n\n")
                    else:
                        f.write(f"Error occurred at {get_dt(for_print=True)}\n\n")
                    f.write(f"Traceback Below:\n{'=' * 80}\n\n")
                    traceback.print_exc(file=f)
                raise e

        return wrapper

    #  "G": OutputStep("fastq", "Guppy Basecalling"),
    #  "A": OutputStep("cat_files", "Filtering Alt. Genomes (currently pretty rough and slow and unimplemented)"),
    #  "T": OutputStep("cat_files", "Trimming TERA-Seq Adapters"),
    #  "M": OutputStep("cat_files", "Minimap2 and SamTools"),
    #  "N": OutputStep("nanopolish", "Nanopolish Index and polyA Calling"),
    #  "F": OutputStep("featureCounts", "FeatureCounts"),
    #  "P": OutputStep("merge_files", "Merging Results w/ Pandas"),
    #  "S": OutputStep("merge_files", "Mapping ENO2 Standards (version 2 assignment method)"),
    #  "L": OutputStep("flair", "Calling Transcripts w/ Flair"),
    #  "X": OutputStep("extras", "Running random eXtra steps"),

    def run_cmd(self, command, command_nickname, save_output_to_file=True):
        """
        Runs a command and logs its output
        :param command: A string of the command to run
        :param command_nickname: A short string to describe the command
        :param save_output_to_file: An optional boolean to save the output to a file
        :return: None
        """
        # set up output file
        output_file = self.log_dir / f"{get_dt(extended_for_file=True)}_{command_nickname}.log"
        with open(output_file, "w") as f:
            # run the command and log/output its output
            log_string = f"Running command ({command_nickname})"
            if save_output_to_file:
                log_string += f", logging to {output_file.parent.name}/{output_file.name}"
            self.logger.info(log_string)
            self.logger.debug(f"Command call: {command}")
            with subprocess.Popen(command, stdout=subprocess.PIPE, shell=True,
                                  bufsize=1, universal_newlines=True) as p:
                for line in p.stdout:
                    print(line, end="")
                    if save_output_to_file:
                        f.write(line)
                p.communicate()

            # check the return code and log the end of the command
            if p.returncode != 0:
                self.logger.error(f"Command ({command_nickname}) failed with return code {p.returncode}")
                raise subprocess.CalledProcessError(p.returncode, p.args)
            else:
                self.logger.info(f"Command ({command_nickname}) completed successfully")
                return p.returncode

    @pipeline_step_decorator
    def guppy_basecaller(self):
        """
        Runs guppy basecaller on the fast5 files in the data_dir.
        Saves and concatenates the fastq files to the fastq_dir.
        """
        guppy_run_tag_file = self.fastq_dir / "guppy_run_tag.txt"
        guppy_flag = self.regenerate or not guppy_run_tag_file.exists()
        if guppy_flag:
            call = rf"""{self.guppy_basecaller_path} -x "cuda:0" """ \
                   rf"""{self.extra_guppy_options}--num_callers 12 """ \
                   rf"""--gpu_runners_per_device 8 """ \
                   rf"""-c {self.guppy_config} -i {self.data_dir}/fast5 -s {self.fastq_dir} """
            self.run_cmd(call, "guppy_basecalling", save_output_to_file=True)
            guppy_run_tag_file.touch()
            guppy_run_tag_file.write_text(f"Last run on {get_dt(for_print=True)}.\n")
            guppy_run_tag_file.write_text(f"Run command: {call}\n")
            self.regenerate = True
        else:
            self.logger.info("Guppy basecalling already completed, skipping")

        concat_fastq = self.cat_files_dir / "cat.fastq"
        concat_flag = self.regenerate or not concat_fastq.exists()
        if concat_flag:
            call = rf"cat {self.fastq_dir}/pass/*.fastq > {self.cat_files_dir}/cat.fastq"
            self.run_cmd(call, "concatenate_fastqs", save_output_to_file=False)
            self.regenerate = True
        else:
            self.logger.info("Fastq concatenation already completed, skipping")

    @pipeline_step_decorator
    def filter_alt_genomes(self):
        raise NotImplementedError

    @pipeline_step_decorator
    def trim_5TERA_adapters(self):
        """
        Trims the 5' TERA-Seq adapters from the fastq file. (Saves a backup of the original fastq file before trimming)
        """
        cat_fastq_path = self.cat_files_dir / "cat.fastq"
        fastq_backup_path = self.cat_files_dir / "cat.untrimmed.fastq"
        if not fastq_backup_path.exists():
            # First backup the fastq file that the real minimap2 call will need:
            call = f"cp {cat_fastq_path} {fastq_backup_path}"
            self.run_cmd(call, "backup_fastq", save_output_to_file=False)
        else:
            self.logger.info(f"Skipping backup of cat.fastq file because it already happened.")

        # Then, we will want to check that adapter trimming hasn't already happened.
        #   This isn't quite as simple as file checking, but the easiest way I can
        #   think of will be to parse the first line of the fastq and check if
        #   'adapter' is in there:
        with open(cat_fastq_path, 'r') as fastq_file:
            first_line = fastq_file.readline()
            cutadapt_was_run = 'TERAADAPTER' in first_line
            self.tera_trimming_ran = cutadapt_was_run
        if not cutadapt_was_run or self.regenerate:
            cutadapt_call = None  # So that we can check if it was set later
            if isinstance(self.tera5adapter, str) and isinstance(self.tera3adapter, str):
                cutadapt_call = f"cutadapt --action=trim -j {self.threads} " \
                                f"-g TERA5={'X' + self.tera5adapter} --overlap 31 --error-rate 0.29 " \
                                "--rename '{id} TERAADAPTER5={adapter_name} {comment}' " \
                                f"{fastq_backup_path} | " \
                                f"cutadapt --action=trim -j {self.threads} " \
                                f"-g TERA3={self.tera3adapter} --overlap 16 --error-rate 0.18 " \
                                "--rename '{id} TERAADAPTER3={adapter_name} {comment}' " \
                                f"- > {cat_fastq_path}"
            elif isinstance(self.tera5adapter, str):
                cutadapt_call = f"cutadapt --action=trim -j {self.threads} " \
                                f"-g TERA5={'X' + self.tera5adapter} --overlap 31 --error-rate 0.29 " \
                                "--rename '{id} TERAADAPTER5={adapter_name} {comment}' " \
                                f"{fastq_backup_path} > {cat_fastq_path}"
            elif isinstance(self.tera3adapter, str):
                cutadapt_call = f"cutadapt --action=trim -j {self.threads} " \
                                f"-a TERA3={self.tera3adapter} --overlap 16 --error-rate 0.18 " \
                                "--rename '{id} TERAADAPTER3={adapter_name} {comment}' " \
                                f"{fastq_backup_path} > {cat_fastq_path}"
            else:
                self.logger.warning(f"Please provide 5TERA and/or TERA3 adapters as strings!! "
                                    f"You passed: {self.tera5adapter} and {self.tera3adapter}")
                self.logger.info("Skipping cutadapt, and moving backup file back to cat.fastq")
                call = f"mv {fastq_backup_path} {cat_fastq_path}"
                self.run_cmd(call, "restore_backup_fastq", save_output_to_file=False)

            if isinstance(cutadapt_call, str):
                print(f"Starting cutadapt at {get_dt(for_print=True)}\nUsing call:\t{cutadapt_call}\n")
                self.run_cmd(cutadapt_call, "tera_adapter_trimming", save_output_to_file=True)
                self.tera_trimming_ran = True
            self.regenerate = True
        else:
            self.logger.info(f"Skipping cutadapt because it has already been run on cat.fastq.")
        self.trim_5TERA_adapters_ran = True

    def _tera_adapter_tagging__(self):
        """
        This method will add the TERA3 and/or TERA5 adapter tags to the bam/sam files.
        :return: 
        """
        import simplesam as ssam
        fastq_path = self.cat_files_dir / "cat.fastq"

        if self.tera_trimming_ran:
            # This open block is just to look at the first line of the fastq
            with open(fastq_path, 'r') as f:
                first_line = f.readline()
                # Check if each of the cutadapt comments were added, save for below.
                tera3_was_run = 'TERAADAPTER3' in first_line
                tera5_was_run = 'TERAADAPTER5' in first_line

            # If neither were run, just skip the rest of this method!
            if not tera3_was_run and not tera5_was_run:
                raise ValueError(f"Neither TERA3 nor TERA5 adapter trimming was run on {fastq_path}!")

        # If the adapter tag IS found in the cat.fastq, we'll add it to the bam/sam files!:
        tagged_fastq_df = FastqFile(fastq_path).df

        # For the two adapters, either extract the info if it's there, or default to false if not.
        if tera3_was_run:  # Parse out TERA3 adapter if it existed
            tagged_fastq_df['t3'] = tagged_fastq_df.comment.str.extract(r'TERAADAPTER3=(\S+)').replace(
                {'no_adapter': '-',
                 'TERA3': '+'})
        else:
            tagged_fastq_df['t3'] = '-'

        if tera5_was_run:  # Parse out TERA5 adapter if it existed
            tagged_fastq_df['t5'] = tagged_fastq_df.comment.str.extract(r'TERAADAPTER5=(\S+)').replace(
                {'no_adapter': '-',
                 'TERA5': '+'})
        else:
            tagged_fastq_df['t5'] = '-'

        # Finally we'll load and iterate through the bam file, creating a new sam file along the
        #   way and adding in the new tags!!
        tagged_fastq_df.set_index('read_id', inplace=True)
        input_bam = self.cat_files_dir / "cat.sorted.bam"
        output_sam = self.cat_files_dir / "cat.sorted.sam"
        self.logger.info(f"Starting sam file tagging with TERA-seq adapter information")
        with ssam.Reader(open(input_bam, 'r')) as in_bam:
            with ssam.Writer(open(output_sam, 'w'), in_bam.header) as out_sam:
                row_iterator = tqdm(in_bam)
                for read in row_iterator:
                    row_iterator.set_description(f"Tagging {read.qname}")
                    read['t5'], read['t3'] = tagged_fastq_df.loc[read.qname, ['t5', 't3']].tolist()
                    out_sam.write(read)
        self.logger.info(f"Finished sam file tagging with TERA-seq adapter information")
        # Finally, we'll overwrite the old bam with the new, tagged sam file, and index it:
        call = f'samtools view -b {output_sam} -o {input_bam}'
        self.run_cmd(call, "tagged_sam_to_bam", save_output_to_file=False)
        
        # We also want to index the new sam/bam files:
        call = f'samtools index {output_sam}'
        self.run_cmd(call, "tagged_sam_to_bam", save_output_to_file=False)
        call = f'samtools index {input_bam}'
        self.run_cmd(call, "tagged_sam_to_bam", save_output_to_file=False)

    @pipeline_step_decorator
    def minimap2_and_samtools(self):
        """
        This method will run minimap2 and samtools to align the cat.fastq file to the genome.
        :return: 
        """
        cat_bam_path = self.cat_files_dir / "cat.bam"
        minimap_flag = self.regenerate or not cat_bam_path.exists()
        if not minimap_flag:
            bam_length = cat_bam_path.stat().st_size
            if bam_length == 0:
                minimap_flag = True  # This is all to catch screwed up runs that have empty bam files!!!
        if minimap_flag:
            genome_fa_file = glob(f"{self.genome_dir}/*allChrs.fa")
            genome_bed_file = glob(f"{self.genome_dir}/*.bed")
            if len(genome_fa_file) != 1:
                raise NotImplementedError(f"Currently this script only supports having genomeDirs "
                                          f"with one fa file that ends with 'allChrs.fa'")
            else:
                genome_fa_file = genome_fa_file[0]
            if len(genome_bed_file) != 1:
                raise NotImplementedError(f"Currently this script only supports having genomeDirs "
                                          f"with one bed file that ends with '.bed'")
            else:
                genome_bed_file = genome_bed_file[0]
            call = f"minimap2 -a {self.minimap_param} {genome_fa_file} {self.cat_files_dir}/cat.fastq " \
                   f"-t {self.threads} --junc-bed {genome_bed_file} | samtools view -b - -o " \
                   f"{cat_bam_path}"
            self.run_cmd(call, "minimap2_mapping", save_output_to_file=True)
            self.regenerate = True
        else:
            self.logger.info(f"Skipping minimap2 because it has already been run on cat.fastq.")

        sorted_cat_bam_path = self.cat_files_dir / "cat.sorted.bam"
        samtools_flag = self.regenerate or not sorted_cat_bam_path.exists()
        if samtools_flag:
            call = f"samtools sort -m 16G -T tmp -o {self.cat_files_dir}/cat.sorted.bam " \
                   f"{self.cat_files_dir}/cat.bam && samtools index " \
                   f"{self.cat_files_dir}/cat.sorted.bam"
            self.run_cmd(call, "samtools_sort_and_index", save_output_to_file=False)
            self.regenerate = True
        else:
            self.logger.info(f"Skipping samtools because it has already been run on cat.bam.")

        if self.tera_trimming_ran:
            # We need to tag BAM file with trimming flags (t5 and t3)
            self._tera_adapter_tagging__()
        else:
            call = f"samtools view {self.cat_files_dir}/cat.sorted.bam " \
                   f"> {self.cat_files_dir}/cat.sorted.sam"
            self.run_cmd(call, "bam_to_sam", save_output_to_file=False)

        bam_file_with_only_mappedAndPrimary = self.output_dir / "cat_files" / "cat.sorted.mappedAndPrimary.bam"
        concat_flag = self.regenerate or not bam_file_with_only_mappedAndPrimary.exists()
        if concat_flag:
            calls = [f"samtools view -b -F 0x904 {sorted_cat_bam_path} > {bam_file_with_only_mappedAndPrimary}",
                     # The above command will build a new bam file w/out reads w/ bit_flags:
                     #    0x004, UNMAP           =   reads whose sequence didn't align to the genome
                     #    0x100, SECONDARY       =   reads that are secondary alignments
                     #    0x800, SUPPLEMENTARY   =   reads that are supplemental alignments
                     f"samtools index {bam_file_with_only_mappedAndPrimary}",
                     f"samtools view {bam_file_with_only_mappedAndPrimary} "
                     f"> {self.cat_files_dir}/cat.sorted.mappedAndPrimary.sam",
                     ]
            self.logger.info(f"Starting SAM/BAM file cleanup")
            for num, call in enumerate(calls):
                self.run_cmd(call, f"bam-sam_cleanup_{num}_of_{len(calls)}", save_output_to_file=False)
            self.regenerate = True
        else:
            self.logger.info(f"Skipping SAM/BAM file cleanup because it has already been run on cat.sorted.bam.")
        self.minimap_ran = True

    @pipeline_step_decorator
    def nanopolish(self):
        """
        This method will run nanopolish on the cat.fastq file to generate a file with polyA tail lengths called.
        :return: 
        """
        # First we have to index the fastq_files!
        nanopolish_index_file = self.cat_files_dir / "cat.fastq.index.readdb"
        nanopolish_index_flag = self.regenerate or not nanopolish_index_file.exists()
        if nanopolish_index_flag:
            sequencing_summary_file = glob(f"{self.data_dir}/sequencing_summary*.txt")
            if len(sequencing_summary_file) > 1:
                raise IndexError(f"Too many matches for 'sequencing_summary' in sequencing directory!!\n"
                                 f"{sequencing_summary_file}")
            if len(sequencing_summary_file) == 0:
                warnings.warn(f"No matches for 'sequencing_summary' in sequencing directory!!\n"
                              f"{sequencing_summary_file}\nGoing to continue w/out any seq summary! THIS WILL BE SLOW!")
                seq_sum_call = ""
            else:
                sequencing_summary_file = sequencing_summary_file[0]
                seq_sum_call = f"--sequencing-summary={sequencing_summary_file} "

            call = f"nanopolish index --directory={self.data_dir}/fast5 " \
                   f"{seq_sum_call}{self.cat_files_dir}/cat.fastq"
            self.logger.warning(f"nanopolish index is very slow and has limited outputs, you've been warned!!")
            self.run_cmd(call, "nanopolish_index", save_output_to_file=False)
            self.regenerate = True
        else:
            self.logger.info(f"Skipping nanopolish index because it has already been run on cat.fastq.")

        nanopolish_polya_output_file = self.nanopolish_dir / "polya.tsv"
        nanopolish_polya_flag = self.regenerate or not nanopolish_polya_output_file.exists()
        if nanopolish_polya_flag:
            genome_fa_file = glob(f"{self.genome_dir}/*allChrs.fa")
            if len(genome_fa_file) != 1:
                raise NotImplementedError(f"Currently this script only supports having genomeDirs "
                                          f"with one fa files that ends with 'allChrs.fa'")
            else:
                genome_fa_file = genome_fa_file[0]
            call = f"nanopolish polya --threads={self.threads} --reads={self.cat_files_dir}/cat.fastq " \
                   f"--bam={self.cat_files_dir}/cat.sorted.mappedAndPrimary.bam --genome={genome_fa_file} " \
                   f"> {nanopolish_polya_output_file}"
            self.run_cmd(call, "nanopolish_polya_calling", save_output_to_file=True)
            self.regenerate = True
        else:
            self.logger.info(f"Skipping nanopolish polya because it has already been run on fast5 files.")

        nanopolish_polya_passed_file = self.nanopolish_dir / "polya.passed.tsv"
        nanopolish_passed_flag = self.regenerate or not nanopolish_polya_passed_file.exists()
        if nanopolish_passed_flag:
            filter_call = f"head -n 1 {nanopolish_polya_output_file} > {nanopolish_polya_passed_file}; " \
                          f"grep PASS {nanopolish_polya_output_file} >> {nanopolish_polya_passed_file}"
            self.run_cmd(filter_call, "nanopolish_polya_filtering", save_output_to_file=False)
            self.regenerate = True
        else:
            self.logger.info(f"Skipping nanopolish polya filtering because it has already been run on polya.tsv.")

        nanopolish_polya_simplified_file = self.nanopolish_dir / "polya.passed.simple.tsv"
        nanopolish_simplified_flag = self.regenerate or not nanopolish_polya_simplified_file.exists()
        if nanopolish_simplified_flag:
            tails_df = pd.read_csv(nanopolish_polya_passed_file, sep="\t")
            tails_df = tails_df[["read_name", "contig", "position", "polya_length", "qc_tag"]]
            tails_df = tails_df.rename(columns={"read_name": "read_name",
                                                "contig": "chr",
                                                "position": "start"})
            tails_df.to_csv(nanopolish_polya_simplified_file, sep="\t", index=False)
            self.logger.info(f"Saved simplified polya file to {nanopolish_polya_simplified_file}")
            self.regenerate = True
        else:
            self.logger.info(f"Skipping nanopolish polya simplifying because "
                             f"it has already been run on polya.passed.tsv.")
        self.nanopolish_ran = True

    @pipeline_step_decorator
    def tailfindr(self):
        """
        This method will run tailfindr on the cat.fastq file to generate a file with polyA tail lengths called.
        :return: 
        """
        tailfindr_output_files = glob(f"{self.tailfindr_dir}/*_tailfindr.parquet")
        # TODO: Add check for previous run!
        # TODO: porechop??!! I don't even know!
        # First, lets make sure we have the fast5 files in out output directory:
        fast5_path = self.fastq_dir / "workspace"
        fast5_files = glob(f"{fast5_path}/*.fast5")
        if not fast5_path.exists() or len(fast5_files) == 0:
            # Raise an error if we don't have any fast5 files!
            raise FileNotFoundError(f"Could not find fast5 files in {fast5_path}! ")
        guppy_run_tag_file = self.fastq_dir / "guppy_run_tag.txt"
        guppy_run_tag_text = guppy_run_tag_file.read_text()
        if "--fast5_out" not in guppy_run_tag_text or "--trim_strategy none" not in guppy_run_tag_text:
            self.logger.warning(f"WARNING: You are not using the recommended guppy options for tailfindr! "
                                f"Tailfindr requires the following options: --fast5_out --trim_strategy none")
            self.logger.warning(f"I'll leave ^this up to you to fix, but I'm going to continue anyway...")
        # Now we can run tailfindr!
        r_command = f"Rscript -e 'library(tailfindr); " \
                    f"library(arrow); " \
                    f"Sys.setenv(HDF5_PLUGIN_PATH = \"/usr/local/hdf5/lib/plugin\"); " \
                    f"r_df <- find_tails(" \
                    f"fast5_dir = \"{fast5_path}\", " \
                    f"num_cores = {self.threads}, " \
                    f"save_dir = \"{self.tailfindr_dir}\"); " \
                    f"write_parquet(r_df, \"{self.tailfindr_dir}/{get_dt()}_tailfindr.parquet\"); " \
                    f"write.csv(r_df, \"{self.tailfindr_dir}/{get_dt()}_tailfindr.csv\")'"
        # This currently doesn't log very well because of the way tailfindr makes progress bars...
        # I'm not sure how to fix this, but it's not a huge deal because it's not a long running process.
        self.run_cmd(r_command, "tailfindr", save_output_to_file=True)
        # nanopore_run_obj.regenerate = True
        self.tailfindr_ran = True

    @pipeline_step_decorator
    def featureCounts(self):
        """
        This method will run featureCounts on the cat.sorted.mappedAndPrimary.bam
        file to identify the gene(s) that each read maps to.
        :return:
        """
        feature_counts_file = self.featureCounts_dir / "cat.sorted.mappedAndPrimary.bam.featureCounts"
        feature_counts_flag = self.regenerate or not feature_counts_file.exists()
        if feature_counts_flag:
            genome_gtf_file = glob(f"{self.genome_dir}/*.gtf")
            if len(genome_gtf_file) != 1:
                raise NotImplementedError(f"Currently this script only supports having genomeDirs "
                                          f"with one gtf files that ends with '.gtf'")
            else:
                genome_gtf_file = genome_gtf_file[0]
            call = f"featureCounts -L -T {self.threads} -R CORE -a {genome_gtf_file} " \
                   f"-o {self.featureCounts_dir}/{get_dt(for_file=True)}_featureCounts " \
                   f"{self.cat_files_dir}/cat.sorted.mappedAndPrimary.bam"
            self.run_cmd(call, "featureCounts", save_output_to_file=True)
            feature_counts_df = pd.read_table(feature_counts_file,
                                              header=None,
                                              names=('read_id', 'qc_tag', 'failed', 'gene_id'))
            feature_counts_df.gene_id.fillna(feature_counts_df.qc_tag, inplace=True)

            names_df = self.get_gene_ids_to_names_df()
            feature_counts_df = feature_counts_df.merge(names_df, on='gene_id', how='left')
            feature_counts_df.gene_name.fillna(feature_counts_df.gene_id, inplace=True)
            feature_counts_counts_df = feature_counts_df.gene_name \
                .value_counts() \
                .to_frame(name='gene_hits') \
                .reset_index(names='gene_name').merge(names_df,
                                                      on='gene_name',
                                                      how='left')
            feature_counts_counts_df = feature_counts_counts_df[['gene_id',
                                                                 'gene_name',
                                                                 'chr',
                                                                 'gene_hits']]
            print(f"Top 10 genes/assignments:")
            print(feature_counts_counts_df.head(10))
            self.logger.info(f"Saving featureCounts output to: "
                             f"{feature_counts_file.parent.name}/{feature_counts_file.name}.parquet")
            feature_counts_df.to_parquet(f"{feature_counts_file}.parquet")
            filter_call = f"grep Assigned {self.featureCounts_dir}/cat.sorted.mappedAndPrimary.bam.featureCounts " \
                          f">> {self.featureCounts_dir}/cat.sorted.mappedAndPrimary.bam.Assigned.featureCounts"
            self.run_cmd(filter_call, "featureCounts_filtering", save_output_to_file=False)
            self.regenerate = True
        else:
            self.logger.info(f"Skipping featureCounts because it has "
                             f"already been run on cat.sorted.mappedAndPrimary.bam.")
        self.featureCounts_ran = True

    @pipeline_step_decorator
    def merge_files(self):
        """
        This method will merge the results from mapping, featureCounts, and nanopolish/tailfindr files into one file.
        :return:
        """
        # Load the sam file and add the strand column:
        sam_df = SamOrBamFile(f"{self.cat_files_dir}/cat.sorted.mappedAndPrimary.sam").df
        sam_df["strand"] = (sam_df.bit_flag & 16).replace(to_replace={16: "-", 0: "+"})
        sam_df = sam_df.astype({'strand': 'category'})
        
        # Load the featureCounts file:
        try:
            featc_df = pd.read_parquet(
                f"{self.featureCounts_dir}/cat.sorted.mappedAndPrimary.bam.featureCounts.parquet")
        except FileNotFoundError:
            featc_df = pd.read_csv(f"{self.featureCounts_dir}/cat.sorted.mappedAndPrimary.bam.Assigned.featureCounts",
                                   sep="\t",
                                   names=["read_id",
                                          "qc_tag_featc",
                                          "qc_pass_featc",
                                          "gene_id"])
        if "gene_name" not in featc_df.columns:
            names_df = self.get_gene_ids_to_names_df()
            featc_df = featc_df.merge(names_df, on="gene_id", how="left")

        if self.tailfindr_ran:
            # Find the tailfindr file and load it into a dataframe:
            tailfinder_parquet_path = glob(f"{self.tailfindr_dir}/*.parquet")
            if len(tailfinder_parquet_path) != 1:
                raise NotImplementedError(f"Currently this script only supports having tailfindrDirs "
                                          f"with one parquet file.")
            else:
                tailfinder_parquet_path = tailfinder_parquet_path[0]
            tail_df = pd.read_parquet(tailfinder_parquet_path)
            # TODO: I have to massage this dataframe a bit to have it similar to nanopolish's output
        elif self.nanopolish_ran:
            tail_df = pd.read_csv(f"{self.nanopolish_dir}/polya.passed.tsv", sep="\t")
            # For some god-awful reason the chr_pos in polyA are -1 to those in the SAM file:
            tail_df["position"] = tail_df["position"] + 1
            tail_df = tail_df.rename(columns={"readname": "read_id",
                                              "qc_tag": "qc_tag_polya",  # b/c featC also has a qc_tag!
                                              "position": "chr_pos",
                                              "contig": "chr_id"})
        raise NotImplementedError

    def get_gene_ids_to_names_df(self):
        """
        This method will return a dataframe with the gene_id and gene_name columns.
        :return: 
        """
        # Find the file in the genome directory that ends with .gtf.parquet:
        gtf_parquet_file = glob(f"{self.genome_dir}/*.gtf.parquet")
        if len(gtf_parquet_file) > 1:
            raise NotImplementedError(f"Currently this script only supports having genomeDirs "
                                      f"with one gtf files that ends with '.gtf.parquet'")
        if len(gtf_parquet_file) == 0:
            self.logger.info(f"Could not find a pre-processed gtf.parquet file in {self.genome_dir}. Making one now.")
            gtf_file = glob(f"{self.genome_dir}/*.gtf")
            if len(gtf_file) != 1:
                raise NotImplementedError(f"Currently this script only supports having genomeDirs "
                                          f"with one gtf files that ends with '.gtf'")
            else:
                gtf_file = gtf_file[0]
            gtf_parquet_file = gtf_to_df(gtf_file)
            gtf_parquet_file.to_parquet(gtf_file + ".parquet")
            names_df = gtf_parquet_file[["gene_name", "gene_id", "chr"]].drop_duplicates(ignore_index=True)
        else:
            gtf_parquet_file = gtf_parquet_file[0]
            names_df = gene_names_to_gene_ids(parquet_path=gtf_parquet_file)
        return names_df

    @pipeline_step_decorator
    def assign_ENO2_standards(self):
        # Used ChatGPT to convert the old method to this one
        # TODO: Determine if this is the correct way to do this!!!
        
        # Check if the merge_files step has already been run
        merge_on_reads_path = self.merge_files_dir / f"{get_dt(for_file=True)}_mergedOnReads.plusStandards.parquet"

        if not self.regenerate and merge_on_reads_path.exists():
            self.logger.info("EN02 standards already assigned, skipping.")
            return

        # Import the necessary classes and modules for standards assignment
        from standardsAlignment.version2_mappingStandardsMethod_classBased import StandardsAlignerENO2

        # Determine the library type based on the guppyConfig
        if self.guppy_config.startswith('rna'):
            library_type = 'dRNA'
        elif self.guppy_config.startswith('dna'):
            library_type = 'cDNA'
        else:
            library_type = None  # Handle unknown library type here if needed

        # Load the merged data if it hasn't been loaded already
        if not hasattr(self, 'merged_df') or self.merged_df is None:
            merge_on_reads_path = find_newest_matching_file(f"{self.merge_files_dir}/*mergedOnReads.parquet")
            self.merged_df = pd.read_parquet(merge_on_reads_path)

        # Perform standards alignment and save the result
        aligner = StandardsAlignerENO2(mjv_compressed_df=self.merged_df.sort_values("chr_id"),
                                       library_type=library_type)
        output_df = aligner.run_alignments()
        out_file_path = self.merge_files_dir / f"{get_dt(for_file=True)}_mergedOnReads.plusStandards.parquet"
        output_df.to_parquet(out_file_path)

        self.logger.info("EN02 standards assigned successfully.")

    @pipeline_step_decorator
    def flair_for_transcripts(self):
        """
        This method will run the FLAIR pipeline on the cat.sorted.mappedAndPrimary.bam file.
        :return: 
        """
        raise NotImplementedError

    @pipeline_step_decorator
    def extra_steps(self):
        raise NotImplementedError

    def run_pipeline(self):
        """
        Runs the pipeline based on the steps_to_run parameter in the settings file
        """
        steps_list = list(self.steps_to_run)
        self.logger.info(f"Starting pipeline. Running steps: {', '.join(steps_list)}")
        # First lets sort the steps to run by their order:
        steps_dict_keys_ordered = sorted(self.steps_dict, key=lambda x: self.steps_dict[x].step_order_number)
        for step_key in steps_dict_keys_ordered:
            step = self.steps_dict[step_key]
            step_function = step.function
            self.logger.info(f"Running step {step_key}: {step_function.__name__}")
            if step_key in steps_list:
                step_function()


if __name__ == '__main__':
    if len(sys.argv) > 1:
        print("Running pipeline with command line arguments...")
        pipeline = NanoporePipeline()
        pipeline.run_pipeline()
    else:
        print("Running test pipeline...")
        pipeline = NanoporePipeline(skip_cli_dict={"settings_file": "/data16/marcus/working/211101_nanoporeSoftLinks"
                                                                    "/211118_nanoporeRun_totalRNA_5108_xrn-1-KD_5TERA"
                                                                    "/211118_totalRNA_5108_xrn-1-KD_settingsFile.txt",
                                                   # "do_not_log": True,
                                                   # "outputDir": f"/tmp/{get_dt()}_temporary_output_dir",
                                                   "stepsToRun": "P",
                                                   "useTailfindr": False,
                                                   "regenerate": True,
                                                   },
                                    )
        pipeline.run_pipeline()
