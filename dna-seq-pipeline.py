#!/usr/bin/env python3
# 0.2rc3 (Build20190401)
# Useage: python gene_panel.py gene_panel.ini
# Required files:
    # One or more pairs of fastq files with suffices _1P and _2P, all belong to different read groups of one same sample

# TODO:
#

import os, sys, glob, re, multiprocessing, shutil,common
import os.path as path
import numpy as np
import pandas as pd
from multiprocessing import Process, Lock, Manager
from common import readable, cprint

CONFIG = common.ini_to_config(sys.argv[1])

if CONFIG.get('PRIMARY', 'case') == '':
    CASE_LST = []
    CASE_AMOUNT_INT = 0
else:
    CASE_LST = re.sub(r'\s+', '', CONFIG.get('PRIMARY', 'case')).split('|')
    CASE_AMOUNT_INT = len(CASE_LST)

if CONFIG.get('PRIMARY', 'control') == '':
    CONTROL_LST = []
    CONTROL_AMOUNT_INT = 0
else:
    CONTROL_LST = re.sub(r'\s+', '', CONFIG.get('PRIMARY', 'control')).split('|')
    CONTROL_AMOUNT_INT = len(CONTROL_LST)

SAMPLE_LST = CASE_LST + CONTROL_LST
SAMPLE_ANOUMT_INT = CASE_AMOUNT_INT + CONTROL_AMOUNT_INT
COHORT_STR = CONFIG.get('PRIMARY', 'cohort_name')
PON_STR = CONFIG.get('SOMATICS', 'pon')
if PON_STR == '':
    PON_STR = COHORT_STR + '.vcf.gz'

TASK_LST = re.sub(r'\s+', '', CONFIG.get('PRIMARY', 'task')).split('|')

REF_FASTA = CONFIG.get('RESOURCES', 'ref_fasta')
CPU_INT = 0 # assign on the fly
BWA = CONFIG.get('TOOLS', 'bwa')
GATK = CONFIG.get('TOOLS', 'gatk')
PICARD = CONFIG.get('TOOLS', 'picard')
SAMTOOLS = CONFIG.get('TOOLS', 'samtools')
VARSCAN = CONFIG.get('TOOLS', 'varscan')
ANNOVAR = CONFIG.get('TOOLS', 'annovar')

CALLING_INTERVALS = CONFIG.get('PRIMARY', 'calling_intervals_file')
VERBOSE = CONFIG.getboolean('PRIMARY', 'verbose', fallback= True)
DBSNP = CONFIG.get('RESOURCES', 'dbsnp')
KNOWN_SITES = CONFIG.get('RESOURCES', 'known_indels_sites')


#===================================================================
########################################################################
class ExternalToolsError(Exception):
    """raise when an enternal tool (such as GATK) has a non-zero returncode"""
    #----------------------------------------------------------------------
    def __init__(self, command:str, message:str):
        """
        Constructor

        Parameters:
            **command**: string
                the external command leading to the error

            **message**: string
                the output from external command, could be from stderr or stdout
        """
        highlighted_error_str = 'An error occured when run the following command in process {pid}:\n\t{command}\n\n'.format(pid= os.getpid(), command= command)
        self.message = message
        tool_match_object = re.search(r'gatk|picard|samtools|bwa|varscan|annovar', command)

        if tool_match_object is None:  # we can't find out which tool you are using
            highlighted_error_str += message + '\n'
            cprint(highlighted_error_str)
            return None

        # we know which tool you are using so we can highlight the error message
        tool_str = tool_match_object.group(0)
        if tool_str == 'picard':
            for line_str in message.split('\n'):
                if line_str.split()[0] == 'Exception' or line_str.split()[0] == 'ERROR' or 'not a valid command' in line_str:
                    highlighted_error_str += line_str + '\n'
                    self.message = line_str

        elif tool_str == 'gatk':
            for line_str in message.split('\n'):
                if 'ERROR' in line_str:
                    highlighted_error_str += line_str + '\n'
                    self.message = line_str

        elif tool_str == 'bwa':
            for line_str in message.split('\n'):
                if 'fail' in line_str:
                    highlighted_error_str += line_str + '\n'
                    self.message = line_str

        elif tool_str == 'varscan':
            for line_str in message.split('\n'):
                if 'Exception' in line_str:
                    highlighted_error_str += line_str + '\n'
                    self.message = line_str

        elif tool_str == 'annovar':
            for line_str in message.split('\n'):
                if line_str.split()[0] == 'Error:':
                    highlighted_error_str += line_str + '\n'
                    self.message = line_str
                    break

        else:
            highlighted_error_str += message + '\n'

        cprint(highlighted_error_str)

        return None

def __get_free_mem() -> int:
    '''
    Return machine free memory in {integer} MB
    '''
    command_str = 'cat /proc/meminfo | grep "MemAvailable:"'
    r = common.run(command_str, False)
    if r[0] == 0:
        free_mem = int(int(r[1].split()[1])/1000)
        return free_mem
    else:
        raise OSError(r[2])


def __get_allele(symbol1:str, symbol2: str) -> str:
    '''
    input the Symbol1 [ATCG]{1} and IUPAC notation [WSMKRY]{1}, output the exact base

    e.g. __get_allele('A', 'W') will return 'T'
    '''

    if symbol2 in ['A', 'T', 'C', 'G']:
        return symbol2

    if symbol1 not in ['A', 'T', 'C', 'G']:
        return ''

    IUPAC_notation_dict = {'W':['A', 'T'],
                           'S':['C', 'G'],
                           'M':['A', 'C'],
                           'K':['G', 'T'],
                           'R':['A', 'G'],
                           'Y':['C', 'T']}

    try:
        return(IUPAC_notation_dict[symbol2][not IUPAC_notation_dict[symbol2].index(symbol1)])
    except:
        return ''

def __varscan_to_variant(*filenames: str) -> pd.DataFrame:
    '''input one or a few varscan output file(s), it will generate a pandas date frame

    Columns: Chr, Start, End, Ref, Alt, VarFreq
    '''

    results_df = pd.DataFrame(data=None, columns=['Chr', 'Start', 'End', 'Ref', 'Alt', 'VarFreq'])
    for s in filenames:
        if not isinstance(s, str):
            next

        filename_str = path.realpath(path.expanduser(s))
        if not readable(filename_str):
            next

        varscan_df = pd.read_csv(filename_str, sep= '\t')
        if len(varscan_df) == 0:
            next

        if sum( varscan_df['Cons'].str.contains(r'\+|-') ) / len(varscan_df) >= 0.5:  #  this is a INDEL file
            for indel_tu in varscan_df.itertuples():
                try:
                    if '+' in indel_tu.Cons:  # this is an insertion
                        chr_str = indel_tu.Chrom
                        start_int = int(indel_tu.Position) + 1
                        end_int = int(indel_tu.Position) + 1
                        ref_str = '-'
                        alt_str = indel_tu.Cons.split('+')[1]
                        ratio_str = indel_tu.VarFreq

                    elif '-' in indel_tu.Cons:  # this is an deletion
                        chr_str = indel_tu.Chrom
                        start_int = int(indel_tu.Position) + 1
                        ref_str = indel_tu.Cons.split('-')[1]
                        alt_str = '-'
                        end_int = start_int + len(ref_str) - 1
                        ratio_str = indel_tu.VarFreq

                    else:
                        next

                except Exception as ex:
                    message = 'Fail to recogonize deletion/inserion in {}:\n{}'.format(s, indel_tu)
                    print(ex)
                    next

                temp_df = pd.DataFrame({'Chr':chr_str, 'Start':start_int, 'End':end_int, 'Ref':ref_str, 'Alt':alt_str, 'VarFreq':ratio_str}, index=[0])
                results_df = pd.concat([results_df, temp_df], ignore_index= True, sort= False)

        else:    #  this is a SNP file
            for snp_tu in varscan_df.itertuples():
                try:
                    chr_str = snp_tu.Chrom
                    start_int = int(snp_tu.Position)
                    end_int = int(snp_tu.Position)
                    ref_str = snp_tu.Ref
                    alt_str = __get_allele(ref_str, snp_tu.Cons)
                    ratio_str = snp_tu.VarFreq
                except Exception as ex:
                    message = 'Fail to recogonize a SNV in {}:\n{}'.format(s, snp_tu)
                    print(ex)
                    next

                temp_df = pd.DataFrame({'Chr':chr_str, 'Start':start_int, 'End':end_int, 'Ref':ref_str, 'Alt':alt_str, 'VarFreq':ratio_str}, index=[0])
                results_df = pd.concat([results_df, temp_df], ignore_index= True, sort= False)


    results_df.drop_duplicates(inplace= True)
    results_df.sort_values(by=['Chr', 'Start', 'End'], inplace = True)
    results_df.index = range(len(results_df))
    try:
        results_df['Start'] = pd.to_numeric(results_df['Start'])
        results_df['End'] = pd.to_numeric(results_df['End'])
    except Exception as ex:
        print(ex)
    return results_df


def __calculate_contamination(tumor_name, normal_name= None):
    '''
    Split the [INTERVALS] calling_intervals_file into N partial interval files
    and put them into {COHORT_STR}/interval_files directory.
    two output files:
    {COHORT_STR}/{tumor_name}_{normal_name}.contamination.table
    {COHORT_STR}/{tumor_name}_{normal_name}.maf_segments.table

    Parameters:
        **tumor_name**: string
            tumor_name

    Returns:
        **value**: type
            its meaning
    '''
    import random
    print('\n==================gatk GetPileupSummaries=========================')
    command_mem = int( (__get_free_mem() - 500) / SAMPLE_ANOUMT_INT )
    contamination_vcf = CONFIG.get('RESOURCES', 'variants_for_contamination')
    tumor_bam = '{tumor_name}.analysis_ready.bam'.format(tumor_name= tumor_name)
    normal_bam = '{normal_name}.analysis_ready.bam'.format(normal_name= normal_name)

    interval = CONFIG.get('PRIMARY', 'calling_intervals_file')
    if interval != '' and readable(interval):
        interval_command = '--interval-set-rule INTERSECTION -L {}'.format(interval)
    else:
        interval_command = ''

    if readable(normal_bam):
        command_str = '{GATK} --java-options "-Xmx{command_mem}m" GetPileupSummaries -I {normal_bam} -V {contamination_vcf} {interval_command} -L {contamination_vcf} -O {COHORT_STR}/{normal_name}_pileups.table'.format(GATK= GATK,
                                                                                                                                                                                                                          command_mem= command_mem,
                                                                                                                                                                                                                          normal_bam= normal_bam,
                                                                                                                                                                                                                          interval_command= interval_command,
                                                                                                                                                                                                                          contamination_vcf= contamination_vcf,
                                                                                                                                                                                                                          COHORT_STR= COHORT_STR,
                                                                                                                                                                                                                          normal_name= normal_name)
        print(command_str)
        print()
        r = common.run(command_str, VERBOSE)
        if r[0] != 0 or not readable('{COHORT_STR}/{normal_name}_pileups.table'.format(COHORT_STR= COHORT_STR, normal_name= normal_name)):
            raise ExternalToolsError(command_str, r[1])


    command_str = '{GATK} --java-options "-Xmx{command_mem}m" GetPileupSummaries -I {tumor_bam} -V {contamination_vcf} {interval_command} -L {contamination_vcf} -O {COHORT_STR}/{tumor_name}_pileups.table'.format(GATK= GATK,
                                                                                                                                                                                                                    command_mem= command_mem,
                                                                                                                                                                                                                    tumor_bam= tumor_bam,
                                                                                                                                                                                                                    interval_command= interval_command,
                                                                                                                                                                                                                    contamination_vcf= contamination_vcf,
                                                                                                                                                                                                                    COHORT_STR= COHORT_STR,
                                                                                                                                                                                                                    tumor_name= tumor_name)
    print(command_str)
    print()
    r = common.run(command_str)
    if r[0] != 0 or not readable('{COHORT_STR}/{tumor_name}_pileups.table'.format(COHORT_STR= COHORT_STR, tumor_name= tumor_name)):
        raise ExternalToolsError(command_str, r[1])

    print('\n==================gatk calculate_contamination=========================')
    normal_pileup = '{COHORT_STR}/{normal_name}_pileups.table'.format( COHORT_STR= COHORT_STR, normal_name= normal_name )
    if readable(normal_pileup):
        normal_pileup_command = '-matched {}'.format(normal_pileup)
    else:
        normal_pileup_command = ''

    command_str =  '{GATK} --java-options "-Xmx{command_mem}m" CalculateContamination -I {COHORT_STR}/{tumor_name}_pileups.table -O {COHORT_STR}/{tumor_name}_{normal_name}.contamination.table --tumor-segmentation {COHORT_STR}/{tumor_name}_{normal_name}.maf_segments.table {normal_pileup_command}'.format(GATK= GATK,
                                                                                                                                                                                                                                                                                                                command_mem= command_mem,
                                                                                                                                                                                                                                                                                                                COHORT_STR= COHORT_STR,
                                                                                                                                                                                                                                                                                                                tumor_name= tumor_name,
                                                                                                                                                                                                                                                                                                                normal_name= normal_name,
                                                                                                                                                                                                                                                                                                                normal_pileup_command= normal_pileup_command)
    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)
    if r[0] != 0:
        raise ExternalToolsError(command_str, r[1])

    return



def __filter_by_orientation_bias(tumor_name, normal_name):
    '''
    Run this function if run_ob_filter is True
    output file : '{COHORT_STR}/{tumor_name}.somatic.ob_filtered.vcf.gz'

    '''
    command_mem = int( (__get_free_mem() - 1000) / SAMPLE_ANOUMT_INT )
    tumor_sequencing_artifact_metrics = CONFIG.get('RESOURCES', 'tumor_sequencing_artifact_metrics')
    tumor_bam = '{tumor_name}.analysis_ready.bam'.format(tumor_name= tumor_name)

    if tumor_sequencing_artifact_metrics == '' or not readable(tumor_sequencing_artifact_metrics): # we have to create our own artifact metrics
        print('\n==================gatk CollectSequencingArtifactMetrics=========================')
        command_str = '{GATK} --java-options "-Xmx{command_mem}m" CollectSequencingArtifactMetrics -I {tumor_name}.analysis_ready.bam -O {COHORT_STR}/{tumor_name} -R {REF_FASTA} -VALIDATION_STRINGENCY LENIENT'.format(GATK= GATK,
                                                                                                                                                                                                                         command_mem= command_mem,
                                                                                                                                                                                                                         tumor_name= tumor_name,
                                                                                                                                                                                                                         COHORT_STR= COHORT_STR,
                                                                                                                                                                                                                         REF_FASTA= REF_FASTA)
        print(command_str)
        print()
        r = common.run(command_str, VERBOSE)
        tumor_sequencing_artifact_metrics = '{COHORT_STR}/{tumor_name}.pre_adapter_detail_metrics'.format(COHORT_STR= COHORT_STR, tumor_name= tumor_name)
        if r[0] != 0 or not readable(tumor_sequencing_artifact_metrics):
            raise ExternalToolsError(command_str, r[1])

    print('\n==================gatk FilterByOrientationBias=========================')
    input_vcf = '{COHORT_STR}/{tumor_name}_{normal_name}.somatic.filtered.vcf.gz'.format(COHORT_STR= COHORT_STR,
                                                                                         tumor_name= tumor_name,
                                                                                         normal_name= normal_name)

    output_vcf = '{COHORT_STR}/{tumor_name}_{normal_name}.somatic.ob_filtered.vcf.gz'.format(COHORT_STR= COHORT_STR,
                                                                                             tumor_name= tumor_name,
                                                                                             normal_name= normal_name)

    artifact_modes_command = ''
    final_artifact_modes = CONFIG.get('SOMATICS', 'artifact_modes')
    if final_artifact_modes == '':
        final_artifact_modes = "'G/T'|'C/T'"
    final_artifact_modes = re.sub(r'\s+', '', final_artifact_modes)

    try:
        for mode in final_artifact_modes.split('|'):
            artifact_modes_command += ' --artifact-modes {}'.format(mode)
        artifact_modes_command = artifact_modes_command.strip()
    except:
        message = 'Can not parse customed artifact modes {}\nCheck configuration file {} [ARGS] section "artifact_modes" entry.'.format(final_artifact_modes, sys.argv[1])
        raise ValueError(message)

    command_str = '{GATK} --java-options "-Xmx{command_mem}m" FilterByOrientationBias -V {input_vcf} {artifact_modes_command} -P {tumor_sequencing_artifact_metrics} -O {output_vcf}'.format(GATK= GATK,
                                                                                                                                                                                             command_mem= command_mem,
                                                                                                                                                                                             input_vcf= input_vcf,
                                                                                                                                                                                             artifact_modes_command= artifact_modes_command,
                                                                                                                                                                                             tumor_sequencing_artifact_metrics= tumor_sequencing_artifact_metrics,
                                                                                                                                                                                             output_vcf= output_vcf)
    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)
    if r[0] != 0 or not readable(output_vcf):
        raise ExternalToolsError(command_str, r[1])

    return



def __filter_alignment_artifacts(tumor_name, normal_name, lock):
    '''

    '''
    command_mem = int( (__get_free_mem() - 500) / SAMPLE_ANOUMT_INT )

    input_vcf1 = '{COHORT_STR}/{tumor_name}_{normal_name}.somatic.ob_filtered.vcf.gz'.format(COHORT_STR= COHORT_STR, tumor_name= tumor_name, normal_name= normal_name) # FilterByOrientationBias output
    input_vcf2 = '{COHORT_STR}/{tumor_name}_{normal_name}.somatic.filtered.vcf.gz'.format(COHORT_STR= COHORT_STR, tumor_name= tumor_name, normal_name= normal_name)   # FilterMutectCalls output
    output_vcf = '{tumor_name}_{normal_name}.final.vcf.gz'.format(tumor_name= tumor_name, normal_name= normal_name)

    if readable(input_vcf1):
        realignment_filter_input = input_vcf1
    elif readable(input_vcf2):
        realignment_filter_input = input_vcf2
    else:
        message = 'Can not find either {} or {}. No input file for task __filter_alignment_artifacts'.format(input_vcf1, input_vcf2)
        raise FileNotFoundError(message)

    lock.acquire()
    realignment_index_bundle = CONFIG.get('RESOURCES', 'realignment_index_bundle')
    if realignment_index_bundle != '' and readable(realignment_index_bundle):
        pass
    else:   # fail to get valid realignment_index_bundle in configuration file.
        ref_image = path.splitext(REF_FASTA)[0] + '.img'
        if readable(ref_image):  # try to use ref image
            realignment_index_bundle = ref_image
        else:  # if not, create one
            print('\n==================gatk BwaMemIndexImageCreator=========================')
            command_str = '{GATK} --java-options "-Xmx{command_mem}m" BwaMemIndexImageCreator -I {REF_FASTA} -O {ref_image}'.format(GATK= GATK,
                                                                                                                                    command_mem= command_mem,
                                                                                                                                    REF_FASTA= REF_FASTA,
                                                                                                                                    ref_image= ref_image)
            print(command_str)
            print()
            r = common.run(command_str, VERBOSE)
            realignment_index_bundle = ref_image
            if r[0] != 0 or not readable(realignment_index_bundle):
                raise ExternalToolsError(command_str, r[1])
    lock.release()

    print('\n==================gatk FilterAlignmentArtifacts=========================')
    realignment_extra_args = CONFIG.get('SOMATICS', 'realignment_extra_args')
    tumor_bam = '{tumor_name}.analysis_ready.bam'.format(tumor_name= tumor_name)
    command_str = '{GATK} --java-options "-Xmx{command_mem}m" FilterAlignmentArtifacts -R {REF_FASTA} -V {realignment_filter_input} -I {tumor_bam} --bwa-mem-index-image {realignment_index_bundle} {realignment_extra_args} -O {output_vcf}'.format(GATK= GATK,
                                                                                                                                                                                                                                                     command_mem= command_mem,
                                                                                                                                                                                                                                                     REF_FASTA= REF_FASTA,
                                                                                                                                                                                                                                                     realignment_filter_input = realignment_filter_input,
                                                                                                                                                                                                                                                     tumor_bam= tumor_bam,
                                                                                                                                                                                                                                                     realignment_index_bundle= realignment_index_bundle,
                                                                                                                                                                                                                                                     realignment_extra_args= realignment_extra_args,
                                                                                                                                                                                                                                                     output_vcf= output_vcf)
    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)
    if r[0] != 0 or not readable(output_vcf):
        raise ExternalToolsError(command_str, r[1])


    return

def Check_ini() -> None:
    '''
    Check all parameters readed from ini file, make sure they are all valid
    '''
    global REF_FASTA

    if CASE_LST == [] and CONTROL_LST == []:
        message = 'No samples provided. Exit.'
        raise OSError(message)

    if CASE_AMOUNT_INT != 0 and CONTROL_AMOUNT_INT != 0 and CASE_AMOUNT_INT != CONTROL_AMOUNT_INT:
        message = 'Case samples {CASE_LST} and control samples {CONTROL_LST} must be equal in amount.'.format(CASE_LST= CASE_LST, CONTROL_LST= CONTROL_LST)
        raise OSError(message)

    if CASE_AMOUNT_INT == 0 and CONTROL_AMOUNT_INT < 2:
        message = 'It seems you only want to create a Pon (Panel of Normal), but too few normal samples.'
        raise OSError(message)

    if COHORT_STR == '':
        message = '[PRIMARY] cohort_name is not optional.'
        raise ValueError(message)

    if TASK_LST == []:
        message = '[PRIMARY] task is not optional.'
        raise ValueError(message)

    for t in TASK_LST:
        if t not in ['preprocessing', 'germline', 'somatic', 'TMB']:
            message = '[PRIMARY] task must be the combination of either "preprocessing", "TMB", "germline" or "somatic".\nCurrent value is {}'.format(TASK_LST)
            raise ValueError(message)

    if REF_FASTA == '' or not readable(REF_FASTA):
        message = 'Reference file {} not found. Check section [RESOURCES] ref_fasta'.format(REF_FASTA)
        raise FileNotFoundError(message)
    else:
        REF_FASTA = path.realpath(path.expanduser(REF_FASTA))

    if path.splitext(REF_FASTA)[1] != '.fa':
        message = '[RESOURCES]ref_fasta must has an extent name of ".fa", currently it is {}'.format(REF_FASTA)
        raise TypeError(message)

    if CALLING_INTERVALS == '' or not readable(CALLING_INTERVALS):
        message = '[PRIMARY] calling_intervals_file {} is empty or not readable. Exit'.format(CALLING_INTERVALS)
        raise FileNotFoundError(message)

    return None

def CreateSequenceGroupingTSV() -> None:

    ref_prefix = path.splitext(REF_FASTA)[0]
    ref_dict = path.realpath(path.expanduser(ref_prefix + '.dict' ))
    if not readable(ref_dict):
        message = '{} does not exist. Use picard CreateSequenceDictionary to create'.format(ref_dict)
        print(message)
        command_str = '{PICARD} CreateSequenceDictionary R={REF_FASTA} O={ref_dict}'.format(PICARD= PICARD,
                                                                                            REF_FASTA= REF_FASTA,
                                                                                            ref_dict= ref_dict)
        print(command_str)
        print()
        r = common.run(command_str, VERBOSE)

        if r[0] != 0 or not readable(ref_dict):
            raise ExternalToolsError(command_str, r[1])

    with open(path.realpath(path.expanduser(ref_dict)), 'r') as ref_dict_f:
        sequence_tuple_lst = []  # to store every reference name and length like [('chr1', 248956422), ('chr2', 242193529), ...]
        longest_length_int = 0
        for s in ref_dict_f:
            if s[:3] == '@SQ':
                line_split_lst = s.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_lst.append( (line_split_lst[1].split('SN:')[1], int(line_split_lst[2].split('LN:')[1])) )

    temp_lst = sequence_tuple_lst
    temp_lst.sort(key=lambda x:x[1], reverse=True)
    longest_length_int = temp_lst[0][1]

    # We are adding this to the intervals because hg38 has contigs named with embedded colons (:) and a bug in
    # some versions of GATK strips off the last element after a colon, so we add this as a sacrificial element.
    #hg38_protection_tag_str = ':1+'
    hg38_protection_tag_str = ''
    tsv_str = sequence_tuple_lst[0][0] + hg38_protection_tag_str # initialize the tsv string with the first sequence, like 'chr1:1+'
    previous_size_int = sequence_tuple_lst[0][1]
    for sequence_tuple in sequence_tuple_lst[1:]:
        if previous_size_int + sequence_tuple[1] <= longest_length_int:
            tsv_str += '\t' + sequence_tuple[0] + hg38_protection_tag_str
        else:
            tsv_str += '\n' + sequence_tuple[0] + hg38_protection_tag_str

        previous_size_int = sequence_tuple[1]

    # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
    with open('./{COHORT_STR}/sequence_grouping.txt'.format(COHORT_STR= COHORT_STR), 'w') as tsv_f:
        tsv_f.write(tsv_str)
        tsv_f.close()

    tsv_str += '\n' + "unmapped"

    with open('./{COHORT_STR}/sequence_grouping_with_unmapped.txt'.format(COHORT_STR= COHORT_STR), 'w') as tsv_file_with_unmapped_f:
        tsv_file_with_unmapped_f.write(tsv_str)
        tsv_file_with_unmapped_f.close()

    return None

def Preprocess(sample_name):
    '''
    GATK preprocessing workflow and generate {sample_name}.analysis_ready.bam and {sample_name}.analysis_ready.bam.bai

    Parameters:
        **sample_name**: string
            the prefix of two input fastq files, forward read: sample_1P, reverse read: sample_2P

    Returns:
        **value**: 0 on success
    '''

    cpu_per_command = CPU_INT / SAMPLE_ANOUMT_INT
    if cpu_per_command < 1:
        cpu_per_command = 1
    else:
        cpu_per_command = int(cpu_per_command)

    ref_prefix = path.splitext(REF_FASTA)[0]
    ref_dict = path.realpath(path.expanduser(ref_prefix + '.dict'))
    #=================================================

    print('\n=============================picard FastqToSam=============================')
    command_str = '{PICARD} FastqToSam FASTQ={sample_name}_1P FASTQ2={sample_name}_2P OUTPUT=./{COHORT_STR}/{sample_name}.unmapped.bam READ_GROUP_NAME={sample_name} SAMPLE_NAME={sample_name} PLATFORM=ILLUMINA'.format(PICARD= PICARD,
                                                                                                                                                                                                                         sample_name= sample_name,
                                                                                                                                                                                                                         COHORT_STR= COHORT_STR)
    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)

    output_file = './{COHORT_STR}/{sample_name}.unmapped.bam'.format(COHORT_STR= COHORT_STR, sample_name= sample_name)
    if r[0] != 0 or not readable(output_file):
        raise ExternalToolsError(command_str, r[1])

    print('\n=============================picard ValidateSamFile=============================')
    command_str = '{PICARD} ValidateSamFile I=./{COHORT_STR}/{sample_name}.unmapped.bam MODE=SUMMARY'.format(PICARD= PICARD,
                                                                                                             COHORT_STR= COHORT_STR,
                                                                                                             sample_name= sample_name)
    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)
    if r[0] != 0:
        raise ExternalToolsError(command_str, r[1])

    print('\n=============================picard SortSam=============================')
    command_str = '{PICARD} SortSam I=./{COHORT_STR}/{sample_name}.unmapped.bam O=./{COHORT_STR}/{sample_name}.sorted.unmapped.bam SORT_ORDER=queryname'.format(PICARD= PICARD,
                                                                                                                                                                COHORT_STR= COHORT_STR,
                                                                                                                                                                sample_name= sample_name)
    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)

    output_file = './{COHORT_STR}/{sample_name}.sorted.unmapped.bam'.format(COHORT_STR= COHORT_STR, sample_name= sample_name)
    if r[0] != 0 or not readable(output_file):
        raise ExternalToolsError(command_str, r[1])

    print('\n==================picard SamToFastq=========================')
    command_str = '{PICARD} SamToFastq INPUT=./{COHORT_STR}/{sample_name}.sorted.unmapped.bam FASTQ=./{COHORT_STR}/{sample_name}.fastq INTERLEAVE=true INCLUDE_NON_PF_READS=true'.format(PICARD= PICARD,
                                                                                                                                                                                         COHORT_STR= COHORT_STR,
                                                                                                                                                                                         sample_name= sample_name)
    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)

    output_file = './{COHORT_STR}/{sample_name}.fastq'.format(COHORT_STR= COHORT_STR, sample_name= sample_name)
    if r[0] != 0 or not readable(output_file):
        raise ExternalToolsError(command_str, r[1])

    print('\n==================bwa mem=========================')
    command_str = '{BWA} mem -K 100000000 -p -t {cpu_per_command} -Y {ref_prefix} ./{COHORT_STR}/{sample_name}.fastq > ./{COHORT_STR}/{sample_name}.sam'.format(BWA= BWA,
                                                                                                                                                                cpu_per_command= cpu_per_command,
                                                                                                                                                                ref_prefix= ref_prefix,
                                                                                                                                                                COHORT_STR= COHORT_STR,
                                                                                                                                                                sample_name= sample_name)
    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)
    if r[0] != 0:
        raise ExternalToolsError(command_str, r[1])

    print('\n==================samtools view=========================')
    command_str = '{SAMTOOLS} view -1 -@ {cpu_per_command} -o ./{COHORT_STR}/{sample_name}.bam ./{COHORT_STR}/{sample_name}.sam'.format(SAMTOOLS= SAMTOOLS,
                                                                                                                                        cpu_per_command= cpu_per_command,
                                                                                                                                        COHORT_STR= COHORT_STR,
                                                                                                                                        sample_name= sample_name)
    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)
    if r[0] != 0:
        raise ExternalToolsError(command_str, r[1])

    print('\n==================picard MergeBamAlignment=========================')
    if not readable(ref_dict):
        message = '{} does not exist. Use picard CreateSequenceDictionary to create'.format(ref_dict)
        print(message)
        command_str = '{PICARD} CreateSequenceDictionary R={REF_FASTA} O={ref_dict}'.format(PICARD= PICARD,
                                                                                            REF_FASTA= REF_FASTA,
                                                                                            ref_dict= ref_dict)
        print(command_str)
        print()
        r = common.run(command_str, VERBOSE)
        if r[0] != 0 or not readable(ref_dict):
            raise ExternalToolsError(command_str, r[1])

    # to get the version of installed bwa
    r = common.run(BWA, False)
    version = r[2].split('Version:')[1].split('\n')[0].strip() # '0.7.12-r1039'
    command_str ='{PICARD} MergeBamAlignment ALIGNED=./{COHORT_STR}/{sample_name}.bam UNMAPPED=./{COHORT_STR}/{sample_name}.sorted.unmapped.bam O=./{COHORT_STR}/{sample_name}_merge.bam R={REF_FASTA} VALIDATION_STRINGENCY=SILENT EXPECTED_ORIENTATIONS=FR ATTRIBUTES_TO_RETAIN=X0 SORT_ORDER=unsorted CLIP_ADAPTERS=false MAX_RECORDS_IN_RAM=2000000 MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant PROGRAM_RECORD_ID=bwamem PROGRAM_GROUP_VERSION="{version}" PROGRAM_GROUP_COMMAND_LINE="{BWA}" PROGRAM_GROUP_NAME=bwamem UNMAPPED_READ_STRATEGY=COPY_TO_TAG ALIGNER_PROPER_PAIR_FLAGS=true UNMAP_CONTAMINANT_READS=true'.format(PICARD= PICARD,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             COHORT_STR= COHORT_STR,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             sample_name= sample_name,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             REF_FASTA= REF_FASTA,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             version= version,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             BWA= BWA)
    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)

    output_file = './{COHORT_STR}/{sample_name}_merge.bam'.format(COHORT_STR= COHORT_STR, sample_name= sample_name)
    if r[0] != 0 or not readable(output_file):
        raise ExternalToolsError(command_str, r[1])

    print('\n==================picard SortSam=========================')
    command_str = '{PICARD} SortSam I=./{COHORT_STR}/{sample_name}_merge.bam O=./{COHORT_STR}/{sample_name}.sort.bam SORT_ORDER="coordinate" CREATE_INDEX=false CREATE_MD5_FILE=false'.format(PICARD= PICARD,
                                                                                                                                                                                              COHORT_STR= COHORT_STR,
                                                                                                                                                                                              sample_name= sample_name)
    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)
    output_file = './{COHORT_STR}/{sample_name}.sort.bam'.format(COHORT_STR= COHORT_STR, sample_name= sample_name)
    if r[0] != 0 or not readable(output_file):
        raise ExternalToolsError(command_str, r[1])

    print('\n==================picard SetNmMdAndUqTags=========================')
    command_str ='{PICARD} SetNmMdAndUqTags I=./{COHORT_STR}/{sample_name}.sort.bam O=./{COHORT_STR}/{sample_name}.set_tags.bam CREATE_INDEX=true R={REF_FASTA}'.format(PICARD= PICARD,
                                                                                                                                                                        COHORT_STR= COHORT_STR,
                                                                                                                                                                        sample_name= sample_name,
                                                                                                                                                                        REF_FASTA= REF_FASTA)
    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)
    output_file = './{COHORT_STR}/{sample_name}.set_tags.bam'.format(COHORT_STR= COHORT_STR, sample_name= sample_name)
    if r[0] != 0 or not readable(output_file):
        raise ExternalToolsError(command_str, r[1])


    is_mark_duplicates = CONFIG.getboolean('PREPROCESS', 'mark_duplicates', fallback= False)
    if is_mark_duplicates:
        print('\n==================picard MarkDuplicates=========================')

    barcode_tag = CONFIG.get('PREPROCESS', 'barcode_tag')
    read_one_barcode_tag = CONFIG.get('PREPROCESS', 'read_one_barcode_tag')
    read_two_barcode_tag = CONFIG.get('PREPROCESS', 'read_two_barcode_tag')
    remove_duplicates = CONFIG.getboolean('PREPROCESS', 'remove_duplicates', fallback= True)

    command_str = '{PICARD} MarkDuplicates TAG_DUPLICATE_SET_MEMBERS=true TAGGING_POLICY=All I=./{COHORT_STR}/{sample_name}.set_tags.bam O=./{COHORT_STR}/{sample_name}.removeDup.bam M=./{COHORT_STR}/{sample_name}.removeDup.txt'.format(COHORT_STR= COHORT_STR,
                                                                                                                                                                                                                                           PICARD= PICARD,
                                                                                                                                                                                                                                           sample_name= sample_name)
    if barcode_tag != '':
        command_str += ' BARCODE_TAG={}'.format(barcode_tag)

    if read_one_barcode_tag != '':
        command_str += ' READ_ONE_BARCODE_TAG={}'.format(read_one_barcode_tag)

    if read_two_barcode_tag != '':
        command_str += ' READ_TWO_BARCODE_TAG={}'.format(read_two_barcode_tag)

    if remove_duplicates:
        command_str += ' REMOVE_DUPLICATES=true'
    else:
        command_str += ' REMOVE_DUPLICATES=false'

    if is_mark_duplicates:
        print(command_str)
        print()
        r = common.run(command_str, VERBOSE)

        output_file = './{COHORT_STR}/{sample_name}.removeDup.bam'.format(COHORT_STR= COHORT_STR, sample_name= sample_name)
        if r[0] != 0 or not readable(output_file):
            raise ExternalToolsError(command_str, r[1])

    command_str = '{SAMTOOLS} sort -@ {cpu_per_command} -O bam -o ./{COHORT_STR}/{sample_name}.removeDup.sort.bam {output_file}'.format(SAMTOOLS= SAMTOOLS,
                                                                                                                                        cpu_per_command= cpu_per_command,
                                                                                                                                        COHORT_STR= COHORT_STR,
                                                                                                                                        sample_name= sample_name,
                                                                                                                                        output_file= output_file)
    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)
    output_file = './{COHORT_STR}/{sample_name}.removeDup.sort.bam'.format(COHORT_STR= COHORT_STR, sample_name= sample_name)
    if r[0] != 0 or not readable(output_file):
        raise ExternalToolsError(command_str, r[1])

    command_str = '{SAMTOOLS} index {output_file}'.format(SAMTOOLS= SAMTOOLS,
                                                          output_file= output_file)
    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)
    output_file = './{COHORT_STR}/{sample_name}.removeDup.sort.bam.bai'.format(COHORT_STR= COHORT_STR, sample_name= sample_name)
    if r[0] != 0 or not readable(output_file):
        raise ExternalToolsError(command_str, r[1])


    print('\n==================gatk BaseRecalibrator=========================')
    known_site_str = ''
    vcf_files_lst = re.sub(r'\s+', '', DBSNP).split('|') + re.sub(r'\s+', '', KNOWN_SITES).split('|')
    for vcf in vcf_files_lst:
        vcf_file = vcf.strip()
        if readable(vcf_file):
            if path.splitext(vcf_file)[1] == '.gz':
                index_file = vcf_file + '.tbi'
            elif path.splitext(vcf_file)[1] == '.vcf':
                index_file = vcf_file + '.idx'
            else:
                message = 'vcf file {} is not a .idx or .gz file. Skip'.format(vcf_file)
                print(message)
                continue

            if readable(index_file):
                pass
            else:
                message = '{} does not exist. Use gatk IndexFeatureFile to create'.format(index_file)
                print(message)
                command_str = '{gatk} IndexFeatureFile -F {vcf_file}'.format(gatk= gatk, vcf_file= vcf_file)
                print(command_str)
                print()
                r = common.run(command_str)
                if r[0] != 0 or not readable(index_file):
                    raise ExternalToolsError(command_str, r[1])

            known_site_str += ' --known-sites {}'.format(vcf_file)

        else:
            message = '{} not found. Skip.'.format(vcf_file)
            print(message)
            continue

    known_site_str = known_site_str.strip()
    if known_site_str == '':
        print(vcf_files_lst)
        message = 'No input known sites file. Exit'
        raise FileNotFoundError(message)

    if is_mark_duplicates:
        input_file = './{COHORT_STR}/{sample_name}.removeDup.sort.bam'.format(COHORT_STR= COHORT_STR, sample_name= sample_name)
    else:
        input_file = './{COHORT_STR}/{sample_name}.set_tags.bam'.format(COHORT_STR= COHORT_STR, sample_name= sample_name)

    proxy_dict = multiprocessing.Manager().dict()
    process_lst = []
    group_f = open('./{COHORT_STR}/sequence_grouping.txt'.format(COHORT_STR= COHORT_STR), 'r')
    for line_str in group_f.readlines():
        intervals_str = ''
        for chrom in line_str.strip().split('\t'):
            intervals_str += ' -L {}'.format(chrom)

        chrom = line_str.strip().split('\t')[0]
        intervals_str = intervals_str.strip()
        command_str = '{GATK} BaseRecalibrator -R {REF_FASTA} -I {input_file} --use-original-qualities -O ./{COHORT_STR}/{sample_name}.{chrom}.recal_data.tsv {known_site_str} {intervals_str}'.format(GATK= GATK,
                                                                                                                                                                                                       REF_FASTA= REF_FASTA,
                                                                                                                                                                                                       input_file= input_file,
                                                                                                                                                                                                       COHORT_STR= COHORT_STR,
                                                                                                                                                                                                       sample_name= sample_name,
                                                                                                                                                                                                       chrom= chrom,
                                                                                                                                                                                                       known_site_str= known_site_str,
                                                                                                                                                                                                       intervals_str= intervals_str)
        print(command_str)
        print()
        p = Process(target=common.run, args=(command_str, VERBOSE, proxy_dict, command_str, ) )
        process_lst.append(p)

    common.multi_run(process_lst, cpu_per_command)
    for command_str, r in proxy_dict.items():
        if r[0] != 0:
            raise ExternalToolsError(command_str, r[1])

    print('\n==================gatk GatherBqsrReports=========================')
    input_str = ''
    group_f = open('./{COHORT_STR}/sequence_grouping.txt'.format(COHORT_STR= COHORT_STR), 'r')
    for line_str in group_f.readlines():
        chrom = line_str.strip().split('\t')[0]
        recal_data_file = './{COHORT_STR}/{sample_name}.{chrom}.recal_data.tsv'.format(sample_name= sample_name,
                                                                                       COHORT_STR= COHORT_STR,
                                                                                       chrom= chrom)
        if readable(recal_data_file):
            input_str += ' -I {}'.format(recal_data_file)
        else:
            message = '{} not found. Skip'.format(recal_data_file)
            print(message)
            continue

    group_f.close()
    input_str = input_str.strip()
    command_str = '{GATK} GatherBQSRReports {input_str} -O ./{COHORT_STR}/{sample_name}.recal_data.tsv'.format(GATK= GATK,
                                                                                                               input_str= input_str,
                                                                                                               COHORT_STR= COHORT_STR,
                                                                                                               sample_name= sample_name)
    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)
    output_file = './{COHORT_STR}/{sample_name}.recal_data.tsv'.format(COHORT_STR= COHORT_STR, sample_name= sample_name)
    if r[0] != 0 or not readable(output_file):
        raise ExternalToolsError(command_str, r[1])

    print('\n==================gatk ApplyBQSR=========================')
    proxy_dict = multiprocessing.Manager().dict()
    process_lst = []
    group_f = open('./{COHORT_STR}/sequence_grouping_with_unmapped.txt'.format(COHORT_STR= COHORT_STR), 'r')
    for line_str in group_f.readlines():
        intervals_str = ''
        for chrom in line_str.strip().split('\t'):
            intervals_str += ' -L {}'.format(chrom)

        chrom = line_str.strip().split('\t')[0]
        intervals_str = intervals_str.strip()
        command_str = '{GATK} ApplyBQSR -R {REF_FASTA} -I ./{COHORT_STR}/{sample_name}.set_tags.bam -O ./{COHORT_STR}/{sample_name}_{chrom}.aligned.duplicates_marked.recalibrated.bam {intervals_str} -bqsr ./{COHORT_STR}/{sample_name}.recal_data.tsv --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 --add-output-sam-program-record --create-output-bam-md5 --use-original-qualities'.format(GATK= GATK,
                                                                                                                                                                                                                                                                                                                                                                                                                                      REF_FASTA= REF_FASTA,
                                                                                                                                                                                                                                                                                                                                                                                                                                      COHORT_STR= COHORT_STR,
                                                                                                                                                                                                                                                                                                                                                                                                                                      sample_name= sample_name,
                                                                                                                                                                                                                                                                                                                                                                                                                                      chrom= chrom,
                                                                                                                                                                                                                                                                                                                                                                                                                                      intervals_str= intervals_str)
        print(command_str)
        p = Process(target=common.run, args=(command_str, VERBOSE, proxy_dict, command_str,) )
        process_lst.append(p)

    group_f.close()
    common.multi_run(process_lst, cpu_per_command)
    for command_str, r in proxy_dict.items():
        if r[0] != 0:
            raise ExternalToolsError(command_str, r[1])

    print('\n==================picard GatherBamFiles=========================')
    input_str = ''
    group_f = open('./{COHORT_STR}/sequence_grouping_with_unmapped.txt'.format(COHORT_STR= COHORT_STR), 'r')
    for line_str in group_f.readlines():
        chrom = line_str.strip().split('\t')[0]
        part_bam_file = './{COHORT_STR}/{sample_name}_{chrom}.aligned.duplicates_marked.recalibrated.bam'.format(COHORT_STR= COHORT_STR,
                                                                                                                 sample_name= sample_name,
                                                                                                                 chrom= chrom )
        if readable(part_bam_file):
            input_str += ' I={}'.format(part_bam_file)
        else:
            message = '{} not found. Skip'.format(part_bam_file)
            print(message)
            continue

    group_f.close()
    input_str = input_str.strip()
    if input_str == '':
        message = 'Can not find any part_bam_file. Exit'
        raise FileNotFoundError(message)

    command_str = '{PICARD} GatherBamFiles {input_str} O=./{COHORT_STR}/{sample_name}.gathered.bam CREATE_INDEX=false CREATE_MD5_FILE=false'.format(PICARD= PICARD,
                                                                                                                                                    input_str= input_str,
                                                                                                                                                    COHORT_STR= COHORT_STR,
                                                                                                                                                    sample_name= sample_name)
    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)
    output_file = './{COHORT_STR}/{sample_name}.gathered.bam'.format(COHORT_STR= COHORT_STR, sample_name= sample_name)
    if r[0] != 0 or not readable(output_file):
        raise ExternalToolsError(command_str, r[1])

    print('\n==================samtools Sort=========================')
    command_str = '{SAMTOOLS} sort -@ {cpu_per_command} -O bam -o {sample_name}.analysis_ready.bam ./{COHORT_STR}/{sample_name}.gathered.bam'.format(SAMTOOLS= SAMTOOLS,
                                                                                                                                                     cpu_per_command= cpu_per_command,
                                                                                                                                                     sample_name= sample_name,
                                                                                                                                                     COHORT_STR= COHORT_STR)
    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)
    output_file = '{sample_name}.analysis_ready.bam'.format(COHORT_STR= COHORT_STR, sample_name= sample_name)
    if r[0] != 0 or not readable(output_file):
        raise ExternalToolsError(command_str, r[1])

    command_str = '{SAMTOOLS} index {sample_name}.analysis_ready.bam'.format(SAMTOOLS= SAMTOOLS,
                                                                             sample_name= sample_name)
    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)
    output_file = '{sample_name}.analysis_ready.bam.bai'.format(COHORT_STR= COHORT_STR, sample_name= sample_name)
    if r[0] != 0 or not readable(output_file):
        raise ExternalToolsError(command_str, r[1])

    print('\n==================picard ValidateSamFile=========================')
    command_str = '{PICARD} ValidateSamFile I={sample_name}.analysis_ready.bam MODE=SUMMARY'.format(PICARD= PICARD,
                                                                                                    sample_name= sample_name)
    print(command_str)
    r = common.run(command_str)
    if r[0] != 0:
        raise ExternalToolsError(command_str, r[1])

    return 0

def CallGermlineVarscan(sample_name):
    '''
    1, Call variants by Varscan from {sample_name}.analysis_ready.bam and {sample_name}.analysis_ready.bam.bai
    2, Concatenate two calling results and write them to ./{COHORT_STR}/{sample_name}.results.txt
    3, Merge the calling results and annotation file and write the annotation into a standard SNPINDEL file {sample_name}.annotated.txt

    Step 1, use varscan to call SNP and INDEL, it will output two files:
    1. {COHORT_STR}/{sample_name}.varscan.snp
    headers:
    Chrom	Position	Ref	Cons	Reads1	Reads2	VarFreq	Strands1	Strands2	Qual1	Qual2	Pvalue	MapQual1	MapQual2	Reads1Plus	Reads1Minus	Reads2Plus	Reads2Minus	VarAllele

    2. {COHORT_STR}/{sample_name}.varscan.indel
    headers:
    Chrom	Position	Ref	Cons	Reads1	Reads2	VarFreq	Strands1	Strands2	Qual1	Qual2	Pvalue	MapQual1	MapQual2	Reads1Plus	Reads1Minus	Reads2Plus	Reads2Minus	VarAllele

    Step 2, merge two results and re-format them into a standard SNPINDEL file {sample_name}.results.txt .
    header-less:
    chr    start    end     ref    alt    (optional columns)

    Step 3, add annotation information in optional columns (gene name, variant name, ratio, clinical meaning... etc.)

    Parameters:
        **sample_name**: string
            the prefix of two input files. {sample_name}.analysis_ready.bam and {sample_name}.analysis_ready.bam.bai must be existed

    Returns:
        **value**: 0 - success
    '''
    import shutil

    print('\n==================varscan=========================')
    # create a BED file (./{COHORT_STR}/{prefix}.bed) from calling_intervals_list
    prefix = path.splitext(path.split(CALLING_INTERVALS)[1])[0]
    output_file = './{COHORT_STR}/{prefix}.bed'.format(COHORT_STR= COHORT_STR, prefix= prefix)

    if path.splitext(CALLING_INTERVALS)[1].upper() == '.BED':
        shutil.copyfile(path.realpath(path.expanduser(CALLING_INTERVALS)), output_file )

    elif path.splitext(CALLING_INTERVALS)[1].lower() == '.interval_list': # picard style
        command_str = '{PICARD} IntervalListToBed I={CALLING_INTERVALS} O={output_file} SORT=true'.format(PICARD= PICARD,
                                                                                                          CALLING_INTERVALS= CALLING_INTERVALS,
                                                                                                          output_file= output_file)
        print(command_str)
        print()
        r = common.run(command_str, VERBOSE)
        if r[0] != 0 or not readable(output_file):
            raise ExternalToolsError(command_str, r[1])

    truncate_depth_int = 1000000  # after test it is assumed that truncate_depth_int should be as large as possible, further test is still welcomed
    min_BQ = CONFIG.getint('GERMLINE', 'min_BQ', fallback= 13)
    min_MQ = CONFIG.getint('GERMLINE', 'min_MQ', fallback= 0)
    mpileup_extra_args = CONFIG.get('GERMLINE', 'mpileup_extra_args', fallback= '')
    command_str = '{SAMTOOLS} mpileup -A --min-BQ {min_BQ} --min-MQ {min_MQ} --max-depth {truncate_depth_int} --positions ./{COHORT_STR}/{prefix}.bed --fasta-ref {REF_FASTA} {mpileup_extra_args} {sample_name}.analysis_ready.bam > {COHORT_STR}/{sample_name}.mpileup'.format(SAMTOOLS= SAMTOOLS,
                                                                                                                                                                                                                                                                                 min_BQ = min_BQ,
                                                                                                                                                                                                                                                                                 min_MQ = min_MQ,
                                                                                                                                                                                                                                                                                 truncate_depth_int= truncate_depth_int,
                                                                                                                                                                                                                                                                                 COHORT_STR= COHORT_STR,
                                                                                                                                                                                                                                                                                 prefix= prefix,
                                                                                                                                                                                                                                                                                 REF_FASTA= REF_FASTA,
                                                                                                                                                                                                                                                                                 mpileup_extra_args = mpileup_extra_args,
                                                                                                                                                                                                                                                                                 sample_name= sample_name)
    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)
    output_file = '{COHORT_STR}/{sample_name}.mpileup'.format(COHORT_STR= COHORT_STR, sample_name= sample_name)
    if r[0] != 0 or not readable(output_file):
        raise ExternalToolsError(command_str, r[1])

    if path.getsize(output_file) == 0:
        message = 'Size of {output_file} is 0. Check BED file {CALLING_INTERVALS}'.format(output_file= output_file, CALLING_INTERVALS= CALLING_INTERVALS)
        cprint(message)
        sys.exit(1)


    min_coverage = CONFIG.getint('GERMLINE', 'min_coverage', fallback= 8)
    min_reads2 = CONFIG.getint('GERMLINE', 'min_read2', fallback= 2)
    min_avg_qual = CONFIG.getint('GERMLINE', 'min_avg_qual', fallback= 15)
    min_var_freq = CONFIG.getfloat('GERMLINE', 'min_var_freq', fallback= 0.01)
    p_value = CONFIG.getfloat('GERMLINE', 'p_value', fallback= 0.99)
    command_str = '{VARSCAN} pileup2snp {COHORT_STR}/{sample_name}.mpileup --min-coverage {min_coverage} --min-reads2 {min_reads2} --min-avg-qual {min_avg_qual} --min-var-freq {min_var_freq} --p-value {p_value} | tee {COHORT_STR}/{sample_name}.varscan.snp'.format(VARSCAN= VARSCAN,
                                                                                                                                                                                                                                                                        COHORT_STR= COHORT_STR,
                                                                                                                                                                                                                                                                        sample_name= sample_name,
                                                                                                                                                                                                                                                                        min_coverage= min_coverage,
                                                                                                                                                                                                                                                                        min_reads2= min_reads2,
                                                                                                                                                                                                                                                                        min_avg_qual= min_avg_qual,
                                                                                                                                                                                                                                                                        min_var_freq= min_var_freq,
                                                                                                                                                                                                                                                                        p_value= p_value)
    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)
    if r[0] == 127:
        message = '[Exception] Program tee not installed.'
        raise ExternalToolsError(command_str, message)
    if 'Exception' in r[1]:
        raise ExternalToolsError(command_str, r[1])

    command_str = '{VARSCAN} pileup2indel {COHORT_STR}/{sample_name}.mpileup --min-coverage {min_coverage} --min-reads2 {min_reads2} --min-avg-qual {min_avg_qual} --min-var-freq {min_var_freq} --p-value {p_value} | tee {COHORT_STR}/{sample_name}.varscan.indel'.format(VARSCAN= VARSCAN,
                                                                                                                                                                                                                                                                            COHORT_STR= COHORT_STR,
                                                                                                                                                                                                                                                                            sample_name= sample_name,
                                                                                                                                                                                                                                                                            min_coverage= min_coverage,
                                                                                                                                                                                                                                                                            min_reads2= min_reads2,
                                                                                                                                                                                                                                                                            min_avg_qual= min_avg_qual,
                                                                                                                                                                                                                                                                            min_var_freq= min_var_freq,
                                                                                                                                                                                                                                                                            p_value= p_value)
    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)
    if r[0] == 127:
        message = '[Exception] Program tee not installed.'
        raise ExternalToolsError(command_str, message)
    if 'Exception' in r[1]:
        raise ExternalToolsError(command_str, r[1])

    varscan_snp_file = '{COHORT_STR}/{sample_name}.varscan.snp'.format(COHORT_STR= COHORT_STR,sample_name= sample_name)
    varscan_indel_file = '{COHORT_STR}/{sample_name}.varscan.indel'.format(COHORT_STR= COHORT_STR,sample_name= sample_name)

    # Step 2, merge two results and re-format them into a standard SNPINDEL file ./{COHORT_STR}/{sample_name}.results.txt
    # header:
    # Chr    Start    End     Ref    Alt    VarFreq (optional columns)
    results_df = __varscan_to_variant(varscan_snp_file, varscan_indel_file)

    # Finish concatenating the results
    if len(results_df) >= 0:
        output_file = './{COHORT_STR}/{sample_name}.results.txt'.format(COHORT_STR= COHORT_STR, sample_name= sample_name)
        print('{} SNPs/INDELs found.\nWrite to {}'.format(len(results_df), output_file))
        results_df.to_csv(output_file, sep= '\t', na_rep = 'NA', index= False)
    else:
        message = 'Varscan failed to call variants.\nCheck {} and {} varscan parameters.'.format(varscan_snp_file, varscan_indel_file)
        raise ExternalToolsError('[Annotate varscan results]', message)

    # merge the calling results and annotation info
    varscan_annotation_db = CONFIG.get('RESOURCES', 'varscan_annotation_db')
    if readable(varscan_annotation_db):
        #annotation_df = pd.read_csv(varscan_annotation_db, sep= '\t', header=None, names=['Chr', 'Start', 'End', 'Ref', 'Alt', 'Gene', 'Name'], dtype={'chr':str, 'start':int, 'end':int, 'ref':str, 'alt':str, 'gene':str, 'name':str})
        annotation_df = pd.read_csv(varscan_annotation_db, sep= '\t')
    else:
        message = '{} is not readable. Annotation failed.'.format(varscan_annotation_db)
        cprint(message)
        sys.exit(1)

    try:
        results_df['Start'] = pd.to_numeric(results_df['Start'])
        results_df['End'] = pd.to_numeric(results_df['End'])
    except:
        annotation_df['Start'] = annotation_df['Start'].astype('object')
        annotation_df['End'] = annotation_df['End'].astype('object')

    merged_df = pd.merge(annotation_df, results_df, how='inner', on = ['Chr', 'Start', 'End', 'Ref', 'Alt'])
    merged_df = merged_df.iloc[:, range(5, len(merged_df.columns))]
    if len(merged_df) >= 2:
        merged_df.sort_values(by= 'VarFreq', ascending= False, inplace = True)
    output_file = '{sample_name}.annotated.txt'.format(sample_name= sample_name)
    merged_df.to_csv(output_file, sep= '\t', na_rep = 'NA', index= False)
    print('{} SNPs/INDELs annotated. Write the results in {}'.format(len(merged_df), output_file))

    return 0

def calc_TMB(sample_name:str) -> float:
    '''
    calculate the tumor mutation burden
    '''
    print('\n==================Calculate TMB=========================')
    min_MQ = CONFIG.getint('TMB', 'min_MQ', fallback= 20)
    mpileup_extra_args = CONFIG.get('TMB', 'mpileup_extra_args', fallback= '')
    command_str = '{SAMTOOLS} mpileup --max-depth 1000000 --min-MQ {min_MQ} --fasta-ref {REF_FASTA} -l {CALLING_INTERVALS} -o {COHORT_STR}/{sample_name}.tmb.mpileup {sample_name}.analysis_ready.bam'.format(SAMTOOLS= SAMTOOLS,
                                                                                                                                                                                                              sample_name= sample_name,
                                                                                                                                                                                                              min_MQ = min_MQ,
                                                                                                                                                                                                              REF_FASTA = REF_FASTA,
                                                                                                                                                                                                              CALLING_INTERVALS = CALLING_INTERVALS,
                                                                                                                                                                                                              COHORT_STR = COHORT_STR
                                                                                                                                                                                                              )
    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)


    min_var_freq = CONFIG.getfloat('TMB', 'min_var_freq', fallback= 0.01)
    min_reads2 = CONFIG.getint('TMB', 'min_read2', fallback= 5)
    min_coverage = CONFIG.getint('GERMLINE', 'min_coverage', fallback= 50)
    command_str = '{VARSCAN} pileup2snp {COHORT_STR}/{sample_name}.tmb.mpileup --min-var-freq {min_var_freq} --min-reads2 {min_reads2} --min-coverage {min_coverage} > {COHORT_STR}/{sample_name}.tmb.snp'.format(VARSCAN = VARSCAN,
                                                                                                                                                                                                              COHORT_STR = COHORT_STR,
                                                                                                                                                                                                              sample_name = sample_name,
                                                                                                                                                                                                              min_var_freq = min_var_freq,
                                                                                                                                                                                                              min_reads2 = min_reads2,
                                                                                                                                                                                                              min_coverage = min_coverage)
    r = common.run(command_str, VERBOSE)
    results_df = __varscan_to_variant('{COHORT_STR}/{sample_name}.tmb.snp'.format(COHORT_STR = COHORT_STR, sample_name = sample_name))
    print('{sample_name} has {length} SNPs'.format(sample_name = sample_name, length = len(results_df)))

    if len(results_df) == 0:
        return -1

    annotation_file_lst = re.sub(r'\s+', '', CONFIG.get('RESOURCES', 'tmb_discard_db')).split('|')
    annotate_lst = [ pd.read_csv(file, sep='\t') for file in annotation_file_lst ]
    for annotate_df in annotate_lst:  # discard the variants already in annotation database
        temp_df = pd.merge(annotate_df, results_df, how='inner', on = ['Chr', 'Start', 'End', 'Ref', 'Alt'], left_index = True)
        results_df = results_df.drop( index = temp_df.index )
        print('Drop {length} SNPs'.format(length = len(temp_df)))

    print('{sample_name} left {length} SNPs'.format(sample_name = sample_name, length = len(results_df)))
    if len(results_df) == 0:
        return -1

    command_str = '{SAMTOOLS} depth -b {CALLING_INTERVALS} -d 9999999 -Q {min_MQ} {sample_name}.analysis_ready.bam > {COHORT_STR}/{sample_name}.tmb.depth'.format(SAMTOOLS = SAMTOOLS,
                                                                                                                                                                  CALLING_INTERVALS = CALLING_INTERVALS,
                                                                                                                                                                  min_MQ = min_MQ,
                                                                                                                                                                  sample_name = sample_name,
                                                                                                                                                                  COHORT_STR = COHORT_STR)
    r = common.run(command_str)

    depth_file = '{COHORT_STR}/{sample_name}.tmb.depth'.format(COHORT_STR  = COHORT_STR, sample_name = sample_name)
    if not readable(depth_file):
        cprint('{depth_file} not readable.'.format(depth_file = depth_file))
        return -1

    length_int = 0
    with open(depth_file) as file_f:
        for line_str in file_f.readlines():
            try:
                if int(line_str.strip().split()[2]) >= 50:
                    length_int += 1
                    if length_int % 100000 == 0:
                        print(length_int, end= '\r')
            except Exception as ex:
                print(ex)
                next
        print(length_int)

    try:
        tmb_float = len(results_df) / length_int * 1000000
        cprint("{sample_name}'s TMB equals to {count}/{length}*1000000 = {tmb_float}".format(sample_name = sample_name,
                                                                                             length = length_int,
                                                                                             count = len(results_df),
                                                                                             tmb_float = round(tmb_float,2)))
    except ZeroDivisionError:
        cprint('{sample_name} sequencing deepth is too low. Skip'.format(sample_name = sample_name))
        return -1

    return tmb_float

def CallGermlineHaplotypeCaller(sample_name):
    '''
    not implemented yet
    '''

    return


def SplitIntervals(N: int) -> None:
    '''
    Split the [INTERVALS] calling_intervals_file into N partial interval files
    and put them into {COHORT_STR}/interval_files directory.

    Parameters:
        **N**: integer
            how many parts do you want to split ?

    Returns:
        **value**: type
            its meaning
    '''

    if isinstance(N, int):
        if N <= 0:
            raise ValueError
    else:
        raise TypeError

    print('\n==================gatk SplitIntervals=========================')
    command_mem = __get_free_mem() - 500
    split_intervals_extra_args = CONFIG.get('SOMATICS', 'split_intervals_extra_args') # optional

    if not path.isdir('{COHORT_STR}/interval_files'.format(COHORT_STR= COHORT_STR)):
        os.mkdir('{COHORT_STR}/interval_files'.format(COHORT_STR= COHORT_STR))

    command_str = '{GATK} --java-options "-Xmx{command_mem}m" SplitIntervals -R {REF_FASTA} -L {CALLING_INTERVALS} -scatter {N} -O {COHORT_STR}/interval_files {split_intervals_extra_args}'.format(GATK= GATK,
                                                                                                                                                                                                    command_mem= command_mem,
                                                                                                                                                                                                    REF_FASTA= REF_FASTA,
                                                                                                                                                                                                    CALLING_INTERVALS= CALLING_INTERVALS,
                                                                                                                                                                                                    N= N,
                                                                                                                                                                                                    COHORT_STR= COHORT_STR,
                                                                                                                                                                                                    split_intervals_extra_args= split_intervals_extra_args)

    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)
    print('r[0] = ', r[0])

    return None

def CallSomaticMutect2(tumor_name, normal_name= None):
    '''
    Call somatic mutation by GATK Mutect2.

    Parameters:
        **tumor_name**: string
            {tumor_name}.analysis_ready.bam is the initial input file

        **normal_name**: string, optional
            {normal_name}.analysis_ready.bam is the initial input file

    Returns:
        **value**: type
            its meaning
    '''

    #===============================LearnReadOrientationModel===============================
    proxy_dict = multiprocessing.Manager().dict()
    command_mem = int( (__get_free_mem() - 1000) / CPU_INT )
    process_1_lst = []
    process_2_lst = []
    for interval in glob.glob('{COHORT_STR}/interval_files/*.intervals'.format(COHORT_STR= COHORT_STR) ):
        num_suffix = path.split(interval)[1].split('-')[0]
        command1_str = '{GATK} --java-options "-Xmx{command_mem}m" CollectF1R2Counts -I {tumor_name}.analysis_ready.bam -R {REF_FASTA} -L {interval} -alt-table {COHORT_STR}/{tumor_name}_alt_{num_suffix}.tsv -ref-hist {COHORT_STR}/{tumor_name}_ref_{num_suffix}.metrics -alt-hist {COHORT_STR}/{tumor_name}_{num_suffix}.metrics'.format(GATK= GATK,
                                                                                                                                                                                                                                                                                                                                             command_mem= command_mem,
                                                                                                                                                                                                                                                                                                                                             tumor_name= tumor_name,
                                                                                                                                                                                                                                                                                                                                             REF_FASTA= REF_FASTA,
                                                                                                                                                                                                                                                                                                                                             interval= interval,
                                                                                                                                                                                                                                                                                                                                             COHORT_STR= COHORT_STR,
                                                                                                                                                                                                                                                                                                                                             num_suffix= num_suffix)
        print(command1_str)
        print()
        p = Process(target= common.run, args=(command1_str, VERBOSE, proxy_dict, ))
        process_1_lst.append(p)

        command2_str = '{GATK} --java-options "-Xmx{command_mem}m" LearnReadOrientationModel -alt-table {COHORT_STR}/{tumor_name}_alt_{num_suffix}.tsv -ref-hist {COHORT_STR}/{tumor_name}_ref_{num_suffix}.metrics -alt-hist {COHORT_STR}/{tumor_name}_{num_suffix}.metrics -O {COHORT_STR}/{tumor_name}_{num_suffix}_artifact_prior_table.tsv'.format(GATK= GATK,
                                                                                                                                                                                                                                                                                                                                             command_mem= command_mem,
                                                                                                                                                                                                                                                                                                                                             COHORT_STR= COHORT_STR,
                                                                                                                                                                                                                                                                                                                                             tumor_name= tumor_name,
                                                                                                                                                                                                                                                                                                                                             num_suffix= num_suffix)
        print(command2_str)
        print()
        p = Process(target= common.run, args=(command2_str, VERBOSE, proxy_dict, ))
        process_2_lst.append(p)

    run_ob_mm_filter = CONFIG.getboolean('SOMATICS', 'run_ob_mm_filter', fallback= True)
    if run_ob_mm_filter:

        print('\n==================gatk CollectF1R2Counts=========================')
        common.multi_run(process_1_lst, CPU_INT)
        for command_str, r in proxy_dict.items():
            if r[0] != 0:
                raise ExternalToolsError(command_str, r[1])

        print('\n==================gatk LearnReadOrientationModel=========================')
        common.multi_run(process_2_lst, CPU_INT)
        for command_str, r in proxy_dict.items():
            if r[0] != 0:
                raise ExternalToolsError(command_str, r[1])

    #===============================Mutect2===============================
    print('\n==================gatk Mutect2=========================')
    proxy_dict = multiprocessing.Manager().dict()
    command_mem = int( (__get_free_mem() - 500) / CPU_INT )
    tumor_command = '-I {tumor_name}.analysis_ready.bam -tumor {tumor_name}'.format(tumor_name= tumor_name)

    if normal_name is None:
        normal_command = ''
    else:
        normal_command = '-I {normal_name}.analysis_ready.bam -normal {normal_name}'.format(normal_name= normal_name)

    if readable(PON_STR):
        pon_command = '--panel-of-normals {PON_STR}'.format(PON_STR= PON_STR)
    else:
        print('PON resource {} not found. Leave --panel-of-normals blank'.format(PON_STR) )
        pon_command = ''

    gnomad = CONFIG.get('RESOURCES', 'germline_resource')
    if gnomad != '' and readable(gnomad): # if germline_resource vcf is available
        gnoman_command = '--germline-resource {}'.format(gnomad)
    else:
        print('Germline resource {} not found. Leave --germline-resource blank'.format(gnomad) )
        gnoman_command = ''

    genotype_given_alleles = CONFIG.get('SOMATICS', 'genotype_given_alleles')
    if genotype_given_alleles != '' and readable(genotype_given_alleles):
        gga_command = '--genotyping-mode GENOTYPE_GIVEN_ALLELES --alleles {}'.format(genotype_given_alleles)
    else:
        print('Genotype alleles resource {} not found. Leave --alleles blank'.format(genotype_given_alleles) )
        gga_command = ''

    m2_extra_args_command = CONFIG.get('SOMATICS', 'm2_extra_args')

    process_lst = []
    for interval in glob.glob('{COHORT_STR}/interval_files/*.intervals'.format(COHORT_STR= COHORT_STR) ):
        num_suffix = path.split(interval)[1].split('-')[0]
        output_vcf = '{COHORT_STR}/{tumor_name}_{normal_name}_{num_suffix}.vcf.gz'.format(COHORT_STR= COHORT_STR,
                                                                                          tumor_name= tumor_name,
                                                                                          normal_name= normal_name,
                                                                                          num_suffix= num_suffix)

        make_bamout = CONFIG.getboolean('SOMATICS', 'make_bamout', fallback= False)
        if make_bamout:
            bamout_command = '--bam-output {COHORT_STR}/{tumor_name}_{normal_name}_{num_suffix}.bamout.bam'.format(COHORT_STR= COHORT_STR,
                                                                                                                   tumor_name= tumor_name,
                                                                                                                   normal_name= normal_name,
                                                                                                                   num_suffix= num_suffix)
        else:
            bamout_command = ''

        # {COHORT_STR}/{tumor_name}_{num_suffix}_artifact_prior_table.tsv
        artifact_prior_table = '{COHORT_STR}/{tumor_name}_{num_suffix}_artifact_prior_table.tsv'.format(COHORT_STR= COHORT_STR,
                                                                                                        tumor_name= tumor_name,
                                                                                                        num_suffix= num_suffix)
        if readable(artifact_prior_table):
            apt_command = '--orientation-bias-artifact-priors {}'.format(artifact_prior_table)
        else:
            apt_command = ''


        command_str = '{GATK} --java-options "-Xmx{command_mem}m" Mutect2 -R {REF_FASTA} {tumor_command} {normal_command} {gnoman_command} {pon_command} -L {interval} {gga_command} -O {output_vcf} {bamout_command} {apt_command} {m2_extra_args_command}'.format(GATK= GATK,
                                                                                                                                                                                                                                                                    command_mem= command_mem,
                                                                                                                                                                                                                                                                    REF_FASTA= REF_FASTA,
                                                                                                                                                                                                                                                                    tumor_command= tumor_command,
                                                                                                                                                                                                                                                                    normal_command= normal_command,
                                                                                                                                                                                                                                                                    gnoman_command= gnoman_command,
                                                                                                                                                                                                                                                                    pon_command= pon_command,
                                                                                                                                                                                                                                                                    interval= interval,
                                                                                                                                                                                                                                                                    gga_command= gga_command,
                                                                                                                                                                                                                                                                    output_vcf= output_vcf,
                                                                                                                                                                                                                                                                    bamout_command= bamout_command,
                                                                                                                                                                                                                                                                    apt_command= apt_command,
                                                                                                                                                                                                                                                                    m2_extra_args_command= m2_extra_args_command)
        print(command_str)
        print()
        p = Process(target= common.run, args=(command_str, VERBOSE, proxy_dict, ))
        process_lst.append(p)

    common.multi_run(process_lst, CPU_INT)
    for command_str, r in proxy_dict.items():
        if r[0] != 0:
            raise ExternalToolsError(command_str, r[1])

    #===============================merge the vcf and bamout bam files===============================
    print('\n==================gatk merge_M2_results=========================')
    # merge vcf files
    command_mem = int( (__get_free_mem() - 500) / CPU_INT )
    input_command = ''
    for interval in glob.glob('{COHORT_STR}/interval_files/*.intervals'.format(COHORT_STR= COHORT_STR) ):
        num_suffix = path.split(interval)[1].split('-')[0]
        input_vcf = '{COHORT_STR}/{tumor_name}_{normal_name}_{num_suffix}.vcf.gz'.format(COHORT_STR= COHORT_STR,
                                                                                         tumor_name= tumor_name,
                                                                                         normal_name= normal_name,
                                                                                         num_suffix= num_suffix)
        if readable(input_vcf):
            input_command += ' -I {input_vcf}'.format(input_vcf= input_vcf)

    input_command = input_command.strip()
    output_vcf = '{COHORT_STR}/{tumor_name}_{normal_name}.somatic.vcf.gz'.format(COHORT_STR= COHORT_STR,
                                                                                 tumor_name= tumor_name,
                                                                                 normal_name= normal_name)
    command_str = '{GATK} --java-options "-Xmx{command_mem}m" MergeVcfs {input_command} -O {output_vcf}'.format(GATK= GATK,
                                                                                                                command_mem= command_mem,
                                                                                                                input_command= input_command,
                                                                                                                output_vcf= output_vcf)
    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)
    if r[0] != 0 or not readable(output_vcf):
        raise ExternalToolsError(command_str, r[1])

    # merge and sort bamout files
    input_command = ''
    for interval in glob.glob('{COHORT_STR}/interval_files/*.intervals'.format(COHORT_STR= COHORT_STR) ):
        num_suffix = path.split(interval)[1].split('-')[0]
        input_bam = '{COHORT_STR}/{tumor_name}_{normal_name}_{num_suffix}.bamout.bam'.format(COHORT_STR= COHORT_STR,
                                                                                             tumor_name= tumor_name,
                                                                                             normal_name= normal_name,
                                                                                             num_suffix= num_suffix)
        if readable(input_bam):
            input_command += ' -I {input_bam}'.format(input_bam= input_bam)
        else:
            message = '{input_bam} not available. Skip.'.format(input_bam= input_bam)
            cprint(message)

    input_command = input_command.strip()
    output_bam = '{COHORT_STR}/{tumor_name}_{normal_name}.unsorted.out.bam'.format(COHORT_STR= COHORT_STR,
                                                                                   tumor_name= tumor_name,
                                                                                   normal_name= normal_name)
    command1_str = '{GATK} --java-options "-Xmx{command_mem}m" GatherBamFiles -R {REF_FASTA} {input_command} -O {output_bam}'.format(GATK= GATK,
                                                                                                                                     command_mem= command_mem,
                                                                                                                                     REF_FASTA= REF_FASTA,
                                                                                                                                     input_command= input_command,
                                                                                                                                     output_bam= output_bam)

    final_output_bam = '{tumor_name}_{normal_name}.bamout.bam'.format(tumor_name= tumor_name, normal_name= normal_name)
    command2_str = '{GATK} --java-options "-Xmx{command_mem}m" SortSam --SORT_ORDER coordinate -VALIDATION_STRINGENCY LENIENT -I {input_bam} -O {final_output_bam}'.format(GATK= GATK,
                                                                                                                                                                           command_mem= command_mem,
                                                                                                                                                                           input_bam= output_bam,
                                                                                                                                                                           final_output_bam= final_output_bam)

    command3_str = '{GATK} --java-options "-Xmx{command_mem}m" BuildBamIndex -I {final_output_bam} -VALIDATION_STRINGENCY LENIENT'.format(GATK= GATK,
                                                                                                                                          command_mem= command_mem,
                                                                                                                                          final_output_bam= final_output_bam)

    if make_bamout:
        print(command1_str)
        print()
        r = common.run(command1_str, VERBOSE)
        if r[0] != 0 or not readable(output_bam):
            raise ExternalToolsError(command_str, r[1])

        print(command2_str)
        print()
        r = common.run(command2_str, VERBOSE)
        if r[0] != 0 or not readable(final_output_bam):
            raise ExternalToolsError(command_str, r[1])

        print(command3_str)
        print()
        r = common.run(command3_str, VERBOSE)
        if r[0] != 0:
            raise ExternalToolsError(command_str, r[1])


    # ===============================calculate_contamination===============================
    contamination_vcf = CONFIG.get('RESOURCES', 'variants_for_contamination')
    if contamination_vcf != '' and readable(contamination_vcf):
        __calculate_contamination(tumor_name, normal_name)
    else:
        print('Contamination vcf file not provided or not found. Check [RESOURCES] section, item "variants_for_contamination" in configuration file.')
        print('Skip task __calculate_contamination.' )

    #===============================FilterMutectCalls===============================
    command_mem = int( (__get_free_mem() - 500) / CPU_INT )
    input_vcf = '{COHORT_STR}/{tumor_name}_{normal_name}.somatic.vcf.gz'.format(COHORT_STR= COHORT_STR,
                                                                                tumor_name= tumor_name,
                                                                                normal_name= normal_name)

    output_vcf = '{COHORT_STR}/{tumor_name}_{normal_name}.somatic.filtered.vcf.gz'.format(COHORT_STR= COHORT_STR,
                                                                                          tumor_name= tumor_name,
                                                                                          normal_name= normal_name)

    contamination_table = '{COHORT_STR}/{tumor_name}_{normal_name}.contamination.table'.format(COHORT_STR= COHORT_STR,
                                                                                               tumor_name= tumor_name,
                                                                                               normal_name= normal_name)

    maf_segments = '{COHORT_STR}/{tumor_name}_{normal_name}.maf_segments.table'.format(COHORT_STR= COHORT_STR,
                                                                                       tumor_name= tumor_name,
                                                                                       normal_name= normal_name)

    if readable(contamination_table) and readable(maf_segments):
        contamination_filter_command = '--contamination-table {contamination_table} --tumor-segmentation {maf_segments}'.format(contamination_table= contamination_table, maf_segments= maf_segments)
    else:
        contamination_filter_command = ''

    m2_extra_filtering_args = CONFIG.get('ARGS', 'm2_extra_filtering_args')
    command_str = '{GATK} --java-options "-Xmx{command_mem}m" FilterMutectCalls {contamination_filter_command} -V {input_vcf} -O {output_vcf} {m2_extra_filtering_args}'.format(GATK= GATK,
                                                                                                                                                                                command_mem= command_mem,
                                                                                                                                                                                contamination_filter_command = contamination_filter_command,
                                                                                                                                                                                input_vcf= input_vcf,
                                                                                                                                                                                output_vcf= output_vcf,
                                                                                                                                                                                m2_extra_filtering_args= m2_extra_filtering_args)

    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)
    if r[0] != 0 or not readable(output_vcf):
        raise ExternalToolsError(command_str, r[1])

    # ob filters and aa filters
    input_vcf1 = '{COHORT_STR}/{tumor_name}_{normal_name}.somatic.ob_filtered.vcf.gz'.format(COHORT_STR= COHORT_STR, tumor_name= tumor_name, normal_name= normal_name) # FilterByOrientationBias output
    input_vcf2 = '{COHORT_STR}/{tumor_name}_{normal_name}.somatic.filtered.vcf.gz'.format(COHORT_STR= COHORT_STR, tumor_name= tumor_name, normal_name= normal_name)   # FilterMutectCalls output
    output_vcf = '{tumor_name}_{normal_name}.final.vcf.gz'.format(tumor_name= tumor_name, normal_name= normal_name)

    run_ob_filter = CONFIG.getboolean('SOMATICS', 'run_ob_filter', fallback= False)
    run_aa_filter = CONFIG.getboolean('SOMATICS', 'run_aa_filter', fallback= True)
    if not run_ob_filter and not run_aa_filter:
        shutil.move(input_vcf2, output_vcf)

    if run_ob_filter:
        __filter_by_orientation_bias(tumor_name, normal_name)
        if not run_aa_filter:
            shutil.move(input_vcf1, output_vcf)

    if run_aa_filter:
        l = Lock()
        __filter_alignment_artifacts(tumor_name, normal_name, l)

    return

def CreatePoN(normal_name_list):
    '''
    Make a panel of normals for use with Mutect2

    Parameters:
        **normal_name_list**: string
            each item represents {normal_name}.analysis_ready.bam

    Returns:
        **value**: type
            its meaning
    input_vcf1 = '{COHORT_STR}/{tumor_name}_{normal_name}.somatic.ob_filtered.vcf.gz'.format(COHORT_STR= COHORT_STR, tumor_name= tumor_name, normal_name= normal_name) # FilterByOrientationBias output
    input_vcf2 = '{COHORT_STR}/{tumor_name}_{normal_name}.somatic.filtered.vcf.gz'.format(COHORT_STR= COHORT_STR, tumor_name= tumor_name, normal_name= normal_name)   # FilterMutectCalls output
    output_vcf = '{tumor_name}_{normal_name}.final.vcf.gz'.format(tumor_name= tumor_name, normal_name= normal_name)
    '''
    process_lst = []
    for normal_name in normal_name_list:
        p = Process(target=CallSomaticMutect2, args=(normal_name, None, ) )
        process_lst.append(p)

    common.multi_run(process_lst, CPU_INT)

    #============================CreateSomaticPanelOfNormals============================
    print('\n==================gatk CreateSomaticPanelOfNormals=========================')
    # set input vcf files -vcfs
    input_vcf_command = ''
    if len(normal_name_list) >= 5:
        pon_list_file_f = open('normals_for_pon_vcf.list'.format(COHORT_STR= COHORT_STR), 'w')
        for normal_name in normal_name_list:
            pon_list_file_f.writelines('{tumor_name}_None.final.vcf.gz'.format(COHORT_STR= COHORT_STR, tumor_name= normal_name) + '\n')
        pon_list_file_f.closed()
        input_vcf_command = '-vcfs normals_for_pon_vcf.list'.format(COHORT_STR= COHORT_STR)
    else:
        for normal_name in normal_name_list:
            input_vcf_command += ' -vcfs {tumor_name}_None.final.vcf.gz'.format(COHORT_STR= COHORT_STR, tumor_name= normal_name)
        input_vcf_command = input_vcf_command.strip()

    # set --duplicate_sample_strategy
    duplicate_sample_strategy = CONFIG.get('SOMATICS', 'duplicate_sample_strategy')
    if duplicate_sample_strategy not in ['THROW_ERROR', 'CHOOSE_FIRST', 'ALLOW_ALL']:
        message = "duplicate_sample_strategy must be one of ['THROW_ERROR', 'CHOOSE_FIRST', 'ALLOW_ALL'], current value is {}. Change to 'CHOOSE_FIRST'.".format(duplicate_sample_strategy)
        print(message)
        duplicate_sample_strategy = 'CHOOSE_FIRST'

    command_str = '{GATK} CreateSomaticPanelOfNormals {input_vcf_command} --duplicate-sample-strategy {duplicate_sample_strategy} -O {PON_STR}'.format(GATK= GATK,
                                                                                                                                                       input_vcf_command= input_vcf_command,
                                                                                                                                                       duplicate_sample_strategy= duplicate_sample_strategy,
                                                                                                                                                       PON_STR= PON_STR)
    print(command_str)
    print()
    r = common.run(command_str, VERBOSE)
    if r[0] != 0 or not readable(PON_STR):
        raise ExternalToolsError(command_str, r[1])

    return

def main(argvList = sys.argv, argv_int = len(sys.argv)):

    #==============
    global CPU_INT
    if not path.isdir(COHORT_STR):
        os.mkdir(COHORT_STR)

    command_str = 'cat /proc/cpuinfo | grep "processor"'
    r = common.run(command_str, False)
    if r[0] == 0:
        CPU_INT = int( len(r[1].split('\n')) * 0.75 )
        CPU_INT = max( [1, CPU_INT] )
        message = '{} CPUs used.'.format(CPU_INT)
    else:
        CPU_INT = 1
        message = 'Fail to find CPU count, default 1'
    print(message)

    #======================check==========================
    Check_ini()
    #======================preprocess=====================
    process_lst = []
    if 'preprocessing' in TASK_LST:

        CreateSequenceGroupingTSV()
        for sample_name in SAMPLE_LST:
            p = Process(target=Preprocess, args=(sample_name, ) )
            process_lst.append(p)

        common.multi_run(process_lst, CPU_INT)

    #======================Call TMB======================
    if 'TMB' in TASK_LST:
        for sample_name in SAMPLE_LST:
            p = Process(target=calc_TMB, args=(sample_name, ) )
            process_lst.append(p)

        common.multi_run(process_lst, CPU_INT)

    #================germline mutation calling============
    process_lst = []
    if 'germline' in TASK_LST:

        for sample_name in SAMPLE_LST:
            p = Process(target=CallGermlineVarscan, args=(sample_name, ) )
            process_lst.append(p)

        common.multi_run(process_lst, CPU_INT)

        # check the calling results {sample_name}.annotated.txt
        for sample_name in SAMPLE_LST:
            temp_str = '{sample_name}.annotated.txt'.format(sample_name= sample_name)
            if readable(temp_str) and path.getsize(temp_str) > 0:
                next
            else:
                message = '{} does not exist.'.format(temp_str)
                cprint(message)



    #===============Somatic calling================
    global SAMPLE_ANOUMT_INT
    scatter_amount = CPU_INT / SAMPLE_ANOUMT_INT
    scatter_amount = int(scatter_amount) if scatter_amount >=1 else 1

    if 'somatic' in TASK_LST:
        if not readable(PON_STR) and CONTROL_AMOUNT_INT >= 2:
            SplitIntervals(scatter_amount)
            CreatePoN(list(set(CONTROL_LST)))

        if CASE_AMOUNT_INT == 0:
            message = 'No case samples provided, create Panel of Normal (PoN) only.\n{PON_STR}'.format(PON_STR= PON_STR)
            print(message)
            sys.exit()

        if CASE_AMOUNT_INT != 0 and CONTROL_AMOUNT_INT == 0:
            SplitIntervals(scatter_amount)
            process_lst = []
            for sample_name in CASE_LST:
                # TODO: check the {case_name}.analysis_ready.bam is available.
                p = Process(target=CallSomaticMutect2, args=(sample_name, ) )
                process_lst.append(p)

            common.multi_run(process_lst, CPU_INT)

        if CASE_AMOUNT_INT != 0 and CONTROL_AMOUNT_INT != 0:
            SAMPLE_ANOUMT_INT = int( SAMPLE_ANOUMT_INT/2 )
            scatter_amount = CPU_INT / SAMPLE_ANOUMT_INT
            scatter_amount = int(scatter_amount) if scatter_amount >=1 else 1
            SplitIntervals(scatter_amount)

            process_lst = []
            for i in range(SAMPLE_ANOUMT_INT):
                # TODO: check the {case_name}.analysis_ready.bam and {control_name}.analysis_ready.bam are both available.
                p = Process( target=CallSomaticMutect2, args=(CASE_LST[i], CONTROL_LST[i], ) )
                process_lst.append(p)

            common.multi_run(process_lst, CPU_INT)

    print()
    print('done')

    return

main()

